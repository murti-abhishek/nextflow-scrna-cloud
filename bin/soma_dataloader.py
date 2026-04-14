#!/usr/bin/env python3
"""
soma_dataloader.py — PyTorch DataLoader streaming from TileDB-SOMA atlas
Produces (input_ids, attention_mask) tensors compatible with LiverTransformer

Replaces: loading precomputed human_liver_tokens.npz
With:     streaming raw counts from SOMA → rank tokenize on the fly
"""

import json
import numpy as np
import torch
from torch.utils.data import DataLoader, IterableDataset
import tiledbsoma as soma

class SOMAIterableDataset(IterableDataset):
    """
    Streams cells from a TileDB-SOMA experiment on S3.
    Applies rank-based tokenization matching LiverTransformer's scheme.
    """

    def __init__(
        self,
        soma_uri: str,
        vocab_path: str,
        seq_len: int = 1024,
        layer: str = "counts",
        chunk_size: int = 512,
        value_filter: str = None,
    ):
        """
        Args:
            soma_uri:     S3 URI of the SOMA experiment
            vocab_path:   path to vocab_v05.json (gene → token_id mapping)
            seq_len:      max sequence length (default 1024, matches LiverTransformer)
            layer:        which X layer to read ('counts' = raw integer counts)
            chunk_size:   cells per chunk streamed from S3
            value_filter: optional SOMA query filter e.g. "tissue == 'blood'"
        """
        self.soma_uri    = soma_uri
        self.seq_len     = seq_len
        self.layer       = layer
        self.chunk_size  = chunk_size
        self.value_filter = value_filter

        # Load vocab: gene_name → token_id
        with open(vocab_path) as f:
            self.vocab = json.load(f)

        # Will be set when we open the experiment
        self._gene_to_token = None

    def _build_gene_index(self, var_df):
        gene_names = None
        for col in ['gene_symbol', 'gene_name', 'feature_name']:
            if col in var_df.columns:
                gene_names = var_df[col].values
                break
        if gene_names is None:
            gene_names = var_df.index.values

        token_ids = np.array([
            self.vocab.get(g, -1) for g in gene_names
        ], dtype=np.int32)
        return token_ids

    def _rank_tokenize(self, counts_row, gene_token_ids):
        """
        Rank-based tokenization matching LiverTransformer's scheme:
        1. Filter to genes in vocab
        2. Rank by expression descending
        3. Take top seq_len
        4. Return token ID sequence + attention mask

        Args:
            counts_row:    1D array of raw counts for one cell (n_genes,)
            gene_token_ids: token ID per gene position (n_genes,)

        Returns:
            input_ids:      (seq_len,) int64 tensor
            attention_mask: (seq_len,) int64 tensor, 1=real 0=pad
        """
        # Only keep genes in vocab
        valid_mask = gene_token_ids >= 0
        valid_counts = counts_row[valid_mask]
        valid_tokens = gene_token_ids[valid_mask]

        # Rank by expression descending (highest expressed first)
        rank_order = np.argsort(-valid_counts)
        ranked_tokens = valid_tokens[rank_order]
        ranked_counts = valid_counts[rank_order]

        # Filter out zero-count genes
        nonzero = ranked_counts > 0
        ranked_tokens = ranked_tokens[nonzero]

        # Truncate to seq_len
        n = min(len(ranked_tokens), self.seq_len)
        input_ids = np.zeros(self.seq_len, dtype=np.int64)
        attention_mask = np.zeros(self.seq_len, dtype=np.int64)

        input_ids[:n] = ranked_tokens[:n]
        attention_mask[:n] = 1

        return (
            torch.tensor(input_ids, dtype=torch.long),
            torch.tensor(attention_mask, dtype=torch.long),
        )

    def __iter__(self):
        with soma.Experiment.open(self.soma_uri) as exp:
            # Build gene → token mapping once
            var_df = exp.ms["RNA"].var.read().concat().to_pandas()
            gene_token_ids = self._build_gene_index(var_df)

            # Query with optional filter
            query = exp.axis_query(
                "RNA",
                obs_query=soma.AxisQuery(value_filter=self.value_filter)
                if self.value_filter else soma.AxisQuery(),
            )

            # Stream X[counts] in chunks
            for chunk in query.X(self.layer).tables():
                # chunk is a PyArrow table with soma_dim_0, soma_dim_1, soma_data
                df = chunk.to_pandas()

                # Pivot: group by cell (soma_dim_0), build dense count vector
                cell_ids = df["soma_dim_0"].unique()

                for cell_id in cell_ids:
                    cell_df = df[df["soma_dim_0"] == cell_id]

                    # Build sparse → dense count vector
                    counts = np.zeros(len(var_df), dtype=np.float32)
                    counts[cell_df["soma_dim_1"].values] = \
                        cell_df["soma_data"].values

                    input_ids, attention_mask = self._rank_tokenize(
                        counts, gene_token_ids
                    )
                    yield input_ids, attention_mask


def make_soma_dataloader(
    soma_uri: str,
    vocab_path: str,
    batch_size: int = 256,
    seq_len: int = 1024,
    layer: str = "counts",
    value_filter: str = None,
    num_workers: int = 0,
) -> DataLoader:
    """
    Drop-in replacement for loading human_liver_tokens.npz.
    Returns a DataLoader yielding (input_ids, attention_mask) batches.

    Usage:
        loader = make_soma_dataloader(
            soma_uri  = "s3://bucket/atlas/pbmc_1k",
            vocab_path = "vocab_v05.json",
            batch_size = 256,
        )
        for input_ids, attention_mask in loader:
            logits = model(input_ids, attention_mask)
    """
    dataset = SOMAIterableDataset(
        soma_uri=soma_uri,
        vocab_path=vocab_path,
        seq_len=seq_len,
        layer=layer,
        value_filter=value_filter,
    )
    return DataLoader(
        dataset,
        batch_size=batch_size,
        num_workers=num_workers,
        pin_memory=torch.cuda.is_available(),
    )

def make_perturb_soma_dataloader(
    soma_uri: str,
    vocab_path: str,
    target_gene: str = None,
    controls_only: bool = False,
    batch_size: int = 256,
    seq_len: int = 1024,
    layer: str = "counts",
    num_workers: int = 0,
) -> DataLoader:
    """
    Perturbation-aware DataLoader. Streams cells from a Perturb-seq SOMA atlas
    filtered by perturbation target — for conditioning models on perturbation identity.

    Usage:
        # All RPL3 knockdown cells
        loader = make_perturb_soma_dataloader(
            soma_uri   = "s3://bucket/atlas/K562_essential",
            vocab_path = "data/vocab_v05.json",
            target_gene = "RPL3",
        )

        # All non-targeting control cells
        loader = make_perturb_soma_dataloader(
            soma_uri      = "s3://bucket/atlas/K562_essential",
            vocab_path    = "data/vocab_v05.json",
            controls_only = True,
        )
    """
    if controls_only:
        value_filter = "is_control == True"
    elif target_gene:
        value_filter = f"target_gene == '{target_gene}'"
    else:
        value_filter = None

    return make_soma_dataloader(
        soma_uri     = soma_uri,
        vocab_path   = vocab_path,
        batch_size   = batch_size,
        seq_len      = seq_len,
        layer        = layer,
        value_filter = value_filter,
        num_workers  = num_workers,
    )

# ── Quick smoke test ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--soma-uri",    required=True)
    parser.add_argument("--vocab-path",  required=True)
    parser.add_argument("--batch-size",  type=int, default=32)
    parser.add_argument("--seq-len",     type=int, default=1024)
    parser.add_argument("--n-batches",   type=int, default=3)
    parser.add_argument("--target-gene", type=str, default=None,
                        help="Filter to specific perturbation target")
    parser.add_argument("--controls-only", action="store_true",
                        help="Only stream non-targeting control cells")
    args = parser.parse_args()

    print(f"Loading SOMA atlas from {args.soma_uri}")

    if args.target_gene or args.controls_only:
        print(f"Filter: {'controls only' if args.controls_only else f'target_gene = {args.target_gene}'}")
        loader = make_perturb_soma_dataloader(
            soma_uri      = args.soma_uri,
            vocab_path    = args.vocab_path,
            target_gene   = args.target_gene,
            controls_only = args.controls_only,
            batch_size    = args.batch_size,
            seq_len       = args.seq_len,
        )
    else:
        loader = make_soma_dataloader(
            soma_uri   = args.soma_uri,
            vocab_path = args.vocab_path,
            batch_size = args.batch_size,
            seq_len    = args.seq_len,
        )

    for i, (input_ids, attention_mask) in enumerate(loader):
        if i >= args.n_batches:
            break
        print(f"Batch {i+1}:")
        print(f"  input_ids shape:      {input_ids.shape}")
        print(f"  attention_mask shape: {attention_mask.shape}")
        print(f"  tokens per cell:      {attention_mask.sum(dim=1).float().mean():.1f} avg")
        print(f"  sample tokens:        {input_ids[0, :8]}")

    print("\nDone.")