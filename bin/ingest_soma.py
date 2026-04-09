#!/usr/bin/env python3
"""
ingest_soma.py — Ingest a cleaned H5AD file into a TileDB-SOMA collection
Input:  clean .h5ad file (output of scanpy_qc.py)
Output: TileDB-SOMA experiment written to S3
"""

import argparse
import logging
import anndata as ad
import tiledbsoma as soma
import tiledbsoma.io as soma_io

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
log = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input',       required=True, help='Input clean .h5ad file')
    parser.add_argument('--uri',         required=True, help='S3 URI for SOMA collection')
    parser.add_argument('--sample-id',   required=True, help='Sample identifier')
    parser.add_argument('--tissue',      default='blood',   help='Tissue type')
    parser.add_argument('--organism',    default='Homo sapiens', help='Organism')
    parser.add_argument('--assay',       default='10x 3v3', help='Assay type')
    return parser.parse_args()

def validate_metadata(adata, sample_id, tissue, organism):
    """Enforce required metadata fields before ingestion."""
    log.info("Validating metadata")

    # Add required fields if missing
    adata.obs['sample_id']  = sample_id
    adata.obs['tissue']     = tissue
    adata.obs['organism']   = organism

    # Ensure string types for categorical fields
    for col in ['sample_id', 'tissue', 'organism']:
        adata.obs[col] = adata.obs[col].astype(str)

    # Validate gene IDs exist
    assert adata.n_vars > 0, "No genes found"
    assert adata.n_obs  > 0, "No cells found"

    log.info(f"Validated: {adata.n_obs} cells × {adata.n_vars} genes")
    log.info(f"  tissue:   {tissue}")
    log.info(f"  organism: {organism}")
    log.info(f"  assay:    {adata.uns.get('assay', 'unknown')}")
    return adata

def main():
    args = parse_args()

    # ── 1. Load ──────────────────────────────────────────────────────────────
    log.info(f"Loading {args.input}")
    adata = ad.read_h5ad(args.input)
    log.info(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    # ── 2. Validate metadata ─────────────────────────────────────────────────
    adata = validate_metadata(adata, args.sample_id, args.tissue, args.organism)

    # ── 3. Set experiment URI inside the collection ──────────────────────────
    experiment_uri = f"{args.uri}/{args.sample_id}"
    log.info(f"Writing SOMA experiment to {experiment_uri}")

    # ── 4. Ingest ────────────────────────────────────────────────────────────
    soma_io.from_anndata(
        experiment_uri,
        adata,
        measurement_name="RNA",
    )
    log.info("Ingestion complete")

    # ── 5. Verify ────────────────────────────────────────────────────────────
    log.info("Verifying SOMA experiment is queryable")
    with soma.Experiment.open(experiment_uri) as exp:
        obs_count = len(exp.obs.read().concat())
        log.info(f"SOMA experiment contains {obs_count} cells — OK")

if __name__ == '__main__':
    main()