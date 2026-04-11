#!/usr/bin/env python3
"""
perturb_qc.py — Perturbation-aware QC and preprocessing for Perturb-seq data
Input:  raw single-cell H5AD with perturbation identity already in obs
Output: clean H5AD with perturbation metadata, ready for SOMA ingestion
"""
import sys
from unittest.mock import MagicMock
sys.modules['numba'] = MagicMock()

import argparse
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
log = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input',                   required=True)
    parser.add_argument('--output',                  required=True)
    parser.add_argument('--sample-id',               required=True)
    parser.add_argument('--n-top-perturbs',          type=int, default=50)
    parser.add_argument('--min-cells-per-perturb',   type=int, default=20)
    parser.add_argument('--min-genes',               type=int, default=200)
    parser.add_argument('--max-genes',               type=int, default=8000)
    parser.add_argument('--max-mito',                type=float, default=25.0)
    parser.add_argument('--min-counts',              type=int, default=500)
    parser.add_argument('--n-hvgs',                  type=int, default=2000)
    return parser.parse_args()

def parse_perturbation_metadata(obs):
    """
    Extract perturbation metadata from Replogle 2022 single-cell obs.
    Guide identity is already assigned — just clean and standardize.
    """
    meta = pd.DataFrame(index=obs.index)
    meta['target_gene']     = obs['gene']
    meta['guide_id']        = obs['transcript'] if 'transcript' in obs.columns else 'unknown'
    meta['perturbation_id'] = obs['gene_transcript'] if 'gene_transcript' in obs.columns else obs['gene']
    meta['is_control']      = obs['gene'] == 'non-targeting'
    return meta

def main():
    args = parse_args()
    sc.settings.verbosity = 2

    # ── 1. Load ──────────────────────────────────────────────────────────────
    log.info(f"Loading {args.input}")
    adata = sc.read_h5ad(args.input)
    adata.var_names_make_unique()
    log.info(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    # ── 2. Parse perturbation metadata ───────────────────────────────────────
    log.info("Parsing perturbation metadata")
    perturb_meta = parse_perturbation_metadata(adata.obs)
    for col in ['target_gene', 'guide_id', 'perturbation_id', 'is_control']:
        adata.obs[col] = perturb_meta[col].values

    n_controls = adata.obs['is_control'].sum()
    n_perturbed = (~adata.obs['is_control']).sum()
    n_targets = adata.obs[~adata.obs['is_control']]['target_gene'].nunique()
    log.info(f"Control cells: {n_controls}")
    log.info(f"Perturbed cells: {n_perturbed} across {n_targets} unique targets")

    # ── 3. Subset to manageable size ─────────────────────────────────────────
    log.info(f"Subsetting to top {args.n_top_perturbs} perturbations + all controls")
    perturb_counts = (
        adata.obs[~adata.obs['is_control']]
        .groupby('target_gene')
        .size()
        .sort_values(ascending=False)
        .head(args.n_top_perturbs)
    )
    top_targets = set(perturb_counts.index)
    keep_mask = adata.obs['is_control'] | adata.obs['target_gene'].isin(top_targets)
    adata = adata[keep_mask].copy()
    log.info(f"After subset: {adata.n_obs} cells")

    # ── 4. QC metrics ────────────────────────────────────────────────────────
    log.info("Computing QC metrics")
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True
    )
    log.info(f"Pre-filter: {adata.n_obs} cells")
    log.info(f"  median genes/cell: {np.median(adata.obs.n_genes_by_counts):.0f}")
    log.info(f"  median UMIs/cell:  {np.median(adata.obs.total_counts):.0f}")
    log.info(f"  median mito %:     {np.median(adata.obs.pct_counts_mt):.1f}%")

    # ── 5. Cell filtering ────────────────────────────────────────────────────
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_cells(adata, max_genes=args.max_genes)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs.pct_counts_mt < args.max_mito].copy()
    adata = adata[adata.obs.total_counts > args.min_counts].copy()
    log.info(f"Post-filter: {adata.n_obs} cells")

    # ── 6. Perturbation-aware QC ─────────────────────────────────────────────
    log.info("Running perturbation-aware QC")
    cells_per_perturb = adata.obs.groupby('target_gene').size()
    valid_targets = cells_per_perturb[
        cells_per_perturb >= args.min_cells_per_perturb
    ].index
    n_dropped = (~adata.obs['target_gene'].isin(valid_targets)).sum()
    adata = adata[adata.obs['target_gene'].isin(valid_targets)].copy()
    log.info(f"Dropped {n_dropped} cells from low-coverage perturbations")
    log.info(f"Remaining targets: {adata.obs['target_gene'].nunique()}")
    log.info(f"Final cell count: {adata.n_obs}")

    # ── 7. Store raw counts ───────────────────────────────────────────────────
    adata.layers['counts'] = adata.X.copy()

    # ── 8. Normalize to median library size ──────────────────────────────────
    median_lib_size = float(np.median(adata.obs.total_counts))
    log.info(f"Normalizing to median library size: {median_lib_size:.0f} UMIs")
    sc.pp.normalize_total(adata, target_sum=median_lib_size)
    sc.pp.log1p(adata)
    adata.layers['log_norm'] = adata.X.copy()
    adata.uns['normalization'] = {
        'method': 'median_library_size',
        'target_sum': median_lib_size,
    }

    # ── 9. HVGs + PCA + UMAP ─────────────────────────────────────────────────
    log.info("HVG selection, PCA, UMAP")
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvgs)
    # adata = adata[:, adata.var['highly_variable']].copy()
    # sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50)
    # sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
    # sc.tl.umap(adata)
    # sc.tl.leiden(adata, resolution=0.5)

    # ── 10. Perturbation summary stats ────────────────────────────────────────
    log.info("Computing per-perturbation summary stats")
    summary = adata.obs.groupby('target_gene').agg(
        n_cells=('target_gene', 'count'),
        mean_genes=('n_genes_by_counts', 'mean'),
        mean_umis=('total_counts', 'mean'),
        mean_mito=('pct_counts_mt', 'mean'),
    ).reset_index()
    adata.uns['perturbation_summary'] = summary.to_dict(orient='list')

    # ── 11. Metadata ──────────────────────────────────────────────────────────
    adata.obs['sample_id'] = args.sample_id
    adata.obs['tissue']    = 'blood'
    adata.obs['organism']  = 'Homo sapiens'
    adata.obs['assay']     = '10x 3v3 Perturb-seq'
    adata.obs['cell_line'] = 'K562'
    adata.obs['dataset']   = 'Replogle_2022_essential'
    adata.uns['qc_params'] = vars(args)

    # ── 12. Save ──────────────────────────────────────────────────────────────
    log.info(f"Saving to {args.output}")
    adata.write_h5ad(args.output, compression='gzip')
    log.info("Done")

    print(f"\n{'='*50}")
    print(f"Perturbation QC complete")
    print(f"  Cells:                {adata.n_obs}")
    print(f"  Genes:                {adata.n_vars}")
    print(f"  Perturbation targets: {adata.obs['target_gene'].nunique()}")
    print(f"  Control cells:        {adata.obs['is_control'].sum()}")
    print(f"  Median lib size:      {median_lib_size:.0f}")
    print(f"  Output:               {args.output}")
    print(f"{'='*50}")

if __name__ == '__main__':
    main()