#!/usr/bin/env python3
"""
scanpy_qc.py — Single-cell QC and preprocessing pipeline
Input:  filtered count matrix (.h5ad)
Output: clean, normalized, annotated .h5ad ready for atlas ingestion
"""

import argparse
import logging
import scanpy as sc
import scrublet as scr
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
log = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input',        required=True,  help='Input .h5ad file')
    parser.add_argument('--output',       required=True,  help='Output .h5ad file')
    parser.add_argument('--sample-id',    required=True,  help='Sample identifier')
    parser.add_argument('--min-genes',    type=int, default=200,   help='Min genes per cell')
    parser.add_argument('--max-genes',    type=int, default=6000,  help='Max genes per cell')
    parser.add_argument('--max-mito',     type=float, default=20.0, help='Max mitochondrial %')
    parser.add_argument('--min-counts',   type=int, default=500,   help='Min UMI counts per cell')
    parser.add_argument('--n-hvgs',       type=int, default=2000,  help='Number of highly variable genes')
    parser.add_argument('--n-pcs',        type=int, default=50,    help='Number of PCs')
    parser.add_argument('--n-neighbors',  type=int, default=15,    help='Neighbors for UMAP')
    return parser.parse_args()

def main():
    args = parse_args()
    sc.settings.verbosity = 2

    # ── 1. Load ──────────────────────────────────────────────────────────────
    log.info(f"Loading {args.input}")
    adata = sc.read_h5ad(args.input)
    adata.var_names_make_unique()
    log.info(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    # ── 2. QC metrics ────────────────────────────────────────────────────────
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True
    )

    # Log pre-filter stats
    log.info(f"Pre-filter: {adata.n_obs} cells")
    log.info(f"  median genes/cell: {np.median(adata.obs.n_genes_by_counts):.0f}")
    log.info(f"  median UMIs/cell:  {np.median(adata.obs.total_counts):.0f}")
    log.info(f"  median mito %:     {np.median(adata.obs.pct_counts_mt):.1f}%")

    # ── 3. Cell filtering ────────────────────────────────────────────────────
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_cells(adata, max_genes=args.max_genes)
    sc.pp.filter_genes(adata, min_cells=3)

    # Mito and count filters
    adata = adata[adata.obs.pct_counts_mt < args.max_mito].copy()
    adata = adata[adata.obs.total_counts > args.min_counts].copy()
    log.info(f"Post-filter: {adata.n_obs} cells remaining")

    # ── 4. Doublet detection (Scrublet) ──────────────────────────────────────
    log.info("Running Scrublet doublet detection")
    scrub = scr.Scrublet(adata.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
    adata.obs['doublet_score'] = doublet_scores
    adata.obs['predicted_doublet'] = predicted_doublets

    n_doublets = predicted_doublets.sum()
    log.info(f"Doublets detected: {n_doublets} ({100*n_doublets/adata.n_obs:.1f}%)")
    adata = adata[~adata.obs['predicted_doublet']].copy()
    log.info(f"Post-doublet removal: {adata.n_obs} cells")

    # ── 5. Normalization ─────────────────────────────────────────────────────
    log.info("Normalizing")
    # Store raw counts before any normalization
    adata.layers['counts'] = adata.X.copy()

    # Normalize and log transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Store log-normalized counts as a layer too
    adata.layers['log_norm'] = adata.X.copy()

    # ── 6. Highly variable genes ─────────────────────────────────────────────
    sc.pp.highly_variable_genes(
    adata, n_top_genes=args.n_hvgs, subset=True
    )
    log.info(f"Selected {adata.n_vars} highly variable genes")

    # ── 7. PCA + neighbors + UMAP ────────────────────────────────────────────
    log.info("Running PCA, neighbors, UMAP")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=args.n_pcs)
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5)

    # ── 8. Basic cell type annotation ────────────────────────────────────────
    log.info("Annotating cell types with marker genes")
    marker_genes = {
        'T cell':      ['CD3D', 'CD3E', 'CD3G'],
        'B cell':      ['CD79A', 'CD79B', 'MS4A1'],
        'NK cell':     ['GNLY', 'NKG7', 'KLRD1'],
        'Monocyte':    ['CD14', 'LYZ', 'CST3'],
        'DC':          ['FCER1A', 'CST3'],
        'Platelet':    ['PPBP', 'PF4'],
    }
    sc.tl.score_genes_cell_cycle(
        adata,
        s_genes=['MCM5','PCNA','TYMS','FEN1','MCM2'],
        g2m_genes=['HMGB2','CDK1','NUSAP1','UBE2C','BIRC5']
    )
    for cell_type, markers in marker_genes.items():
        valid = [g for g in markers if g in adata.var_names]
        if valid:
            sc.tl.score_genes(adata, valid, score_name=f'score_{cell_type.replace(" ", "_")}')

    # ── 9. Metadata ──────────────────────────────────────────────────────────
    adata.obs['sample_id'] = args.sample_id
    adata.uns['qc_params'] = vars(args)

    # ── 10. Save ─────────────────────────────────────────────────────────────
    log.info(f"Saving to {args.output}")
    adata.write_h5ad(args.output, compression='gzip')
    log.info("Done")

if __name__ == '__main__':
    main()