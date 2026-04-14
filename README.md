# nextflow-scrna-cloud

Cloud-scale single-cell and Perturb-seq data engineering pipeline with Nextflow DSL2, Docker, and AWS Batch. Produces a TileDB-SOMA atlas queryable by cell type, tissue, and perturbation target — directly feeding PyTorch training pipelines for biological foundation models.

Built as a portfolio project bridging HPC-based bioinformatics (UCSF SLURM) to cloud-native data engineering — a core requirement for virtual cell and computational biology roles.

---

## What this pipeline does

```
FASTQs (scRNA-seq or Perturb-seq)
    │
    ▼  nf-core/scrnaseq (STARsolo, AWS Batch)
    │
Count matrix (.h5ad)
    │
    ▼  Scanpy QC module (AWS Batch)
    │  ambient RNA · doublets · normalization · HVGs · PCA
    │
Clean annotated .h5ad
    │
    ▼  SOMA ingestion module (AWS Batch)
    │  schema validation · ontology · TileDB-SOMA write
    │
TileDB-SOMA atlas (S3)
    │
    ▼  SOMAIterableDataset (PyTorch)
       rank tokenization · batch streaming · perturbation filtering
       → LiverTransformer / any foundation model training loop
```

For Perturb-seq data, a dedicated `PERTURB_QC` module handles perturbation-aware cell filtering and ingests the atlas with `target_gene` metadata, enabling queries like "give me all RPL3 knockdown cells."

---

## Architecture

```
Your Mac (control plane only)
    │
    ├── nextflow run main.nf -profile aws
    ├── Seqera Platform monitoring (Tower)
    │
    ▼
AWS Batch
    │   nextflow-batch-queue
    │   nextflow-compute-env (spot · x86 · max 256 vCPU)
    │
    ├──▶ EC2 Spot — nf-core/scrnaseq     (r6i.4xlarge, STARsolo)
    ├──▶ EC2 Spot — SCANPY_QC            (m6i, scanpy_qc.py)
    ├──▶ EC2 Spot — PERTURB_QC           (m6i, perturb_qc.py)
    └──▶ EC2 Spot — INGEST_SOMA          (m6i, ingest_soma.py)
         └── all write to S3 via Fusion filesystem

S3 (nextflow-scrna-abhishek)
    ├── raw/              input FASTQs + H5ADs
    ├── work/             Nextflow intermediate files
    ├── results/          pipeline outputs
    ├── indexes/          cached STAR genome index (GRCh38)
    └── atlas/
        ├── pbmc_1k/      PBMC 1k v3 — 1,092 cells, blood
        └── K562_essential/  Replogle 2022 Perturb-seq — 39K cells, 51 perturbation targets
```

---

## Stack

| Layer | Tool |
|---|---|
| Workflow orchestration | Nextflow DSL2 |
| Alignment | nf-core/scrnaseq, STARsolo |
| Single-cell QC | Scanpy, Scrublet |
| Perturb-seq QC | Custom perturbation-aware pipeline |
| Cloud compute | AWS Batch + EC2 Spot |
| Storage | Amazon S3 |
| S3 filesystem | Fusion (Wave/Seqera) |
| Container registry | Docker Hub, AWS ECR |
| Atlas format | TileDB-SOMA |
| ML DataLoader | PyTorch IterableDataset |
| Monitoring | Seqera Platform (Tower) |

---

## Repository structure

```
.
├── main.nf                      # workflow — chains all modules
├── nextflow.config              # params, local + aws profiles
├── conf/
│   ├── aws_batch.config         # Batch queue, Fusion, Wave, Tower
│   └── local.config
├── modules/
│   ├── fastqc/main.nf           # FastQC process
│   ├── scanpy_qc/main.nf        # Scanpy QC process
│   ├── perturb_qc/main.nf       # Perturb-seq QC process
│   └── ingest_soma/main.nf      # TileDB-SOMA ingestion process
├── bin/
│   ├── scanpy_qc.py             # QC script: H5AD → clean H5AD
│   ├── perturb_qc.py            # Perturb-seq QC: guide parsing, perturbation-aware filtering
│   ├── ingest_soma.py           # SOMA ingestion: H5AD → TileDB-SOMA on S3
│   └── soma_dataloader.py       # PyTorch DataLoader streaming from SOMA atlas
├── docker/
│   ├── scanpy/Dockerfile        # murtiabhishek/scanpy-qc:1.0.1
│   ├── soma/Dockerfile          # murtiabhishek/soma-ingest:1.0.1
│   └── perturb/Dockerfile       # murtiabhishek/perturb-qc:1.0.0
└── data/
    ├── samplesheet.csv          # scRNA-seq H5AD inputs
    ├── nfcore_samplesheet.csv   # FASTQ inputs for nf-core/scrnaseq
    ├── perturb_samplesheet.csv  # Perturb-seq H5AD inputs
    └── vocab_v05.json           # 36K gene vocabulary for tokenization
```

---

## Quickstart

### Prerequisites

- Nextflow 25.x (`curl -s https://get.nextflow.io | bash`)
- Docker Desktop
- AWS CLI configured (`aws configure`)
- AWS Batch compute environment and job queue (see Setup)
- Seqera Platform account (free at seqera.io)

### Run on AWS Batch

```bash
export TOWER_ACCESS_TOKEN=your_seqera_token

# scRNA-seq QC + atlas ingestion
nextflow run main.nf -profile aws --input data/samplesheet.csv

# Perturb-seq QC + atlas ingestion
nextflow run main.nf -profile aws --input data/perturb_samplesheet.csv

# Full alignment from FASTQs (nf-core/scrnaseq)
nextflow run nf-core/scrnaseq -profile aws -r 2.7.1 \
  --input data/nfcore_samplesheet.csv \
  --star_index s3://nextflow-scrna-abhishek/indexes/star_GRCh38/ \
  --aligner star --protocol 10XV3 \
  --outdir s3://nextflow-scrna-abhishek/results/nfcore_scrnaseq
```

Monitor runs at: https://cloud.seqera.io

### Query the atlas

```python
import tiledbsoma as soma

# Query PBMC atlas by leiden cluster
with soma.Experiment.open("s3://nextflow-scrna-abhishek/atlas/pbmc_1k") as exp:
    obs = exp.obs.read(value_filter="leiden == '0'").concat().to_pandas()

# Query Perturb-seq atlas by perturbation target
with soma.Experiment.open("s3://nextflow-scrna-abhishek/atlas/K562_essential") as exp:
    rpl3_cells = exp.obs.read(value_filter="target_gene == 'RPL3'").concat().to_pandas()
```

### Stream to PyTorch

```python
from bin.soma_dataloader import make_soma_dataloader, make_perturb_soma_dataloader

# Standard DataLoader (PBMC)
loader = make_soma_dataloader(
    soma_uri   = "s3://nextflow-scrna-abhishek/atlas/pbmc_1k",
    vocab_path = "data/vocab_v05.json",
    batch_size = 256,
)

# Perturbation-aware DataLoader (K562)
loader = make_perturb_soma_dataloader(
    soma_uri    = "s3://nextflow-scrna-abhishek/atlas/K562_essential",
    vocab_path  = "data/vocab_v05.json",
    target_gene = "RPL3",   # or controls_only=True
    batch_size  = 256,
)

for input_ids, attention_mask in loader:
    # input_ids: (batch, seq_len=1024) — rank-tokenized gene IDs
    # attention_mask: (batch, seq_len) — 1 for real tokens, 0 for padding
    logits = model(input_ids, attention_mask)
```

---

## Setup

### AWS IAM

Create IAM user `nextflow-dev` with AdministratorAccess for development.

Create instance role `nextflow-batch-instance-role` with:
- `AmazonEC2ContainerServiceforEC2Role`
- `AmazonS3FullAccess`
- `AmazonEC2ContainerRegistryReadOnly`

Trust policy must include `ec2.amazonaws.com` and `ecs-tasks.amazonaws.com`.

### AWS Batch

- Compute environment: Managed, Spot, `default_x86_64`, max 256 vCPU, `SPOT_CAPACITY_OPTIMIZED`
- Job queue: `nextflow-batch-queue`, priority 1

### S3

```bash
aws s3 mb s3://nextflow-scrna-abhishek --region us-east-1
```

---

## Samplesheet formats

**scRNA-seq (H5AD input):**
```csv
sample_id,h5ad
pbmc_1k,s3://nextflow-scrna-abhishek/results/scanpy_qc/pbmc_1k/pbmc_1k_filtered_matrix.h5ad
```

**Perturb-seq (H5AD input):**
```csv
sample_id,perturb_h5ad
K562_essential,s3://nextflow-scrna-abhishek/raw/replogle/K562_essential_subset_50k.h5ad
```

**FASTQ input (nf-core/scrnaseq):**
```csv
sample,fastq_1,fastq_2,expected_cells
pbmc_1k,s3://bucket/pbmc_1k_L001_R1.fastq.gz,s3://bucket/pbmc_1k_L001_R2.fastq.gz,1000
```

---

## HPC vs Cloud: lessons learned

| | UCSF HPC (SLURM) | AWS Batch |
|---|---|---|
| Job submission | `sbatch script.sh` | Nextflow submits automatically |
| Storage | Scratch NFS filesystem | S3 object storage |
| Environments | Conda / environment modules | Docker containers |
| Scaling | Fixed cluster size | Elastic, scales to zero |
| Cost visibility | Invisible (institutional) | Per-job, per-second billing |
| Reproducibility | Environment drift risk | Containers = exact reproducibility |
| Data locality | Data must be on scratch | S3 mounted via Fusion as local filesystem |
| Spot interruptions | N/A | Automatic retry with `maxRetries = 2` |

The biggest mental shift: on HPC you go to the data; on cloud the data comes to the compute. Fusion makes this transparent — containers see S3 paths as if they were local files, with no AWS CLI required inside the container.

---

## Datasets

**PBMC 1k v3** — 10x Genomics, publicly available. 2 lanes, paired R1/R2, ~5GB FASTQs. Used for scRNA-seq pipeline validation.

**Replogle et al. 2022** — K562 essential gene Perturb-seq screen. 39K cells (subset), 51 perturbation targets + non-targeting controls. Source: [figshare](https://plus.figshare.com/articles/dataset/20029387).

---

## Benchmarks

### nf-core/scrnaseq — PBMC 1k v3 (GRCh38, STARsolo)

| Metric | Value |
|---|---|
| Total reads | 66.6M |
| Uniquely mapped | 87.74% |
| Alignment throughput | 322M reads/hour |
| Alignment time | 12m 24s |
| Total pipeline duration | 58m 52s |
| Total CPU hours | 10.7 |
| Instance type | r6i.4xlarge (spot) |
| Estimated cost | ~$4–6 |

### Perturb-seq pipeline — K562 essential subset

| Metric | Value |
|---|---|
| Input cells | 39,002 |
| Perturbation targets | 51 |
| Control cells | 10,691 |
| Pipeline duration | 5m 2s |
| Total CPU hours | 0.2 |
| Atlas size on S3 | ~500MB |

STAR genome index cached at `s3://nextflow-scrna-abhishek/indexes/star_GRCh38/` — eliminates 20-minute index generation on subsequent runs.

---

## Milestones

- [x] Phase 0 — Nextflow DSL2 + Docker + AWS Batch + Seqera monitoring
- [x] Phase 0 — nf-core/scrnaseq on real PBMC data, 87.74% mapping rate
- [x] Phase 1 — Scanpy QC module on AWS Batch
- [x] Phase 2 — TileDB-SOMA atlas ingestion on AWS Batch
- [x] Phase 3 — PyTorch SOMAIterableDataset, perturbation-aware DataLoader
- [x] Phase 4 — Perturb-seq QC + SOMA atlas queryable by perturbation target

## Planned

- [ ] Perturb-seq processing stack as full cloud-native Nextflow pipeline
- [ ] Multiome (scRNA-seq + scATAC-seq) processing module
- [ ] SOMA DataLoader integration with LiverTransformer training loop