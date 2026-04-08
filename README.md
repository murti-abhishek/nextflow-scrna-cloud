# nextflow-scrna-cloud

Cloud-scale scRNA-seq pipeline with Nextflow DSL2, Docker, and AWS Batch.

Built as a portfolio project to bridge HPC-based bioinformatics (UCSF SLURM) to cloud-native pipeline orchestration — a core requirement for data engineering roles in virtual cell and computational biology.

---

## Architecture

```
Your Mac (control plane)
    │
    ├── nextflow run main.nf -profile aws
    │
    ▼
AWS Batch (job scheduler)
    │   nextflow-batch-queue
    │   nextflow-compute-env (spot, x86, max 256 vCPU)
    │
    ├──▶ EC2 Spot Instance — sample 1
    │        └── Docker container (FastQC)
    │             ├── reads from  s3://bucket/raw/
    │             └── writes to   s3://bucket/work/
    │
    └──▶ EC2 Spot Instance — sample 2  (parallel)
             └── Docker container (FastQC)
                  ├── reads from  s3://bucket/raw/
                  └── writes to   s3://bucket/work/
                  
S3 (nextflow-scrna-abhishek)
    ├── raw/        input FASTQs
    ├── work/       Nextflow intermediate files
    └── results/    final outputs (HTML/ZIP reports)
```

Key components:
- **Nextflow DSL2** — workflow orchestrator, runs as a lightweight coordinator on your Mac
- **AWS Batch** — managed job scheduler, provisions EC2 instances per task
- **EC2 Spot** — actual compute, ~70% cheaper than on-demand, terminates after each job
- **Docker** — each process runs in an isolated container pulled from a registry
- **Wave (Seqera)** — augments containers with Fusion filesystem support
- **Fusion** — mounts S3 directly into containers, eliminating the need for AWS CLI in images
- **Seqera Platform** — run monitoring, task logs, resource usage dashboard

---

## Stack

| Layer | Tool |
|---|---|
| Workflow orchestration | Nextflow DSL2 |
| Containerization | Docker |
| Cloud compute | AWS Batch + EC2 Spot |
| Storage | Amazon S3 |
| Container registry | AWS ECR / quay.io biocontainers |
| S3 filesystem | Fusion (Wave/Seqera) |
| Monitoring | Seqera Platform (Tower) |

---

## Repository structure

```
.
├── main.nf                  # main workflow
├── nextflow.config          # base config with local and aws profiles
├── conf/
│   ├── local.config         # local execution profile
│   └── aws_batch.config     # AWS Batch + Wave + Fusion config
├── modules/
│   ├── cellranger/          # (upcoming) Cell Ranger alignment module
│   └── scanpy_qc/           # (upcoming) Scanpy QC module
├── docker/
│   ├── cellranger/          # (upcoming) Cell Ranger Dockerfile
│   └── scanpy/              # (upcoming) Scanpy Dockerfile
└── data/
    └── samplesheet.csv      # sample manifest (points to S3 paths)
```

---

## Quickstart

### Prerequisites

- Nextflow 25.x (`curl -s https://get.nextflow.io | bash`)
- Docker Desktop
- AWS CLI configured (`aws configure`)
- AWS Batch compute environment and job queue (see Setup below)

### Run locally

```bash
nextflow run main.nf -profile local
```

### Run on AWS Batch

```bash
export TOWER_ACCESS_TOKEN=your_seqera_token
nextflow run main.nf -profile aws
```

Monitor at: https://cloud.seqera.io

---

## Setup

### 1. AWS IAM

Create an IAM user (`nextflow-dev`) with AdministratorAccess for development. In production, scope down to S3, Batch, ECR, and CloudWatch only.

Create an instance role (`nextflow-batch-instance-role`) with:
- `AmazonEC2ContainerServiceforEC2Role`
- `AmazonS3FullAccess`
- `AmazonEC2ContainerRegistryReadOnly`

Trust policy must include both `ec2.amazonaws.com` and `ecs-tasks.amazonaws.com`.

### 2. S3 bucket

```bash
aws s3 mb s3://nextflow-scrna-abhishek --region us-east-1
```

### 3. AWS Batch

- Compute environment: Managed, Spot, `default_x86_64`, max 256 vCPU
- Job queue: connected to the compute environment, priority 1

### 4. Samplesheet

`data/samplesheet.csv` format:

```csv
sample_id,fastq_r1,fastq_r2
pbmc_1k_L001,s3://bucket/raw/pbmc_1k_v3_S1_L001_R1_001.fastq.gz,s3://bucket/raw/pbmc_1k_v3_S1_L001_R2_001.fastq.gz
pbmc_1k_L002,s3://bucket/raw/pbmc_1k_v3_S1_L002_R1_001.fastq.gz,s3://bucket/raw/pbmc_1k_v3_S1_L002_R2_001.fastq.gz
```

---

## HPC vs Cloud: lessons learned

Coming from SLURM-based pipelines on UCSF's HPC, the key differences:

| | UCSF HPC (SLURM) | AWS Batch |
|---|---|---|
| Job submission | `sbatch script.sh` | Nextflow submits automatically |
| Storage | Scratch filesystem (shared NFS) | S3 object storage |
| Environments | Conda/modules | Docker containers |
| Scaling | Fixed cluster size | Elastic, scales to zero |
| Cost visibility | Invisible (institutional) | Per-job, per-second billing |
| Reproducibility | Environment drift risk | Containers = exact reproducibility |
| Data locality | Data must be on scratch | Data lives in S3, pulled on demand |

The biggest mental shift: **on HPC you go to the data; on cloud the data comes to the compute.** Fusion makes this seamless by mounting S3 directly into the container as if it were a local filesystem.

---

## Milestones completed

- [x] Milestone 1 — Containerize a pipeline step (FastQC in Docker)
- [x] Milestone 2 — Nextflow DSL2 workflow running locally
- [x] Milestone 3 — Deployed to AWS Batch with Fusion filesystem
- [x] Milestone 4 — Seqera Platform monitoring and run dashboard
- [x] Milestone 5 — Multi-sample parallel batch run on real PBMC data

## Upcoming

- [ ] Phase 1 — Cell Ranger alignment + Scanpy QC pipeline
- [ ] Phase 2 — TileDB-SOMA data lake and atlas building
- [ ] Phase 3 — PyTorch DataLoader integration for model training
- [ ] Phase 4 — Perturb-seq processing pipeline

---

## Dataset

10x Genomics PBMC 1k v3 — publicly available at [10xgenomics.com](https://www.10xgenomics.com/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-1-standard-3-0-0)

2 lanes × paired R1/R2 reads, ~5GB total. Used for benchmarking pipeline parallelism across samples.

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

Pipeline ran on AWS Batch with Fusion filesystem and spot instances.
STAR genome index cached to S3 for subsequent runs.