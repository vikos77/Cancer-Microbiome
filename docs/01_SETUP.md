# CLAUDE Pipeline — Part 01: Setup & Architecture

> **Navigation:** 01_SETUP → [02_DATA_ACQUISITION](02_DATA_ACQUISITION.md) → [Phase 1-2 SOPs](sops/phase_1_sop.md) → [Phase 3-4 SOPs](sops/phase_3_sop.md) → [Phase 5 SOP](sops/phase_5_sop.md) → [Phase 6 SOP](sops/phase_6_sop.md) → [Phase 7 SOP](sops/phase_7_sop.md) → [Phase 8 SOP](sops/phase_8_sop.md) → [Phase 9 SOP](sops/phase_9_sop.md) → [Phase 10-11 SOPs](sops/phase_10_sop.md)

---

## Pipeline Name

**CLAUDE: Cancer-Linked Analysis of Underlying DNA Elements**

**Version:** 0.1.0-dev
**Date:** 2026-02-24
**License:** MIT

---

## 1. Pipeline Philosophy

Three non-negotiable principles derived from the TCGA cancer microbiome controversy
(Gihawi et al. 2023; Ge et al. 2025; Dohlman et al. 2025):

**Principle 1 — Decontaminate first, ask questions later.**
Every signal must survive multi-layer decontamination before it is considered real.
False positives in low-biomass settings (tumor WGS: ~0.35% microbial reads) are
far more damaging than false negatives.

**Principle 2 — Consensus over any single tool.**
No single phage detection tool is sufficient. We require agreement from ≥2 independent
methods with orthogonal algorithms before calling a contig as viral.

**Principle 3 — Every step must be auditable.**
Each phase produces a QC checkpoint file. If a checkpoint fails, the pipeline halts
with a diagnostic message. No silent failures.

---

## 2. Architecture Overview

```
┌─────────────────────────────────────────────────────────────────────┐
│                    CLAUDE Pipeline Architecture                      │
├─────────────────────────────────────────────────────────────────────┤
│                                                                      │
│  RAW READS (TCGA/Hartwig unmapped BAM/CRAM)                         │
│       │                                                              │
│       ▼                                                              │
│  ┌─────────────────────────────────────────┐                         │
│  │ PHASE 1: Data Acquisition               │  02_DATA_ACQUISITION    │
│  │   BAM → FASTQ conversion                │  ── QC-1: Read counts   │
│  └──────────────┬──────────────────────────┘                         │
│                 ▼                                                     │
│  ┌─────────────────────────────────────────┐                         │
│  │ PHASE 2: Host Depletion (3-pass)        │  03_HOST_DEPLETION      │
│  │   Pass 1: Hostile + T2T-CHM13           │                         │
│  │   Pass 2: HPRC pangenome sweep          │  ── QC-2: Depletion %   │
│  │   Pass 3: Kraken2 human re-check        │                         │
│  └──────────────┬──────────────────────────┘                         │
│                 ▼                                                     │
│  ┌─────────────────────────────────────────┐                         │
│  │ PHASE 3-4: Read QC + Profiling          │  04_READ_QC_PROFILING   │
│  │   Adapter trim, quality filter,         │  ── QC-3: Quality dist  │
│  │   KrakenUniq + Kraken2/Bracken          │  ── QC-4: Known taxa    │
│  └──────────────┬──────────────────────────┘                         │
│                 ▼                                                     │
│  ┌─────────────────────────────────────────┐                         │
│  │ PHASE 5: Assembly                       │  05_ASSEMBLY             │
│  │   metaSPAdes + MEGAHIT + metaviralSPAdes│  ── QC-5: Assembly stats│
│  │   → merge & deduplicate contigs         │                         │
│  └──────────────┬──────────────────────────┘                         │
│                 ▼                                                     │
│  ┌─────────────────────────────────────────┐                         │
│  │ PHASE 6-7: Decontam + MAG Binning       │  06_DECONTAM_MAGS       │
│  │   BLASTn vs T2T, UniVec screen          │  ── QC-6: Contam. rate  │
│  │   MetaBAT2/DAS_Tool → CheckM2 → GTDB   │  ── QC-7: MAG quality   │
│  └──────────┬──────────────┬───────────────┘                         │
│             ▼              ▼                                          │
│  ┌──────────────────────────────────────────┐                        │
│  │ PHASE 8-9: Viral ID + QC                 │  07_VIRAL_ID            │
│  │  geNomad + DVF + VirSorter2 (≥2/3)      │  ── QC-8: Agreement     │
│  │  CheckV → vRhyme binning                 │  ── QC-9: Completeness  │
│  └──────────────────┬───────────────────────┘                        │
│                     ▼                                                 │
│  ┌──────────────────────────────────────────┐                        │
│  │ PHASE 10-11: Prophage + Annotation        │  08_PROPHAGE_ANNOTATION│
│  │  geNomad + PHASTEST + PropagAtE           │  ── QC-10: Concordance │
│  │  iPHoP host + PHROG function              │  ── QC-11: Annotation  │
│  └──────────────────┬───────────────────────┘                        │
│                     ▼                                                 │
│  ┌──────────────────────────────────────────┐                        │
│  │ PHASE 12: Mobilome                        │  09_MOBILOME            │
│  │  ISEScan + IntegronFinder + mobileOG-db   │  ── QC-12: MGE count   │
│  └──────────────────┬───────────────────────┘                        │
│                     ▼                                                 │
│  ┌──────────────────────────────────────────┐                        │
│  │ PHASE 13: Integration                     │  10_INTEGRATION         │
│  │  Epitope prediction, batch correction,    │  ── QC-13: Validation  │
│  │  clinical linkage, QC dashboard           │                        │
│  └──────────────────────────────────────────┘                        │
└─────────────────────────────────────────────────────────────────────┘
```

---

## 3. Compute Requirements

| Resource        | Minimum        | Recommended     | Notes                           |
|-----------------|----------------|-----------------|---------------------------------|
| CPU cores       | 16             | 64              | Assembly is CPU-bound           |
| RAM             | 64 GB          | 256 GB          | metaSPAdes peak usage           |
| Storage         | 2 TB           | 10 TB           | Raw TCGA WGS is large           |
| GPU             | Optional       | NVIDIA T4/A100  | VAMB binning, ViraLM            |
| OS              | Ubuntu 22.04+  | Ubuntu 24.04    |                                 |

---

## 4. Conda Environment Setup

```bash
# Create base environment
mamba create -n claude_pipeline python=3.10 -y
mamba activate claude_pipeline

# --- Core bioinformatics ---
mamba install -c bioconda -c conda-forge \
    samtools=1.20 \
    bedtools=2.31 \
    bowtie2=2.5.4 \
    minimap2=2.28 \
    fastp=0.23.4 \
    fastqc=0.12.1 \
    multiqc=1.22 \
    seqkit=2.8 \
    bbmap=39.08 \
    pigz=2.8 \
    -y

# --- Host depletion ---
pip install hostile

# --- Taxonomic profiling ---
mamba install -c bioconda \
    kraken2=2.1.3 \
    krakenuniq=1.0.4 \
    bracken=2.9 \
    -y

# --- Assembly ---
mamba install -c bioconda \
    megahit=1.2.9 \
    spades=4.0.0 \
    -y

# --- Binning ---
mamba install -c bioconda \
    metabat2=2.17 \
    maxbin2=2.2.7 \
    das_tool=1.1.7 \
    -y
pip install vamb

# --- Viral detection ---
pip install genomad
# NOTE: hostile is installed via conda, not pip (pip is system-locked on Ubuntu 24.04)
mamba install -c bioconda \
    virsorter=2.2.4 \
    checkv=1.0.3 \
    -y
git clone https://github.com/jessieren/DeepVirFinder.git

# --- Bacterial QC & taxonomy ---
mamba install -c bioconda \
    checkm2=1.0.2 \
    gtdbtk=2.4.0 \
    -y

# --- Phage annotation ---
pip install iphop
git clone https://github.com/KennthShang/PhaBOX.git

# --- Mobilome ---
mamba install -c bioconda \
    isescan=1.7.2.3 \
    integron_finder=2.0.5 \
    -y

# --- Prophage ---
git clone https://github.com/AnantharamanLab/PropagAtE.git

# --- Viral binning ---
pip install vrhyme

# --- Utilities ---
mamba install -c bioconda \
    prodigal=2.6.3 \
    diamond=2.1.9 \
    blast=2.16.0 \
    cd-hit=4.8.1 \
    coverm=0.7.0 \
    -y
```

---

## 5. Database Downloads

```bash
export CLAUDE_DB="/path/to/claude_databases"
mkdir -p ${CLAUDE_DB}/{hostile,kraken2,krakenuniq,genomad,checkv,gtdbtk,iphop,univec}

# --- Hostile: T2T-CHM13 + HLA masked index ---
# hostile 2.x caches indexes at ~/.local/share/hostile/ (not a custom path)
# Use phage-masked index: prevents phage reads from being removed as "human"
# This supersedes the original human-t2t-hla-argos985 recommendation
hostile index fetch --name human-t2t-hla.rs-viral-202401_ml-phage-202401 --bowtie2
# Symlink cache into project database dir for reference:
# ln -s ~/.local/share/hostile ${CLAUDE_DB}/hostile

# --- Kraken2: Standard + human + viral ---
kraken2-build --standard --db ${CLAUDE_DB}/kraken2/standard_plus_human
kraken2-build --download-library human --db ${CLAUDE_DB}/kraken2/standard_plus_human
kraken2-build --build --db ${CLAUDE_DB}/kraken2/standard_plus_human

# --- KrakenUniq: Curated complete genomes only ---
# Use Ge et al. 2025 approach: 50,651 genomes, 30,355 species
# Build from NCBI RefSeq complete genomes ONLY (no draft/contig-level)
krakenuniq-build --db ${CLAUDE_DB}/krakenuniq/complete_genomes \
    --threads 32

# --- geNomad database ---
genomad download-database ${CLAUDE_DB}/genomad

# --- CheckV database ---
checkv download_database ${CLAUDE_DB}/checkv

# --- GTDB-Tk (R220) ---
download-db.sh ${CLAUDE_DB}/gtdbtk

# --- iPHoP database ---
iphop download --db_dir ${CLAUDE_DB}/iphop

# --- UniVec for vector screening ---
wget -P ${CLAUDE_DB}/univec \
    https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec
makeblastdb -in ${CLAUDE_DB}/univec/UniVec -dbtype nucl \
    -out ${CLAUDE_DB}/univec/UniVec_db

# --- T2T-CHM13 v2.0 for contig decontamination ---
wget -P ${CLAUDE_DB}/hostile \
    https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
gunzip ${CLAUDE_DB}/hostile/chm13v2.0.fa.gz
makeblastdb -in ${CLAUDE_DB}/hostile/chm13v2.0.fa -dbtype nucl \
    -out ${CLAUDE_DB}/hostile/chm13v2_blastdb

# --- mobileOG-db ---
wget https://mobileogdb.flsi.cloud.vt.edu/entries/mobileOG-db_beatrix-1.6.All.faa \
    -P ${CLAUDE_DB}/
diamond makedb --in ${CLAUDE_DB}/mobileOG-db_beatrix-1.6.All.faa \
    --db ${CLAUDE_DB}/mobileOG-db
```

---

## 6. Directory Structure

```bash
PROJECT_DIR="/path/to/project"
mkdir -p ${PROJECT_DIR}/{
    00_raw,
    01_fastq,
    02_host_depleted,
    03_qc_filtered,
    04_read_profiles,
    05_assembly,
    06_decontam_contigs,
    07_mags,
    08_viral_contigs,
    09_viral_qc,
    10_prophages,
    11_phage_annotation,
    12_mobilome,
    13_integration,
    qc_checkpoints,
    logs,
    scripts
}
```

---

> **Next:** [02_DATA_ACQUISITION.md](02_DATA_ACQUISITION.md) — BAM/CRAM extraction and QC-1
