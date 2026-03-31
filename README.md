# CLAUDE Pipeline

**Cancer-Linked Analysis of Underlying DNA Elements**

A 13-phase bioinformatics pipeline for extracting and characterising bacteriophages and mobile genetic elements from tumour whole-genome sequencing (WGS) data. The pipeline was developed to address the growing evidence that the tumour microbiome, particularly phages and mobilome elements, plays a role in colorectal cancer (CRC) and other malignancies.

---

## Background

Microbial reads in tumour WGS data typically make up around 0.35% of total sequencing output. Extracting meaningful biological signal from this fraction requires careful host depletion, stringent quality control, and a multi-tool approach at every phase to minimise false positives. This pipeline was built with that problem in mind, drawing on lessons from the TCGA microbiome controversy (Gihawi et al. 2023) and following the recommendations of Ge et al. 2025 for low-biomass metagenomics.

The full pipeline covers everything from raw read acquisition through to an integrated master table linking viral contigs, host predictions, prophage activity, antibiotic resistance genes, and mobile elements.

---

## Design Principles

Three principles underpin every decision in this pipeline, derived directly from the TCGA microbiome controversy:

**Principle 1 — Decontaminate first, ask questions later.**
Every signal must survive multi-layer decontamination before it is considered real. False positives in low-biomass settings (tumour WGS: ~0.35% microbial reads) are far more damaging than false negatives.

**Principle 2 — Consensus over any single tool.**
No single phage detection tool is sufficient. We require agreement from at least two independent methods with orthogonal algorithms before calling a contig as viral.

**Principle 3 — Every step must be auditable.**
Each phase produces a QC checkpoint file. If a checkpoint fails, the pipeline halts with a diagnostic message. No silent failures.

---

## Architecture

```
┌─────────────────────────────────────────────────────────────────────┐
│                    CLAUDE Pipeline Architecture                      │
├─────────────────────────────────────────────────────────────────────┤
│                                                                      │
│  RAW READS (TCGA / ENA / Hartwig unmapped BAM/CRAM)                 │
│       │                                                              │
│       ▼                                                              │
│  ┌─────────────────────────────────────────┐                        │
│  │ PHASE 1: Data Acquisition               │                        │
│  │   fasterq-dump / BAM → FASTQ            │  ── QC-1: Read counts  │
│  └──────────────┬──────────────────────────┘                        │
│                 ▼                                                    │
│  ┌─────────────────────────────────────────┐                        │
│  │ PHASE 2: Host Depletion (3-pass)        │                        │
│  │   Pass 1: BWA-MEM2 (hg38)              │                        │
│  │   Pass 2: Bowtie2 (T2T-CHM13)          │  ── QC-2: Depletion %  │
│  │   Pass 3: Hostile (phage-masked index)  │                        │
│  └──────────────┬──────────────────────────┘                        │
│                 ▼                                                    │
│  ┌─────────────────────────────────────────┐                        │
│  │ PHASE 3: Read QC                        │                        │
│  │   fastp adapter trim + quality filter   │  ── QC-3: Quality dist │
│  └──────────────┬──────────────────────────┘                        │
│                 ▼                                                    │
│  ┌─────────────────────────────────────────┐                        │
│  │ PHASE 4: Taxonomic Profiling            │                        │
│  │   Kraken2 + Bracken (abundance)         │  ── QC-4: Known taxa   │
│  │   KrakenUniq (false-positive control)   │                        │
│  └──────────────┬──────────────────────────┘                        │
│                 ▼                                                    │
│  ┌─────────────────────────────────────────┐                        │
│  │ PHASE 5: Assembly                       │                        │
│  │   metaSPAdes + MEGAHIT + metaviralSPAdes│  ── QC-5: Assembly     │
│  │   → merge & deduplicate at 99% ANI      │       stats            │
│  └──────────────┬──────────────────────────┘                        │
│                 ▼                                                    │
│  ┌─────────────────────────────────────────┐                        │
│  │ PHASE 6: Decontam + MAG Binning         │                        │
│  │   BLASTn vs T2T + UniVec screen         │  ── QC-6: Contam rate  │
│  │   MetaBAT2 → DAS Tool → CheckM2         │  ── QC-7: MAG quality  │
│  └──────────┬──────────────┬───────────────┘                        │
│             ▼              ▼                                         │
│  ┌──────────────────────────────────────────┐                       │
│  │ PHASE 6B: MAG Annotation                 │                       │
│  │   Bakta → mobileOG-db                    │  ── QC-6B: CDS annot  │
│  └──────────────────┬───────────────────────┘                       │
│                     ▼                                                │
│  ┌──────────────────────────────────────────┐                       │
│  │ PHASE 7: Viral Identification            │                       │
│  │   geNomad + DeepVirFinder + VirSorter2   │  ── QC-8: Agreement   │
│  │   ≥2/3 consensus → CheckV QC            │  ── QC-9: Completeness│
│  └──────────────────┬───────────────────────┘                       │
│                     ▼                                                │
│  ┌──────────────────────────────────────────┐                       │
│  │ PHASE 8: Phage Annotation                │                       │
│  │   Pharokka (PHROG/CARD/VFDB/INPHARED)    │  ── QC-10: Prophage   │
│  │   PropagAtE (activity) + iPHoP (hosts)   │  ── QC-11: Host pred  │
│  └──────────────────┬───────────────────────┘                       │
│                     ▼                                                │
│  ┌──────────────────────────────────────────┐                       │
│  │ PHASE 9: Mobilome                        │                       │
│  │   ISEScan + IntegronFinder + mobileOG-db │  ── QC-12: MGE count  │
│  └──────────────────┬───────────────────────┘                       │
│                     ▼                                                │
│  ┌──────────────────────────────────────────┐                       │
│  │ PHASE 10: Integration                    │                       │
│  │   Master table, QC dashboard,            │  ── QC-13: Validation │
│  │   biological summary                     │                       │
│  └──────────────────┬───────────────────────┘                       │
│                     ▼                                                │
│  ┌──────────────────────────────────────────┐                       │
│  │ PHASE 11: Defence Systems                │                       │
│  │   DefenseFinder + PADLOC                 │  ── QC-14: Defence    │
│  │   + CRISPRCasFinder                      │       systems count   │
│  └──────────────────────────────────────────┘                       │
└─────────────────────────────────────────────────────────────────────┘
```

---

## Pipeline Overview

The pipeline runs in 13 phases, grouped into four broad stages:

**Stage 1: Data Preparation (Phases 1 to 3)**

Phase 1 handles data acquisition from public repositories (ENA/SRA) or local sources. Phase 2 runs a three-pass host depletion strategy: BWA-MEM2 against hg38, then Bowtie2 against the T2T-CHM13 genome, and finally Hostile using a phage-masked T2T+HLA index. Using a phage-masked index in the final pass is important because it prevents genuine phage reads from being discarded when they share sequence similarity with human genomic regions. Phase 3 applies quality filtering and adapter trimming with fastp, followed by FastQC and MultiQC reporting.

**Stage 2: Community Profiling and Assembly (Phases 4 to 6)**

Phase 4 runs taxonomic profiling using two complementary approaches: Kraken2 with Bracken for broad bacterial/viral classification and abundance estimation, and KrakenUniq with a complete-genomes-only viral database for false-positive control via unique k-mer counting. Phase 5 runs three assemblers in parallel: metaSPAdes, MEGAHIT, and metaviralSPAdes. Each recovers different parts of the metagenome, and the outputs are merged and deduplicated at 99% ANI using CD-HIT-EST. Phase 6 removes residual human contigs using BLASTn against T2T-CHM13 and UniVec, then bins contigs into metagenome-assembled genomes (MAGs) using MetaBAT2, with quality assessment by CheckM2 and DAS Tool. MAGs are annotated with Bakta and mobile elements identified via mobileOG-db.

**Stage 3: Viral Identification and Annotation (Phases 7 to 9)**

Phase 7 identifies viral contigs using a consensus approach across geNomad, DeepVirFinder, and VirSorter2, requiring agreement from at least two tools. CheckV assesses completeness and quality. Phase 8 annotates phage genomes using Pharokka (PHROG, CARD, VFDB, INPHARED databases), predicts prophage activity using PropagAtE, and predicts phage-host relationships using iPHoP against the full Jun 2025 database. Phase 9 characterises the mobilome: plasmid detection via geNomad, insertion sequences via ISEScan, integrons via IntegronFinder, and mobile element annotation via mobileOG-db.

**Stage 4: Defence Systems and Integration (Phases 10 to 11)**

Phase 10 consolidates all results into a master integration table and QC dashboard linking viral contigs, host predictions, prophage activity, ARGs, and mobile elements. Phase 11 profiles bacterial defence systems across all MAGs using DefenseFinder, PADLOC, and CRISPRCasFinder.

---

## Requirements

### Compute

| Resource   | Minimum   | Recommended  | Notes                         |
|------------|-----------|--------------|-------------------------------|
| CPU cores  | 16        | 64           | Assembly is CPU-bound         |
| RAM        | 64 GB     | 256 GB       | metaSPAdes and iPHoP peak RAM |
| Storage    | 2 TB      | 10 TB        | Raw WGS and databases         |
| GPU        | Optional  | NVIDIA T4    | VAMB binning                  |
| OS         | Ubuntu 22.04+ | Ubuntu 24.04 |                           |

### Conda Environments

The pipeline uses 9 separate conda environments to manage conflicting dependencies. Environment YAML files are in the `envs/` directory. To recreate any environment:

```bash
mamba env create -f envs/claude_pipeline.yml
```

| Environment | Key Tools |
|-------------|-----------|
| claude_pipeline | samtools, bowtie2, fastp, kraken2, megahit, spades, metabat2, blast, seqkit, ISEScan, IntegronFinder, PropagAtE |
| dvf_env | DeepVirFinder (Python 3.6, keras 2.2.4, theano) |
| genomad_env | geNomad 1.9.7 |
| checkv_env | CheckV 1.0.3 |
| virsorter2_env | VirSorter2 2.2.4 |
| checkm2_env | CheckM2 1.0.2 |
| das_tool_env | DAS Tool 1.1.7 |
| pharokka_env | Pharokka 1.9.1 |
| iphop_env | iPHoP 1.4.2 |

For fresh installation commands, see [`docs/01_SETUP.md`](docs/01_SETUP.md).

### Databases

The following databases need to be downloaded separately before running the pipeline.

| Database | Tool | Approx. Size |
|----------|------|--------------|
| Hostile phage-masked T2T+HLA index | Hostile 2.x | ~4 GB |
| k2_standard_08gb | Kraken2 + Bracken | 8 GB |
| RefSeq complete viral genomes | KrakenUniq | ~21 GB |
| T2T-CHM13v2.0 BLAST db | BLASTn | ~3 GB |
| UniVec | BLASTn | minimal |
| genomad_db | geNomad | ~3 GB |
| checkv-db-v1.5 | CheckV | ~2 GB |
| combined.hmm | VirSorter2 | ~8 GB |
| CheckM2_database | CheckM2 | ~3 GB |
| PHROG v4 + CARD + VFDB + INPHARED | Pharokka | ~5 GB |
| iPHoP Jun 2025 full database | iPHoP | ~280 GB |

Download commands for all databases are in [`docs/01_SETUP.md`](docs/01_SETUP.md).

### Directory Structure

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
    07b_bakta_annotation,
    08_viral_contigs,
    09_viral_qc,
    10_prophages,
    11_phage_annotation,
    12_mobilome,
    13_integration,
    14_defence_systems,
    qc_checkpoints,
    logs,
    scripts
}
```

---

## Usage

Each phase has a dedicated SOP file in `docs/sops/` with step-by-step commands, expected outputs, QC thresholds, and troubleshooting notes. Read the relevant SOP before executing each phase.

```bash
# Example: run Phase 4 taxonomic profiling
# Read docs/sops/phase_4_sop.md first, then:

conda run -n claude_pipeline bash -c "
    kraken2 \
        --db databases/kraken2 \
        --paired \
        --gzip-compressed \
        --output pipeline_run/04_read_profiles/sample_kraken2.out \
        --report pipeline_run/04_read_profiles/sample_kraken2.report \
        --confidence 0.1 \
        --threads 8 \
        pipeline_run/03_qc_filtered/sample_R1.final.fastq.gz \
        pipeline_run/03_qc_filtered/sample_R2.final.fastq.gz
"
```

---

## Validation

The pipeline was validated end-to-end on a VLP-enriched gut virome dataset and on a pair of stool metagenomes from the Wirbel et al. 2019 Nature Medicine CRC cohort: one stage III CRC patient and one healthy control. QC checkpoints for both validation runs are in `pipeline_run/qc_checkpoints_public/`. Sample accession numbers have been replaced with aliases (CRC_sample and control_sample) in all public-facing files.

Key validation results:

* CRC sample: 335 viral contigs (1 complete 108 kb phage genome), 242 host predictions (8/8 iPHoP tools including RaFAH), 1 active prophage, 1 ARG (tet(O))
* Control sample: 542 viral contigs (17 complete), 352 host predictions, 5 active prophages, 0 ARGs
* Phase 11 defence profiling: 208 defence systems across 30 CRC MAGs vs 371 across 92 control MAGs
* Primary biological finding: 6 Fusobacterium-associated phages exclusive to the CRC sample

Full results are in `pipeline_run/13_integration/Wirbel_CRC_public/`.

---

## Data Availability

Raw sequencing data is not included in this repository. The pipeline scripts, SOPs, environment YAMLs, and anonymised QC outputs are provided for reproducibility. Researchers wishing to reproduce the full analysis should download the relevant samples from ENA/SRA and follow the phase SOPs.

---

## References

* Wirbel et al. 2019. Meta-analysis of fecal metagenomes reveals global microbial signatures that are specific for colorectal cancer. Nature Medicine 25:679-689.
* Gihawi et al. 2023. Major data analysis errors invalidate cancer microbiome findings. mBio.
* Ge et al. 2025. Recommendations for low-biomass metagenomics.
* Dohlman et al. 2025. Identifying and correcting false-positive cancer microbiome findings.
* Roux et al. 2019. Minimum information about an uncultivated virus genome (MIUViG). Nature Biotechnology.

---

## Notes

This pipeline was developed and tested on Ubuntu 24.04 with an AMD Ryzen 5 3550H CPU and 30 GB RAM. Some steps, particularly iPHoP with the full Jun 2025 database and the RaFAH random forest model, push close to the 30 GB limit. Running on a machine with 64 GB RAM or more is strongly recommended for production use. The `TROUBLESHOOTING_public.md` file documents issues encountered during development and the fixes applied.
