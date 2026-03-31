# CLAUDE Pipeline

**Cancer-Linked Analysis of Underlying DNA Elements**

A 13-phase bioinformatics pipeline for extracting and characterising bacteriophages and mobile genetic elements from tumour whole-genome sequencing (WGS) data. The pipeline was developed to address the growing evidence that the tumour microbiome, particularly phages and mobilome elements, plays a role in colorectal cancer (CRC) and other malignancies.

---

## Background

Microbial reads in tumour WGS data typically make up around 0.35% of total sequencing output. Extracting meaningful biological signal from this fraction requires careful host depletion, stringent quality control, and a multi-tool approach at every phase to minimise false positives. This pipeline was built with that problem in mind, drawing on lessons from the TCGA microbiome controversy (Gihawi et al. 2023) and following the recommendations of Ge et al. 2025 for low-biomass metagenomics.

The full pipeline covers everything from raw read acquisition through to an integrated master table linking viral contigs, host predictions, prophage activity, antibiotic resistance genes, and mobile elements.

---

## Pipeline Overview

The pipeline runs in 13 phases, grouped into four broad stages:

**Stage 1: Data Preparation (Phases 1 to 3)**

Phase 1 handles data acquisition from public repositories (ENA/SRA) or local sources. Phase 2 runs a three-pass host depletion strategy: BWA-MEM2 against hg38, then Bowtie2 against the T2T-CHM13 genome, and finally Hostile using a phage-masked T2T+HLA index. Using a phage-masked index in the final pass is important because it prevents genuine phage reads from being discarded when they share sequence similarity with human genomic regions. Phase 3 applies quality filtering and adapter trimming with fastp, followed by FastQC and MultiQC reporting.

**Stage 2: Community Profiling and Assembly (Phases 4 to 6)**

Phase 4 runs taxonomic profiling using two complementary approaches: Kraken2 with Bracken for broad bacterial/viral classification and abundance estimation, and KrakenUniq with a complete-genomes-only viral database for false-positive control via unique k-mer counting. Phase 5 runs three assemblers in parallel: metaSPAdes, MEGAHIT, and metaviralSPAdes. Each recovers different parts of the metagenome, and the outputs are merged and deduplicated at 99% ANI using CD-HIT-EST. Phase 6 removes residual human contigs using BLASTn against T2T-CHM13 and UniVec, then bins contigs into metagenome-assembled genomes (MAGs) using MetaBAT2, with quality assessment by CheckM2 and DAS Tool.

**Stage 3: Viral Identification and Annotation (Phases 7 to 9)**

Phase 7 identifies viral contigs using a consensus approach across geNomad, DeepVirFinder, and VirSorter2, requiring agreement from at least two tools. CheckV assesses completeness and quality. Phase 8 annotates phage genomes using Pharokka (PHROG, CARD, VFDB, INPHARED databases), predicts prophage activity using PropagAtE, and predicts phage-host relationships using iPHoP against the full Jun 2025 database. Phase 9 characterises the mobilome: plasmid detection via geNomad, insertion sequences via ISEScan, integrons via IntegronFinder, and mobile element annotation via mobileOG-db.

**Stage 4: Defence Systems and Integration (Phases 10 to 13)**

Phase 10 focuses on bacterial MAG annotation using Bakta for comprehensive gene prediction and functional annotation. Phase 11 profiles bacterial defence systems using DefenseFinder, PADLOC, and CRISPRCasFinder across all MAGs. Phase 12 consolidates all results into a master integration table and QC dashboard. Phase 13 provides biological interpretation of the full dataset.

---

## Requirements

### Conda Environments

The pipeline uses 9 separate conda environments to manage conflicting dependencies. All environments should be created before running any phase.

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

Environment YAML files for each are provided in the `envs/` directory.

### Databases

The following databases need to be downloaded separately. Approximate sizes are listed to help with storage planning.

| Database | Tool | Size |
|----------|------|------|
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

Database setup instructions are in `01_SETUP.md`.

---

## Usage

Each phase has a dedicated SOP file (phase_1_sop.md through phase_11_sop.md) with step-by-step commands, expected outputs, QC thresholds, and troubleshooting notes. The recommended workflow is to read the relevant SOP before executing each phase.

The current pipeline state (which phases are complete, input/output locations, key metrics) is tracked in `PIPELINE_STATE.md`.

```
# Example: run Phase 4 taxonomic profiling
# Read phase_4_sop.md first, then:

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

The pipeline was validated end-to-end on a VLP-enriched gut virome dataset (Wahida et al.) and on a pair of stool metagenomes from the Wirbel et al. 2019 Nature Medicine CRC cohort: one stage III CRC patient and one healthy control. QC checkpoints for both validation runs are in `pipeline_run/qc_checkpoints_public/`. Sample accession numbers have been replaced with aliases (CRC_sample and control_sample) in all public-facing files.

Key validation results:

* CRC sample: 335 viral contigs (1 complete 108 kb phage genome), 242 host predictions (8/8 iPHoP tools including RaFAH), 1 active prophage, 1 ARG (tet(O))
* Control sample: 542 viral contigs (17 complete), 352 host predictions, 5 active prophages, 0 ARGs
* Phase 11 defence profiling: 208 defence systems in CRC vs 371 in control across 30 vs 92 MAGs

---

## Data Availability

Raw sequencing data is not included in this repository. The pipeline scripts, SOPs, and anonymised QC outputs are provided for reproducibility. Researchers wishing to reproduce the full analysis should download the relevant samples from ENA/SRA and follow the phase SOPs.

---

## References

* Wirbel et al. 2019. Meta-analysis of fecal metagenomes reveals global microbial signatures that are specific for colorectal cancer. Nature Medicine 25:679-689.
* Gihawi et al. 2023. Major data analysis errors invalidate cancer microbiome findings. mBio.
* Ge et al. 2025. Recommendations for low-biomass metagenomics.
* Roux et al. 2019. Minimum information about an uncultivated virus genome (MIUViG). Nature Biotechnology.

---

## Notes

This pipeline was developed and tested on Ubuntu 24.04 with an AMD Ryzen 5 3550H CPU and 30 GB RAM. Some steps (particularly iPHoP with the full database and RaFAH random forest model) require close to the full 30 GB. Running on a machine with 64 GB RAM or more is strongly recommended for production use. The `TROUBLESHOOTING.md` file documents every issue encountered during development along with the fixes applied.
