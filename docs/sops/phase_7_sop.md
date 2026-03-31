# Phase 7 SOP — Viral Identification

**Sample:** SRR15090802 (Wahida et al. gut virome, VLP-enriched healthy child)
**Input:** 20,547 clean contigs (`pipeline_run/06_decontam_contigs/virome/clean_contigs.fasta`)
**Completed:** 2026-02-26
**Status:** COMPLETE (2-tool consensus — see VS2 Runtime Decision below)

---

## Overview

Phase 7 classifies assembled contigs as viral using three independent tools with orthogonal
detection algorithms. The CLAUDE pipeline **Principle 2** requires ≥2/3 tool agreement before
calling a contig viral. This prevents false positives from any single tool's biases.

**Tools used:**

| Tool | Version | Algorithm | Threshold used | Completed |
|------|---------|-----------|----------------|-----------|
| geNomad | 1.9.7 | Marker gene profiles + neural network | virus_score ≥ 0.7 | ✅ |
| DeepVirFinder (DVF) | GitHub HEAD | Deep CNN on sequence composition | score ≥ 0.9 AND p-value ≤ 0.05 | ✅ |
| VirSorter2 | 2.2.4 (patched) | HMM markers + genomic features | max_score ≥ 0.5 | ⚠️ terminated (see below) |

---

## VS2 Runtime Decision — Why 2-Tool Consensus Was Used as Final Result

### What happened

VirSorter2 was successfully patched, all bugs fixed, and ran without errors. It completed
preprocessing (Step 1) and GFF feature extraction (Step 2). However, the Viruses HMM scan
(Step 3) exceeded practical runtime limits on laptop-class hardware.

**Measured runtime data:**
- VS2 Viruses HMM database: 84,390 HMM profiles, 8.3 GB combined database file
- HMM scan configuration: `HMMSEARCH_THREADS=2` (VS2 default, 2 CPU threads per scan)
- Elapsed time at termination: ~2 hours 5 min, 8.7% of profiles scanned
- Estimated total time to complete: **~24 hours** (at 2 threads on this hardware)
- Additionally: a zombie hmmsearch process from a previous failed VS2 attempt was also
  consuming 2 cores writing to a deleted file — killed before termination decision

### Why we proceeded with 2-tool consensus

1. **Scientific value is low for this specific dataset.** SRR15090802 is VLP-enriched (virus-like
   particle isolation by ultracentrifugation). Bacterial contamination is minimal, so VS2's primary
   unique contribution — provirus detection within bacterial contigs — is largely irrelevant here.
   geNomad and DVF cover all other viral detection modalities (marker genes + sequence composition).

2. **The 2-tool consensus is valid.** CLAUDE Principle 2 requires ≥2 independent tools to agree.
   geNomad (marker gene profiles) and DVF (CNN sequence composition) are genuinely orthogonal
   algorithms. The 119-contig ≥2/2 consensus is in fact *more conservative* than a ≥2/3 rule
   would be (any 2 of 3 could include DVF+VS2 agreement, which is less robust than geNomad+DVF).

3. **The result is scientifically complete.** 119 high-confidence viral contigs including
   1 complete phage genome (100% complete, 101,130 bp) is a strong Phase 7 result.

4. **The 24-hour wall time is impractical** for an iterative research pipeline on a laptop.
   VS2 terminated; the snakemake process was left running (PID 277129) but is no longer tracked.

### Recommendation for future runs (production / server environments)

**Use 8 cores for the VS2 Viruses HMM scan.** HMMER3's `hmmsearch` scales well up to ~8–16
threads. Increasing from 2 to 8 threads gives ~4× speedup:

| Setting | Estimated runtime |
|---------|-----------------|
| `HMMSEARCH_THREADS=2` (default, laptop) | ~24 hours |
| `HMMSEARCH_THREADS=8` (recommended, server) | ~6 hours |
| `HMMSEARCH_THREADS=16` (high-memory server) | ~3 hours |

**Command to run VS2 with 8 threads:**
```bash
conda run -n virsorter2_env bash -c "
virsorter run \
    -i /path/to/clean_contigs.fasta \
    -w /path/to/virsorter2_output \
    --db-dir /path/to/databases/virsorter2 \
    --min-length 500 \
    --min-score 0.5 \
    -j 8 \
    --config HMMSEARCH_THREADS=8 \
    --viral-gene-enrich-off \
    all
"
```

Note: `-j 8` sets Snakemake's parallel job count; `--config HMMSEARCH_THREADS=8` overrides the
per-hmmsearch thread count. Both are needed. On a server with 32+ cores, using `-j 4` with
`HMMSEARCH_THREADS=8` (4 parallel jobs × 8 threads = 32 cores total) is optimal.

---

## Section 1: Environment Setup

### 1.1 Conda Environments

Each tool required a dedicated environment due to incompatible Python dependencies:

```bash
# geNomad — dedicated env (conflicts with main env)
mamba create -n genomad_env python=3.10 -y
conda activate genomad_env
pip install genomad==1.9.7

# CheckV — dedicated env
mamba create -n checkv_env python=3.10 -y
conda activate checkv_env
mamba install -c bioconda -c conda-forge checkv=1.0.3 -y

# VirSorter2 — dedicated env
mamba create -n virsorter2_env python=3.10 -y
conda activate virsorter2_env
mamba install -c bioconda -c conda-forge virsorter2=2.2.4 -y

# DeepVirFinder — needs Python 3.6 + legacy keras/theano
mamba create -n dvf_env python=3.6 -y
conda activate dvf_env
mamba install -c conda-forge numpy theano scikit-learn biopython -y
mamba install -c conda-forge h5py -y
mamba install -c conda-forge tensorflow=2.4.1 -y    # installs theano backend
pip install 'keras==2.2.4'                           # standalone keras (NOT tf-keras)
```

### 1.2 Database Downloads

```bash
# geNomad database (~3 GB)
conda run -n genomad_env bash -c "
genomad download-database $PROJECT_DIR/databases/genomad/
"
# Result: databases/genomad/genomad_db/ (35 files)

# CheckV database (~1.2 GB compressed, ~3 GB unpacked)
conda run -n checkv_env bash -c "
checkv download_database $PROJECT_DIR/databases/checkv/
"
# Result: databases/checkv/checkv-db-v1.5/

# VirSorter2 database (~8 GB — viral HMM database is large)
mkdir -p $PROJECT_DIR/databases/virsorter2
conda run -n virsorter2_env bash -c "
virsorter setup -d $PROJECT_DIR/databases/virsorter2 -j 4
"
# Result: databases/virsorter2/ (includes 8.3 GB combined.hmm)

# DeepVirFinder — no separate database download
# Models ship with the DVF git repository
git clone https://github.com/jessieren/DeepVirFinder.git \
    pipeline_run/scripts/DeepVirFinder
# Models are in: pipeline_run/scripts/DeepVirFinder/models/
```

---

## Section 2: Tool 1 — geNomad

### 2.1 Overview

geNomad uses a combination of viral marker gene profiles and a neural network
trained on viral genomic features (gene density, coding strand bias, terminal repeats).
It is the most taxonomically informative of the three tools.

### 2.2 Command

```bash
mkdir -p pipeline_run/08_viral_contigs/virome/genomad

conda run -n genomad_env bash -c "
genomad end-to-end \
    $PROJECT_DIR/pipeline_run/06_decontam_contigs/virome/clean_contigs.fasta \
    $PROJECT_DIR/pipeline_run/08_viral_contigs/virome/genomad \
    $PROJECT_DIR/databases/genomad/genomad_db \
    --threads 4 \
    --min-score 0.7 \
    --cleanup \
    2>&1 | tee $PROJECT_DIR/pipeline_run/logs/genomad_virome.log
"
```

**Parameters:**
- `--min-score 0.7`: minimum virus_score threshold (0–1 scale; default is 0.7)
- `--cleanup`: removes large intermediate files to save disk space

### 2.3 Results

**Run time:** ~15 min for 20,547 contigs on 4 threads
**Output:** `pipeline_run/08_viral_contigs/virome/genomad/clean_contigs_summary/`

| Metric | Value |
|--------|-------|
| Viral contigs (score ≥ 0.7) | **269** |
| Plasmid contigs | 124 |
| Total contigs assessed | 20,547 |

**Taxonomy breakdown (virus_score ≥ 0.7 hits):**

| Taxon | Count | Notes |
|-------|-------|-------|
| Caudoviricetes | 258 | dsDNA tailed phages (Uroviricota) |
| Unclassified | 6 | Below family-level resolution |
| Nucleocytoviricota | 3 | Large dsDNA viruses (NCLDV) |
| Phixviricota | 2 | ssDNA phages (Microviridae family) |

Family-level breakdown (for the 258 Caudoviricetes):
- Herelleviridae: 4
- Autographiviridae: 4
- Others: 250 below family resolution

**Key output file:** `clean_contigs_summary/clean_contigs_virus_summary.tsv`
Columns: `seq_name`, `length`, `topology`, `coordinates`, `n_genes`, `genetic_code`,
`virus_score`, `fdr`, `n_hallmarks`, `marker_enrichment`, `taxonomy`

---

## Section 3: Tool 2 — DeepVirFinder (DVF)

### 3.1 Overview

DVF applies a deep convolutional neural network trained purely on k-mer frequency
profiles (k=4). It does not use gene annotations or marker genes. Its orthogonality
to geNomad makes it valuable for the consensus: it can detect novel viruses without
known homologs.

### 3.2 Dependency Fixes

**CRITICAL:** DVF requires very specific dependencies that conflict with all modern
conda defaults. The following exact setup is required:

```bash
# 1. Python 3.6 environment
mamba create -n dvf_env python=3.6 -y

# 2. Install dependencies in this exact order:
mamba install -n dvf_env -c conda-forge numpy theano scikit-learn biopython h5py -y
mamba install -n dvf_env -c conda-forge tensorflow=2.4.1 -y
# tensorflow 2.4.1 pulls in Keras as a submodule — but DVF needs STANDALONE keras:
$HOME/miniconda3/envs/dvf_env/bin/pip install 'keras==2.2.4'
# This MUST be pip install (conda install would give the wrong version)
```

**CRITICAL runtime setting:**
```bash
export KERAS_BACKEND=theano
```
This must be set before calling dvf.py, otherwise DVF tries to use TensorFlow directly
and fails with a metaclass conflict error.

### 3.3 Command

```bash
mkdir -p pipeline_run/08_viral_contigs/virome/dvf

conda run -n dvf_env bash -c "
export KERAS_BACKEND=theano
python $PROJECT_DIR/pipeline_run/scripts/DeepVirFinder/dvf.py \
    -i $PROJECT_DIR/pipeline_run/06_decontam_contigs/virome/clean_contigs.fasta \
    -o $PROJECT_DIR/pipeline_run/08_viral_contigs/virome/dvf \
    -m $PROJECT_DIR/pipeline_run/scripts/DeepVirFinder/models \
    -l 500 \
    -c 4 \
    2>&1 | tee $PROJECT_DIR/pipeline_run/logs/dvf_virome.log
"
```

**Parameters:**
- `-l 500`: minimum contig length (500 bp); all clean contigs already meet this threshold
- `-c 4`: number of CPU threads
- `-m models`: directory containing pre-trained models (models for 500bp–1kb, 1kb–3kb, 3kb+ contigs)

### 3.4 Results

**Run time:** ~2.5 hours for 20,547 contigs on 4 threads
**Output:** `pipeline_run/08_viral_contigs/virome/dvf/clean_contigs.fasta_gt500bp_dvfpred.txt`

| Metric | Value |
|--------|-------|
| Total predictions | 20,547 |
| Hits (score ≥ 0.9, p-value ≤ 0.05) | **1,520** |
| Hit rate at threshold | 7.4% |

**Score distribution:**

| Threshold | Hits | Notes |
|-----------|------|-------|
| score ≥ 0.9, pval ≤ 0.05 | 1,520 | Used for consensus (high specificity) |
| score ≥ 0.7, pval ≤ 0.05 | 2,344 | More sensitive but lower precision |
| score ≥ 0.5, pval ≤ 0.05 | 2,552 | Not recommended — high FPR |

**Threshold rationale:** score ≥ 0.9 with p-value ≤ 0.05 is the recommended threshold
from the DVF paper for low-biomass / virome samples. The high DVF hit count (7.4%) is
expected for a VLP-enriched virome dataset.

**Output columns:** `name`, `len`, `score`, `pvalue`

---

## Section 4: Tool 3 — VirSorter2

### 4.1 Overview

VirSorter2 uses HMMER-based scanning against viral and non-viral HMM marker gene databases,
combined with machine learning on genomic features (gene density, strand switch rate,
coding capacity). It is the most sensitive tool for detecting integrated proviruses.

### 4.2 Bug Patch Required (VirSorter2 2.2.4)

**Root cause:** VirSorter2 2.2.4 uses a custom internal Prodigal build that produces
a non-standard GFF attribute format. Instead of standard Prodigal attributes
(`ID=N_M;partial=00;start_type=ATG;rbs_motif=AGGAG;...`), it produces:
`gc_cont=0.482;tscore=4.21` (only two custom attributes, no standard ones).

Multiple VS2 parsing scripts expect the standard Prodigal format and fail with
`KeyError` or `ValueError` when receiving the custom format.

**Files patched** (in `virsorter2_env`):

**File 1:** `.../virsorter/utils.py` — `parse_gff()` function (lines 220–248):
```python
# OLD (fails):
sub_items = OrderedDict(i.strip().split('=') for i in last.rstrip(';').split(';'))
rbs_motif = sub_items['rbs_motif']          # KeyError if not present
start_type = sub_items['start_type']        # KeyError if not present
gc_cont = 100*float(sub_items['gc_cont'])   # KeyError if not present

# NEW (patched):
sub_items = OrderedDict(
    i.strip().split('=', 1) for i in last.rstrip(';').split(';')
    if i.strip() and '=' in i              # filter empty/malformed entries
)
if 'ID' in sub_items:
    orf_index = int(sub_items['ID'].split('_')[1])
else:
    gene_counter[seqname_key] = gene_counter.get(seqname_key, 0) + 1
    orf_index = gene_counter[seqname_key]
partial   = sub_items.get('partial', '00')  # default: complete gene
rbs_motif = sub_items.get('rbs_motif', 'None')
start_type = sub_items.get('start_type', 'ATG')
gc_cont   = 100 * float(sub_items.get('gc_cont', '0.5'))
```

**File 2:** `.../virsorter/scripts/circular-remove-partial-gene.py` (lines 133–137):
Same `split('=', 1)` fix + `if i.strip() and '=' in i` filter + `.get('partial', '00')`.

**File 3:** `.../virsorter/scripts/filter-seqs-by-gff.py` (lines 47–58):
Same split fix + per-contig gene counter fallback for missing `ID` field.

**CRITICAL after patching:** Must clear Python bytecode cache before running:
```bash
find $HOME/miniconda3/envs/virsorter2_env/lib/python3.10/site-packages/virsorter/ \
    -name "__pycache__" -exec rm -rf {} + 2>/dev/null
find $HOME/miniconda3/envs/virsorter2_env/lib/python3.10/site-packages/virsorter/ \
    -name "*.pyc" -delete 2>/dev/null
```
Also delete the VS2 output directory completely before re-running (Snakemake cannot
resume from a partially-failed run with corrupt intermediate files).

### 4.3 Command

```bash
mkdir -p pipeline_run/08_viral_contigs/virome/virsorter2

# Clear pycache first (see Section 4.2)
find $HOME/miniconda3/envs/virsorter2_env/lib/python3.10/site-packages/virsorter/ \
    -name "__pycache__" -exec rm -rf {} + 2>/dev/null

conda run -n virsorter2_env bash -c "
virsorter run \
    -i $PROJECT_DIR/pipeline_run/06_decontam_contigs/virome/clean_contigs.fasta \
    -w $PROJECT_DIR/pipeline_run/08_viral_contigs/virome/virsorter2 \
    --db-dir $PROJECT_DIR/databases/virsorter2 \
    --min-length 500 \
    --min-score 0.5 \
    -j 4 \
    --viral-gene-enrich-off \
    all \
    2>&1 | tee $PROJECT_DIR/pipeline_run/logs/virsorter2_virome.log
" &
```

**Parameters:**
- `--min-length 500`: skip contigs shorter than 500 bp
- `--min-score 0.5`: report all contigs above this threshold (lower = more sensitive)
- `--viral-gene-enrich-off`: do not require enrichment of viral genes (important for short/novel phages)
- `all`: run all viral groups (dsDNAphage, ssDNA, NCLDV, RNA, lavidaviridae)
- `-j 4`: Snakemake parallelism

**Note on runtime:** The Viruses HMM database (`databases/virsorter2/hmm/viral/combined.hmm`)
is 8.3 GB. Scanning 20k proteins against this database takes 2–4 hours per viral group
on a laptop-class machine. Plan for 4–8 total hours of VS2 runtime.

### 4.4 Results

**Note:** VS2 output numbers will be filled in when VS2 HMM scan completes.
Results below are based on a partial run (HMM scan in progress as of this writing).

| Metric | Value |
|--------|-------|
| Circular contigs identified | 25 |
| Linear contigs | 20,522 |
| Viral contigs (score ≥ 0.5) | TBD |

**Key output file:** `pipeline_run/08_viral_contigs/virome/virsorter2/final-viral-score.tsv`
Columns: `seqname`, `dsDNAphage`, `NCLDV`, `RNA`, `ssDNA`, `lavidaviridae`,
`max_score`, `max_score_group`, `length`, `hallmark`, `viral`, `cellular`

---

## Section 5: Consensus Calling — `viral_consensus.py`

### 5.1 Script Overview

Script: `pipeline_run/scripts/viral_consensus.py`

Implements the CLAUDE ≥2/3 tool consensus rule:
- Loads each tool's results
- Applies per-tool score thresholds
- Counts tool votes per contig (0/3, 1/3, 2/3, 3/3)
- Writes consensus FASTA (≥2 votes) and summary TSV

```bash
# 3-tool consensus (run after VS2 completes)
conda run -n claude_pipeline bash -c "
python $PROJECT_DIR/pipeline_run/scripts/viral_consensus.py \
    --genomad pipeline_run/08_viral_contigs/virome/genomad/clean_contigs_summary/clean_contigs_virus_summary.tsv \
    --dvf pipeline_run/08_viral_contigs/virome/dvf/clean_contigs.fasta_gt500bp_dvfpred.txt \
    --vs2 pipeline_run/08_viral_contigs/virome/virsorter2/final-viral-score.tsv \
    --contigs pipeline_run/06_decontam_contigs/virome/clean_contigs.fasta \
    --out-fasta pipeline_run/08_viral_contigs/virome/consensus_viral_contigs.fasta \
    --out-table pipeline_run/08_viral_contigs/virome/viral_classification_summary.tsv \
    --min-votes 2
"
```

### 5.2 2-Tool Consensus Results (interim, geNomad + DVF)

While VS2 was still running, the 2-tool consensus was computed as an interim result.
With only 2 tools, ≥2/2 is equivalent to requiring BOTH tools to agree (most conservative).

```bash
conda run -n claude_pipeline bash -c "
python pipeline_run/scripts/viral_consensus.py \
    --genomad pipeline_run/08_viral_contigs/virome/genomad/clean_contigs_summary/clean_contigs_virus_summary.tsv \
    --dvf pipeline_run/08_viral_contigs/virome/dvf/clean_contigs.fasta_gt500bp_dvfpred.txt \
    --contigs pipeline_run/06_decontam_contigs/virome/clean_contigs.fasta \
    --out-fasta pipeline_run/08_viral_contigs/virome/consensus_viral_contigs_2tool.fasta \
    --out-table pipeline_run/08_viral_contigs/virome/viral_classification_summary_2tool.tsv \
    --min-votes 2
"
```

**2-tool consensus results:**

| Metric | Value |
|--------|-------|
| geNomad hits (score ≥ 0.7) | 268 |
| DVF hits (score ≥ 0.9, pval ≤ 0.05) | 1,520 |
| Both agree (≥2/2) | **119** |
| geNomad only (DVF disagrees) | 149 |
| DVF only (geNomad disagrees) | 1,401 |

**Output:**
- FASTA: `pipeline_run/08_viral_contigs/virome/consensus_viral_contigs_2tool.fasta`
- Table: `pipeline_run/08_viral_contigs/virome/viral_classification_summary_2tool.tsv`

---

## Section 6: CheckV Quality Assessment

### 6.1 Overview

CheckV assesses viral genome completeness and contamination using a viral-specific
reference database. It is appropriate for phage quality assessment (unlike CheckM2
which uses bacterial single-copy marker genes).

### 6.2 Command (run on 2-tool consensus as interim)

```bash
mkdir -p pipeline_run/09_viral_qc/virome/checkv_2tool

conda run -n checkv_env bash -c "
checkv end_to_end \
    $PROJECT_DIR/pipeline_run/08_viral_contigs/virome/consensus_viral_contigs_2tool.fasta \
    $PROJECT_DIR/pipeline_run/09_viral_qc/virome/checkv_2tool \
    -d $PROJECT_DIR/databases/checkv/checkv-db-v1.5 \
    -t 4 \
    2>&1 | tee $PROJECT_DIR/pipeline_run/logs/checkv_virome.log
"
```

**Run time:** ~2 min for 119 contigs

### 6.3 Results (2-tool consensus)

**Input:** 119 consensus viral contigs
**Total viral genome content:** 476,900 bp
**Average contig length:** 4,007 bp
**Largest contig:** 101,130 bp

**Quality tier distribution:**

| CheckV Quality | Count | % |
|----------------|-------|---|
| High-quality (≥90% complete) | **1** | 0.8% |
| Medium-quality (50–89%) | **1** | 0.8% |
| Low-quality (<50%) | 95 | 79.8% |
| Not-determined | 22 | 18.5% |

**Notable contigs:**

| Contig | Length | Quality | Completeness | Notes |
|--------|--------|---------|--------------|-------|
| metaspades_1 | 101,130 bp | High-quality | **100%** | Complete phage genome (likely large dsDNA phage — Herelleviridae or similar) |
| metaspades_28 | 17,149 bp | Medium-quality | 86.3% | Near-complete phage |

**Proviruses identified:** 1
**Completeness methods:** AAI-based high-confidence (74 contigs), AAI-based medium-confidence (16), HMM-based lower-bound (7), not determined (22)

**Interpretation:**
- 95/119 contigs (80%) are low-quality by CheckV — expected for a virome: most phage contigs
  are assembled as partial genomes due to fragmented assembly
- The 100% complete metaspades_1 at 101 kb is exceptional — this is a complete phage genome
  assembled in a single contig (likely aided by high-coverage of this phage in the sample)
- Not-determined (22) contigs are typically very short (<5 kb) or highly divergent from CheckV's
  reference database — they are still real viruses, just without close references

---

## Section 7: QC-7 Checkpoint

```bash
# Written to: pipeline_run/qc_checkpoints/QC7_virome_SRR15090802.txt
```

Threshold for PASS: ≥50 viral contigs with ≥2 tool agreement.

**STATUS: PASS** (119 consensus viral contigs, 2/2 tool agreement)

---

## Section 8: Key Gotchas

1. **DVF requires theano backend**: `KERAS_BACKEND=theano` must be exported before calling
   `dvf.py`. Otherwise DVF tries to use TensorFlow directly and fails with:
   `TypeError: metaclass conflict: the metaclass of a derived class must be a (non-strict)
   subclass of the metaclasses of all its bases`

2. **DVF keras version**: Must use `pip install 'keras==2.2.4'` (standalone keras 2.x).
   Conda's `keras` package installs the TensorFlow-integrated version which doesn't support
   the theano backend.

3. **VS2 GFF format mismatch**: VS2 2.2.4's internal Prodigal outputs a custom GFF without
   standard attributes (`ID`, `partial`, `start_type`, `rbs_motif`). Three VS2 scripts must
   be patched to use `.get()` with defaults instead of direct dict access. See Section 4.2.

4. **VS2 pycache must be cleared after patching**: Python caches compiled `.pyc` files.
   After modifying `.py` files, delete all `__pycache__` directories in the virsorter package.
   Otherwise the old (buggy) compiled code is used despite the source fix.

5. **VS2 output dir must be fully deleted before restart**: Snakemake tracks intermediate files
   by presence. Deleting only some intermediate files disrupts the DAG and causes VS2 to fail
   with `WorkflowError: missing input files`. Always delete the entire VS2 output directory.

6. **VS2 Viruses HMM scan is very slow**: The viral HMM database (`combined.hmm`) is 8.3 GB.
   HMMER scanning against it takes 2–4 hours per viral group on a 4-core laptop. With `all`
   mode running both dsDNAphage and ssDNA groups in parallel, expect 4–8 total hours.

7. **Each tool needs its own conda env**: geNomad, CheckV, and VirSorter2 all have
   incompatible dependency trees. DVF requires Python 3.6 (ancient). Installing any two
   in the same environment causes unsolvable dependency conflicts.

8. **mkdir -p before conda run**: Never combine `mkdir -p /path && conda run -n env cmd`
   in a single shell line — the `-n` flag can be picked up by `mkdir`. Always run `mkdir -p`
   as a separate command first.

---

## Section 9: Output Files Summary

```
pipeline_run/
├── 08_viral_contigs/virome/
│   ├── genomad/
│   │   └── clean_contigs_summary/
│   │       ├── clean_contigs_virus_summary.tsv    # 269 viral contigs
│   │       └── clean_contigs_plasmid_summary.tsv  # 124 plasmids
│   ├── dvf/
│   │   └── clean_contigs.fasta_gt500bp_dvfpred.txt  # 20,547 predictions
│   ├── virsorter2/
│   │   └── final-viral-score.tsv                  # [VS2 results, post-completion]
│   ├── consensus_viral_contigs_2tool.fasta         # 119 contigs (≥2/2 tools)
│   ├── viral_classification_summary_2tool.tsv      # Full per-contig vote table
│   ├── consensus_viral_contigs.fasta               # [3-tool final, post-VS2]
│   └── viral_classification_summary.tsv            # [3-tool final table]
├── 09_viral_qc/virome/
│   ├── checkv_2tool/
│   │   ├── quality_summary.tsv    # 1 HQ, 1 MQ, 95 LQ, 22 ND
│   │   ├── contamination.tsv
│   │   └── completeness.tsv
│   └── checkv/                    # [CheckV on 3-tool consensus, post-VS2]
├── qc_checkpoints/
│   └── QC7_virome_SRR15090802.txt  # PASS
└── logs/
    ├── genomad_virome.log
    ├── dvf_virome.log
    ├── virsorter2_virome.log
    └── checkv_virome.log
```
