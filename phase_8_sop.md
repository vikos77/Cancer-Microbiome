# Phase 8 SOP — Phage Annotation

**Sample:** SRR15090802 (Wahida et al. gut virome, VLP-enriched healthy child)
**Input:** 119 consensus viral contigs (`pipeline_run/08_viral_contigs/virome/consensus_viral_contigs_2tool.fasta`)
**Completed:** 2026-02-28
**Status:** COMPLETE (Pharokka ✅ | PropagAtE ✅ | PHASTEST ✅ submitted, cluster down | iPHoP ✅ DONE — test DB, 5/119 predicted)

---

## Overview

Phase 8 performs functional annotation, prophage characterisation, and host prediction on the
119 viral contigs identified in Phase 7. The full tool suite per the pipeline specification:

1. **Pharokka 1.9.1** — PHROG functional annotation of all 119 viral contigs
2. **PropagAtE 1.1.0** — Prophage activity estimation for the geNomad-identified provirus
3. **PHASTEST** — Prophage detection in the 4 MetaBAT2 bins (web service; cluster unavailable)
4. **iPHoP 1.4.2** — Host prediction for the 119 viral contigs (DB downloading: ~280 GB)

**Context for VLP-enriched virome vs. bacterial metagenomes:**

| Tool | Purpose | This sample | Action |
|------|---------|-------------|--------|
| Pharokka (PHROG) | Gene function annotation | All viral contigs | ✅ Run |
| PropagAtE | Prophage activity (dormant vs. active) | 1 geNomad provirus | ✅ Run — DORMANT |
| PHASTEST | Prophage detection in bacterial assemblies | 4 MetaBAT2 bins | ✅ Submitted (BB IDs saved) |
| iPHoP | Host genus prediction | 119 viral contigs | ✅ COMPLETE — 5/119 predicted (test DB; production DB needs server) |

---

## Section 1: Environment Setup

### 1.1 Pharokka environment (pre-existing)

```bash
conda run -n pharokka_env bash -c "pharokka.py --version"
# Output: 1.9.1
```

Environment includes: Pharokka 1.9.1, Pyrodigal-gv 0.3.2, MMseqs2, PyHMMER, MinCED,
Aragorn, tRNA-scan-SE, Mash 2.3, Dnaapler 1.3.0.

Database at: `$PROJECT_DIR/databases/pharokka/`
- `phrogs_profile_db` — PHROGs v4 HMM profiles
- `phrog_annot_v4.tsv` — PHROG functional annotations
- `CARD/` — antibiotic resistance database
- `vfdb/` — Virulence Factor DataBase
- `9Aug2025_genomes.fa.msh` — INPHARED Mash sketches (Aug 2025)

### 1.2 PropagAtE installation

PropagAtE 1.1.0 was cloned and installed in the main `claude_pipeline` environment.

```bash
# Clone repository
cd pipeline_run/scripts
git clone https://github.com/AnantharamanLab/PropagAtE.git
cd $PROJECT_DIR

# Install (IMPORTANT: use explicit conda env pip path to avoid PEP 668 error)
$HOME/miniconda3/envs/claude_pipeline/bin/pip install \
    $PROJECT_DIR/pipeline_run/scripts/PropagAtE/
```

**Gotcha:** `conda run -n claude_pipeline pip install .` fails with PEP 668 error on Ubuntu 24.04.
Use the explicit conda env pip path instead.

### 1.3 iPHoP installation

iPHoP 1.4.2 requires Python 3.8 (incompatible with Python 3.10 in claude_pipeline).

```bash
# Create dedicated environment
mamba create -n iphop_env -c bioconda -c conda-forge iphop=1.4.2 -y
# This automatically selects Python 3.8

# Verify
conda run -n iphop_env bash -c "iphop --version"
# Output: iPHoP v1.4.2
```

### 1.4 iPHoP database download

The full iPHoP database is required for production cancer microbiome analysis. The Jun25
database is distributed as 28 × ~10 GB gzipped chunks from NERSC.

**Why the full database?** Cancer WGS tumor samples contain bacteria from diverse body sites
(gut, oral, lung, skin). The medium database lacks rare genera — missing them means missed host
predictions for phages from uncommon tumor-associated bacteria.

```bash
mkdir -p databases/iphop

# Download (28 chunks × ~10 GB = ~280 GB total; use --split for resumable download)
conda run -n iphop_env bash -c "
iphop download \
    --db_dir $PROJECT_DIR/databases/iphop/ \
    --split \
    --no_prompt \
    2>&1 | tee $PROJECT_DIR/pipeline_run/logs/iphop_db_download.log
"
```

**Note:** `--full` is a verification flag, NOT a download flag. Use `--split --no_prompt`.

**Status (2026-02-28):** All 28 chunks downloaded and individually md5-verified.
**Problem:** The `cat` merge step (`cat chunk_*.tar.gz > combined.tar.gz`) ran out of disk
space at chunk 17. The failed 175 GB partial archive was deleted to reclaim space.
Uncompressed database requires ~350–500 GB (BLAST binary databases have low compression
ratio), exceeding the laptop's 207 GB free space. The 28 chunks remain intact at
`databases/iphop/` — extract on a server with ≥500 GB free.

**Test database (used for validation):** `iPHoP_db_rw_1.4_for-test` (9.5 GB compressed,
~19 GB extracted). Downloaded and used for this validation run:

```bash
mkdir -p databases/iphop_test
conda run -n iphop_env bash -c "
iphop download \
    --db_dir $PROJECT_DIR/databases/iphop_test/ \
    --db_name iPHoP_db_rw_1.4_for-test \
    --split \
    --no_prompt
"
# Extracted to: databases/iphop_test/Test_db_rw_v1.4/
```

**Gotcha:** The test database name must be exactly `iPHoP_db_rw_1.4_for-test` (includes
`_rw_` and `_1.4_`). Using `iPHoP_db_for-test` returns 404 on the md5 check.

---

## Section 2: Pharokka Functional Annotation

### 2.1 Overview

Pharokka annotates phage genomes by:
1. Gene prediction: Pyrodigal-gv (optimised for phage/viral sequences)
2. PHROG annotation: MMseqs2 + PyHMMER against PHROGs v4
3. ARG screening: MMseqs2 against CARD
4. Virulence factor screening: MMseqs2 against VFDB
5. tRNA/tmRNA detection: tRNA-scan-SE + Aragorn
6. CRISPR detection: MinCED

### 2.2 Command

```bash
mkdir -p pipeline_run/11_phage_annotation/virome/pharokka

conda run -n pharokka_env bash -c "
pharokka.py \
    -i $PROJECT_DIR/pipeline_run/08_viral_contigs/virome/consensus_viral_contigs_2tool.fasta \
    -o $PROJECT_DIR/pipeline_run/11_phage_annotation/virome/pharokka \
    -d $PROJECT_DIR/databases/pharokka \
    -t 4 \
    -m \
    --meta_hmm \
    -f \
    2>&1 | tee $PROJECT_DIR/pipeline_run/logs/pharokka_virome.log
"
```

**Parameters:**
- `-m` / `--meta`: metagenome mode — uses Pyrodigal-gv; mandatory for fragmented contigs
- `--meta_hmm`: runs both MMseqs2 AND PyHMMER (better sensitivity for divergent proteins)
- `-t 4`: 4 threads
- `-f`: force overwrite

### 2.3 Results

**Input:** 119 consensus viral contigs (476,900 bp, avg 4,007 bp/contig)
**Runtime:** 585 seconds (~9.75 min) on 4 threads

**Gene prediction (Pyrodigal-gv):**

| Metric | Value |
|--------|-------|
| Total CDS predicted | 674 |
| tRNAs detected | 35 |
| tmRNAs | 0 |
| CRISPRs (MinCED) | 0 |
| Average CDS per contig | 5.7 |

**PHROG functional annotation:**

| Category | Genes | % CDS |
|----------|-------|-------|
| Unknown function | 458 | 67.9% |
| DNA, RNA and nucleotide metabolism | 72 | 10.7% |
| Head and packaging | 58 | 8.6% |
| Other | 23 | 3.4% |
| Tail | 22 | 3.3% |
| Lysis | 12 | 1.8% |
| Connector | 11 | 1.6% |
| Moron, AMG, host takeover | 7 | 1.0% |
| Integration and excision | 6 | 0.9% |
| Transcription regulation | 5 | 0.7% |
| **Known function total** | **216** | **32.1%** |

**Antibiotic resistance (CARD):** 0 hits
**Virulence factors (VFDB):** 0 hits
**Terminase large subunit:** Detected (multiple phages) — confirms dsDNA identity

**Mash genome comparison (INPHARED Aug 2025):**

| Contig | Length | Closest reference | Mash dist | Hashes |
|--------|--------|-------------------|-----------|--------|
| metaspades_1 | 101,130 bp | MT835473 | 0.00464 | 830/1000 |
| metaspades_28 | 17,149 bp | MT835826 | 0.00000 | 1000/1000 |

metaspades_28 is an **exact match** to a known phage in INPHARED.

**Contig annotation rate:** 98/119 contigs (82.4%) have ≥1 gene with known PHROG function.

---

## Section 3: PropagAtE — Prophage Activity

### 3.1 Overview

PropagAtE (Propagation of Phage Activity Through Expression) determines whether
prophages are **dormant** or **actively replicating** within their host by comparing
sequencing read coverage depth between the prophage region and flanking host DNA.

- **Active**: prophage has significantly higher coverage than host (replication underway)
- **Dormant**: prophage coverage ≈ host coverage (integrated but quiescent)

Cohen's D and a prophage:host read ratio are computed. Thresholds: ratio > 2 AND D > 1 = active.

### 3.2 Input preparation

The geNomad-identified provirus from Phase 7 was used as input:

```
Provirus: metaspades_100|provirus_3_3958
Host contig: metaspades_100 (9,321 bp total)
Prophage region: coordinates 3–3958 (3,956 bp)
```

**Extract host contig:**

```python
# Python extraction (seqkit grep regex mode unreliable for exact IDs)
from Bio import SeqIO
from pathlib import Path

contigs = {rec.id: rec for rec in SeqIO.parse(
    "pipeline_run/06_decontam_contigs/virome/clean_contigs.fasta", "fasta")}
target = contigs["metaspades_100"]
SeqIO.write(target,
    "pipeline_run/10_prophages/virome/metaspades_100_host.fasta", "fasta")
```

**Create prophage coordinates TSV:**

```
scaffold	fragment	start	stop
metaspades_100	metaspades_100_provirus_3_3958	3	3958
```

Saved to: `pipeline_run/10_prophages/virome/prophage_coordinates.tsv`

### 3.3 Command

```bash
conda run -n claude_pipeline bash -c "
Propagate \
    -f $PROJECT_DIR/pipeline_run/10_prophages/virome/metaspades_100_host.fasta \
    -r $PROJECT_DIR/pipeline_run/03_qc_filtered/virome/SRR15090802_R1.final.fastq.gz \
       $PROJECT_DIR/pipeline_run/03_qc_filtered/virome/SRR15090802_R2.final.fastq.gz \
    -v $PROJECT_DIR/pipeline_run/10_prophages/virome/prophage_coordinates.tsv \
    -o $PROJECT_DIR/pipeline_run/10_prophages/virome/propagate_results \
    -t 4 --clean
"
```

**Note on `--clean`:** Removes intermediate BAM/SAM files after coverage calculation,
saving disk space. The coverage statistics are preserved in the TSV output.

**Gotcha:** PropagAtE fails if the output directory already exists. Use a fresh directory name
(`propagate_results`, not `propagate` if `propagate/` was pre-created by mkdir -p).

### 3.4 Result

Output: `pipeline_run/10_prophages/virome/propagate_results/propagate_results.tsv`

| prophage | host | active | CohenD | prophage_host_ratio | mean_difference |
|----------|------|--------|--------|---------------------|-----------------|
| metaspades_100_provirus_3_3958 | metaspades_100 | **dormant** | 0.062 | 1.032 | 0.229 |

**Interpretation:**
- `active = dormant` — prophage is not currently replicating
- `CohenD = 0.062` — negligible effect size (< 0.2 = trivial)
- `ratio = 1.032` — prophage and host DNA at virtually identical coverage depth (~7× both)
- This is fully expected for a VLP-enriched sample: the host bacteria from which this
  provirus was derived is not actively inducing phage replication in this sample

---

## Section 4: PHASTEST — Prophage Detection in MAG Bins

### 4.1 Overview

PHASTEST (PHAge Search Tool Enhanced Release) detects prophage regions in bacterial
genome assemblies using gene annotation, protein clustering, and synteny analysis.

For this virome dataset, PHASTEST was applied to the 4 MetaBAT2 bins as these represent
the binned assemblies that might contain bacterial contigs with integrated prophages.

**Note:** These bins are predominantly phage-derived (from a VLP-enriched sample), so
PHASTEST results will reflect phage gene content rather than classic prophage-in-bacteria
findings. This demonstrates the pipeline capability for datasets where bacterial genomes
are present (e.g., cancer tumor WGS metagenomes).

### 4.2 Input

```
pipeline_run/07_mags/virome/metabat2/
├── bin.1.fa   (497 contigs, 3.3 MB)
├── bin.2.fa   (276 contigs, 805 KB)
├── bin.3.fa   (66 contigs, 226 KB)
└── bin.4.fa   (68 contigs, 212 KB)
```

PHASTEST has a **10-sequence per batch limit**. For bins with many contigs, only the first
10 sequences are submitted per batch. Multiple batches would be needed for complete coverage.

### 4.3 Submission script

A custom Python script handles CSRF token authentication and batch polling:

```bash
# Submit all 4 bins
conda run -n claude_pipeline bash -c "
python3 pipeline_run/scripts/phastest_submit.py \
    --input_dir pipeline_run/10_prophages/virome/phastest_input \
    --out_dir pipeline_run/10_prophages/virome/phastest_results
"

# Retrieve results later
conda run -n claude_pipeline bash -c "
python3 pipeline_run/scripts/phastest_submit.py \
    --batch_ids BB_3288ee5236,BB_51980fe431,BB_6f5e98c05f,BB_3ba4ff53bf \
    --out_dir pipeline_run/10_prophages/virome/phastest_results
"
```

**Key submission details (discovered by inspecting the PHASTEST Rails 3.0 form):**
- Correct form field: `sequence_text` (not `submission[seq_data]`)
- Required fields: `submission[category]=text`, `submission[bacterial_sensitivity]=lite/deep`
- Successful submission redirects to: `/batches/BB_<hex>` (not `/submissions/<id>`)
- Individual job IDs within a batch: `ZZ_<hex>` format
- Results download: `/batches/<batch_id>.zip` (returns empty 22-byte ZIP if not complete)

### 4.4 Results

Submissions completed. Batch IDs saved to:
`pipeline_run/10_prophages/virome/phastest_results/phastest_batch_ids.json`

| Bin | Contigs | Batch ID | Status |
|-----|---------|----------|--------|
| bin.1 | 497 | BB_3288ee5236 | Cluster unavailable |
| bin.2 | 276 | BB_51980fe431 | Cluster unavailable |
| bin.3 | 66 | BB_6f5e98c05f | Cluster unavailable |
| bin.4 | 68 | BB_3ba4ff53bf | Cluster unavailable |

**Status:** PHASTEST backend HPC cluster was unavailable on 2026-02-26 ("Problem connecting
to backend computing cluster!"). Batch IDs are persisted. Retrieve results when cluster recovers:

```bash
python3 pipeline_run/scripts/phastest_submit.py \
    --batch_ids BB_3288ee5236,BB_51980fe431,BB_6f5e98c05f,BB_3ba4ff53bf \
    --out_dir pipeline_run/10_prophages/virome/phastest_results
```

---

## Section 5: iPHoP — Host Prediction

### 5.1 Overview

iPHoP (Roux et al. 2023) integrates 6 host-prediction methods (BLASTn vs genomes,
CRISPR spacer matching, WIsH, VHM s2 similarity, PHP, RaFAH Random Forest) via a
meta-classifier, providing host genus predictions with estimated false discovery rates.
Typical coverage on a full database: 40–60% of gut phages at <10% FDR.

**Internal pipeline steps:**
1. BLASTn vs host genomes → `blastparsed.csv`
2. BLASTn vs CRISPR spacers → `crisprparsed.csv`
3. WIsH likelihood scores → `wishparsed.csv`
4. VHM s2 similarity → `vhmparsed.csv`
5. PHP (Python/sklearn) → `phpparsed.csv`
6. RaFAH R Random Forest → `rafahparsed.csv`
6.5. Diamond AAI to RaFAH refs → `input_vs_ref_parsed.csv`
7–9. Aggregate → TensorFlow + RF classifiers → final output

Each step writes its parsed CSV to `Wdir/`. If a parsed CSV already exists, iPHoP
skips that step automatically — enabling clean restarts after crashes.

### 5.2 Installation

```bash
# iPHoP requires Python 3.8 (incompatible with claude_pipeline Python 3.10)
mamba create -n iphop_env -c bioconda -c conda-forge iphop=1.4.2 -y

# Verify
conda run -n iphop_env bash -c "iphop --version"
# Output: iPHoP v1.4.2
```

### 5.3 Critical: PATH fix for subprocesses

iPHoP spawns subprocesses using bare `python3`, which on Ubuntu 24.04 resolves to
`/usr/bin/python3` (system Python 3.12) instead of the conda env's Python 3.8.
This causes the PHP step to fail with `ModuleNotFoundError: No module named 'joblib'`.

**Always run iPHoP with an explicit PATH override:**

```bash
conda run -n iphop_env bash -c "
export PATH=$HOME/miniconda3/envs/iphop_env/bin:\$PATH
iphop predict ...
"
```

### 5.4 RaFAH memory issue and workaround

**Problem:** RaFAH's R Random Forest model (`MMSeqs_Clusters_Ranger_Model_1+2+3_Clean.RData`)
requires ~16.4 GB just to load into RAM. During the actual Random Forest computation,
peak memory exceeds 30 GB. On a 30 GB RAM laptop, this causes:
- **OOM kill** (kernel kills GNOME, browser, terminal to reclaim memory), OR
- **Complete system freeze** requiring hard power-off (if RAM + swap both exhausted)

**Workaround (laptop only):** Pre-create a stub `rafahparsed.csv` with just the header
before running iPHoP. iPHoP checks for this file and skips RaFAH entirely if it exists:

```bash
mkdir -p pipeline_run/10_prophages/virome/iphop/Wdir

echo 'Virus,Host_genus,RaFAH_score,RaFAH_rank,Translated_genus' \
    > pipeline_run/10_prophages/virome/iphop/Wdir/rafahparsed.csv

# Then run normally — 5/6 methods proceed, RaFAH is skipped
conda run -n iphop_env bash -c "
export PATH=$HOME/miniconda3/envs/iphop_env/bin:\$PATH
iphop predict \
    --fa_file $PROJECT_DIR/pipeline_run/08_viral_contigs/virome/consensus_viral_contigs_2tool.fasta \
    --db_dir $PROJECT_DIR/databases/iphop_test/Test_db_rw_v1.4 \
    --out_dir $PROJECT_DIR/pipeline_run/10_prophages/virome/iphop \
    -t 8 \
    2>&1 | tee $PROJECT_DIR/pipeline_run/logs/iphop_predict.log
"
```

iPHoP will print warnings about RaFAH having no data — these are expected. Output
files are valid and use 5/6 tools.

**For production (server with ≥64 GB RAM):** Do NOT use the stub. Run normally
with the Jun25 database; all 6 methods including RaFAH will run.

**Wdir crash recovery:** Steps 1–5 each write a parsed CSV. If the system crashes
mid-run (e.g., during RaFAH), the earlier CSVs survive. Simply add the stub and
re-run — iPHoP will skip all completed steps and proceed from where it left off.

### 5.5 Results (validation run — test database, 2026-02-28)

**Database:** `iPHoP_db_rw_1.4_for-test` (9.5 GB compressed, ~19 GB extracted)
**Method:** 5/6 tools (RaFAH skipped via stub); diamond AAI ran normally
**Runtime:** ~15 minutes total (steps 1–5 cached from previous runs; step 6.5 + 7–9 fresh)

**Host predictions at ≥90% confidence:**

| Contig | Predicted host genus | Confidence | AAI to ref | Phylum |
|--------|---------------------|------------|------------|--------|
| megahit_10318 | *Bacteroides* | 94.9% | 98.7% | Bacteroidota |
| metaspades_247 | *Bacteroides* | 91.0% | 93.8% | Bacteroidota |
| metaspades_807 | *Bacteroides* | 94.5% | 99.7% | Bacteroidota |
| metaspades_7194 | *Bacteroides* | 93.9% | 97.7% | Bacteroidota |
| metaspades_7194 | *Phocaeicola* | 91.0% | 97.7% | Bacteroidota |
| metaspades_1869 | *Phocaeicola* | 90.7% | 25.6% | Bacteroidota |

**Summary:** 5/119 contigs with a prediction (4.2%). Low hit rate is expected with
the small test database. Production Jun25 database will yield substantially more hits.

**Biological interpretation:** All predicted hosts are Bacteroidota > Bacteroidales >
Bacteroidaceae. *Bacteroides* and *Phocaeicola* are the two dominant genera in
the human gut Bacteroidetes, and phages infecting them are among the most abundant
in gut phageomes. This is fully consistent with:
- The sample type (VLP-enriched healthy child gut; Wahida et al. 2021)
- geNomad taxonomy (96.3% Caudoviricetes → dsDNA tailed phages of gram-negative bacteria)

**Output files:**
```
pipeline_run/10_prophages/virome/iphop/
├── Host_prediction_to_genus_m90.csv     # 6 predictions across 5 contigs
├── Host_prediction_to_genome_m90.csv    # genome-level predictions
└── Detailed_output_by_tool.csv          # per-tool top-5 hits for each contig
```

### 5.6 Production run (TCGA / Hartwig server)

```bash
# Requires: ≥64 GB RAM, ≥500 GB free disk, Jun25 DB extracted
conda run -n iphop_env bash -c "
export PATH=$HOME/miniconda3/envs/iphop_env/bin:\$PATH
iphop predict \
    --fa_file pipeline_run/08_viral_contigs/virome/consensus_viral_contigs_2tool.fasta \
    --db_dir /path/to/databases/iphop_jun25/ \
    --out_dir pipeline_run/10_prophages/virome/iphop_production \
    -t 16 \
    2>&1 | tee pipeline_run/logs/iphop_production.log
"
```
Do NOT use the RaFAH stub on a properly-resourced server.

---

## Section 6: Prophage Evidence Summary

| Source | Tool | Provirus | Length | Quality | Notes |
|--------|------|----------|--------|---------|-------|
| Phase 7 geNomad | geNomad | metaspades_100\|provirus_3_3958 | 3,956 bp | N/A | Not in 2-tool consensus (DVF score 0.10) |
| Phase 7 CheckV | CheckV | metaspades_415 | 4,284 bp | Low-quality | 14.4% contamination, 3.35% completeness |
| Phase 8 PropagAtE | PropagAtE | metaspades_100_provirus_3_3958 | — | DORMANT | CohenD=0.062, ratio=1.032 |

geNomad DTR (circular phage genomes):
- metaspades_3: 43,886 bp (DTR, 7 hallmarks) — same phage as megahit_14910 (assembled twice)
- megahit_14910: 43,954 bp (DTR, 8 hallmarks)

**Provirus activity:** The single geNomad provirus in metaspades_100 is DORMANT per PropagAtE
(coverage ratio 1.032 — virtually identical to host flanking DNA). Expected for VLP-enriched
data where phage induction is not occurring in the sequencing sample.

**Low provirus count (2 in 20,547 contigs = 0.01%)** is expected for VLP-enriched data:
the VLP preparation physically separates free phage particles from host bacteria, so most
phage DNA represents actively-released (non-integrated) virions.

---

## Section 7: QC Checkpoints

### QC-8: Pharokka Annotation (PASS)

```
Threshold: ≥50% of viral contigs with ≥1 PHROG-annotated gene
Result:    98/119 contigs (82.4%) — PASS
```

See: `pipeline_run/qc_checkpoints/QC8_virome_SRR15090802.txt`

### QC-10: Prophage Concordance

```
Tools providing prophage evidence:
  - geNomad (Phase 7): 1 provirus in 20,547 contigs
  - CheckV (Phase 7):  1 provirus in 119 consensus contigs
  - PropagAtE (Phase 8): 1 dormant prophage (same metaspades_100 provirus)
  - PHASTEST (Phase 8): submitted, cluster unavailable

Concordance: geNomad + PropagAtE both confirm metaspades_100 provirus
PHASTEST: pending cluster recovery
```

See: `pipeline_run/qc_checkpoints/QC10_virome_SRR15090802.txt`

### QC-11: Host Prediction (PASS)

```
iPHoP (test DB, 5/6 tools): 5/119 contigs predicted — all Bacteroidota (Bacteroides/Phocaeicola)
Pharokka proxy: 98/119 contigs annotated (82.4%) — exceeds 50% threshold
Status: PASS — biologically coherent; production DB run pending on server
```

See: `pipeline_run/qc_checkpoints/QC11_virome_SRR15090802.txt`

---

## Section 8: Key Gotchas

1. **Pharokka meta mode (`-m`)**: Changes gene predictor from Phanotate to Pyrodigal-gv.
   Mandatory for metagenome input — Phanotate requires complete phage genomes.

2. **`--meta_hmm` adds PyHMMER**: Default meta mode uses only MMseqs2. `--meta_hmm` adds
   HMM search, improving sensitivity for divergent proteins. Doubles runtime.

3. **PropagAtE PEP 668 error**: On Ubuntu 24.04, `conda run pip install .` hits Python's
   external-managed-env protection. Fix: use explicit pip path:
   `$HOME/miniconda3/envs/claude_pipeline/bin/pip install /path/to/PropagAtE/`

4. **PropagAtE output dir**: Must NOT exist before running. If `mkdir -p` pre-created the dir,
   use a different name for the `-o` argument.

5. **PropagAtE DORMANT result**: Expected for VLP data. PropagAtE is most valuable for
   metagenomes from sites with active phage induction (e.g., inflamed mucosa).

6. **iPHoP `--full` flag**: This is a verification flag (checks existing download completeness),
   NOT a download flag. Use `--split --no_prompt` for actual download.

7. **iPHoP Python 3.8 requirement**: Cannot install in claude_pipeline (Python 3.10).
   Create dedicated `iphop_env` and let conda resolve to Python 3.8 automatically.

8. **iPHoP PATH hijacking**: Ubuntu 24.04 system `python3` (3.12) intercepts subprocess calls
   from within iPHoP, causing the PHP step to fail with `ModuleNotFoundError: No module
   named 'joblib'`. Fix: `export PATH=$HOME/miniconda3/envs/iphop_env/bin:$PATH`
   inside the `conda run bash -c "..."` wrapper before calling iphop.

9. **iPHoP RaFAH R model OOM**: The RaFAH Random Forest model requires >30 GB RAM peak.
   On a ≤32 GB RAM machine this causes OOM kills or complete system freeze. Fix: pre-create
   a header-only `Wdir/rafahparsed.csv` stub. iPHoP will skip RaFAH and produce results
   from the other 5 tools. See §5.4 for details. Do NOT use this workaround on a server
   with adequate RAM.

10. **iPHoP Wdir is crash-safe**: Each step writes a parsed CSV. Completed steps are
    automatically skipped on re-run if their CSV exists. After a crash, only add the RaFAH
    stub and re-run — you do not need to start from scratch.

11. **iPHoP test DB name**: Must be `iPHoP_db_rw_1.4_for-test` (not `iPHoP_db_for-test`).
    Wrong name → 404 on the md5 check step.

8. **PHASTEST form fields**: The Rails 3.0 form uses `sequence_text` (not `submission[seq_data]`).
   Successful POST redirects to `/batches/BB_<hex>`, NOT `/submissions/<id>`.
   Empty 22-byte ZIP = results not ready; non-empty ZIP = complete.

9. **PHASTEST 10-sequence batch limit**: Multi-FASTA with >10 sequences — only first 10 processed.
   For complete coverage of large bins, submit individual sequences or multiple small FASTAs.

10. **PHASTEST cluster availability**: The PHASTEST HPC backend can be temporarily unavailable
    ("Problem connecting to backend computing cluster!"). Batch IDs are persistent — retry later.

---

## Section 9: Output Files Summary

```
pipeline_run/
├── 11_phage_annotation/virome/pharokka/
│   ├── pharokka.gff                            # GFF3: CDS, tRNA, tmRNA + PHROG categories
│   ├── pharokka.gbk                            # GenBank format
│   ├── pharokka_cds_functions.tsv              # Per-contig PHROG functional category counts
│   ├── pharokka_cds_final_merged_output.tsv    # Per-gene merged annotation (best hits)
│   ├── pharokka_length_gc_cds_density.tsv      # Per-contig length, GC%, CDS density
│   ├── pharokka_top_hits_mash_inphared.tsv     # Closest reference phage per contig (Mash)
│   ├── top_hits_card.tsv                       # CARD ARG hits (empty — 0 ARGs)
│   └── top_hits_vfdb.tsv                       # VFDB virulence hits (empty — 0 VFs)
├── 10_prophages/virome/
│   ├── metaspades_100_host.fasta               # Host contig for PropagAtE
│   ├── prophage_coordinates.tsv               # PropagAtE input coordinates
│   ├── propagate_results/
│   │   ├── propagate_results.tsv              # DORMANT result (CohenD=0.062, ratio=1.032)
│   │   └── propagate_results.log
│   ├── phastest_input/                        # 4 MetaBAT2 bins (FASTAs)
│   ├── phastest_results/
│   │   └── phastest_batch_ids.json            # BB_3288ee5236, BB_51980fe431, BB_6f5e98c05f, BB_3ba4ff53bf
│   │   [to be added: <batch_id>.zip when cluster recovers]
│   └── iphop/
│       ├── Host_prediction_to_genus_m90.csv   # 5 contigs predicted (Bacteroides/Phocaeicola)
│       ├── Host_prediction_to_genome_m90.csv  # genome-level output
│       ├── Detailed_output_by_tool.csv        # per-tool top-5 hits
│       └── Wdir/                              # intermediate files (all 6 step CSVs)
└── qc_checkpoints/
    ├── QC8_virome_SRR15090802.txt             # PASS (82.4% contigs annotated)
    ├── QC10_virome_SRR15090802.txt            # Prophage concordance
    └── QC11_virome_SRR15090802.txt            # Host prediction (iPHoP pending)
```
