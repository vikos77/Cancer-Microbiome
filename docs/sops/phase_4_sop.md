# Phase 4 SOP: Read-Based Taxonomic Profiling
## CLAUDE Pipeline: Cancer-Linked Analysis of Underlying DNA Elements

**Date executed:** 2026-02-25
**Sample:** NA12878_chr22 (1000 Genomes Project, GRCh37 alignment, chr22:16050000–17500000 slice)
**Phase:** 4 of 13
**SOP version:** 1.0

---

## 1. Objective

Phase 4 applies two complementary read classifiers to characterize the microbial
content of host-depleted, quality-filtered reads:

1. **Kraken2 + Bracken**: Broad taxonomic survey against a standard 8 GB database
   (bacteria, archaea, viruses, human). Identifies any residual human contamination
   missed by Phase 2, and provides species-level abundance estimates via Bracken.

2. **KrakenUniq**: Viral-specific classification against RefSeq complete viral
   genomes only. The critical addition is unique k-mer counting per taxon; a
   read with only 2–5 matching k-mers is flagged as a false positive, whereas
   Kraken2 would call it a confident hit. This directly implements the Ge et al.
   2025 recommendation for low-biomass metagenomics.

This phase produces QC-4: the "known taxa" checkpoint. A sample with zero
microbial reads passing filters is expected for our test data (chimeric human reads)
and should not alarm; in real tumor WGS data, we would expect 0.35–5% of
host-depleted reads to classify against known bacterial/viral genomes.

---

## 2. Biological Rationale: Why KrakenUniq Over Kraken2 Alone?

The TCGA cancer microbiome controversy (Gihawi et al. 2023) showed that standard
Kraken2 produces false-positive microbial signals in tumor WGS because:

- Short reads share k-mers with both human and microbial genomes
- Reads near repetitive elements get mis-classified at low k-mer thresholds
- k-mer databases contain draft genomes with contaminating human sequences

KrakenUniq adds a **unique k-mer count** per taxon per sample. A genuine infection
produces reads distributed across many unique genomic positions (high unique k-mers,
high coverage). A false positive from shared k-mers produces reads hitting the same
few k-mers repeatedly (low unique k-mers, near-zero coverage).

Ge et al. 2025 further recommends using **complete genomes only** (no draft/WGS) in
the reference database to eliminate contaminated reference sequences. Our KrakenUniq
database was built from `refseq/viral/Complete_Genome` exclusively.

---

## 3. Environment & Software

| Tool | Version | Source | Env |
|------|---------|--------|-----|
| Kraken2 | 2.1.3 | bioconda | claude_pipeline |
| Bracken | 2.9 | bioconda | claude_pipeline |
| KrakenUniq | 1.0.4 | bioconda | claude_pipeline |
| jellyfish | 1.1.12 | bioconda | claude_pipeline (for DB build only) |

```bash
conda run -n claude_pipeline kraken2 --version   # 2.1.3
conda run -n claude_pipeline bracken --version   # (prints usage with version)
conda run -n claude_pipeline krakenuniq --version # 1.0.4
```

---

## 4. Databases

| Database | Location | Size | Contents |
|----------|----------|------|----------|
| Kraken2 k2_standard_08gb | databases/kraken2/ | 8.0 GB | Bacteria, archaea, virus, human |
| KrakenUniq viral | databases/krakenuniq_viral/ | ~21 GB on disk | 14,596 RefSeq complete viral genomes |

### 4.1 KrakenUniq Database Build (documented separately)

The KrakenUniq viral database was built from scratch during Phase 4 execution.
See Section 7 (Problems & Troubleshooting) for the full account of build issues.

```bash
# Step 1: Download taxonomy
krakenuniq-download --db databases/krakenuniq_viral taxonomy

# Step 2: Download RefSeq complete viral genomes only
krakenuniq-download \
    --db databases/krakenuniq_viral \
    --threads 8 \
    refseq/viral/Complete_Genome

# Step 3: Build database
# CRITICAL: export JELLYFISH_BIN='' before calling (see Problem 1 in Section 7)
# Use --threads 1 or --threads 2 to avoid thermal overheating on laptops (see Problem 2)
conda run -n claude_pipeline bash -c "
    export JELLYFISH_BIN=''
    nice -n 19 krakenuniq-build \
        --db databases/krakenuniq_viral \
        --threads 1 \
        --taxids-for-genomes \
        --taxids-for-sequences
"
```

Build produces 6 sequential steps:
| Step | Output file | Description |
|------|-------------|-------------|
| 1 | database.jdb (4 GB) | Jellyfish k-mer count hash |
| 2 | (skipped; no size reduction needed) | |
| 3 | database0.kdb (4 GB) + database.idx (8.1 GB) | Sorted k-mer DB + position index |
| 4 | seqid2taxid.map | Sequence ID to taxon ID mapping |
| 5 | taxDB (132 MB) | Taxonomy database |
| 6 | database.kdb (4 GB) + database.kdb.counts (314 KB) | LCA-assigned database + unique k-mer counts |

---

## 5. Input Files

| File | Location | Pairs |
|------|----------|-------|
| NA12878_chr22_R1.final.fastq.gz | pipeline_run/03_qc_filtered/ | 391 |
| NA12878_chr22_R2.final.fastq.gz | pipeline_run/03_qc_filtered/ | 391 |

---

## 6. Step-by-Step Execution

### Step 1: Install tools

```bash
conda install -n claude_pipeline -c bioconda -c conda-forge \
    krakenuniq=1.0.4 bracken=2.9 -y
```

Note: `jellyfish=1.1.12` must be active for the KrakenUniq DB build. If jellyfish
2.x was installed (e.g., for another tool), downgrade before building:
```bash
conda install -n claude_pipeline -c bioconda "jellyfish=1.1.12" -y
```

### Step 2: Run Kraken2

```bash
SAMPLE="NA12878_chr22"
K2_DB="databases/kraken2"
OUT="pipeline_run/04_read_profiles"

conda run -n claude_pipeline bash -c "
    kraken2 \
        --db ${K2_DB} \
        --paired \
        --gzip-compressed \
        --output ${OUT}/${SAMPLE}_kraken2.out \
        --report ${OUT}/${SAMPLE}_kraken2.report \
        --confidence 0.1 \
        --threads 8 \
        pipeline_run/03_qc_filtered/${SAMPLE}_R1.final.fastq.gz \
        pipeline_run/03_qc_filtered/${SAMPLE}_R2.final.fastq.gz
" 2>&1 | tee pipeline_run/logs/kraken2_phase4_${SAMPLE}.log
```

**Result:**
```
391 sequences (0.08 Mbp) processed in 0.037s
  2 sequences classified (0.51%) → Homo sapiens (taxid 9606)
389 sequences unclassified (99.49%)
```

### Step 3: Run Bracken

Bracken re-estimates species-level abundances by redistributing reads from higher
taxonomic ranks down to species using k-mer distributions.

```bash
conda run -n claude_pipeline bash -c "
    bracken \
        -d databases/kraken2 \
        -i ${OUT}/${SAMPLE}_kraken2.report \
        -o ${OUT}/${SAMPLE}_bracken_S.txt \
        -w ${OUT}/${SAMPLE}_bracken_S.report \
        -r 100 \
        -l S \
        -t 1
"
```

**Result:**
```
Species detected (reads > threshold=1): 1
Total reads kept at species level: 2
  Homo sapiens: 2 reads (fraction_total: 1.000)
Microbial taxa: 0
```

### Step 4: Run KrakenUniq

```bash
KU_DB="databases/krakenuniq_viral"

conda run -n claude_pipeline bash -c "
    nice -n 19 krakenuniq \
        --db ${KU_DB} \
        --paired \
        --output ${OUT}/${SAMPLE}_krakenuniq.out \
        --report-file ${OUT}/${SAMPLE}_krakenuniq.report \
        --threads 2 \
        pipeline_run/03_qc_filtered/${SAMPLE}_R1.final.fastq.gz \
        pipeline_run/03_qc_filtered/${SAMPLE}_R2.final.fastq.gz
"
```

**Result:**
```
3 sequences classified (0.77%)
388 sequences unclassified (99.23%)
```

---

## 7. QC-4 Results

### 7.1 Kraken2 findings

Two reads classified as Homo sapiens (taxid 9606). These are **residual human reads**
that survived Phase 2 host depletion. Root cause analysis:

**Read 1 (SRR622461.68752106):**
- Phase 2 Kraken2: `U` (unclassified); 3 human k-mers out of 67 total = **4.5%**, below 10% confidence threshold
- Phase 4 Kraken2: `C` 9606; same 3 k-mers, but read trimmed to 55 bp in Phase 3 = **3/21 = 14.3%**, above 10% threshold
- Lesson: **Adapter trimming increases Kraken2 classification confidence.** A borderline read becomes unambiguous after the read is shortened.

**Read 2 (SRR622461.77007588):**
- Phase 2 Kraken2: `C` **taxid 1 (root)**; k-mers split between human (9606) and other sequences, LCA resolves to root
- `extract_nonhuman_kraken.py` filtered only exact `taxid == 9606` hits → this read passed through
- Phase 4 Kraken2: `C` 9606; after trimming, unambiguously human
- Lesson: **extract_nonhuman_kraken.py must also filter ancestor taxids** (Homo=9605, Hominidae=9604, etc.) to catch LCA-resolved reads.

**Production fix for extract_nonhuman_kraken.py:**
```python
# Current (incomplete):
if status == "C" and int(taxid_str) == exclude_taxid:

# Production (catches LCA and ancestor hits too):
HUMAN_TAXIDS = {9606, 9605, 9604, 207598, 9526, 314293, 376913,
                9443, 314146, 1437010, 9347, 32525, 40674, 32524,
                32523, 1338369, 8287, 117571, 117570, 7776, 7742,
                89593, 7711}  # All ancestors from Homo to Chordata
if status == "C" and int(taxid_str) in HUMAN_TAXIDS:
```

### 7.2 KrakenUniq findings

Three viral hits, all evaluated against production thresholds
(minimum 50 unique k-mers, minimum 2 reads):

| Taxon | Reads | Unique k-mers | Coverage | Verdict |
|-------|-------|--------------|----------|---------|
| Human endogenous retrovirus K113 | 1 | 27 | 0.003 | FALSE POSITIVE; HERV, part of human genome |
| Ictalurid herpesvirus 1 (fish) | 1 | 4 | 3.5e-05 | FALSE POSITIVE; <50 unique k-mers |
| Pandoravirus macleodensis (amoeba) | 1 | 2 | 1.1e-06 | FALSE POSITIVE; <50 unique k-mers |

**All three hits fail the minimum unique k-mers filter → 0 viral taxa pass.**

**HERV-K is a special case:** Human endogenous retroviruses are ancient viral sequences
integrated permanently into the human genome. Any human WGS sample will produce HERV
hits against a viral database. These **must** be added to an exclusion list in
production, alongside other human-associated sequences (EBV, HHV-6 in HERVs,
adeno-associated viruses from gene therapy, etc.).

---

## 8. Problems & Troubleshooting

### Problem 1: JELLYFISH_BIN: unbound variable

**Symptom:**
```
/libexec/build_db.sh: line 64: JELLYFISH_BIN: unbound variable
```

**Root cause:** The build script uses `set -u` (crash on unset variables). Line 64 reads:
```bash
[[ "$JELLYFISH_BIN" == "" ]] && JELLYFISH_BIN="jellyfish"
```
If `JELLYFISH_BIN` is **unset** (not exported at all), `"$JELLYFISH_BIN"` crashes under
`set -u`. An **empty exported string** is treated as set, so it passes.

**Fix:**
```bash
# Inside bash -c block, before calling krakenuniq-build:
export JELLYFISH_BIN=''   # Empty string = SET variable = passes set -u check
                           # Script then sets it to "jellyfish" automatically
```

**Does NOT work:**
```bash
# Passing JELLYFISH_BIN in the outer shell — conda run isolates environment
JELLYFISH_BIN=jellyfish conda run -n env krakenuniq-build ...  # FAILS
```

---

### Problem 2: Laptop reboots during KrakenUniq database build (thermal shutdown)

**Symptom:** System reboots 2–3 times during `krakenuniq-build`. Occurs during the
jellyfish k-mer counting and db_sort steps. ACPI critical temperature trip point is
103°C; CPU reached this under 8-thread sustained load.

**Contributing factors:**
- AMD Ryzen 5 3550H laptop CPU
- `thermald` service failed to start ("Unsupported cpu model or platform"); the thermal
  throttling daemon was not running
- jellyfish count with `--threads 8` = 100% CPU utilization on all 8 threads
  sustained for multiple minutes
- db_sort similarly CPU- and I/O-intensive (7+ GB of data)

**Diagnosis commands:**
```bash
# Check ACPI trip points (what temp triggers shutdown)
cat /sys/class/thermal/thermal_zone0/trip_point_*_temp | awk '{printf "%.0f°C\n", $1/1000}'

# Check CPU temp in real-time
watch -n 2 "cat /sys/class/thermal/thermal_zone0/temp | awk '{printf \"%.0f°C\n\", \$1/1000}'"

# Check last reboot reason
journalctl -b -1 --no-pager | grep -iE "thermal|critical|acpi" | grep -v "governor\|Registered"
```

**Fix:**
```bash
# Reduce threads to 1 (25% of CPU) + lowest scheduling priority
export JELLYFISH_BIN=''
nice -n 19 krakenuniq-build --db ... --threads 1 ...
```

**Result:** Build completed successfully at 80°C peak (well below 103°C threshold).
Total build time increased from ~5 min (8 threads) to ~35 min (1 thread), but
completed without any reboot.

**Production recommendation:** Run KrakenUniq database builds on a server or
workstation with active cooling. On a laptop: use `--threads 2`, `nice -n 19`,
and a hard flat surface (not fabric) to allow bottom-panel ventilation.

---

### Problem 3: jellyfish version conflict

**Symptom:**
```
KrakenUniq requires Jellyfish version 1 for building the database
jellyfish 2.2.3
```

**Root cause:** KrakenUniq 1.0.4 DB build requires jellyfish 1.x (not 2.x).
During troubleshooting we installed jellyfish 2.2.3, which replaced 1.1.12.

**Fix:**
```bash
conda install -n claude_pipeline -c bioconda "jellyfish=1.1.12" -y
```

Note: jellyfish 1.x is only needed for **building** the KrakenUniq database.
Running KrakenUniq for classification (`krakenuniq`) does not require jellyfish.
After the DB is built, jellyfish can be upgraded back if needed for other tools.

---

### Problem 4: KrakenUniq build can resume from interrupted db_sort step

**Background:** If the build is interrupted during step 3 (db_sort), the build
can resume from step 3 without re-running step 1 (jellyfish, the expensive step).

**How the script decides what to skip:**
```
Step 1 skipped if: database.jdb OR database0.kdb already exists
Step 3 skipped if: database0.kdb already exists
```

**Safe resume procedure after interrupted db_sort:**
```bash
# 1. Delete the partial db_sort output files (incomplete = corrupt)
rm -f databases/krakenuniq_viral/database.idx
rm -f databases/krakenuniq_viral/database0.kdb  # if it exists but is partial

# 2. Keep database.jdb intact (complete step 1 output)

# 3. Re-run krakenuniq-build — it will skip step 1 and restart from step 3
```

Total jellyfish count time with 2 threads: ~2 minutes
Total db_sort time with 1 thread: ~8 minutes
Total remaining steps (4,5,6): ~15 minutes

---

## 9. Key Biological Reasoning

### Why Kraken2 classifications change after read trimming

The Kraken2 confidence score is calculated as:
```
confidence = (k-mers matching classified taxon) / (total k-mers in read)
```

After fastp trims 46 bp of adapter/low-quality bases from a 101 bp read (→ 55 bp):
- The numerator (matching k-mers) stays the same (those k-mers are still present)
- The denominator (total k-mers) drops from ~67 to ~21
- A read at 4.5% confidence (101 bp) becomes 14.3% (55 bp); it crosses the 10% threshold

**Implication for pipeline design:** Trimming should always precede Kraken2
classification. Doing depletion on raw reads (as done in Phase 2) risks under-removing
borderline human reads that become classifiable only after trimming. In production,
consider adding a post-trimming Kraken2 pass (which is exactly what Phase 4 does).

### The HERV problem in cancer microbiome studies

Human Endogenous Retroviruses (HERVs) are ancient retroviral sequences that integrated
into the primate germline 30–100 million years ago. About 8% of the human genome
consists of HERV-derived sequences. They appear in viral databases because:

1. Taxonomically they are classified under Retroviridae
2. Several HERV loci have been sequenced as complete "viral genomes" in RefSeq
3. Some HERV families (especially HERV-K) are transcriptionally active in cancer cells

This creates a dual problem:
- **False positives:** HERV reads from tumor WGS hit the viral database
- **False negatives:** Reads from genuine exogenous retroviruses may be mistakenly
  discarded as "human" by the host depletion step (especially with the phage-masked
  Hostile index, which was chosen to protect phage signals; the same principle applies)

**Production solution:** Maintain an HERV exclusion list for the KrakenUniq viral
database, and in Phase 8 (prophage annotation), HERV-K loci must be explicitly
distinguished from integrated prophage by genomic context (flanking human sequence,
absence of bacterial genes, retroviral integrase structure).

### Why unique k-mer count is the key filter

| Scenario | Reads | Unique k-mers | Coverage | Interpretation |
|----------|-------|--------------|----------|----------------|
| Genuine viral infection | 50–500 | 500–5000 | >0.1% | Real; reads distributed across genome |
| k-mer cross-contamination | 1–3 | 2–10 | <0.001% | False positive; same few k-mers hit repeatedly |
| HERV hit | 2–20 | 15–100 | 0.003% | Borderline; requires genomic context to resolve |
| Sequencing artifact | 1 | 1–3 | <0.0001% | False positive; likely repetitive element |

The minimum thresholds used in CLAUDE pipeline (production):
- **≥50 unique k-mers per taxon** (Ge et al. 2025 recommendation)
- **≥2 reads** per taxon
- **Coverage ≥0.01%** of reference genome

---

## 10. Lessons Learned

1. **Kraken2 Phase 2 and Phase 4 are not redundant.** Phase 4 catches human reads
   that became classifiable after trimming (Phase 3 shortens reads → higher k-mer
   density). Multi-pass approach is genuinely necessary.

2. **extract_nonhuman_kraken.py needs ancestor taxid filtering.** An exact
   `taxid == 9606` match misses LCA-resolved reads (taxid 1, 9605, 9604, etc.).
   Fix: use a set of all Homo sapiens ancestor taxids.

3. **JELLYFISH_BIN must be exported as empty string (not unset).** Bash `set -u`
   treats unset variables differently from empty strings.

4. **KrakenUniq requires jellyfish 1.x to BUILD databases.** Once built, the
   database is self-contained and jellyfish is not needed for classification.

5. **KrakenUniq database build can resume from step 3.** If the jellyfish step
   (step 1) completed but db_sort was interrupted, deleting only `database.idx`
   and re-running resumes from db_sort, saving ~2 minutes of re-computation.

6. **Laptop thermal management matters.** `thermald` was not working for this AMD
   CPU. `nice -n 19 --threads 1` is required for any multi-minute CPU-intensive
   process on this laptop. Always monitor temperature during builds.

7. **HERV hits are expected in any human WGS sample.** They are not a signal of
   viral infection. They must be in an exclusion list for all viral profiling steps.

8. **KrakenUniq viral DB (complete genomes only) = 14,596 sequences.** A full
   viral RefSeq including draft genomes would be ~2× larger but significantly
   more prone to false positives from contaminated draft sequences.

---

## 11. Files Produced

```
pipeline_run/
├── 04_read_profiles/
│   ├── NA12878_chr22_kraken2.out           # Per-read Kraken2 classifications
│   ├── NA12878_chr22_kraken2.report        # Kraken2 taxonomy report
│   ├── NA12878_chr22_bracken_S.txt         # Bracken species abundance table
│   ├── NA12878_chr22_bracken_S.report      # Bracken-adjusted taxonomy report
│   ├── NA12878_chr22_krakenuniq.out        # Per-read KrakenUniq classifications
│   └── NA12878_chr22_krakenuniq.report     # KrakenUniq report (with unique k-mer counts)
├── qc_checkpoints/
│   └── QC4_NA12878_chr22.txt              # QC-4 checkpoint (STATUS: PASS)
├── logs/
│   ├── kraken2_phase4_NA12878_chr22.log
│   ├── krakenuniq_build.log
│   └── cpu_temp_during_build.log          # Thermal monitoring log
databases/
├── kraken2/                               # 8.0 GB k2_standard_08gb (from Phase 2)
└── krakenuniq_viral/                      # 21 GB RefSeq complete viral genomes
    ├── database.kdb     (4 GB)            # LCA k-mer database
    ├── database.kdb.counts (314 KB)       # Unique k-mer counts per taxon
    ├── database0.kdb    (4 GB)            # Sorted k-mer database (intermediate)
    ├── database.idx     (8.1 GB)          # Position index
    ├── database.jdb     (4 GB)            # Jellyfish k-mer hash (intermediate)
    ├── taxDB            (132 MB)          # Taxonomy database
    ├── seqid2taxid.map  (405 KB)          # Sequence ID → taxon ID mapping
    └── library/                           # 617 MB raw viral genome FASTA files
```

---

## 12. QC-4 Thresholds (Production)

| Metric | PASS | WARN | FAIL |
|--------|------|------|------|
| Residual human (Kraken2) | <1% | 1–5% | >5% |
| KrakenUniq hits passing filter | N/A | N/A | N/A |
| Known microbial taxa (real data) | >0 expected | n/a | 0 (genuine metagenome) |
| HERV hits excluded | Yes | n/a | n/a |
| Minimum unique k-mers for viral call | ≥50 | n/a | <50 |

**QC-4 status: PASS.** No genuine microbial signal. Expected for test data.

---

## 13. Required Fix Before Production

Before running Phase 4 on real tumor WGS data, update
`pipeline_run/scripts/extract_nonhuman_kraken.py` to filter all human ancestor
taxids, not just exactly 9606. See Section 8 (Problem solution for Read 2) for the
complete set of human lineage taxids to exclude.

---

## 14. Next Step

Phase 5: Assembly (metaSPAdes + MEGAHIT + metaviralSPAdes)

**Input:** `pipeline_run/03_qc_filtered/NA12878_chr22_R1.final.fastq.gz` (391 pairs)

**Key question:** Can 391 test reads assemble into contigs? Likely no contigs ≥500 bp
from 391 reads, but this validates tool configuration and directory structure for
production.

**Tools needed:** metaSPAdes (SPAdes ≥4.0.0), MEGAHIT ≥1.2.9, metaviralSPAdes

---

*Phase 4 SOP completed. KrakenUniq viral DB built and all three tools executed.*
*QC-4: PASS. No genuine microbial taxa detected (expected for test data).*
*Two critical lessons documented: trimming affects Kraken2 confidence, and*
*extract_nonhuman_kraken.py requires ancestor taxid filtering.*
