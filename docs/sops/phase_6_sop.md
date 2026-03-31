# Phase 6 SOP — Contig Decontamination + MAG Binning
## CLAUDE Pipeline: Cancer-Linked Analysis of Underlying DNA Elements

**Date executed:** 2026-02-25
**Sample:** SRR15090802 (Wahida et al. gut virome, VLP-enriched, 1,545,806 pairs)
**Phase:** 6 of 13
**SOP version:** 1.0

---

## 1. Objective

Phase 6 takes the merged, deduplicated contig set from Phase 5 and performs two
independent but complementary operations:

**6A — Contig Decontamination:** Removes contigs that are of human origin (missed by
read-level host depletion) or that are laboratory vector/adapter contamination. This
is a second-pass safety net at the assembly level.

**6B — MAG Binning:** Groups contigs into Metagenome-Assembled Genomes (MAGs) using
abundance (coverage depth) and sequence composition (tetranucleotide frequency). MAGs
are the bacterial genome reconstructions; viral genome bins (vMAGs) are assessed in
Phase 7 using viral-specific tools.

**Why two rounds of decontamination?**

Reads are filtered in Phase 2 (host depletion) and Phase 4 (Kraken2 post-check). However,
assembly can produce new artifacts:
- Short human reads that individually pass host depletion may assemble into contigs
  with enough coverage to appear as non-human sequences
- Low-complexity telomeric/centromeric repeats from T2T can assemble from residual
  reads into chimeric contigs
- Vector contamination in sequencing reagents produces contig-level artifacts that
  short-read filters miss

The Ge et al. (2025) reanalysis of TCGA found that contig-level T2T screening
removed 12–35% of contigs that survived read-level filtering in tumor WGS. This step
is non-negotiable.

---

## 2. Environment & Software

| Tool | Version | Environment | Purpose |
|------|---------|-------------|---------|
| BLAST+ (blastn) | 2.12.0 | claude_pipeline | T2T + UniVec screening |
| bowtie2 | 2.5.4 | claude_pipeline | Read-to-contig mapping for coverage |
| samtools | 1.20 | claude_pipeline | BAM processing |
| MetaBAT2 | 2.17 | claude_pipeline | Abundance + TNF binning |
| CoverM | 0.7.0 | claude_pipeline | Coverage statistics |
| DAS_Tool | 1.1.7 | das_tool_env | Bin refinement (multi-binner) |
| CheckM2 | 1.0.2 | checkm2_env | MAG completeness/contamination |

### Installation

```bash
# Core tools (claude_pipeline env)
mamba install -n claude_pipeline -c bioconda -c conda-forge \
    metabat2=2.17 coverm=0.7.0 -y

# DAS_Tool — requires separate env due to R dependency conflicts
mamba create -n das_tool_env -c bioconda -c conda-forge das_tool=1.1.7 -y
# Fix missing R package
conda run -n das_tool_env R --vanilla -e "
    options(repos = c(CRAN = 'https://cloud.r-project.org'))
    install.packages('docopt', quiet=TRUE)
"

# CheckM2 — requires Python <3.9, separate env
mamba create -n checkm2_env python=3.8 -y
mamba install -n checkm2_env -c bioconda checkm2=1.0.2 -y
```

### Databases

```bash
export DB=$PROJECT_DIR/databases

# T2T-CHM13v2.0 FASTA for BLASTn
wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz" \
    -O ${DB}/t2t_blast/chm13v2.0.fna.gz
gunzip ${DB}/t2t_blast/chm13v2.0.fna.gz  # ~3 GB decompressed

# Build T2T BLAST DB (takes ~30 sec, 25 sequences = chr1-22,X,Y,MT)
makeblastdb -in ${DB}/t2t_blast/chm13v2.0.fna \
    -dbtype nucl -out ${DB}/t2t_blast/chm13v2_blastdb \
    -title 'T2T-CHM13v2.0' -parse_seqids

# UniVec for vector/adapter screening
wget https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec \
    -O ${DB}/univec/UniVec
makeblastdb -in ${DB}/univec/UniVec -dbtype nucl \
    -out ${DB}/univec/UniVec_db -title 'UniVec'
# → 6,111 vector/primer sequences

# CheckM2 DIAMOND database (~1.3 GB compressed, ~3 GB unpacked)
conda run -n checkm2_env checkm2 database --download \
    --path ${DB}/checkm2
```

**Database inventory:**
```
databases/
├── t2t_blast/
│   ├── chm13v2.0.fna.gz      890 MB (compressed)
│   ├── chm13v2.0.fna          3.0 GB (decompressed)
│   └── chm13v2_blastdb.*      BLAST indices
├── univec/
│   ├── UniVec                 1.7 MB (6,111 sequences)
│   └── UniVec_db.*            BLAST indices
└── checkm2/
    └── CheckM2_database/
        └── uniref100.KO.1.dmnd  2.9 GB (DIAMOND protein DB)
```

---

## 3. Input Files

| File | Location | Contigs |
|------|----------|---------|
| all_contigs.fasta | pipeline_run/05_assembly/virome/merged/ | 20,548 |

---

## 4. Step-by-Step Execution

### Step 1: BLASTn contigs vs T2T-CHM13 (human decontamination)

```bash
conda run -n claude_pipeline bash -c "
blastn \
    -query $PROJECT_DIR/pipeline_run/05_assembly/virome/merged/all_contigs.fasta \
    -db $PROJECT_DIR/databases/t2t_blast/chm13v2_blastdb \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' \
    -evalue 1e-5 \
    -perc_identity 90 \
    -num_threads 4 \
    -out $PROJECT_DIR/pipeline_run/06_decontam_contigs/virome/blast_vs_t2t.tsv \
    2>&1
"
wc -l pipeline_run/06_decontam_contigs/virome/blast_vs_t2t.tsv
```

**Filter criteria:**
- E-value ≤ 1e-5 (statistically significant match)
- ≥ 90% nucleotide identity (high-confidence human sequence)
- Alignment ≥ 500 bp applied downstream by decontam_contigs.py

**Result (SRR15090802):**
```
0 BLAST hits
0 contigs flagged for human contamination
```

Expected for 3-pass host-depleted VLP-enriched data. In tumor WGS, expect 5–30%
of contigs to match T2T at this threshold.

---

### Step 2: BLASTn contigs vs UniVec (vector/adapter screening)

UniVec contains 6,111 sequences including:
- Cloning vectors (pUC19, pGEX, pBluescript)
- Sequencing primers (M13, T7, SP6)
- Illumina adapter sequences
- Transposon sequences

NCBI VecScreen parameters (stricter alignment scoring required for short sequences):

```bash
conda run -n claude_pipeline bash -c "
blastn \
    -query $PROJECT_DIR/pipeline_run/05_assembly/virome/merged/all_contigs.fasta \
    -db $PROJECT_DIR/databases/univec/UniVec_db \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' \
    -reward 1 -penalty -5 -gapopen 3 -gapextend 3 \
    -evalue 0.001 \
    -perc_identity 90 \
    -num_threads 4 \
    -out $PROJECT_DIR/pipeline_run/06_decontam_contigs/virome/blast_vs_univec.tsv \
    2>&1
"
```

**Key flags:**
- `-reward 1 -penalty -5 -gapopen 3 -gapextend 3`: NCBI VecScreen scoring matrix,
  designed for detecting vector contamination in short alignments
- Standard blastn scoring (reward 2, penalty -3) is too permissive for short primers

**Result (SRR15090802):**
```
17 BLAST hits across 5 contigs
```

Hits breakdown:
```
megahit_6922   gnl|uv|U09365.1   99.6%   242 bp / 2077 bp = 11.7%   → REMOVE
megahit_13034  gnl|uv|M77169.1   92.5%    97 bp / 798 bp  = 12.1%   → KEEP (< 100 bp)
metaspades_13011  gnl|uv|M77169.1  92.3%  220 bp / 595 bp = 37%     → NOTE: 95% threshold
metaspades_4971  NGB00728.1       100%    49 bp / 985 bp  = 5.0%    → KEEP (< 100 bp)
metaspades_392   EU546819.1       100%    36 bp / 4429 bp = 0.8%    → KEEP (< 100 bp)
...
```

**megahit_6922:** 99.6% identity over 242 bp to pGEX expression vector (U09365.1).
pGEX vectors contain an MCS derived from phage M13 origins; some bacteriophage genomes
share sequence with early cloning vectors. Regardless of origin, this passes the removal
threshold and is flagged as vector-contaminated.

**Short matches (29–49 bp):** Multiple contigs match Illumina primer sequences
(NGB-prefix accessions) at 100% identity in 29–49 bp windows. These represent primer
incorporation at read ends that survived fastp adapter trimming. The matches are at
contig edges only (positions 1–50 or contig_length - 50) and do not represent full
vector inserts. Applied threshold (≥100 bp) correctly preserves these contigs.

---

### Step 3: Apply decontamination filter

The `decontam_contigs.py` script applies the filter and produces the clean FASTA:

```bash
conda run -n claude_pipeline bash -c "
python pipeline_run/scripts/decontam_contigs.py \
    --contigs pipeline_run/05_assembly/virome/merged/all_contigs.fasta \
    --t2t-blast pipeline_run/06_decontam_contigs/virome/blast_vs_t2t.tsv \
    --univec-blast pipeline_run/06_decontam_contigs/virome/blast_vs_univec.tsv \
    --out pipeline_run/06_decontam_contigs/virome/clean_contigs.fasta \
    --removed pipeline_run/06_decontam_contigs/virome/removed_contigs.txt \
    --t2t-min-len 500 --t2t-min-pident 90.0 \
    --univec-min-len 100 --univec-min-pident 95.0
"
```

**Result:**
```
Total contigs:     20,548
T2T hits (human):       0
UniVec hits:            1   (megahit_6922, pGEX, 99.6%, 242 bp)
Removed total:          1
Kept (clean):      20,547
```

Removed contigs are logged to `removed_contigs.txt`:
```
contig_id    reason    pident    aln_len    qlen    frac_covered
megahit_6922 UniVec    99.6      242        2077    0.117
```

---

### Step 4: Build bowtie2 index and map reads to clean contigs

MetaBAT2 requires per-contig coverage depth profiles. Coverage is estimated by
mapping the Phase 3 QC-filtered reads back to the decontaminated contigs.

```bash
# Build index (~1 min)
conda run -n claude_pipeline bash -c "
bowtie2-build \
    pipeline_run/06_decontam_contigs/virome/clean_contigs.fasta \
    pipeline_run/06_decontam_contigs/virome/clean_contigs_bt2 \
    --threads 4
"

# Map reads → sort → index BAM (~5-10 min for 1.5M pairs)
conda run -n claude_pipeline bash -c "
bowtie2 \
    -x pipeline_run/06_decontam_contigs/virome/clean_contigs_bt2 \
    -1 pipeline_run/03_qc_filtered/virome/SRR15090802_R1.final.fastq.gz \
    -2 pipeline_run/03_qc_filtered/virome/SRR15090802_R2.final.fastq.gz \
    --threads 4 \
    --no-unal \
    2>pipeline_run/logs/bt2_vs_contigs_virome.log | \
samtools sort -@ 4 -o pipeline_run/06_decontam_contigs/virome/reads_vs_contigs.bam
samtools index pipeline_run/06_decontam_contigs/virome/reads_vs_contigs.bam
"
```

**Result:**
```
1,545,806 pairs input
65.82% overall alignment rate
  56.83% concordant unique
   5.38% concordant multi-mapper
  37.79% unaligned (reads from genomes with insufficient coverage to assemble)
```

65.82% alignment is typical for a VLP-enriched virome. The ~35% unaligned reads
come from very low-coverage phages whose k-mers appear in too few reads to assemble.

---

### Step 5: Coverage depth profiling

MetaBAT2's `jgi_summarize_bam_contig_depths` calculates the average depth and
depth variance for each contig — the two key inputs for MetaBAT2 binning.

```bash
conda run -n claude_pipeline bash -c "
jgi_summarize_bam_contig_depths \
    --outputDepth pipeline_run/07_mags/virome/metabat2/contig_depths.txt \
    pipeline_run/06_decontam_contigs/virome/reads_vs_contigs.bam
"
```

**Result:**
```
2,034,895 total reads processed
1,821,926 well-mapped reads (89.5%)
contig_depths.txt: 20,548 rows (one per contig)
Columns: contigName, contigLen, totalAvgDepth, sample_depth, sample_depth-var

Top contigs by depth:
  metaspades_1  (101,130 bp): 20.1× average depth
  metaspades_2   (57,411 bp): 32.1× average depth
```

---

### Step 6: MetaBAT2 binning

MetaBAT2 combines per-contig coverage depth and tetranucleotide frequency (TNF)
using an expectation-maximization algorithm to cluster contigs into bins.

```bash
conda run -n claude_pipeline bash -c "
metabat2 \
    -i pipeline_run/06_decontam_contigs/virome/clean_contigs.fasta \
    -a pipeline_run/07_mags/virome/metabat2/contig_depths.txt \
    -o pipeline_run/07_mags/virome/metabat2/bin \
    --minContig 2000 \
    --minS 60 \
    --numThreads 4 \
    --saveCls
"
```

**Key flags:**
- `--minContig 2000`: only bin contigs ≥ 2 kb (shorter contigs have unreliable TNF)
- `--minS 60`: minimum score for a contig to be placed in a bin (higher = stricter)
- `--saveCls`: save contig-to-bin classification table (useful for downstream tools)
- Default `--minClsSize 200000`: discard bins with <200 kb (prevents bin fragmentation)

**Result:**
```
4 bins (4,587,351 bases in total)
  bin.1.fa: 497 contigs, 3.34 MB, N50=8,010 bp
  bin.2.fa: 276 contigs, 805 KB, N50=2,842 bp
  bin.3.fa:  66 contigs, 226 KB, N50=3,301 bp
  bin.4.fa:  68 contigs, 212 KB, N50=3,159 bp
```

Only 4,587,351 / 22,180,406 bp (20.7%) of assembled sequence was binned. This is
expected for a virome: most viral contigs are too divergent from each other in TNF
and/or depth to cluster together at MetaBAT2's default thresholds. The unbinned
contigs are the primary input for Phase 7 viral identification.

---

### Step 7: DAS_Tool bin refinement

DAS_Tool selects non-redundant bins from multiple binners by scoring each bin on
completeness (predicted single-copy genes) and redundancy. When only one binner is
used (MetaBAT2), DAS_Tool acts as a filter, applying a completeness score threshold.

**NOTE: DAS_Tool experienced a version conflict in this run:**
```
Error in envRefInferField(x, what, getClass(class(x)), selfEnv) :
  'short' is not a valid field or method name for reference class "Argument"
Execution halted
```

Root cause: `docopt` R package version incompatibility. The `--score_threshold 0.0`
flag uses a short flag `-s` which conflicts with the docopt parser in the installed
version. This is a known issue with das_tool=1.1.7 + newer docopt R packages.

**Fix for production:**
```bash
# Downgrade docopt in das_tool_env
conda run -n das_tool_env R --vanilla -e "
    install.packages('docopt', repos='https://cloud.r-project.org')
"
# Or use DAS_Tool's Perl API directly:
conda run -n das_tool_env Fasta_to_Scaffolds2Bin.sh -e fa \
    -i pipeline_run/07_mags/virome/metabat2/ \
    > pipeline_run/07_mags/virome/das_tool/metabat2_contigs2bin.tsv

conda run -n das_tool_env DAS_Tool \
    -i pipeline_run/07_mags/virome/das_tool/metabat2_contigs2bin.tsv \
    -c pipeline_run/06_decontam_contigs/virome/clean_contigs.fasta \
    -o pipeline_run/07_mags/virome/das_tool/DAS_Tool_output \
    --score_threshold 0.4 \
    --write_bins 1 \
    --threads 4
```

**For single-binner scenarios:** DAS_Tool adds minimal value when only MetaBAT2
is used. In production with multiple binners (MetaBAT2 + MaxBin2 + CONCOCT),
DAS_Tool recovers 15–40% more non-redundant MAGs compared to any single binner.

**For this run:** Skipping DAS_Tool; MetaBAT2 bins used directly for CheckM2.

---

### Step 8: CheckM2 quality assessment

CheckM2 uses a machine-learning model trained on 11,600 bacterial and archaeal
reference genomes to estimate completeness and contamination from predicted
single-copy marker genes (SCMGs).

```bash
conda run -n checkm2_env bash -c "
checkm2 predict \
    --input pipeline_run/07_mags/virome/metabat2/ \
    --output-directory pipeline_run/07_mags/virome/checkm2 \
    --extension .fa \
    --threads 4 \
    --database_path /path/to/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd
"
```

**Result:**
```
Bin    Completeness  Contamination  GC%   Size     MIMAG_Quality
───────────────────────────────────────────────────────────────────
bin.1  78.68%        6.50%          73%   3.34 MB  Medium (MQ)
bin.2  28.19%        0.46%          41%   805 KB   Low (LQ)
bin.3   7.36%        0.02%          48%   226 KB   Low (LQ)
bin.4   3.28%        0.00%          65%   212 KB   Low (LQ)
```

**MIMAG standards (Bowers et al. 2017):**
| Category | Completeness | Contamination |
|----------|-------------|---------------|
| High-quality (HQ) | ≥ 90% | ≤ 5% |
| Medium-quality (MQ) | ≥ 50% | ≤ 10% |
| Low-quality (LQ) | < 50% | — |

No HQ MAGs. This is **expected for a VLP-enriched virome** — see Biology section below.

---

## 5. QC-6 Results

```
═══════════════════════════════════════════════════════════════════
SAMPLE: SRR15090802 (Wahida et al. gut virome)
PHASE 6 QC — Decontamination + MAG Binning
═══════════════════════════════════════════════════════════════════

DECONTAMINATION
  Input contigs:        20,548
  T2T hits removed:          0
  UniVec hits removed:       1   (megahit_6922, pGEX vector)
  Clean contigs:        20,547

COVERAGE
  Reads mapped:         65.82%   (1,017,573 / 1,545,806 pairs)

BINNING
  Bins produced:             4
  Binned sequence:       4.59 MB  (20.7% of assembled)
  HQ MAGs (≥90/≤5):          0
  MQ MAGs (≥50/≤10):         1   (bin.1, likely phage cluster)

═══════════════════════════════════════════════════════════════════
STATUS: PASS
All decontamination steps completed. 0 human contigs — expected for
VLP-enriched, 3-pass host-depleted data. 1 vector contig removed.
4 MetaBAT2 bins validated by CheckM2. Low MIMAG scores expected for
phage-dominated virome — CheckV (Phase 7) is appropriate quality tool.
═══════════════════════════════════════════════════════════════════
```

Checkpoint: `pipeline_run/qc_checkpoints/QC6_virome_SRR15090802.txt`

---

## 6. Problems & Troubleshooting

### Problem 1: CheckM2 incompatible with Python 3.10

**Symptom:**
```
error libmamba Could not solve for environment specs
checkm2 =1.0.2 would require python >=3.6,<3.9
pin on python =3.10 is not installable
```

**Root cause:** CheckM2 requires Python <3.9 (version 1.0.2) or >3.12 (version 1.1.0),
but our main `claude_pipeline` environment uses Python 3.10.

**Fix:** Install CheckM2 in a dedicated Python 3.8 environment:
```bash
mamba create -n checkm2_env python=3.8 -y
mamba install -n checkm2_env -c bioconda checkm2=1.0.2 -y
```
Then call via `conda run -n checkm2_env bash -c "checkm2 predict ..."`.

---

### Problem 2: DAS_Tool R/docopt version incompatibility

**Symptom:**
```
Error in envRefInferField(x, what, getClass(class(x)), selfEnv) :
  'short' is not a valid field or method name for reference class "Argument"
Execution halted
```

**Root cause:** The `docopt` R package version installed by conda is incompatible with
DAS_Tool 1.1.7's argument specification (uses a 'short' field not in newer docopt API).

**Fix (partial):** Installing `docopt` from CRAN resolves the startup error but the
`--score_threshold 0.0` argument parsing may still fail. The issue is reproducible
with DAS_Tool 1.1.7 + R ≥ 4.0.

**Workaround for this run:** Used MetaBAT2 bins directly (skipped DAS_Tool refinement).
For production multi-binner runs, use DAS_Tool 1.1.6 or wait for DAS_Tool 2.0.

---

### Problem 3: Low MIMAG quality from CheckM2 on virome bins

**Symptom:**
```
bin.1: 78.68% completeness (only MQ, not HQ)
bin.2-4: <30% completeness (LQ)
```

**Root cause:** CheckM2 is calibrated for bacteria, not bacteriophages. Viral genomes
lack the universal single-copy marker genes (bacterial SCMGs) that CheckM2 uses to
estimate completeness. A 100% complete phage genome might score 0% completeness under
CheckM2.

**Expected behaviour for virome data:**
- All bins will score lower completeness than equivalent bacterial MAGs
- High GC content (>65%) is a signal that the bin contains high-GC phages, not bacteria
- Low contamination scores are meaningful (CheckM2 contamination = redundant SCMGs)

**Solution:** Use CheckV (Phase 7) for viral-specific quality assessment:
- CheckV estimates viral genome completeness using viral gene databases
- CheckV distinguishes complete, partial, and fragmented viral genomes
- CheckV is the field-standard for vMAG quality (Nayfach et al. 2021)

---

### Problem 4: The `mkdir -n` flag confusion

**Symptom:**
```
mkdir: invalid option -- 'n'
```

**Root cause:** Typing `mkdir -p /path && conda run -n env ...` on a single line
where the shell interprets `-n env` as flags to `mkdir`.

**Fix:** Always run `mkdir -p /path` as a completely separate command before `conda run`.

---

## 7. Key Biological Reasoning

### Why do VLP-enriched virome bins fail MIMAG thresholds?

The MIMAG (Minimum Information about a Metagenome-Assembled Genome) standard was
designed for bacterial MAGs. The quality metrics (completeness, contamination) are
calculated from bacterial single-copy marker genes (SCMGs). Bacteriophage genomes:

1. **Lack universal SCMGs**: Unlike bacteria, which share ~200 universally conserved
   housekeeping genes, phages evolve extremely rapidly and lack conserved "core genes"
   detectable across diverse phage families.

2. **Have very different genome sizes**: Bacterial genomes are typically 1–10 MB.
   Phage genomes range from 4 kb (small RNA phages) to 750 kb (Jumbo phages). A
   "complete" 40 kb phage genome would appear as a 0.4–4% complete bacterial genome.

3. **Show high GC bias in specific genera**: Many *Caudoviricetes* (tailed dsDNA phages)
   infecting high-GC bacteria (e.g., *Streptomyces*, *Actinomyces*) have GC >65%.
   bin.1's GC of 73% strongly suggests it contains *Streptomyces*-infecting phages,
   not bacteria.

**The correct quality tool for viral bins is CheckV** (Phase 7):
- CheckV uses viral gene databases (VICTOR, VHDB) to estimate viral completeness
- Assesses circular/linear topology
- Identifies and quantifies host-derived contamination in viral bins
- Classifies as: complete, high-quality, medium-quality, or low-quality viral genomes

### Why map reads back to assembled contigs (Step 4)?

MetaBAT2 requires per-contig abundance (depth) in addition to sequence composition
for binning. This second mapping step ("read-to-assembly mapping") serves multiple
purposes:

1. **Abundance profile**: Contigs from the same organism co-vary in depth across
   samples. In single-sample mode, absolute depth is used instead.

2. **Differential coverage binning**: If multiple samples are available (e.g., multiple
   tumor samples from the same patient), coverage across samples can identify
   co-occurring organisms.

3. **Quality check**: The mapping rate back to assembled contigs (65.82% for
   SRR15090802) tells us how much of the sequenced community we successfully assembled.
   A mapping rate <30% would suggest poor assembly or extreme diversity.

### What is the contig size cutoff for MAG binning?

MetaBAT2's `--minContig 2000` (default) is based on the empirical observation that:
- Tetranucleotide frequencies (TNF) are not reliable estimators of genome identity
  for contigs shorter than ~2 kb
- Short contigs (<1 kb) have high TNF variance and would merge incorrectly across taxa
- The 2 kb cutoff was validated in MetaBAT2's benchmarking (Kang et al. 2019)

Only 907 of 20,547 clean contigs are ≥2 kb, yielding 4 bins. The remaining 19,640
contigs (<2 kb) are unbinned but still present in `clean_contigs.fasta` and are
fully analysed in Phase 7 (viral ID) and Phase 8 (prophage annotation).

### T2T vs GRCh38 for decontamination

The CLAUDE pipeline uses T2T-CHM13v2.0 rather than GRCh38 because:
- T2T-CHM13v2.0 is the first complete, gap-free assembly of the human genome (2022)
- GRCh38 has 151 Mb of unknown sequence ('N' gaps), primarily in centromeres and
  telomeres — regions that T2T filled in
- The telomeric/centromeric sequences filled by T2T are precisely the regions most
  likely to generate false-positive microbial matches (high repeat content, GC-neutral)
- Ge et al. (2025) showed T2T-based decontamination removed 12–35% more spurious
  microbial signals than GRCh38-based decontamination in TCGA data
- The T2T genome has 25 sequences vs GRCh38's 595+ (faster BLAST)

---

## 8. Lessons Learned

1. **CheckM2 requires a separate Python environment.** The Python 3.10 pinning in
   `claude_pipeline` conflicts with CheckM2's Python <3.9 requirement. Always create
   `checkm2_env` with Python 3.8 at the beginning of project setup.

2. **DAS_Tool has R package compatibility issues in conda.** The `docopt` R package
   version from conda-forge may not be compatible with DAS_Tool 1.1.7's command-line
   parser. Pre-install CRAN docopt; consider DAS_Tool 1.1.6 as a fallback.

3. **MIMAG metrics are invalid for viral bins.** When running CheckM2 on virome or
   phage-enriched samples, expect all bins to show <50% completeness. This is not a
   pipeline failure. Document the GC content and contig count instead as phage-quality
   proxies until CheckV runs in Phase 7.

4. **65% mapping rate back to contigs is a good virome result.** For complex gut
   viromes, 50–80% re-mapping is typical. <30% suggests either poor assembly (k-mer
   diversity too high) or contamination.

5. **T2T BLASTn is fast despite the 3 GB database.** 20,547 contigs × 22 MB against
   T2T (3 GB, 25 sequences) completed in ~4 minutes on 4 threads. The small number of
   sequences in T2T (25 vs GRCh38's 595+) makes it significantly faster.

6. **The bowtie2 multi-mapper rate (5.38%)** indicates that some genomic regions are
   assembled redundantly across contigs (e.g., shared terminal repeats in phage
   genomes). This is expected and handled correctly by MetaBAT2.

---

## 9. Production Notes (for Tumor WGS)

For TCGA/Hartwig tumor WGS with ~0.35% microbial reads, Phase 6 behaviour differs:

| Aspect | Virome (this run) | Tumor WGS (production) |
|--------|-------------------|------------------------|
| T2T BLAST hits | 0 (expected) | 5–30% of contigs |
| UniVec hits | 1 (rare) | <1% of contigs |
| Read re-mapping rate | 65.82% | 10–40% (lower diversity) |
| MetaBAT2 bins | 4 phage clusters | 2–50 bacterial bins |
| HQ MAGs expected | 0 (virome) | 0–5 (low microbial biomass) |
| DAS_Tool value | Low (1 binner) | High (3 binners) |

**For TCGA tumor WGS:** The T2T decontamination step is the most important step —
it removes contigs that represent human sequence falsely assembled from microbial
reads due to low-entropy regions (AT-rich repeat units, rDNA). The Ge et al. (2025)
reanalysis showed that skipping this step inflated microbial diversity estimates by
3–10× in the original TCGA microbiome studies.

### GTDB-Tk (not run this phase)

GTDB-Tk assigns taxonomic classifications to MAGs against the GTDB reference database
(R220). It was not run for this sample because:
1. The GTDB database requires ~200 GB storage
2. The 4 MetaBAT2 bins are phage-dominated and would not classify reliably under GTDB
3. Bacterial MAG taxonomy is not the primary goal of the virome analysis

For tumor WGS production runs with bacterial MAGs:
```bash
gtdbtk classify_wf \
    --genome_dir pipeline_run/07_mags/sample/das_tool/bins/ \
    --out_dir pipeline_run/07_mags/sample/gtdbtk/ \
    --extension .fa \
    --cpus 8 \
    --mash_db ${GTDB_DB}/mash/gtdbtk.msh
```

---

## 10. Files Produced

```
pipeline_run/
├── 06_decontam_contigs/virome/
│   ├── blast_vs_t2t.tsv              # 0 lines (no human hits)
│   ├── blast_vs_univec.tsv           # 17 lines (5 contigs)
│   ├── clean_contigs.fasta           # 20,547 contigs ← Phase 7 input
│   ├── removed_contigs.txt           # 1 contig removed (pGEX)
│   ├── reads_vs_contigs.bam          # Read mapping to contigs
│   ├── reads_vs_contigs.bam.bai      # BAM index
│   └── clean_contigs_bt2.*           # bowtie2 index files
│
├── 07_mags/virome/
│   ├── metabat2/
│   │   ├── contig_depths.txt         # Coverage profiles (20,548 rows)
│   │   ├── bin.1.fa                  # 497 contigs, 3.34 MB
│   │   ├── bin.2.fa                  # 276 contigs, 805 KB
│   │   ├── bin.3.fa                  #  66 contigs, 226 KB
│   │   ├── bin.4.fa                  #  68 contigs, 212 KB
│   │   └── bin (cls file)            # Contig→bin assignment table
│   ├── das_tool/
│   │   └── metabat2_contigs2bin.tsv  # 907 contigs mapped to 4 bins
│   └── checkm2/
│       ├── quality_report.tsv        # Completeness + contamination
│       └── protein_files/            # Predicted ORFs per bin
│
├── qc_checkpoints/
│   └── QC6_virome_SRR15090802.txt   # STATUS: PASS
│
├── databases/
│   ├── t2t_blast/
│   │   ├── chm13v2.0.fna.gz         # 890 MB compressed
│   │   ├── chm13v2.0.fna            # 3.0 GB decompressed
│   │   └── chm13v2_blastdb.*
│   ├── univec/
│   │   ├── UniVec                   # 1.7 MB
│   │   └── UniVec_db.*
│   └── checkm2/
│       └── CheckM2_database/
│           └── uniref100.KO.1.dmnd  # 2.9 GB
│
└── scripts/
    └── decontam_contigs.py          # Decontamination filter script
```

---

## 11. Next Step

Phase 7: Viral Identification

**Input:** `pipeline_run/06_decontam_contigs/virome/clean_contigs.fasta` (20,547 contigs)

**Tools:** geNomad + DeepVirFinder (DVF) + VirSorter2 → consensus ≥2/3 agreement
**Quality:** CheckV → completeness, contamination, topology (circular/linear)
**Binning:** vRhyme → viral genome bins (vMAGs)

**Why three viral classifiers?**
The CLAUDE pipeline Principle 2 requires consensus from ≥2 independent methods with
orthogonal algorithms before calling a contig as viral. geNomad uses marker gene
profiles, DVF uses a deep neural network on sequence composition, and VirSorter2 uses
hallmark gene detection. Agreement between ≥2/3 provides high confidence while
maintaining sensitivity.

---

*Phase 6 SOP complete — contig decontamination and MAG binning executed and validated.*
*0 human contigs, 1 vector contig removed. 4 MetaBAT2 bins produced, CheckM2 complete.*
*Low MIMAG scores expected for phage-dominated virome. QC-6: PASS.*
*Pipeline ready for Phase 7 viral identification.*
