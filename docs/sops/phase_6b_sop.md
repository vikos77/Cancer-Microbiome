# Phase 6B SOP: Bacterial MAG Annotation (Bakta)
## CLAUDE Pipeline: Cancer-Linked Analysis of Underlying DNA Elements

**Date created:** 2026-03-02
**Phase:** 6B of 13 (new phase; inserted between Phase 6 and Phase 7)
**SOP version:** 1.0

---

## 1. Objective

Phase 6B takes the bacterial MAGs produced by Phase 6 (MetaBAT2/SemiBin2 binning) and annotates them with Bakta. This converts raw nucleotide FASTA files into fully annotated genomes with:
- Gene calls (CDS, tRNA, rRNA, tmRNA, ncRNA, CRISPR arrays, sORFs, pseudogenes)
- Functional annotations (COG categories, GO terms, EC numbers)
- Virulence factor detection (VFDB expert system)
- AMR gene detection (AMRFinderPlus expert system)
- Origins of replication and transfer (oriC, oriV, oriT)
- Database cross-references (RefSeq, UniRef100, UniParc)

**Why this phase is essential:**

DefenseFinder and PADLOC (Phase 11) require annotated protein sequences as input; they cannot run on raw FASTA assemblies. The Acinetobacter pipeline consumed pre-annotated RefSeq genomes. MAGs from metagenomes have no annotations. Bakta bridges this gap.

Bakta also provides virulence factor and AMR gene detection that no other phase in the pipeline currently performs. For CRC metagenomes, this means automated identification of FadA (F. nucleatum), BFT/fragilysin (B. fragilis), colibactin biosynthesis genes (pks island in E. coli), and other cancer-relevant virulence determinants.

Additionally, Bakta's protein predictions (Pyrodigal-based + sORF detection) replace the independent Prodigal run in Phase 9 (mobilome), ensuring consistent gene calls across phases.

**Rationale for placement between Phase 6 and Phase 7:**

Phase 7 (viral identification) operates on the full clean contig set, not on MAG bins. Phase 6B operates specifically on MAG bins. These are independent operations that can run in parallel. However, Phase 6B's output is consumed by:
- Phase 9 (mobilome detection; uses Bakta FAA instead of Prodigal)
- Phase 11 (defence system profiling; requires Bakta GFF3/FAA)
- Phase 10 (integration; merges virulence/AMR/defence data)

---

## 2. Environment & Software

| Tool | Version | Environment | Purpose |
|------|---------|-------------|---------|
| Bakta | ≥1.10.0 | bakta_env | Bacterial genome/MAG annotation |
| AMRFinderPlus | (bundled) | bakta_env | AMR gene detection (Bakta dependency) |

### 2.1 Installation

```bash
# Create dedicated environment (Bakta has many dependencies)
mamba create -n bakta_env -c conda-forge -c bioconda bakta -y

# Activate and verify
conda activate bakta_env
bakta --version
# Expected: bakta 1.10.x or higher
```

### 2.2 Database Download

Bakta requires a mandatory annotation database. Two options:

**Full database (~30 GB compressed, ~84 GB unpacked, recommended):**
```bash
# Create database directory
mkdir -p $PROJECT_DIR/databases/bakta

# Download full database
bakta_db download --output $PROJECT_DIR/databases/bakta --type full

# This downloads the latest version with:
# - UniRef protein clusters (UniRef100, UniRef90, UniRef50)
# - RefSeq protein sequences
# - VFDB virulence factors
# - AMRFinderPlus database
# - Rfam, IS elements, oriC/oriV/oriT references
# - Pfam-A domain models
```

**Light database (~1.5 GB compressed, ~3.5 GB unpacked, fallback option):**
```bash
# Use if storage is severely constrained
bakta_db download --output $PROJECT_DIR/databases/bakta --type light

# Light DB omits UniRef50 clusters — results in more "hypothetical protein" annotations
# Acceptable for pipeline validation but suboptimal for publication-quality analysis
```

**AMRFinderPlus internal database setup:**
```bash
# Required after Bakta DB download — Bakta's own download procedure handles this
# But if needed manually:
amrfinder_update --force_update \
    --database $PROJECT_DIR/databases/bakta/db/amrfinderplus-db/
```

**Set environment variable (optional convenience):**
```bash
export BAKTA_DB=$PROJECT_DIR/databases/bakta/db
# Or pass --db flag every time
```

---

## 3. Input

| File | Source | Description |
|------|--------|-------------|
| `pipeline_run/07_mags/{sample}/metabat2/bin.*.fa` | Phase 6 | MetaBAT2 MAG bins (FASTA) |
| `pipeline_run/07_mags/{sample}/checkm2/quality_report.tsv` | Phase 6 | MAG quality metrics |

**Which bins to annotate:**

Annotate ALL bins from Phase 6 regardless of CheckM2 quality scores. CheckM2 MIMAG metrics (completeness, contamination) are calibrated for bacteria. For CRC stool metagenomes, even medium-quality MAGs (≥50% completeness, <10% contamination) can contain biologically significant contigs with virulence factors or defence systems. High-quality bins (≥90% completeness, <5% contamination) will produce the most reliable annotations, but we do not discard lower-quality bins at this stage; we flag them in Phase 10 integration.

**Exception:** Bins that Phase 6 QC identifies as purely viral (phage-dominated bins from VLP data) should still be annotated. Bakta will produce minimal bacterial gene calls, which is itself informative and confirms the bin's viral nature.

---

## 4. Procedure

### 4.1 Create Output Directory

```bash
mkdir -p pipeline_run/07b_bakta_annotation/{sample}/
```

### 4.2 Run Bakta on Each MAG Bin

```bash
# Set variables
SAMPLE="CRC_sample"  # or "ctrl_control_sample"
DB_PATH="$PROJECT_DIR/databases/bakta/db"
BIN_DIR="pipeline_run/07_mags/${SAMPLE}/metabat2"
OUT_BASE="pipeline_run/07b_bakta_annotation/${SAMPLE}"

# Loop over all bins
for BIN in ${BIN_DIR}/bin.*.fa; do
    BIN_NAME=$(basename ${BIN} .fa)
    OUT_DIR="${OUT_BASE}/${BIN_NAME}"
    mkdir -p ${OUT_DIR}

    conda run -n bakta_env --no-banner \
        bakta \
        --db ${DB_PATH} \
        --output ${OUT_DIR} \
        --prefix ${BIN_NAME} \
        --meta \
        --threads 4 \
        --min-contig-length 200 \
        --verbose \
        ${BIN} \
        2>&1 | tee ${OUT_DIR}/${BIN_NAME}_bakta.log
done
```

**Key parameters explained:**

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `--meta` | (flag) | Metagenome mode: uses MetaProdigal for CDS prediction, optimised for fragmented/mixed assemblies |
| `--threads 4` | 4 | Balance between speed and resource usage; adjust based on system |
| `--min-contig-length 200` | 200 bp | Include short contigs that may carry sORFs or partial genes |
| `--verbose` | (flag) | Detailed logging for debugging |

### 4.3 Expected Runtime

| MAG size | Approximate runtime | Notes |
|----------|-------------------|-------|
| 1–2 MB | 2–5 min | Small MAG, few contigs |
| 2–5 MB | 5–10 min | Typical bacterial MAG |
| 5–10 MB | 10–20 min | Large MAG or high contamination |

For a typical CRC stool metagenome producing 10–30 MAG bins, total annotation time is 1–4 hours on 4 threads.

---

## 5. Output Files

Bakta produces multiple output files per bin. The critical ones for downstream phases:

| File | Extension | Consumed by | Purpose |
|------|-----------|-------------|---------|
| `{bin}.gff3` | .gff3 | Phase 11 (DefenseFinder/PADLOC) | Gene coordinates + annotations |
| `{bin}.faa` | .faa | Phase 9 (mobileOG-db), Phase 11 | Protein sequences |
| `{bin}.ffn` | .ffn | Phase 11 (CRISPRCasFinder) | Nucleotide sequences of genes |
| `{bin}.gbff` | .gbff | Visualisation, manual inspection | GenBank format |
| `{bin}.tsv` | .tsv | Phase 10 (integration) | Tab-separated feature table |
| `{bin}.json` | .json | Programmatic access | Complete annotation data |
| `{bin}.txt` | .txt | QC | Summary statistics |
| `{bin}.embl` | .embl | Submission (ENA) | INSDC-compliant format |
| `{bin}.hypotheticals.faa` | .faa | Optional further analysis | Hypothetical proteins only |
| `{bin}.hypotheticals.tsv` | .tsv | Optional further analysis | Hypothetical protein details |
| `{bin}.png` | .png | QC, visualisation | Circular genome plot |

### 5.1 Consolidate Key Outputs

After all bins are annotated, consolidate for downstream consumption:

```bash
SAMPLE="CRC_sample"
OUT_BASE="pipeline_run/07b_bakta_annotation/${SAMPLE}"

# Concatenate all protein predictions for Phase 9 mobileOG-db
cat ${OUT_BASE}/*/bin.*.faa > ${OUT_BASE}/all_bins_proteins.faa

# Count total CDS, tRNA, etc. across all bins
echo "=== Bakta Annotation Summary ===" > ${OUT_BASE}/bakta_summary.txt
for BIN_DIR in ${OUT_BASE}/bin.*; do
    BIN_NAME=$(basename ${BIN_DIR})
    echo "--- ${BIN_NAME} ---" >> ${OUT_BASE}/bakta_summary.txt
    cat ${BIN_DIR}/${BIN_NAME}.txt >> ${OUT_BASE}/bakta_summary.txt
    echo "" >> ${OUT_BASE}/bakta_summary.txt
done

# Extract virulence factor hits
echo -e "bin\tfeature_type\tstart\tstop\tstrand\tlocus_tag\tgene\tproduct\tdbxrefs" \
    > ${OUT_BASE}/virulence_factors.tsv
for BIN_DIR in ${OUT_BASE}/bin.*; do
    BIN_NAME=$(basename ${BIN_DIR})
    # Parse GFF3 for VFDB annotations
    grep -i "VFDB" ${BIN_DIR}/${BIN_NAME}.gff3 2>/dev/null | \
        awk -v bin="${BIN_NAME}" 'BEGIN{OFS="\t"}{print bin,$3,$4,$5,$7,$9}' \
        >> ${OUT_BASE}/virulence_factors.tsv
done

# Extract AMR gene hits
echo -e "bin\tfeature_type\tstart\tstop\tstrand\tlocus_tag\tgene\tproduct\tdbxrefs" \
    > ${OUT_BASE}/amr_genes.tsv
for BIN_DIR in ${OUT_BASE}/bin.*; do
    BIN_NAME=$(basename ${BIN_DIR})
    grep -i "AMRFinder" ${BIN_DIR}/${BIN_NAME}.gff3 2>/dev/null | \
        awk -v bin="${BIN_NAME}" 'BEGIN{OFS="\t"}{print bin,$3,$4,$5,$7,$9}' \
        >> ${OUT_BASE}/amr_genes.tsv
done

# Extract oriT hits (HGT potential)
echo -e "bin\tstart\tstop\tstrand\tannotation" > ${OUT_BASE}/oriT_hits.tsv
for BIN_DIR in ${OUT_BASE}/bin.*; do
    BIN_NAME=$(basename ${BIN_DIR})
    grep "oriT\|oriC\|oriV" ${BIN_DIR}/${BIN_NAME}.gff3 2>/dev/null | \
        awk -v bin="${BIN_NAME}" 'BEGIN{OFS="\t"}{print bin,$4,$5,$7,$9}' \
        >> ${OUT_BASE}/oriT_hits.tsv
done
```

---

## 6. QC Checkpoint: QC-6B

```bash
SAMPLE="CRC_sample"
OUT_BASE="pipeline_run/07b_bakta_annotation/${SAMPLE}"
QC_FILE="pipeline_run/qc_checkpoints/QC6B_${SAMPLE}.txt"

echo "=== QC-6B: Bakta MAG Annotation ===" > ${QC_FILE}
echo "Date: $(date)" >> ${QC_FILE}
echo "Sample: ${SAMPLE}" >> ${QC_FILE}
echo "" >> ${QC_FILE}

# Count bins annotated
N_BINS=$(ls -d ${OUT_BASE}/bin.* 2>/dev/null | wc -l)
echo "Bins annotated: ${N_BINS}" >> ${QC_FILE}

# Per-bin statistics
for BIN_DIR in ${OUT_BASE}/bin.*; do
    BIN_NAME=$(basename ${BIN_DIR})
    GFF="${BIN_DIR}/${BIN_NAME}.gff3"
    FAA="${BIN_DIR}/${BIN_NAME}.faa"

    if [ -f "${GFF}" ] && [ -f "${FAA}" ]; then
        N_CDS=$(grep -c "^[^#].*\tCDS\t" ${GFF} 2>/dev/null || echo 0)
        N_TRNA=$(grep -c "^[^#].*\ttRNA\t" ${GFF} 2>/dev/null || echo 0)
        N_RRNA=$(grep -c "^[^#].*\trRNA\t" ${GFF} 2>/dev/null || echo 0)
        N_CRISPR=$(grep -c "^[^#].*\tCRISPR\t" ${GFF} 2>/dev/null || echo 0)
        N_PROTEINS=$(grep -c "^>" ${FAA} 2>/dev/null || echo 0)
        N_HYPO=$(grep -c "hypothetical protein" ${FAA} 2>/dev/null || echo 0)

        if [ ${N_PROTEINS} -gt 0 ]; then
            PCT_HYPO=$(echo "scale=1; ${N_HYPO} * 100 / ${N_PROTEINS}" | bc)
        else
            PCT_HYPO="N/A"
        fi

        echo "${BIN_NAME}: ${N_CDS} CDS, ${N_TRNA} tRNA, ${N_RRNA} rRNA, ${N_CRISPR} CRISPR, ${PCT_HYPO}% hypothetical" >> ${QC_FILE}
    else
        echo "${BIN_NAME}: FAILED — output files missing" >> ${QC_FILE}
    fi
done

echo "" >> ${QC_FILE}

# Virulence factor count
N_VF=$(tail -n +2 ${OUT_BASE}/virulence_factors.tsv 2>/dev/null | wc -l)
echo "Total virulence factor hits: ${N_VF}" >> ${QC_FILE}

# AMR gene count
N_AMR=$(tail -n +2 ${OUT_BASE}/amr_genes.tsv 2>/dev/null | wc -l)
echo "Total AMR gene hits: ${N_AMR}" >> ${QC_FILE}

# oriT count
N_ORI=$(tail -n +2 ${OUT_BASE}/oriT_hits.tsv 2>/dev/null | wc -l)
echo "Total oriC/oriV/oriT hits: ${N_ORI}" >> ${QC_FILE}

echo "" >> ${QC_FILE}

# QC verdict
if [ ${N_BINS} -gt 0 ]; then
    echo "QC-6B VERDICT: PASS" >> ${QC_FILE}
else
    echo "QC-6B VERDICT: FAIL — no bins annotated" >> ${QC_FILE}
fi

cat ${QC_FILE}
```

**QC-6B thresholds:**

| Metric | PASS | WARN | FAIL |
|--------|------|------|------|
| Bins annotated | ≥1 | n/a | 0 |
| CDS per bin (stool metagenome) | >100 | 10–100 | <10 |
| % hypothetical proteins | <60% | 60–80% | >80% |
| GFF3 + FAA files present | All bins | n/a | Any missing |

**Interpreting high hypothetical rates:**

For well-studied organisms (E. coli, B. fragilis), expect <40% hypothetical. For novel/uncultured MAGs, 50–70% is normal. >80% suggests the MAG may be from an organism poorly represented in databases; flag but do not discard.

---

## 7. Gotchas

1. **`--meta` flag is critical for MAGs.** Without it, Bakta uses standard Prodigal mode which assumes a single complete genome. MetaProdigal (enabled by `--meta`) handles mixed/fragmented assemblies correctly.

2. **Bakta requires its own conda environment.** It depends on tRNAscan-SE, Aragorn, Infernal, Pyrodigal, DIAMOND, BLAST+, hmmer, and AMRFinderPlus; version conflicts with `claude_pipeline` are likely.

3. **AMRFinderPlus database must be initialised.** If Bakta fails with AMRFinder errors, run `amrfinder_update` manually (see Section 2.2).

4. **Database path must be passed explicitly** unless `BAKTA_DB` environment variable is set. `conda run` does not preserve environment variables from the outer shell; either set the variable inside the `conda run` command or use `--db` flag.

5. **Large MAGs (>10 MB) may indicate high contamination.** A 10 MB "bacterial" MAG likely contains contigs from multiple organisms. Bakta will annotate everything, but defence system analysis in Phase 11 should interpret results from such bins cautiously.

6. **CRISPR detection overlap with Phase 11.** Bakta detects CRISPR arrays using its own method. CRISPRCasFinder in Phase 11 provides more detailed subtyping (I-A, I-B, II-A, etc.). Use Bakta's CRISPR count as QC and CRISPRCasFinder's output for analysis.

7. **sORF detection.** Bakta finds small proteins (<30 aa) that Prodigal misses. These can include small antimicrobial peptides and regulatory proteins. In the cancer context, some bacteriocins and small secreted effectors fall into this category.

---

## 8. Connection to Downstream Phases

```
Phase 6B outputs → consumed by:

{bin}.faa ──────→ Phase 9 (mobileOG-db DIAMOND search, replaces Prodigal)
{bin}.gff3 + .faa → Phase 11 (DefenseFinder + PADLOC)
{bin}.gff3 ─────→ Phase 11 (CRISPRCasFinder)
virulence_factors.tsv → Phase 10 (integration: virulence × defence × mobilome)
amr_genes.tsv ──→ Phase 10 (integration: AMR × defence correlation)
oriT_hits.tsv ──→ Phase 10 (integration: HGT potential × defence trade-off)
bakta_summary.txt → Phase 10 (QC dashboard)
```
