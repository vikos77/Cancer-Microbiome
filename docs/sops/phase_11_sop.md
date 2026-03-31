# Phase 11 SOP: Bacterial Defence System Profiling
## CLAUDE Pipeline: Cancer-Linked Analysis of Underlying DNA Elements

**Date created:** 2026-03-02
**Phase:** 11 of 13 (new phase; inserted between Phase 10 integration draft and final integration)
**SOP version:** 1.0

---

## 1. Objective

Phase 11 profiles bacterial defence systems in annotated MAGs from Phase 6B. This is the direct application of the Acinetobacter defence system analytical framework to cancer-associated bacteria.

Three tools are used, each detecting different but overlapping sets of defence systems:

| Tool | Scope | Systems detected |
|------|-------|-----------------|
| DefenseFinder | 150+ defence system types | CRISPR-Cas, R-M, BREX, DISARM, Gabija, Druantia, Thoeris, Zorya, Lamassu, retrons, abortive infection, etc. |
| PADLOC | 150+ defence system types | Overlapping but independently curated; different HMM models catch systems DefenseFinder may miss and vice versa |
| CRISPRCasFinder | CRISPR-Cas only | Detailed CRISPR array structure: spacer sequences, repeat types, Cas subtyping (I-A, I-B, II-A, etc.) |

**Why two general tools (DefenseFinder + PADLOC)?**

This is the same dual-tool approach validated in the Acinetobacter pipeline. Neither tool detects all defence systems. DefenseFinder uses MacSyFinder-based HMM models; PADLOC uses its own curated HMM database. In the Acinetobacter analysis, the intersection of both tools provided the most reliable catalogue, while tool-specific detections flagged candidates requiring manual verification. Fisher's exact test co-occurrence matrices from both tools were compared side-by-side.

**What biological questions this phase answers:**

1. What is the defence system repertoire of each cancer-associated bacterial MAG?
2. Do CRC-enriched bacteria (F. nucleatum, pks+ E. coli, B. fragilis) carry CRISPR-Cas arrays, and if so, do their spacers match phages identified in Phase 7?
3. Is there a correlation between defence system burden and mobile genetic element load (Phase 9)?
4. Do bacteria with fewer defence systems carry more virulence factors (Phase 6B Bakta VFDB hits)?
5. Is the defence-vs-HGT trade-off observed in Acinetobacter replicated in cancer-associated gut bacteria?

---

## 2. Environment & Software

### 2.1 DefenseFinder

```bash
# DefenseFinder requires its own environment (MacSyFinder dependency conflicts)
mamba create -n defensefinder_env -c conda-forge -c bioconda defense-finder -y

# Verify
conda activate defensefinder_env
defense-finder --version
# Expected: defense-finder 2.x

# Update models (required on first run)
defense-finder update
```

### 2.2 PADLOC

```bash
# PADLOC can share the defensefinder_env or use its own
mamba create -n padloc_env -c conda-forge -c bioconda padloc -y

# Verify
conda activate padloc_env
padloc --version

# Download PADLOC database
padloc --db-update
```

### 2.3 CRISPRCasFinder

```bash
# CRISPRCasFinder is best installed via Docker/Singularity
# or use the standalone Perl installation
# If using conda:
mamba create -n crispr_env -c conda-forge -c bioconda ccfinder -y

# Alternative: use Bakta's CRISPR calls as baseline and
# CRISPRCasFinder for detailed subtyping when needed
```

---

## 3. Input

| File | Source | Description |
|------|--------|-------------|
| `pipeline_run/07b_bakta_annotation/{sample}/bin.*/bin.*.faa` | Phase 6B | Protein sequences per MAG |
| `pipeline_run/07b_bakta_annotation/{sample}/bin.*/bin.*.gff3` | Phase 6B | Gene annotations per MAG |
| `pipeline_run/07_mags/{sample}/metabat2/bin.*.fa` | Phase 6 | Raw MAG FASTA (for CRISPRCasFinder) |

---

## 4. Procedure

### 4.1 Create Output Directory

```bash
SAMPLE="CRC_sample"
mkdir -p pipeline_run/14_defence_systems/${SAMPLE}/{defensefinder,padloc,crisprcasfinder}
```

### 4.2 Run DefenseFinder on Each MAG

DefenseFinder takes protein FASTA (`.faa`) as input:

```bash
SAMPLE="CRC_sample"
BAKTA_DIR="pipeline_run/07b_bakta_annotation/${SAMPLE}"
OUT_DIR="pipeline_run/14_defence_systems/${SAMPLE}/defensefinder"

for BIN_DIR in ${BAKTA_DIR}/bin.*; do
    BIN_NAME=$(basename ${BIN_DIR})
    FAA="${BIN_DIR}/${BIN_NAME}.faa"

    if [ -f "${FAA}" ] && [ -s "${FAA}" ]; then
        mkdir -p ${OUT_DIR}/${BIN_NAME}

        conda run -n defensefinder_env --no-banner \
            defense-finder run \
            --models-dir $(defense-finder update --models-dir 2>/dev/null || echo "$HOME/.macsyfinder/models") \
            -o ${OUT_DIR}/${BIN_NAME} \
            -w 4 \
            ${FAA} \
            2>&1 | tee ${OUT_DIR}/${BIN_NAME}/defensefinder.log
    fi
done
```

**DefenseFinder output files:**

| File | Description |
|------|-------------|
| `defense_finder_systems.tsv` | System-level summary (type, subtype, genes) |
| `defense_finder_genes.tsv` | Gene-level detail (each protein in each system) |
| `defense_finder_hmmer.tsv` | Raw HMM hit details |

### 4.3 Run PADLOC on Each MAG

PADLOC takes protein FASTA (`.faa`) and GFF3 (`.gff3`) as input:

```bash
SAMPLE="CRC_sample"
BAKTA_DIR="pipeline_run/07b_bakta_annotation/${SAMPLE}"
OUT_DIR="pipeline_run/14_defence_systems/${SAMPLE}/padloc"

for BIN_DIR in ${BAKTA_DIR}/bin.*; do
    BIN_NAME=$(basename ${BIN_DIR})
    FAA="${BIN_DIR}/${BIN_NAME}.faa"
    GFF="${BIN_DIR}/${BIN_NAME}.gff3"

    if [ -f "${FAA}" ] && [ -s "${FAA}" ] && [ -f "${GFF}" ]; then
        mkdir -p ${OUT_DIR}/${BIN_NAME}

        conda run -n padloc_env --no-banner \
            padloc \
            --faa ${FAA} \
            --gff ${GFF} \
            --outdir ${OUT_DIR}/${BIN_NAME} \
            --cpu 4 \
            2>&1 | tee ${OUT_DIR}/${BIN_NAME}/padloc.log
    fi
done
```

**PADLOC output files:**

| File | Description |
|------|-------------|
| `*_padloc.csv` | System-level summary (system, type, subtype) |
| `*_padloc_genes.csv` | Gene-level detail |

### 4.4 Run CRISPRCasFinder on Each MAG (Optional: Detailed Subtyping)

CRISPRCasFinder takes nucleotide FASTA as input:

```bash
SAMPLE="CRC_sample"
MAG_DIR="pipeline_run/07_mags/${SAMPLE}/metabat2"
OUT_DIR="pipeline_run/14_defence_systems/${SAMPLE}/crisprcasfinder"

for BIN in ${MAG_DIR}/bin.*.fa; do
    BIN_NAME=$(basename ${BIN} .fa)
    mkdir -p ${OUT_DIR}/${BIN_NAME}

    # If using Docker:
    # docker run -v $(pwd):/data -w /data \
    #     crisprfinder \
    #     -i ${BIN} \
    #     -o ${OUT_DIR}/${BIN_NAME}

    # If using conda ccfinder:
    conda run -n crispr_env --no-banner \
        CRISPRCasFinder.pl \
        -i ${BIN} \
        -o ${OUT_DIR}/${BIN_NAME} \
        -cas \
        -meta \
        2>&1 | tee ${OUT_DIR}/${BIN_NAME}/crisprcasfinder.log
done
```

CRISPRCasFinder provides what DefenseFinder/PADLOC do not:
- **Spacer sequences**: these can be BLASTed against Phase 7 viral contigs to find direct phage-host CRISPR targeting evidence
- **Repeat consensus sequences**: for CRISPR array classification
- **Cas protein subtyping**: detailed I-A, I-B, II-A, III-A classification

### 4.5 Consolidate Results

```bash
SAMPLE="CRC_sample"
BASE="pipeline_run/14_defence_systems/${SAMPLE}"

# --- DefenseFinder consolidation ---
echo -e "bin\tsystem_type\tsystem_subtype\tgenes\tprotein_ids" \
    > ${BASE}/defensefinder_consolidated.tsv

for BIN_DIR in ${BASE}/defensefinder/bin.*; do
    BIN_NAME=$(basename ${BIN_DIR})
    SYSTEMS_FILE="${BIN_DIR}/defense_finder_systems.tsv"
    if [ -f "${SYSTEMS_FILE}" ]; then
        tail -n +2 ${SYSTEMS_FILE} | \
            awk -v bin="${BIN_NAME}" 'BEGIN{OFS="\t"}{print bin,$0}' \
            >> ${BASE}/defensefinder_consolidated.tsv
    fi
done

# --- PADLOC consolidation ---
echo -e "bin\tsystem\tsubtype\tgenes" \
    > ${BASE}/padloc_consolidated.tsv

for BIN_DIR in ${BASE}/padloc/bin.*; do
    BIN_NAME=$(basename ${BIN_DIR})
    PADLOC_FILE=$(ls ${BIN_DIR}/*_padloc.csv 2>/dev/null | head -1)
    if [ -f "${PADLOC_FILE}" ]; then
        tail -n +2 ${PADLOC_FILE} | \
            awk -v bin="${BIN_NAME}" -F',' 'BEGIN{OFS="\t"}{print bin,$0}' \
            >> ${BASE}/padloc_consolidated.tsv
    fi
done

# --- Summary counts ---
echo "=== Defence System Summary ===" > ${BASE}/defence_summary.txt
echo "Sample: ${SAMPLE}" >> ${BASE}/defence_summary.txt
echo "" >> ${BASE}/defence_summary.txt

echo "--- DefenseFinder ---" >> ${BASE}/defence_summary.txt
echo "Total systems: $(tail -n +2 ${BASE}/defensefinder_consolidated.tsv | wc -l)" >> ${BASE}/defence_summary.txt
echo "Unique system types:" >> ${BASE}/defence_summary.txt
tail -n +2 ${BASE}/defensefinder_consolidated.tsv | cut -f2 | sort | uniq -c | sort -rn >> ${BASE}/defence_summary.txt
echo "" >> ${BASE}/defence_summary.txt

echo "--- PADLOC ---" >> ${BASE}/defence_summary.txt
echo "Total systems: $(tail -n +2 ${BASE}/padloc_consolidated.tsv | wc -l)" >> ${BASE}/defence_summary.txt
echo "Unique system types:" >> ${BASE}/defence_summary.txt
tail -n +2 ${BASE}/padloc_consolidated.tsv | cut -f2 | sort | uniq -c | sort -rn >> ${BASE}/defence_summary.txt

cat ${BASE}/defence_summary.txt
```

---

## 5. Cross-Phase Analysis: CRISPR Spacers vs Viral Contigs

This is the critical analysis that connects Phase 7 (viral identification) to Phase 11 (defence systems). If a bacterial MAG carries CRISPR spacers that match a viral contig from the same sample, this is direct evidence of a phage-host interaction.

```bash
SAMPLE="CRC_sample"

# Extract spacer sequences from CRISPRCasFinder output
# (adjust path based on CRISPRCasFinder output format)
SPACER_DIR="pipeline_run/14_defence_systems/${SAMPLE}/crisprcasfinder"
VIRAL_FASTA="pipeline_run/08_viral_contigs/${SAMPLE}/consensus_viral_contigs_2tool.fasta"

# Combine all spacers into one FASTA
cat ${SPACER_DIR}/*/TSV/Crisprs_REPORT.tsv 2>/dev/null | \
    grep -v "^#" | \
    awk -F'\t' '{print ">"$1"_spacer_"NR"\n"$NF}' \
    > ${SPACER_DIR}/all_spacers.fasta

# BLAST spacers against viral contigs
if [ -s "${SPACER_DIR}/all_spacers.fasta" ] && [ -f "${VIRAL_FASTA}" ]; then
    makeblastdb -in ${VIRAL_FASTA} -dbtype nucl -out ${SPACER_DIR}/viral_db

    blastn \
        -query ${SPACER_DIR}/all_spacers.fasta \
        -db ${SPACER_DIR}/viral_db \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
        -evalue 1e-5 \
        -word_size 7 \
        -num_threads 4 \
        -out ${SPACER_DIR}/spacer_vs_viral_hits.tsv

    echo "CRISPR spacer-viral contig matches:" >> pipeline_run/14_defence_systems/${SAMPLE}/defence_summary.txt
    wc -l < ${SPACER_DIR}/spacer_vs_viral_hits.tsv >> pipeline_run/14_defence_systems/${SAMPLE}/defence_summary.txt
fi
```

A spacer match to a viral contig means: "This bacterium has previously encountered this phage (or a close relative) and incorporated a spacer targeting it." This is direct molecular evidence of phage predation pressure on cancer-associated bacteria.

---

## 6. QC Checkpoint: QC-14

```bash
SAMPLE="CRC_sample"
BASE="pipeline_run/14_defence_systems/${SAMPLE}"
QC_FILE="pipeline_run/qc_checkpoints/QC14_${SAMPLE}.txt"

echo "=== QC-14: Defence System Profiling ===" > ${QC_FILE}
echo "Date: $(date)" >> ${QC_FILE}
echo "Sample: ${SAMPLE}" >> ${QC_FILE}
echo "" >> ${QC_FILE}

# DefenseFinder results
N_DF_SYSTEMS=$(tail -n +2 ${BASE}/defensefinder_consolidated.tsv 2>/dev/null | wc -l)
N_DF_TYPES=$(tail -n +2 ${BASE}/defensefinder_consolidated.tsv 2>/dev/null | cut -f2 | sort -u | wc -l)
echo "DefenseFinder: ${N_DF_SYSTEMS} systems, ${N_DF_TYPES} unique types" >> ${QC_FILE}

# PADLOC results
N_PL_SYSTEMS=$(tail -n +2 ${BASE}/padloc_consolidated.tsv 2>/dev/null | wc -l)
N_PL_TYPES=$(tail -n +2 ${BASE}/padloc_consolidated.tsv 2>/dev/null | cut -f2 | sort -u | wc -l)
echo "PADLOC: ${N_PL_SYSTEMS} systems, ${N_PL_TYPES} unique types" >> ${QC_FILE}

# CRISPR spacer matches
N_SPACER_HITS=$(wc -l < ${BASE}/crisprcasfinder/spacer_vs_viral_hits.tsv 2>/dev/null || echo 0)
echo "CRISPR spacer-viral matches: ${N_SPACER_HITS}" >> ${QC_FILE}

echo "" >> ${QC_FILE}

# Verdict
if [ ${N_DF_SYSTEMS} -gt 0 ] || [ ${N_PL_SYSTEMS} -gt 0 ]; then
    echo "QC-14 VERDICT: PASS" >> ${QC_FILE}
else
    echo "QC-14 VERDICT: WARN — 0 defence systems detected (possible for viral bins or sparse MAGs)" >> ${QC_FILE}
fi

cat ${QC_FILE}
```

**QC-14 thresholds:**

| Metric | Expected (CRC stool MAGs) | Notes |
|--------|--------------------------|-------|
| Defence systems per bacterial MAG | 3–15 | Typical for gut bacteria |
| Common types | R-M, CRISPR-Cas, Abi | Restriction-modification is nearly universal |
| CRISPR arrays per MAG | 0–5 | Not all bacteria carry CRISPR |
| Spacer-viral matches | 0–few | Depends on phage coverage in Phase 7 |

**0 defence systems in a MAG** can mean: (a) sparse/fragmented MAG; defence loci are typically 5–20 kb operons that require contiguous assembly; (b) the organism genuinely has few defences (rare but possible for HGT-permissive organisms); (c) novel/divergent defence systems not in current databases.

---

## 7. Expected Results for CRC Samples

For a CRC stool metagenome from the Wirbel cohort (CRC_sample), the expected defence landscape:

**F. nucleatum MAG (if assembled):**
- Likely sparse CRISPR-Cas (Fusobacterium genomes typically carry 0–2 CRISPR arrays)
- R-M systems (Type I or Type II, common in Fusobacteriaceae)
- Possible abortive infection systems
- Low overall defence burden; consistent with F. nucleatum's role as a hub organism that acquires diverse MGEs

**E. coli MAG (if assembled):**
- CRISPR-Cas Type I-E or I-F (common in E. coli)
- Multiple R-M systems (E. coli carries 2–6 R-M systems typically)
- Possibly BREX, Gabija, or other recently discovered systems
- Higher defence burden than F. nucleatum; but pks+ strains specifically may show reduced defences (enabling pks island acquisition)

**B. fragilis MAG (if assembled):**
- Typically high R-M system diversity
- CRISPR-Cas presence variable
- Shufflon-based phase variation of defence systems

**CRC vs Control comparison:**
- The key question is not whether CRC bacteria have more/fewer defences overall, but whether defence architecture differs for CRC-enriched species versus their counterparts in healthy controls
- Hypothesis: CRC-enriched strains may have reduced defence systems enabling virulence factor acquisition via HGT (paralleling the Acinetobacter finding)

---

## 8. Gotchas

1. **DefenseFinder models must be updated** before first run. `defense-finder update` downloads the latest MacSyFinder models. Without this, DefenseFinder may use outdated or empty model sets.

2. **PADLOC GFF3 format sensitivity.** PADLOC requires GFF3 files with specific formatting. Bakta's GFF3 output is PADLOC-compatible. If you encounter parsing errors, verify the GFF3 has standard column structure (seqid, source, type, start, end, score, strand, phase, attributes).

3. **DefenseFinder and PADLOC count the same system differently.** DefenseFinder may call "RM_Type_II" while PADLOC calls "restriction_modification_type_II". Reconciliation requires manual mapping; the Acinetobacter pipeline's R script (`defense_distribution.R`) already handles this.

4. **CRISPRCasFinder is the most computationally demanding** of the three tools. For large MAGs (>5 MB), expect 10–30 min per bin. If runtime is prohibitive, use Bakta's CRISPR calls as the primary source and reserve CRISPRCasFinder for bins where detailed spacer analysis is needed.

5. **Spacer BLAST parameters matter.** CRISPR spacers are short (25–45 bp). Use `word_size 7` and relaxed e-value (`1e-5`) to capture divergent matches. Mismatches in spacers can indicate historical rather than current targeting.

6. **Defence system operons span multiple contigs in fragmented MAGs.** If a CRISPR-Cas system requires both cas genes and the CRISPR array, and these are on different contigs in a MAG, DefenseFinder may miss it. This is an inherent limitation of short-read metagenome assembly.

---

## 9. Connection to Phase 10 Integration

Phase 11 outputs feed directly into the revised Phase 10 integration script:

```
Phase 10 integration now merges:
├── Phase 7: viral contigs (master_viral_table.tsv)
├── Phase 8: phage annotation + host prediction
├── Phase 9: mobilome (IS elements, integrons, plasmids, mobileOG-db)
├── Phase 6B: Bakta annotation (virulence factors, AMR genes, oriT)
├── Phase 11: defence systems (DefenseFinder + PADLOC)
│
└── NEW correlation analyses:
    ├── Defence system count vs MGE count per MAG (Fisher's exact test)
    ├── Defence system type vs virulence factor presence (chi-square)
    ├── CRISPR spacers vs viral contigs (direct phage-host evidence)
    ├── Defence burden: CRC-enriched vs control-enriched MAGs
    └── Co-occurrence matrices (DefenseFinder × PADLOC, as in Acinetobacter)
```

The integration script (`integrate_pipeline.py`) will need to be updated to consume Phase 6B and Phase 11 outputs. This is Phase 10 SOP territory.
