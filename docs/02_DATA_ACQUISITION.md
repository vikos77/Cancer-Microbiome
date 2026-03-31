# CLAUDE Pipeline — Part 02: Data Acquisition

> **Navigation:** [01_SETUP](01_SETUP.md) → 02_DATA_ACQUISITION → [03_HOST_DEPLETION](03_HOST_DEPLETION.md)

---

## PHASE 1 — Data Acquisition & Pre-processing

### Rationale

TCGA and Hartwig store aligned data as BAM/CRAM. We extract unmapped read pairs
(SAM flag 12: both mates unmapped) plus reads where only one mate mapped (chimeric
pairs — the unmapped mate may be microbial). Ge et al. (2025) started from TCGA
unmapped reads; Battaglia et al. (2024) used SAM flag 12 on Hartwig CRAMs.

**Critical detail:** Some TCGA BAMs were aligned to GRCh37, others GRCh38. The
unmapped fraction varies by alignment reference and aligner settings. Always record
which reference was used upstream.

### Commands

```bash
SAMPLE="TCGA-XX-XXXX"
INPUT_BAM="${PROJECT_DIR}/00_raw/${SAMPLE}.bam"

# Extract unmapped pairs (both mates unmapped)
samtools view -b -f 12 -F 256 ${INPUT_BAM} \
    -@ 16 -o ${PROJECT_DIR}/01_fastq/${SAMPLE}.unmapped.bam

# Also extract reads where one mate is mapped, one is not
# (chimeric reads — microbial mate may be here)
samtools view -b -f 4 -F 264 ${INPUT_BAM} -@ 16 \
    -o ${PROJECT_DIR}/01_fastq/${SAMPLE}.mate_unmapped_R1.bam
samtools view -b -f 8 -F 260 ${INPUT_BAM} -@ 16 \
    -o ${PROJECT_DIR}/01_fastq/${SAMPLE}.mate_unmapped_R2.bam

# Merge all unmapped reads
samtools merge -@ 16 ${PROJECT_DIR}/01_fastq/${SAMPLE}.all_unmapped.bam \
    ${PROJECT_DIR}/01_fastq/${SAMPLE}.unmapped.bam \
    ${PROJECT_DIR}/01_fastq/${SAMPLE}.mate_unmapped_R1.bam \
    ${PROJECT_DIR}/01_fastq/${SAMPLE}.mate_unmapped_R2.bam

# Sort by name (required for FASTQ conversion)
samtools sort -n -@ 16 \
    ${PROJECT_DIR}/01_fastq/${SAMPLE}.all_unmapped.bam \
    -o ${PROJECT_DIR}/01_fastq/${SAMPLE}.sorted.bam

# Convert to FASTQ
samtools fastq -@ 16 \
    -1 ${PROJECT_DIR}/01_fastq/${SAMPLE}_R1.fastq.gz \
    -2 ${PROJECT_DIR}/01_fastq/${SAMPLE}_R2.fastq.gz \
    -0 ${PROJECT_DIR}/01_fastq/${SAMPLE}_other.fastq.gz \
    -s ${PROJECT_DIR}/01_fastq/${SAMPLE}_singleton.fastq.gz \
    ${PROJECT_DIR}/01_fastq/${SAMPLE}.sorted.bam

# Clean up intermediates
rm ${PROJECT_DIR}/01_fastq/${SAMPLE}.unmapped.bam \
   ${PROJECT_DIR}/01_fastq/${SAMPLE}.mate_unmapped_R1.bam \
   ${PROJECT_DIR}/01_fastq/${SAMPLE}.mate_unmapped_R2.bam \
   ${PROJECT_DIR}/01_fastq/${SAMPLE}.all_unmapped.bam \
   ${PROJECT_DIR}/01_fastq/${SAMPLE}.sorted.bam
```

### CRAM Input (Hartwig / newer TCGA data)

```bash
# CRAMs require the reference genome used during alignment
# Hartwig typically uses GRCh38; check the CRAM header to confirm:
#   samtools view -H ${SAMPLE}.cram | grep "@SQ" | head -3

REFERENCE="/path/to/GRCh38.fa"   # must match the CRAM's alignment reference

samtools view -b -f 12 -F 256 \
    --reference ${REFERENCE} \
    ${PROJECT_DIR}/00_raw/${SAMPLE}.cram \
    -@ 16 -o ${PROJECT_DIR}/01_fastq/${SAMPLE}.unmapped.bam

samtools view -b -f 4 -F 264 \
    --reference ${REFERENCE} \
    ${PROJECT_DIR}/00_raw/${SAMPLE}.cram -@ 16 \
    -o ${PROJECT_DIR}/01_fastq/${SAMPLE}.mate_unmapped_R1.bam

samtools view -b -f 8 -F 260 \
    --reference ${REFERENCE} \
    ${PROJECT_DIR}/00_raw/${SAMPLE}.cram -@ 16 \
    -o ${PROJECT_DIR}/01_fastq/${SAMPLE}.mate_unmapped_R2.bam

# Then proceed with merge → sort → fastq → cleanup as above
```

### QC Checkpoint 1 — Read Extraction

```bash
# ═══════════════════════════════════════════════════════════
# QC-1: Verify read extraction completeness
# ═══════════════════════════════════════════════════════════

R1_COUNT=$(zcat ${PROJECT_DIR}/01_fastq/${SAMPLE}_R1.fastq.gz | wc -l | awk '{print $1/4}')
R2_COUNT=$(zcat ${PROJECT_DIR}/01_fastq/${SAMPLE}_R2.fastq.gz | wc -l | awk '{print $1/4}')
TOTAL_BAM=$(samtools view -c ${INPUT_BAM})

echo "SAMPLE: ${SAMPLE}" > ${PROJECT_DIR}/qc_checkpoints/QC1_${SAMPLE}.txt
echo "Total reads in BAM: ${TOTAL_BAM}" >> ${PROJECT_DIR}/qc_checkpoints/QC1_${SAMPLE}.txt
echo "R1 extracted: ${R1_COUNT}" >> ${PROJECT_DIR}/qc_checkpoints/QC1_${SAMPLE}.txt
echo "R2 extracted: ${R2_COUNT}" >> ${PROJECT_DIR}/qc_checkpoints/QC1_${SAMPLE}.txt
echo "Unmapped fraction: $(echo "scale=4; (${R1_COUNT}+${R2_COUNT})/${TOTAL_BAM}" | bc)" \
    >> ${PROJECT_DIR}/qc_checkpoints/QC1_${SAMPLE}.txt
echo "Pairs balanced: $([ "${R1_COUNT}" = "${R2_COUNT}" ] && echo "YES" || echo "NO — INVESTIGATE")" \
    >> ${PROJECT_DIR}/qc_checkpoints/QC1_${SAMPLE}.txt
```

### QC-1 Thresholds

| Metric | PASS | WARN | FAIL |
|--------|------|------|------|
| Unmapped fraction | 0.1–3% | 3–5% | <0.01% or >5% |
| R1 == R2 count | Yes | — | No (broken pairs) |

**Red flags:**
- Unmapped fraction >5%: possible alignment issue or unfiltered data
- Unmapped fraction <0.01%: data may already be pre-filtered of unmapped reads
- R1 ≠ R2: broken pairs — re-sort by name and re-extract

### Batch Processing

For cohort-scale runs (hundreds of TCGA/Hartwig samples), use GNU parallel to process
multiple samples simultaneously without oversubscribing CPU:

```bash
#!/bin/bash
# batch_extract.sh
# Usage: bash batch_extract.sh sample_list.txt /path/to/project [bam|cram] [/path/to/ref.fa]

SAMPLE_LIST="$1"
PROJECT_DIR="$2"
INPUT_TYPE="${3:-bam}"          # bam or cram
REFERENCE="${4:-}"              # required if INPUT_TYPE=cram
THREADS_PER_SAMPLE=8
MAX_PARALLEL=4                  # adjust to available CPUs / THREADS_PER_SAMPLE

extract_sample() {
    local SAMPLE="$1"
    local PROJECT_DIR="$2"
    local INPUT_TYPE="$3"
    local REFERENCE="$4"
    local THREADS="$5"

    local INPUT_FILE="${PROJECT_DIR}/00_raw/${SAMPLE}.${INPUT_TYPE}"

    if [ ! -f "${INPUT_FILE}" ]; then
        echo "SKIP: ${SAMPLE} — file not found" >&2
        return 1
    fi

    local REF_FLAG=""
    [ "${INPUT_TYPE}" = "cram" ] && REF_FLAG="--reference ${REFERENCE}"

    samtools view -b -f 12 -F 256 ${REF_FLAG} ${INPUT_FILE} \
        -@ ${THREADS} -o ${PROJECT_DIR}/01_fastq/${SAMPLE}.unmapped.bam
    samtools view -b -f 4 -F 264 ${REF_FLAG} ${INPUT_FILE} -@ ${THREADS} \
        -o ${PROJECT_DIR}/01_fastq/${SAMPLE}.mate_unmapped_R1.bam
    samtools view -b -f 8 -F 260 ${REF_FLAG} ${INPUT_FILE} -@ ${THREADS} \
        -o ${PROJECT_DIR}/01_fastq/${SAMPLE}.mate_unmapped_R2.bam

    samtools merge -@ ${THREADS} ${PROJECT_DIR}/01_fastq/${SAMPLE}.all_unmapped.bam \
        ${PROJECT_DIR}/01_fastq/${SAMPLE}.unmapped.bam \
        ${PROJECT_DIR}/01_fastq/${SAMPLE}.mate_unmapped_R1.bam \
        ${PROJECT_DIR}/01_fastq/${SAMPLE}.mate_unmapped_R2.bam

    samtools sort -n -@ ${THREADS} \
        ${PROJECT_DIR}/01_fastq/${SAMPLE}.all_unmapped.bam \
        -o ${PROJECT_DIR}/01_fastq/${SAMPLE}.sorted.bam

    samtools fastq -@ ${THREADS} \
        -1 ${PROJECT_DIR}/01_fastq/${SAMPLE}_R1.fastq.gz \
        -2 ${PROJECT_DIR}/01_fastq/${SAMPLE}_R2.fastq.gz \
        -0 ${PROJECT_DIR}/01_fastq/${SAMPLE}_other.fastq.gz \
        -s ${PROJECT_DIR}/01_fastq/${SAMPLE}_singleton.fastq.gz \
        ${PROJECT_DIR}/01_fastq/${SAMPLE}.sorted.bam

    rm ${PROJECT_DIR}/01_fastq/${SAMPLE}.unmapped.bam \
       ${PROJECT_DIR}/01_fastq/${SAMPLE}.mate_unmapped_R1.bam \
       ${PROJECT_DIR}/01_fastq/${SAMPLE}.mate_unmapped_R2.bam \
       ${PROJECT_DIR}/01_fastq/${SAMPLE}.all_unmapped.bam \
       ${PROJECT_DIR}/01_fastq/${SAMPLE}.sorted.bam

    # QC-1
    local R1=$(zcat ${PROJECT_DIR}/01_fastq/${SAMPLE}_R1.fastq.gz | wc -l | awk '{print $1/4}')
    local R2=$(zcat ${PROJECT_DIR}/01_fastq/${SAMPLE}_R2.fastq.gz | wc -l | awk '{print $1/4}')
    local TOTAL=$(samtools view -c ${INPUT_FILE})
    local FRAC=$(echo "scale=4; (${R1}+${R2})/${TOTAL}" | bc)
    local BALANCED=$([ "${R1}" = "${R2}" ] && echo "YES" || echo "NO — INVESTIGATE")

    cat > ${PROJECT_DIR}/qc_checkpoints/QC1_${SAMPLE}.txt << EOF
SAMPLE: ${SAMPLE}
Total reads in BAM/CRAM: ${TOTAL}
R1 extracted:            ${R1}
R2 extracted:            ${R2}
Unmapped fraction:       ${FRAC}
Pairs balanced:          ${BALANCED}
EOF
    echo "[$(date)] DONE: ${SAMPLE} | unmapped_frac=${FRAC} | balanced=${BALANCED}"
}

export -f extract_sample

# Run up to MAX_PARALLEL samples at once
parallel -j ${MAX_PARALLEL} \
    extract_sample {} "${PROJECT_DIR}" "${INPUT_TYPE}" "${REFERENCE}" "${THREADS_PER_SAMPLE}" \
    :::: "${SAMPLE_LIST}"
```

---

> **Next:** [03_HOST_DEPLETION.md](03_HOST_DEPLETION.md) — 3-pass host removal (most critical step)
