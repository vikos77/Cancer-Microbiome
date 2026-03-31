#!/usr/bin/env bash
set -euo pipefail

PROJECT=/home/vicky/Microbiome_cancer/pipeline_run
LOGS=${PROJECT}/logs

for SAMPLE in CRC_ERR2726414 control_ERR2726517; do
    # Map sample name to ERR accession
    if [[ "$SAMPLE" == "CRC_ERR2726414" ]]; then
        ACC="ERR2726414"
    else
        ACC="ERR2726517"
    fi

    IN_R1="${PROJECT}/02_host_depleted/${SAMPLE}/${ACC}_R1.clean.fastq.gz"
    IN_R2="${PROJECT}/02_host_depleted/${SAMPLE}/${ACC}_R2.clean.fastq.gz"
    OUTDIR="${PROJECT}/03_qc_filtered/${SAMPLE}"

    echo "[$(date)] === Phase 3 QC: ${SAMPLE} ==="

    # Step 1: fastp
    echo "[$(date)] Running fastp on ${SAMPLE}..."
    conda run -n claude_pipeline bash -c "
        fastp \
            -i ${IN_R1} -I ${IN_R2} \
            -o ${OUTDIR}/${ACC}_R1.qc.fastq.gz \
            -O ${OUTDIR}/${ACC}_R2.qc.fastq.gz \
            --qualified_quality_phred 20 \
            --length_required 60 \
            --low_complexity_filter \
            --complexity_threshold 30 \
            --correction \
            --dedup \
            --thread 4 \
            --json ${OUTDIR}/${ACC}_fastp.json \
            --html ${OUTDIR}/${ACC}_fastp.html \
            2>&1
    " | tee ${LOGS}/fastp_${SAMPLE}.log
    echo "[$(date)] fastp done for ${SAMPLE}"

    # Step 2: bbduk entropy filter
    echo "[$(date)] Running bbduk on ${SAMPLE}..."
    /usr/bin/bbduk.sh \
        in1=${OUTDIR}/${ACC}_R1.qc.fastq.gz \
        in2=${OUTDIR}/${ACC}_R2.qc.fastq.gz \
        out1=${OUTDIR}/${ACC}_R1.final.fastq.gz \
        out2=${OUTDIR}/${ACC}_R2.final.fastq.gz \
        entropy=0.5 entropywindow=50 threads=4 \
        stats=${OUTDIR}/${ACC}_bbduk_stats.txt \
        2>&1 | tee ${LOGS}/bbduk_${SAMPLE}.log
    echo "[$(date)] bbduk done for ${SAMPLE}"

    # Step 3: FastQC
    echo "[$(date)] Running FastQC on ${SAMPLE}..."
    conda run -n claude_pipeline bash -c "
        fastqc -t 4 \
            -o ${OUTDIR}/ \
            ${OUTDIR}/${ACC}_R1.final.fastq.gz \
            ${OUTDIR}/${ACC}_R2.final.fastq.gz \
            2>&1
    " | tee ${LOGS}/fastqc_${SAMPLE}.log
    echo "[$(date)] FastQC done for ${SAMPLE}"

    # Step 4: parse metrics and write QC-3 checkpoint
    echo "[$(date)] Writing QC-3 checkpoint for ${SAMPLE}..."
    conda run -n claude_pipeline python3 << PYEOF
import json, gzip, re, datetime

sample = "${SAMPLE}"
acc    = "${ACC}"
outdir = "${OUTDIR}"
qc_dir = "${PROJECT}/qc_checkpoints"

# Parse fastp JSON
with open(f"{outdir}/{acc}_fastp.json") as f:
    fp = json.load(f)

reads_in        = fp["summary"]["before_filtering"]["total_reads"] // 2
reads_passed    = fp["filtering_result"]["passed_filter_reads"] // 2
low_qual        = fp["filtering_result"]["low_quality_reads"] // 2
too_short       = fp["filtering_result"]["too_short_reads"] // 2
low_complex     = fp["filtering_result"]["low_complexity_reads"] // 2
dup_rate        = fp["duplication"]["rate"] * 100
q20_post        = fp["summary"]["after_filtering"]["q20_rate"] * 100
q30_post        = fp["summary"]["after_filtering"]["q30_rate"] * 100
pct_passed_fp   = reads_passed / reads_in * 100

# Parse bbduk log for input/output counts
bbduk_in = bbduk_out = bbduk_removed = 0
with open(f"{outdir}/{acc}_bbduk_stats.txt") as f:
    for line in f:
        m = re.match(r"#Total\s+(\d+)", line)
        if m: bbduk_in = int(m.group(1)) // 2
with open("${LOGS}/bbduk_${SAMPLE}.log") as f:
    content = f.read()
m = re.search(r"Result:\s+(\d+) reads", content)
if m: bbduk_out = int(m.group(1)) // 2
m = re.search(r"Total Removed:\s+(\d+) reads", content)
if m: bbduk_removed = int(m.group(1)) // 2
else: bbduk_removed = bbduk_in - bbduk_out if bbduk_out else 0

final_pairs = bbduk_out if bbduk_out else reads_passed
pct_survive = final_pairs / reads_in * 100
status = "PASS" if pct_survive >= 50 else ("WARN" if pct_survive >= 30 else "FAIL")

chk = f"""QC-3 Read-Level QC Checkpoint
Sample: {sample}
Date: {datetime.date.today()}
Tools: fastp 0.23.4, bbduk (system), FastQC 0.12.1
Input files: {acc}_R1/R2.clean.fastq.gz
═══════════════════════════════════════════════════════════
INPUT (from Phase 2):           {reads_in:>12,} pairs

--- fastp (Q20, len≥60, dedup, complexity≥30) ---
Passed filter:                  {reads_passed:>12,} pairs  ({pct_passed_fp:.1f}%)
Low quality removed:            {low_qual:>12,} ({low_qual/reads_in*100:.1f}%)
Too short removed:              {too_short:>12,} ({too_short/reads_in*100:.1f}%)
Low complexity removed:         {low_complex:>12,} ({low_complex/reads_in*100:.1f}%)
Duplication rate:               {dup_rate:>11.1f}%
Q20 rate (post-filter):         {q20_post:>11.1f}%
Q30 rate (post-filter):         {q30_post:>11.1f}%

--- bbduk (entropy≥0.5, window=50bp) ---
Input:                          {bbduk_in:>12,} pairs
Low entropy removed:            {bbduk_removed:>12,} ({bbduk_removed/bbduk_in*100:.2f}% if bbduk_in else 0)
Output (FINAL):                 {final_pairs:>12,} pairs

FINAL OUTPUT:                   {final_pairs:>12,} pairs  ({pct_survive:.1f}% of input) → Phase 4 input
QC-3 STATUS: {status}
═══════════════════════════════════════════════════════════
OUTPUT FILES:
  Final R1: pipeline_run/03_qc_filtered/{sample}/{acc}_R1.final.fastq.gz
  Final R2: pipeline_run/03_qc_filtered/{sample}/{acc}_R2.final.fastq.gz
  fastp report: pipeline_run/03_qc_filtered/{sample}/{acc}_fastp.html
  FastQC: pipeline_run/03_qc_filtered/{sample}/{acc}_R1.final_fastqc.html
═══════════════════════════════════════════════════════════
"""
with open(f"{qc_dir}/QC3_{sample}.txt", "w") as f:
    f.write(chk)
print(chk)
PYEOF

    echo "[$(date)] === ${SAMPLE} Phase 3 COMPLETE ==="
done

echo "[$(date)] Phase 3 complete for all samples."
