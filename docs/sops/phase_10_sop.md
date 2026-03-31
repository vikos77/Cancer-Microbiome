# Phase 10 SOP: Pipeline Integration

**Sample:** SRR15090802 (Wahida et al. gut virome, VLP-enriched healthy child)
**Input:** All Phase 1–9 outputs
**Completed:** 2026-02-26
**Status:** COMPLETE

---

## Overview

Phase 10 integrates all upstream results into three deliverables:

1. **`master_viral_table.tsv`**: per-viral-contig summary (119 rows × 43 columns)
2. **`qc_dashboard.txt`**: aggregated pass/fail status for all QC1–QC12 checkpoints
3. **`pipeline_summary.md`**: biological narrative tying all findings together
4. **`mobilome_summary.tsv`**: IS elements + top plasmid contigs

All outputs go to `pipeline_run/13_integration/virome/`.

---

## Section 1: Environment Setup

Integration requires only Python 3.10 + standard library (csv, os, sys, collections).
No additional conda packages needed.

**Critical:** Use explicit conda Python path to avoid Ubuntu 24.04 PATH issue:
```bash
PYTHON=$HOME/miniconda3/envs/claude_pipeline/bin/python3.10
```

---

## Section 2: Run Command

```bash
conda run -n claude_pipeline \
    $HOME/miniconda3/envs/claude_pipeline/bin/python3.10 \
    pipeline_run/scripts/integrate_pipeline.py \
    --base_dir $PROJECT_DIR \
    --out_dir pipeline_run/13_integration/virome \
    2>&1 | tee pipeline_run/logs/integration_virome.log
```

**Runtime:** ~5 seconds (pure Python, no external tools)

---

## Section 3: Master Viral Table

### Input files joined

| Column group | Source file |
|-------------|-------------|
| Contig quality | `09_viral_qc/virome/checkv_2tool/quality_summary.tsv` |
| Tool votes | `08_viral_contigs/virome/viral_classification_summary_2tool.tsv` |
| Length, GC, CDS density | `11_phage_annotation/virome/pharokka/pharokka_length_gc_cds_density.tsv` |
| PHROG functional categories | `11_phage_annotation/virome/pharokka/pharokka_cds_functions.tsv` |
| geNomad taxonomy | `08_viral_contigs/virome/genomad/clean_contigs_summary/clean_contigs_virus_summary.tsv` |
| INPHARED closest hit | `11_phage_annotation/virome/pharokka/pharokka_top_hits_mash_inphared.tsv` |
| PropagAtE activity | `10_prophages/virome/propagate_results/propagate_results.tsv` |

### Backbone: CheckV quality_summary.tsv
The CheckV file contains exactly the 119 consensus viral contigs and is used as the join
backbone. All other tables are left-joined onto it; missing values become 'NA'.

### Column schema (43 columns)

```
contig_id, length, gc_perc, cds_coding_density
checkv_quality, checkv_completeness, checkv_contamination, provirus
checkv_gene_count, checkv_viral_genes, checkv_host_genes
total_votes, genomad_hit, genomad_score, dvf_hit, dvf_score
genomad_taxonomy, genomad_topology, genomad_hallmarks
cds_total, cds_known, cds_unknown, cds_known_pct
phrog_head_packaging, phrog_tail, phrog_lysis, phrog_dna_metabolism
phrog_connector, phrog_integration_excision, phrog_moron_amg
phrog_other, phrog_transcription, trnas, tmrnas, crisprs
inphared_accession, inphared_mash_dist, inphared_hashes
inphared_description, inphared_host, inphared_family, inphared_ref_length
propagate_activity
```

### Sorting
Rows are sorted by CheckV quality tier (HQ → MQ → LQ → ND), then by completeness
descending within each tier.

### Key stats (SRR15090802)

| Metric | Value |
|--------|-------|
| Total rows | 119 |
| HQ / MQ / LQ / ND | 1 / 1 / 95 / 22 |
| Proviruses | 1 (metaspades_415, LQ, 3.35%) |
| Contigs with ≥1 known CDS | 98 (82.4%) |
| Total CDS | 674 |
| Known-function CDS | 216 (32.1%) |
| tRNAs | 32 |
| Exact INPHARED matches (dist=0) | 2 |

---

## Section 4: QC Dashboard

`qc_dashboard.txt` aggregates all QC checkpoint files.

### Status extraction logic
The script scans each QC file for a line containing `STATUS` + one of
`PASS / FAIL / PARTIAL / COMPLETE / WARN`. Two legacy files (QC1, QC2) use
a box-drawing format without an explicit STATUS line; these are handled by
a hardcoded fallback dict (`_IMPLICIT_STATUS`).

### QC files scanned

```
QC1_NA12878_chr22.txt         QC7_virome_SRR15090802.txt
QC2_NA12878_chr22.txt         QC8_virome_SRR15090802.txt
QC3_NA12878_chr22.txt         QC10_virome_SRR15090802.txt
QC4_NA12878_chr22.txt         QC11_virome_SRR15090802.txt
QC5_NA12878_chr22.txt         QC12_virome_SRR15090802.txt
QC5_virome_SRR15090802.txt
QC6_virome_SRR15090802.txt
```

Note: QC9 does not exist (there is no Phase 9 QC in this numbering; Phases 8/9 = QC8–12).

### Result (SRR15090802): 10 PASS | 2 PARTIAL | 0 FAIL | 0 MISSING

---

## Section 5: Mobilome Summary Table

`mobilome_summary.tsv` has one row per IS element (from ISEScan .sum files)
plus the top 20 plasmid contigs by geNomad score.

Columns: `bin, element_type, contig_id, family, n_copies, pct_genome, bp_covered, contig_len, notes`

**Note on ISEScan .sum paths:** ISEScan writes output under `isescan/<bin>/metabat2/<bin>.fa.sum`
(it mirrors the input directory structure). bin.3 has 0 IS elements; no .sum file is written.

---

## Section 6: Pipeline Summary Report

`pipeline_summary.md` is a narrative document covering:
- Data acquisition and host depletion
- Assembly statistics
- Viral contig identification and quality
- Phage annotation highlights (complete phage, exact matches)
- Prophage analysis
- Host prediction (pending iPHoP)
- Mobilome findings
- Biosafety summary (0 ARGs, 0 virulence)
- Pipeline validation assessment
- Readiness statement for TCGA/Hartwig

---

## Section 7: Pending Items

Three tools are pending server/database availability; none block QC-13 PASS:

### PHASTEST (cluster down on 2026-02-26)
```bash
python3 pipeline_run/scripts/phastest_submit.py \
    --batch_ids BB_3288ee5236,BB_51980fe431,BB_6f5e98c05f,BB_3ba4ff53bf \
    --out_dir pipeline_run/10_prophages/virome/phastest_results
```
After results available: update QC-10 with PHASTEST prophage annotations.

### iPHoP (~280 GB database, downloading)
```bash
# Monitor download:
tail -f pipeline_run/logs/iphop_db_download.log

# Run after completion:
conda run -n iphop_env bash -c "
iphop predict \
    --fa_file $PROJECT_DIR/pipeline_run/08_viral_contigs/virome/consensus_viral_contigs_2tool.fasta \
    --db_dir $PROJECT_DIR/databases/iphop/ \
    --out_dir $PROJECT_DIR/pipeline_run/10_prophages/virome/iphop \
    -t 4
" 2>&1 | tee pipeline_run/logs/iphop_predict.log
```
After results: update QC-11 from PARTIAL PASS to PASS.

### mobileOG-db (server 504)
```bash
# When server recovers:
wget https://mobileogdb.flsi.cloud.vt.edu/entries/mobileOG-db_beatrix-1.6.All.faa \
    -P databases/mobileog/

conda run -n claude_pipeline bash -c "
diamond makedb \
    --in $PROJECT_DIR/databases/mobileog/mobileOG-db_beatrix-1.6.All.faa \
    --db $PROJECT_DIR/databases/mobileog/mobileOG-db \
    --threads 4

diamond blastp \
    --query $PROJECT_DIR/pipeline_run/12_mobilome/virome/mobileog/all_bin_proteins.faa \
    --db $PROJECT_DIR/databases/mobileog/mobileOG-db \
    --out $PROJECT_DIR/pipeline_run/12_mobilome/virome/mobileog/mobileog_all_hits.tsv \
    --outfmt 6 --evalue 1e-10 \
    --query-cover 60 --subject-cover 60 \
    --threads 4
"
```
After results: update QC-12 from COMPLETE to COMPLETE ✅.

---

## Section 8: Output Files Summary

```
pipeline_run/
├── 13_integration/virome/
│   ├── master_viral_table.tsv      119 rows × 43 cols (all viral contig annotations)
│   ├── qc_dashboard.txt            QC1–QC12 status dashboard
│   ├── mobilome_summary.tsv        IS elements + top plasmids
│   └── pipeline_summary.md         Full biological narrative
├── qc_checkpoints/
│   └── QC13_virome_SRR15090802.txt STATUS: PASS
└── logs/
    └── integration_virome.log      Script stdout/stderr
```

---

## Section 9: Gotchas

1. **System Python 3.12 in PATH**: `conda run -n claude_pipeline python3` resolves to
   `/usr/bin/python3` (Python 3.12, system pandas has NumPy 2.x conflict).
   Fix: always use explicit path `bin/python3.10` from conda env.

2. **pharokka_cds_functions.tsv format**: Three columns: `Description, Count, contig`
   (not the standard two-column wide format). The script pivots this table manually
   (defaultdict pivot, not pandas) to avoid the system Python issue.

3. **ISEScan .sum path**: Output is `isescan/<bin>/metabat2/<bin>.fa.sum`, not
   `isescan/<bin>/<bin>.fa.sum`. The `metabat2/` subdir mirrors the input filename
   structure (`pipeline_run/07_mags/virome/metabat2/bin.N.fa`).

4. **QC1/QC2 STATUS**: These files (written early in the project) use box-drawing
   decorators instead of a `STATUS:` line. The script uses a hardcoded fallback dict.

5. **geNomad virus_summary has 269 rows, not 119**: The geNomad summary covers all
   geNomad-called viral contigs (before 2-tool consensus filtering). The 119 consensus
   contigs are a subset; 21 consensus contigs are DVF-only (not in geNomad summary),
   so `genomad_taxonomy` will be NA for those rows.

---

## Section 10: Adapting to TCGA/Hartwig Tumour WGS

When applying Phase 10 to production data, update `integrate_pipeline.py` to:

1. Add iPHoP host prediction column (`Host_prediction_to_genus_m90.csv`)
2. Add PHASTEST result column (prophage classification per bin)
3. Add mobileOG-db hit column (MGE protein count per bin protein)
4. Add clinical metadata columns (tumour type, patient ID, MSI status, etc.)
5. Add GTDB-Tk taxonomy for bacterial MAG bins (when present)
6. Add cancer-MGE positive control hits (cag PAI, pks island, FadA, BFT)
7. Consider a pandas-based rewrite when running on cluster (no PATH issue)
