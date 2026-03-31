# Phase 9 SOP — Mobilome Detection

**Sample:** SRR15090802 (Wahida et al. gut virome, VLP-enriched healthy child)
**Input:** 4 MetaBAT2 bins (`pipeline_run/07_mags/virome/metabat2/bin.{1,2,3,4}.fa`)
**Completed:** 2026-02-26
**Status:** COMPLETE (mobileOG-db pending server recovery)

---

## Overview

Phase 9 systematically profiles mobile genetic elements (MGEs) in the MAG bins from Phase 6.
For cancer microbiome research, the key MGEs are:

| MGE type | Tool | Cancer relevance |
|----------|------|-----------------|
| Plasmids | geNomad (Phase 7 reuse) | Oncotoxin delivery (pks island), ARG transfer |
| Insertion sequences | ISEScan 1.7.2.3 | Genome plasticity, IS-mediated gene disruption |
| Integrons | IntegronFinder 2.0.5 | ARG cassette acquisition, virulence gene capture |
| All MGE proteins | mobileOG-db + DIAMOND | Comprehensive MGE gene catalogue |

**Note on VLP virome context:** The 4 MetaBAT2 bins from this VLP-enriched sample are
phage-dominated (not bacterial MAGs). Results reflect phage/co-purified DNA rather than
a classic bacterial mobilome. For production TCGA/Hartwig tumor WGS, the MAG bins would
be bacterial, enabling detection of cancer-relevant MGEs (cag PAI, pks island, etc.).

---

## Section 1: Environment Setup

All mobilome tools installed into the existing `claude_pipeline` conda environment:

```bash
mamba install -n claude_pipeline -c bioconda -c conda-forge \
    isescan=1.7.2.3 \
    integron_finder=2.0.5 \
    diamond=2.1.11 \
    -y
```

**Versions installed:**
- ISEScan 1.7.2.3 (bioconda)
- IntegronFinder 2.0.5 (bioconda)
- DIAMOND 2.1.11 (bioconda)
- Prodigal 2.6.3 (already present in claude_pipeline)

### ISEScan PATH issue (Ubuntu 24.04)

ISEScan fails if called as `isescan.py` directly because `/usr/bin/python3` (system Python 3.12)
precedes the conda env's Python 3.10 in the user's PATH, and ISEScan has `#!/usr/bin/env python3`
as its shebang.

**Fix:** Use the explicit conda Python path:
```bash
# Required wrapper
PYTHON=$HOME/miniconda3/envs/claude_pipeline/bin/python3.10
ISESCAN=$HOME/miniconda3/envs/claude_pipeline/bin/isescan.py
export LD_LIBRARY_PATH=$HOME/miniconda3/envs/claude_pipeline/lib:$LD_LIBRARY_PATH
$PYTHON $ISESCAN --seqfile <input.fa> --output <outdir> --nthread 4
```

---

## Section 2: Output Directory Setup

```bash
mkdir -p pipeline_run/12_mobilome/virome/{plasmids,isescan,integrons,mobileog}
```

---

## Section 3: Plasmid Detection (geNomad — reuse from Phase 7)

geNomad was run in Phase 7 on all 20,547 clean contigs. The plasmid summary is directly
reused — no additional computation required.

```bash
cp pipeline_run/08_viral_contigs/virome/genomad/clean_contigs_summary/clean_contigs_plasmid_summary.tsv \
    pipeline_run/12_mobilome/virome/plasmids/
cp pipeline_run/08_viral_contigs/virome/genomad/clean_contigs_summary/clean_contigs_plasmid.fna \
    pipeline_run/12_mobilome/virome/plasmids/
```

### Results

| Metric | Value |
|--------|-------|
| Total plasmid contigs | 124 |
| With conjugation genes | 39 |
| With AMR genes | 0 |
| DTR (circular) topology | 5 |
| Size range | 514–29,043 bp |

Notable contigs with conjugation machinery:
- metaspades_392 (4,429 bp, MOBV)
- megahit_1934 (3,472 bp, MOBP1)
- metaspades_103 (9,069 bp, F_traW/U/V/F/H/G — F-type conjugation)
- metaspades_338 (4,838 bp, F_traV/K/E/L)

**Interpretation:** 39 plasmids carry MOB relaxases or F-type tra genes suggesting
conjugative transfer capacity. In a cancer microbiome dataset, these would represent
potential vectors for horizontal gene transfer of virulence or resistance determinants.

---

## Section 4: Insertion Sequences (ISEScan)

### Command

```bash
PYTHON=$HOME/miniconda3/envs/claude_pipeline/bin/python3.10
ISESCAN=$HOME/miniconda3/envs/claude_pipeline/bin/isescan.py
MAGS=$PROJECT_DIR/pipeline_run/07_mags/virome/metabat2

conda run -n claude_pipeline bash -c "
export LD_LIBRARY_PATH=$HOME/miniconda3/envs/claude_pipeline/lib:\$LD_LIBRARY_PATH
PYTHON=$HOME/miniconda3/envs/claude_pipeline/bin/python3.10
ISESCAN=$HOME/miniconda3/envs/claude_pipeline/bin/isescan.py
MAGS=$PROJECT_DIR/pipeline_run/07_mags/virome/metabat2
OUTBASE=$PROJECT_DIR/pipeline_run/12_mobilome/virome/isescan

for i in 1 2 3 4; do
    \$PYTHON \$ISESCAN \
        --seqfile \${MAGS}/bin.\${i}.fa \
        --output \${OUTBASE}/bin.\${i} \
        --nthread 4
done
" 2>&1 | tee pipeline_run/logs/isescan_virome.log
```

**Output structure:** `isescan/bin.{N}/metabat2/bin.{N}.fa.{sum,tsv,csv,is.fna}`
Note: If 0 IS elements found, no output files are written (only intermediate proteome/hmm dirs).

### Results

| Bin | Contigs | IS elements | IS families |
|-----|---------|-------------|-------------|
| bin.1 | 497 | 9 | IS256(3), IS3(3), IS110, IS30, IS481 |
| bin.2 | 276 | 3 | IS4(2), IS91 |
| bin.3 | 66 | 0 | (none) |
| bin.4 | 68 | 4 | IS630(2), IS481, IS256 |
| **Total** | **907** | **16** | 7 families |

IS family breakdown:
- IS256 ×4, IS3 ×3, IS4 ×2, IS630 ×2, IS481 ×2, IS110 ×1, IS30 ×1, IS91 ×1

Notable observations:
- IS elements found on metaviralspades contigs (metaviralspades_16, metaviralspades_46)
  — possibly bacterial DNA co-assembled with viral sequences
- IS91 in bin.2 (metaspades_1370): 81% of the 2,003 bp contig is transposase
  — likely a contig composed almost entirely of a transposase fragment
- IS630 belongs to the Tc1/mariner superfamily — active in many organisms

---

## Section 5: Integrons (IntegronFinder)

### Command

```bash
conda run -n claude_pipeline bash -c "
MAGS=$PROJECT_DIR/pipeline_run/07_mags/virome/metabat2
OUTBASE=$PROJECT_DIR/pipeline_run/12_mobilome/virome/integrons

for i in 1 2 3 4; do
    integron_finder \${MAGS}/bin.\${i}.fa \
        --outdir \${OUTBASE}/bin.\${i} \
        --cpu 4 \
        --local-max \
        --linear
done
" 2>&1 | tee pipeline_run/logs/integron_finder_virome.log
```

**Note on `--linear` flag:** Use `--linear` for individual phage/metagenomic contigs.
The default assumes circular topology; `--linear` is more appropriate for fragmented
assemblies where circular topology cannot be assumed.

**Gotcha:** `--no-func-annot` does NOT exist in IntegronFinder 2.0.5.
Functional annotation is disabled by default (must explicitly pass `--func-annot` to enable).

**Output structure:** `integrons/bin.{N}/Results_Integron_Finder_bin.{N}/bin.{N}.{summary,integrons}`

### Results

| Bin | Complete integrons | CALIN | In0 (degenerate) |
|-----|--------------------|-------|-----------------|
| bin.1 | 0 | 0 | 0 |
| bin.2 | 0 | 0 | 0 |
| bin.3 | 0 | 0 | 0 |
| bin.4 | 0 | 0 | 0 |

**Result: 0 integrons across all bins.** Expected for phage-dominated bins.

Integron categories:
- **Complete integron**: has integrase (IntI) + attI site + ≥1 attC site
- **CALIN**: Class 1 integron fragment — attC sites without adjacent integrase
- **In0**: degenerate integron — integrase without attC sites

---

## Section 6: mobileOG-db (DIAMOND)

### Database download

```bash
mkdir -p databases/mobileog

# Primary source (use when server is available)
wget https://mobileogdb.flsi.cloud.vt.edu/entries/mobileOG-db_beatrix-1.6.All.faa \
    -P databases/mobileog/

# Build DIAMOND database
conda run -n claude_pipeline bash -c "
diamond makedb \
    --in $PROJECT_DIR/databases/mobileog/mobileOG-db_beatrix-1.6.All.faa \
    --db $PROJECT_DIR/databases/mobileog/mobileOG-db \
    --threads 4
"
```

**Status:** Server returning 504 Gateway Timeout on 2026-02-26. Retry when available.

### Protein prediction (pre-computed, ready for DIAMOND)

```bash
conda run -n claude_pipeline bash -c "
MAGS=$PROJECT_DIR/pipeline_run/07_mags/virome/metabat2
OUTBASE=$PROJECT_DIR/pipeline_run/12_mobilome/virome/mobileog

for i in 1 2 3 4; do
    prodigal -i \${MAGS}/bin.\${i}.fa \
        -a \${OUTBASE}/bin.\${i}_proteins.faa \
        -p meta -q
done
cat \${OUTBASE}/bin.*_proteins.faa > \${OUTBASE}/all_bin_proteins.faa
"
```
Result: **4,848 proteins** in `pipeline_run/12_mobilome/virome/mobileog/all_bin_proteins.faa`

### DIAMOND search (run when database available)

```bash
conda run -n claude_pipeline bash -c "
diamond blastp \
    --query $PROJECT_DIR/pipeline_run/12_mobilome/virome/mobileog/all_bin_proteins.faa \
    --db $PROJECT_DIR/databases/mobileog/mobileOG-db \
    --out $PROJECT_DIR/pipeline_run/12_mobilome/virome/mobileog/mobileog_all_hits.tsv \
    --outfmt 6 \
    --evalue 1e-10 \
    --query-cover 60 \
    --subject-cover 60 \
    --threads 4 \
    2>&1 | tee $PROJECT_DIR/pipeline_run/logs/mobileog_diamond.log
"
```

---

## Section 7: Cancer-Associated MGE Cross-Reference

For production cancer microbiome analysis, add these positive controls when the relevant
bacteria are present in MAGs:

| Organism | MGE | Relevance | Detection approach |
|----------|-----|-----------|-------------------|
| H. pylori | cag PAI (type IV secretion island) | CagA oncoprotein → gastric cancer | BLAST CagA/VirB4 vs MAG proteins |
| E. coli B2 | pks genomic island (colibactin) | CRC mutational signature | BLAST ClbA-S vs MAG proteins |
| F. nucleatum | FadA adhesin locus | CRC invasion | BLAST FadA vs MAG proteins |
| B. fragilis | bft toxin locus | CRC tumour promotion | BLAST BFT (fragilysin) vs MAG proteins |

```bash
# Example: check for colibactin biosynthesis genes in E. coli MAGs
# First get pks island protein sequences from GenBank (NC_011741.1)
# Then BLAST against all MAG proteins
conda run -n claude_pipeline bash -c "
blastp \
    -query databases/cancer_mge/pks_island_proteins.faa \
    -db pipeline_run/12_mobilome/virome/mobileog/all_bin_proteins.faa \
    -outfmt 6 -evalue 1e-10 \
    -out pipeline_run/12_mobilome/virome/cancer_mge_hits.tsv \
    -num_threads 4
"
```

For this VLP virome test dataset: no bacteria known to carry cancer-associated MGEs are
present — this cross-reference would be performed on TCGA/Hartwig tumor WGS MAGs.

---

## Section 8: Key Gotchas

1. **ISEScan PATH (Ubuntu 24.04):** `python3` resolves to system Python 3.12, not conda Python 3.10.
   `isescan.py` imports `libssw.so` (SSW aligner) from the conda env's lib. Fix: use explicit
   paths `python3.10 isescan.py` + `LD_LIBRARY_PATH=/conda/env/lib:...`.

2. **ISEScan 0-result behaviour:** When no IS elements are found, ISEScan writes no output
   files (only intermediate proteome/ and hmm/ dirs remain). Do not mistake absence of .sum
   file for a crashed run — check the log for "ISEScan ends at" message.

3. **IntegronFinder `--no-func-annot` does not exist:** The flag in the plan template is wrong.
   Functional annotation is OFF by default; `--func-annot` enables it. Use `--linear` for
   metagenome contigs (do not assume circular topology).

4. **mobileOG-db server availability:** `mobileogdb.flsi.cloud.vt.edu` is a VT cloud service
   that experiences periodic outages. Pre-compute proteins with Prodigal; run DIAMOND when
   the server comes back.

5. **Prodigal output format:** Without `-f gff`, prodigal writes GenBank-format `.gbk` to
   stdout. Use `-q` (quiet) to suppress non-FAA output, and redirect stdout or use `-o`.

6. **VLP virome mobilome interpretation:** These bins are phage-dominated. IS elements found on
   phage contigs (especially metaviralspades) likely represent co-assembled bacterial DNA,
   not phage IS elements per se. Zero integrons is the expected result.

---

## Section 9: Output Files Summary

```
pipeline_run/
├── 12_mobilome/virome/
│   ├── plasmids/
│   │   ├── clean_contigs_plasmid_summary.tsv   # 124 plasmid contigs (geNomad)
│   │   └── clean_contigs_plasmid.fna           # Plasmid sequences
│   ├── isescan/
│   │   ├── bin.1/metabat2/bin.1.fa.{sum,tsv,csv,is.fna}   # 9 IS elements
│   │   ├── bin.2/metabat2/bin.2.fa.{sum,tsv,csv,is.fna}   # 3 IS elements
│   │   ├── bin.3/                               # 0 IS (no output files)
│   │   └── bin.4/metabat2/bin.4.fa.{sum,tsv,csv,is.fna}   # 4 IS elements
│   ├── integrons/
│   │   ├── bin.1/Results_Integron_Finder_bin.1/ # 0 integrons
│   │   ├── bin.2/Results_Integron_Finder_bin.2/ # 0 integrons
│   │   ├── bin.3/Results_Integron_Finder_bin.3/ # 0 integrons
│   │   └── bin.4/Results_Integron_Finder_bin.4/ # 0 integrons
│   └── mobileog/
│       ├── bin.{1,2,3,4}_proteins.faa           # Per-bin Prodigal proteins
│       └── all_bin_proteins.faa                 # Combined (4,848 proteins; ready for DIAMOND)
│       [to add: mobileog_all_hits.tsv when mobileOG-db server recovers]
└── qc_checkpoints/
    └── QC12_virome_SRR15090802.txt              # Complete
```
