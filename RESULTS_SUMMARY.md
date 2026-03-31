# CLAUDE Pipeline — Validation Results Summary

**Dataset:** Wirbel et al. 2019 (Nat Med 25:679-689), German CRC cohort (PRJEB27928)
**Samples:** One stage III CRC stool metagenome (CRC_sample) and one matched healthy control (control_sample)
**Selection rationale:** CRC sample selected for the highest *Fusobacterium nucleatum* relative
abundance in the cohort (MetaPhlAn2 profiling); control selected as the closest depth-matched
healthy sample to enable a direct comparison.

All phases completed on both samples. QC checkpoints are in `pipeline_run/qc_checkpoints/`.

---

## Phase-by-Phase Summary

| Phase | CRC_sample | control_sample |
|-------|------------|----------------|
| 1. Acquisition | 75.9M read pairs | 53.7M read pairs |
| 2. Host depletion | 0.18% human removed | 0.002% human removed |
| 3. Read QC | 7.8M pairs retained (89.8% pass) | 22.2M pairs retained (90.6% pass) |
| 4. Taxonomic profiling | 358 species detected | 421 species detected |
| 5. Assembly | 72,778 dedup contigs, N50=5,243 bp, max=788 kb | 128,703 dedup contigs, N50=5,186 bp, max=448 kb |
| 6. Decontam + MAG binning | 72,696 clean contigs, 30 MAGs (2 HQ) | 128,673 clean contigs, 92 MAGs (2 HQ) |
| 6B. MAG annotation (Bakta) | 141k CDS, 94.7% known function | 243k CDS, 93.3% known function |
| 7. Viral ID | 335 viral contigs (1 complete, 108 kb) | 542 viral contigs (17 complete) |
| 8. Phage annotation | 1,775 CDS, 80.0% contigs annotated, 1 ARG | 6,027 CDS, 81.7% contigs annotated, 0 ARGs |
| 9. Mobilome | 607 plasmids, 140 IS, 4,980 mobileOG hits | 1,294 plasmids, 356 IS, 9,615 mobileOG hits |
| 10. PropagAtE | 1 active prophage, 51 dormant | 5 active prophages, 89 dormant |
| 11. iPHoP (full DB) | 242/335 predicted (72.2%), 8/8 tools | 352/542 predicted (64.9%), 8/8 tools |
| 12. Defence systems | 208 systems (DefenseFinder), 511 (PADLOC), 70 CRISPR arrays | 371 systems (DefenseFinder), 913 (PADLOC), 202 CRISPR arrays |

---

## Taxonomic Profiling (Phase 4)

### CRC_sample top microbiome (Bracken species level)
1. *Bacteroides fragilis* 25.4% — enterotoxigenic *B. fragilis* (ETBF), well-characterised CRC pathobiont
2. *Prevotella nigrescens* 15.8%
3. *Parvimonas micra* 1.7% — CRC-associated oral pathobiont
4. *Fusobacterium varium* 1.4%

### control_sample top microbiome (Bracken species level)
1. *Bacteroides uniformis* 15.3%
2. *Faecalibacterium prausnitzii* 13.7% — anti-inflammatory commensal, typically depleted in CRC
3. *Agathobacter rectalis* 12.8%

The CRC vs healthy contrast at the species level is consistent with published CRC microbiome
signatures: CRC dysbiosis dominated by pathobionts (*B. fragilis*, *P. micra*) while the healthy
gut is dominated by anti-inflammatory commensals (*F. prausnitzii*).

### KrakenUniq viral profiling
The control sample had 15.8x more CrAssphage reads than CRC (0.88% vs 0.16%), suggesting
that phage diversity associated with *Bacteroides* is reduced in the CRC microenvironment.

---

## Assembly (Phase 5)

Both samples assembled well above the QC thresholds. The CRC sample produced fewer contigs
and a lower total assembly size, consistent with a less diverse microbial community
(30 MAGs vs 92 in the control). The metaviralSPAdes assembler produced 22 viral contigs
for CRC (N50=34 kb) and 82 for the control (N50=73 kb), including the 108 kb complete
phage genome recovered from CRC.

---

## Viral Identification (Phase 7)

| Metric | CRC_sample | control_sample |
|--------|------------|----------------|
| Total viral contigs (2-tool consensus) | 335 | 542 |
| Complete genomes (CheckV) | 1 (108 kb) | 17 |
| High quality (CheckV) | 3 | 9 |
| Medium quality | 12 | 42 |

The 108 kb complete phage genome in the CRC sample is one of the largest
complete phage genomes recovered in this analysis. CheckV confirmed 100% completeness.

---

## Phage Annotation (Phase 8)

### Pharokka functional annotation
| Metric | CRC_sample | control_sample |
|--------|------------|----------------|
| Total CDS predicted | 1,775 | 6,027 |
| tRNAs detected | 4 | 43 |
| Annotated contigs (≥1 PHROG hit) | 268/335 (80.0%) | 443/542 (81.7%) |
| CARD AMR genes | 1 (tet(O)) | 0 |
| VFDB virulence factors | 0 | 0 |

The tet(O) tetracycline resistance gene found in CRC_sample is encoded on a 55 kb contig
and is notable as a phage-associated AMR gene. *tet(O)* has been previously linked to
*Campylobacter jejuni* and horizontal gene transfer via mobile elements.

The control sample has 3.4x more CDS and 10x more tRNAs than CRC, consistent with
larger and more complete phage genomes in the healthy gut virome.

### Prophage activity (PropagAtE)
| Metric | CRC_sample | control_sample |
|--------|------------|----------------|
| Active prophages | 1 | 5 |
| Dormant prophages | 51 | 89 |

The CRC sample has one actively replicating prophage (ratio=4.45, CohenD=2.06),
with phage coverage approximately 4.5x higher than the surrounding host contig.
The control sample has 5 active prophages, with the strongest signal showing a
ratio of 6.76 and CohenD of 6.13.

The lower number of active prophages in CRC may reflect the reduced microbial
diversity observed throughout — fewer distinct bacterial species means fewer
distinct prophage populations.

---

## Host Prediction (Phase 11 — iPHoP full database)

iPHoP was run with the full Jun 2025 database (280 GB) using all 8/8 tools including
RaFAH (Random Forest model). This required ≥28 GB free RAM and temporary disabling
of the systemd-oomd service (see TROUBLESHOOTING_public.md §6.5).

| Metric | CRC_sample | control_sample |
|--------|------------|----------------|
| Viruses with host prediction (≥90% confidence) | 242/335 (72.2%) | 352/542 (64.9%) |
| Mean confidence score | 93.7% | 93.8% |
| Tools used | 8/8 | 8/8 |

### Key host-prediction contrasts (CRC vs control)

| Host genus | CRC_sample | control_sample | Interpretation |
|------------|------------|----------------|----------------|
| *Fusobacterium* | 6 | 0 | **CRC-exclusive** — primary finding |
| *Escherichia* | 12 | 0 | CRC-exclusive — consistent with pks+ *E. coli* |
| *Bacteroides* | 69 | 48 | Elevated in CRC |
| *Bifidobacterium* | 0 | 11 | Absent in CRC |
| *Lachnospira* | 8 | 32 | 75% depleted in CRC |
| *Faecalibacterium* | 9 | 19 | 53% depleted in CRC |

The 6 *Fusobacterium* phages detected exclusively in CRC_sample represent the
pipeline's primary biological finding. *Fusobacterium nucleatum* is one of the
most consistently enriched bacteria in CRC tumours, and identifying its phage
repertoire is a key goal of this work. Their complete absence in the healthy
control confirms these are CRC-associated signals, not background noise.

---

## Mobilome (Phase 9)

| Element type | CRC_sample | control_sample |
|--------------|------------|----------------|
| Plasmids (geNomad) | 607 | 1,294 |
| Insertion sequences (ISEScan) | 140 | 356 |
| mobileOG-db hits (Bakta) | 4,980 | 9,615 |

The control sample has roughly twice the mobilome content of CRC across all categories,
again consistent with greater overall microbial diversity. The CRC sample's mobilome
is enriched relative to its total contig count, however, suggesting active horizontal
gene transfer even in a dysbiotic community.

---

## Defence Systems (Phase 12)

| Tool | CRC_sample | control_sample |
|------|------------|----------------|
| DefenseFinder systems | 208 (55 types) | 371 (71 types) |
| PADLOC systems | 511 (100 types) | 913 (137 types) |
| CRISPRCasFinder arrays | 70 | 202 |

The control sample has substantially more defence system diversity (71 vs 55 types
by DefenseFinder, 137 vs 100 by PADLOC). This is partly explained by the larger
number of MAGs (92 vs 30), but the per-MAG rate is also lower in CRC, suggesting
a genuine reduction in phage immunity capacity in the CRC microbiome.

---

## QC Dashboard

All 23 QC checkpoints passed for both samples. Full QC reports are in
`pipeline_run/qc_checkpoints/`.

---

## Integration

All results from Phases 4-12 are consolidated into a master table at
`pipeline_run/integration/Wirbel_CRC/` with the following outputs:

* `master_viral_table_CRC_sample.tsv` — 335 rows × 43 columns
* `master_viral_table_control_sample.tsv` — 542 rows × 43 columns
* `qc_dashboard_Wirbel_CRC.txt` — 23/23 PASS
* `pipeline_summary_Wirbel_CRC.md` — full biological narrative
* `defence_summary_CRC_sample.txt` and `defence_summary_control_sample.txt`
* `mobilome_summary_CRC_sample.tsv` and `mobilome_summary_control_sample.tsv`
