# Artefact or Ecology? Leveraging bacteria–phage defence and mobilome signatures to distinguish human-associated MAGs from probable contaminants in cancer metagenomes

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24060765.svg)](https://doi.org/10.5281/zenodo.24060765)

Analysis code for our study of bacterial defence systems, anti-defence systems, insertion
sequences and phages in metagenome-assembled genomes (MAGs) recovered from metastatic tumour
biopsies of the Hartwig Medical Foundation cohort (Battaglia et al., 2024), with an independent
colorectal cancer (CRC) stool cohort (Zeller et al., 2014; Wirbel et al., 2019) and complete RefSeq reference genomes
for comparison.

The scripts start from the compiled outputs of the annotation tools (in `data/`) and reproduce
Figures 2 to 6, S1 and S2, Tables 3 and S4 to S6, and the MAG classifications behind Tables S1
to S3. Figure 1 is a schematic of the workflow and is not generated from data.

## Repository layout

```
data/
  hartwig/            514 Hartwig MAGs: GTDB-Tk/CheckM2, DefenseFinder, AntiDefenseFinder,
                      ISEScan, geNomad, CheckV, iPHoP and Pharokka outputs
  wirbel/             785 Wirbel CRC MAGs: the same tools, plus genome size and GC
  reference_genomes/  DefenseFinder and ISEScan results for complete RefSeq genomes
                      (chromosomes only), with each genome's NCBI BioSample source
scripts/              analysis scripts, run in the order below
output/               created when the scripts run (not tracked)
```

## Running the analysis

Run every script from the repository root.

```sh
Rscript scripts/00_hartwig_mag_classification.R
Rscript scripts/01_wirbel_mag_classification.R
Rscript scripts/02_genus_composition.R
Rscript scripts/03_defence_antidefence.R
Rscript scripts/04_insertion_sequences.R
Rscript scripts/05_mobile_genetic_elements.R
Rscript scripts/06_reference_genome_validation.R
Rscript scripts/07_phage_host_prediction.R
```

Scripts 00 and 01 classify the MAGs and must run first. Scripts 02 to 07 read their output
and can run in any order. The scripts check key results against the values reported in the
manuscript and stop if any of them differ.

| Script | Produces |
|---|---|
| `00_hartwig_mag_classification.R` | Hartwig MAG classification, human-associated or probable contaminant (Table S1); screening allowlists (Table S3) |
| `01_wirbel_mag_classification.R` | Wirbel CRC MAG classification (Table S2) |
| `02_genus_composition.R` | Figure 2 |
| `03_defence_antidefence.R` | Figure 3, Table 3, Tables S4 and S5 |
| `04_insertion_sequences.R` | Figure 4, Table S6 |
| `05_mobile_genetic_elements.R` | Figure 5 |
| `06_reference_genome_validation.R` | Figure 6, Figure S2 |
| `07_phage_host_prediction.R` | Figure S1 |

Figures are written to `output/figures/` (PDF and 600 dpi PNG), the data behind each panel to
`output/figure_data/`, tables to `output/tables/` and the MAG classifications to
`output/mag_classification/`. Script 06 also draws the reference-genome comparison for the
probable-contaminant genera (`output/figures/figure_S3/`), which is not part of the paper.

Scripts 00 and 01 can also compare their classifications with the published Tables 1, S1, S2
and S3. To run these checks, save those tables as tab-separated files in a `published_tables/`
folder at the repository root, named as in the scripts. Without them the comparison is skipped.

### Requirements

R 4.5.3 with tidyverse 2.0.0, ggplot2 4.0.2, patchwork 1.3.2 and scales 1.4.0.

```r
install.packages(c("tidyverse", "patchwork", "scales"))
```

Running all scripts takes a few minutes on a laptop.

## Data sources

- **Hartwig MAGs**: the 514 MAGs and their GTDB-Tk v2.3.0 (R214) taxonomy and CheckM2 quality
  scores from Battaglia et al. (2024), Zenodo
  [10.5281/zenodo.10777510](https://doi.org/10.5281/zenodo.10777510).
- **Wirbel CRC MAGs**: 785 MAGs we assembled from 55 CRC stool metagenomes in ENA projects
  PRJEB6070 and PRJEB27928 (Zeller et al., 2014; Wirbel et al., 2019).
- **Reference genomes**: complete RefSeq genomes for the genera compared in Figure 6 and
  Figure S2, downloaded from NCBI, with source (human or non-human) assigned from BioSample
  metadata.

The annotation tools and versions are listed in the manuscript Methods. The tables in `data/`
are those tools' outputs, compiled across MAGs.

## Citation

If you use this code, please cite the paper (citation to follow) and this repository,
archived on Zenodo at [10.5281/zenodo.24060765](https://doi.org/10.5281/zenodo.24060765)
(see `CITATION.cff`).

Battaglia TW, Mimpen IL, Traets JJH, et al. A pan-cancer analysis of the microbiome in
metastatic cancer. *Cell* 2024;187(9):2324-2335.e19. https://doi.org/10.1016/j.cell.2024.03.021

Wirbel J, Pyl PT, Kartal E, et al. Meta-analysis of fecal metagenomes reveals global microbial
signatures that are specific for colorectal cancer. *Nat Med* 2019;25(4):679-689.
https://doi.org/10.1038/s41591-019-0406-6

Zeller G, Tap J, Voigt AY, et al. Potential of fecal microbiota for early-stage detection of
colorectal cancer. *Mol Syst Biol* 2014;10(11):766. https://doi.org/10.15252/msb.20145645

## Licence

The code is released under the MIT licence (see `LICENSE`). The Hartwig MAG tables in
`data/hartwig/` are derived from the Battaglia et al. (2024) Zenodo deposit and remain under its
CC BY 4.0 licence.
