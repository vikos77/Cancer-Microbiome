# Hartwig MAG classification: human-associated vs probable contaminant.
#
# Starts from all 514 Hartwig MAGs and classifies each one in three steps:
#   1. Taxonomic screen: a MAG passes if its GTDB-Tk phylum or genus is on the allowlist.
#   2. iPHoP host ecology: corroborating evidence for every flagged MAG.
#   3. Literature review: decides each flagged MAG's final call.
#
# Output (output/mag_classification/):
#   hartwig_mag_classification.tsv  read by every analysis script
#   screening_allowlists.tsv         also read by 01_wirbel_mag_classification.R
# Run this script first.

library(tidyverse)

# Run from the repository root: Rscript scripts/00_hartwig_mag_classification.R
stopifnot("Run this script from the repository root" = dir.exists("data"))
DATA   <- file.path("data", "hartwig")
OUTDIR <- file.path("output", "mag_classification")
REFDIR <- "published_tables"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Allowlists
# -----------------------------------------------------------------------------

SCREEN_PHYLA <- c("Firmicutes", "Bacteroidota", "Fusobacteriota", "Actinomycetota",
                  "Campylobacterota", "Spirochaetota", "Synergistota", "Verrucomicrobiota")

SCREEN_GENERA <- c("Escherichia", "Klebsiella", "Haemophilus", "Haemophilus_D", "Haemophilus_A",
                   "Neisseria", "Acinetobacter", "Moraxella_A", "Citrobacter", "Citrobacter_A",
                   "Enterobacter", "Aggregatibacter", "Salmonella", "Pasteurella",
                   "Pseudenterobacter", "Bibersteinia", "Kingella_B", "Desulfovibrio",
                   "Bilophila", "Methanobrevibacter_A", "Sutterella", "Parasutterella", "Eikenella")

# iPHoP reports GTDB phylum names, so Firmicutes appears as the Bacillota phyla
IPHOP_PHYLA <- c("Bacillota", "Bacillota_A", "Bacillota_B", "Bacillota_C", "Bacteroidota",
                 "Fusobacteriota", "Actinomycetota", "Campylobacterota", "Spirochaetota",
                 "Verrucomicrobiota", "Synergistota", "Elusimicrobiota")

# Gut lineages confirmed by our literature review (Table 2) are on the iPHoP list but not the
# screen list, so their MAGs are still flagged and reviewed. iPHoP often predicts a MAG's own
# lineage as its phages' host, so leaving them off would count those hosts as environmental.
IPHOP_GENERA <- c(SCREEN_GENERA,
                  "51-20", "Scatacola_A", "Scatocola", "CAG-267", "HGM16780", "CAZU01", "MGBC133411",
                  "Duodenibacillus", "Aphodousia", "Mesosutterella",
                  "CAG-196", "Limenecus", "Scatousia", "Stercorousia", "CALUVH01", "Scatenecus", "SIG32")

allowlists <- bind_rows(
  tibble(list = "Taxonomic screen",          rank = "Phylum", taxon = SCREEN_PHYLA),
  tibble(list = "Taxonomic screen",          rank = "Genus",  taxon = SCREEN_GENERA),
  tibble(list = "iPHoP host classification", rank = "Phylum", taxon = IPHOP_PHYLA),
  tibble(list = "iPHoP host classification", rank = "Genus",  taxon = IPHOP_GENERA)
)
stopifnot(!anyDuplicated(allowlists))
write_tsv(allowlists, file.path(OUTDIR, "screening_allowlists.tsv"))

# -----------------------------------------------------------------------------
# The 514 Hartwig MAGs
# -----------------------------------------------------------------------------

gtdb <- read_csv(file.path(DATA, "gtdbtk_checkm2_hartwig.csv"),
                  show_col_types = FALSE) %>%
  rename(mag_id = bin_id) %>%
  mutate(
    phylum_raw   = str_remove(Phylum, "^p__"),
    phylum_clean = case_when(
      str_detect(phylum_raw, "^Bacillota|^Firmicutes")          ~ "Firmicutes",
      str_detect(phylum_raw, "^Bacteroidota")                   ~ "Bacteroidota",
      str_detect(phylum_raw, "^Pseudomonadota|^Proteobacteria") ~ "Proteobacteria",
      str_detect(phylum_raw, "^Fusobacteriota")                 ~ "Fusobacteriota",
      str_detect(phylum_raw, "^Actinomycetota")                 ~ "Actinomycetota",
      str_detect(phylum_raw, "^Campylobacterota")               ~ "Campylobacterota",
      str_detect(phylum_raw, "^Spirochaetota")                  ~ "Spirochaetota",
      str_detect(phylum_raw, "^Synergistota")                   ~ "Synergistota",
      TRUE                                                       ~ phylum_raw
    )
  ) %>%
  mutate(
    genome_length_mbp = round(Length / 1e6, 3),
    gc_pct            = round(`GC(%)` * 100, 2)   # source column is a fraction
  ) %>%
  select(mag_id, phylum_clean, Order, Genus, Completeness, Contamination,
         genome_length_mbp, gc_pct)

stopifnot(nrow(gtdb) == 514)

is_counts <- read_tsv(file.path(DATA, "isescan_hartwig.tsv"),
                      show_col_types = FALSE) %>%
  count(mag_id, name = "n_IS_elements")

# -----------------------------------------------------------------------------
# Step 1: taxonomic screen
# -----------------------------------------------------------------------------

# Genera with a documented history as reagent or environmental contaminants
KNOWN_CONTAMINANT_GENERA <- "^Ralstonia|^Bradyrhizobium|^Pelomonas|^Novosphingobium|^Sphingomonas|^Pseudomonas_E|^Curvibacter|^Stagnimonas"

# -----------------------------------------------------------------------------
# Step 2: iPHoP host ecology per MAG
# -----------------------------------------------------------------------------

iphop_pct_by_mag <- read_csv(file.path(DATA, "iphop_hartwig.csv"),
                             show_col_types = FALSE) %>%
  filter(mag_id %in% gtdb$mag_id) %>%
  transmute(mag_id,
            host_gtdb   = `Host genus`,
            host_genus  = str_extract(host_gtdb, "g__([^;,]+)") %>% str_remove("g__"),
            is_gut_host = str_detect(host_gtdb, paste(IPHOP_PHYLA, collapse = "|")) |
                          host_genus %in% IPHOP_GENERA) %>%
  group_by(mag_id) %>%
  summarise(n_iphop_predictions       = n(),
            pct_environmental_host    = round(100 * sum(!is_gut_host, na.rm = TRUE) / n(), 1),
            pct_human_associated_host = round(100 * sum(is_gut_host,  na.rm = TRUE) / n(), 1),
            .groups = "drop")

# -----------------------------------------------------------------------------
# Step 3: literature review of flagged MAGs (summarised in Table 1)
# -----------------------------------------------------------------------------

# Flagged genera the literature confirms as environmental lineages
ENVIRONMENTAL_GENERA <- c("Ralstonia", "Bradyrhizobium", "Pelomonas", "Novosphingobium",
                          "Sphingomonas", "Pseudomonas_E", "Curvibacter", "Stagnimonas",
                          "Brachymonas")

# Two Caulobacterales MAGs with no GTDB genus, called environmental at order level
ENVIRONMENTAL_MAGS <- c("S508C470", "S3891C470")

# The literature does not place the other two flagged MAGs (S79C1114, S3591C1867) in an
# environmental lineage, so they are retained.

# -----------------------------------------------------------------------------
# Final classification
# -----------------------------------------------------------------------------

# A flagged MAG is a probable contaminant only when the literature confirms an
# environmental lineage. Literature takes precedence over iPHoP where they disagree.
mag_classification <- gtdb %>%
  left_join(iphop_pct_by_mag, by = "mag_id") %>%
  left_join(is_counts, by = "mag_id") %>%
  mutate(
    across(c(n_iphop_predictions, n_IS_elements), ~ replace_na(., 0L)),
    taxonomic_screen = case_when(
      str_detect(Genus, regex(KNOWN_CONTAMINANT_GENERA, ignore_case = TRUE)) %in% TRUE
                                                                ~ "Flagged: known contaminant genus",
      phylum_clean %in% SCREEN_PHYLA | Genus %in% SCREEN_GENERA ~ "Passed",
      TRUE                                                      ~ "Flagged: not on allowlist"
    ),
    iphop_evidence = case_when(
      taxonomic_screen == "Passed"  ~ "Not assessed (passed taxonomic screen)",
      n_iphop_predictions == 0      ~ "No iPHoP host predictions",
      pct_environmental_host >= 75 & pct_human_associated_host == 0 ~ "Environmental host profile",
      TRUE                          ~ "Human-associated host signal"
    ),
    classification = if_else(
      taxonomic_screen != "Passed" &
        (Genus %in% ENVIRONMENTAL_GENERA | mag_id %in% ENVIRONMENTAL_MAGS),
      "Probable contaminant", "Human-associated"),
    Genus = replace_na(Genus, "Unclassified")
  ) %>%
  select(mag_id, phylum_clean, Order, Genus, Completeness, Contamination,
         taxonomic_screen, n_iphop_predictions, pct_environmental_host,
         pct_human_associated_host, iphop_evidence, classification,
         genome_length_mbp, gc_pct, n_IS_elements) %>%
  arrange(desc(classification == "Probable contaminant"), Genus, mag_id)

# -----------------------------------------------------------------------------
# Checks against the counts reported in the manuscript
# -----------------------------------------------------------------------------

print(count(mag_classification, taxonomic_screen, classification))
flagged_ha <- mag_classification %>%
  filter(taxonomic_screen != "Passed", classification == "Human-associated") %>% pull(mag_id)
stopifnot(sum(mag_classification$classification == "Human-associated") == 383,
          sum(mag_classification$classification == "Probable contaminant") == 131,
          sum(mag_classification$taxonomic_screen == "Passed") == 381,
          sum(mag_classification$taxonomic_screen == "Flagged: known contaminant genus") == 128,
          sum(mag_classification$taxonomic_screen == "Flagged: not on allowlist") == 5,
          setequal(flagged_ha, c("S79C1114", "S3591C1867")))
cat("Confirmed: 514 MAGs = 383 human-associated + 131 probable contaminants.\n")

evidence_counts <- mag_classification %>%
  filter(classification == "Probable contaminant") %>%
  count(iphop_evidence)
print(evidence_counts)
n_ev <- function(x) evidence_counts$n[evidence_counts$iphop_evidence == x]
stopifnot(n_ev("Environmental host profile") == 116,
          n_ev("No iPHoP host predictions") == 14,
          n_ev("Human-associated host signal") == 1)
cat("Confirmed iPHoP evidence for the 131 probable contaminants: 116 environmental / 14 none / 1 human-associated.\n")

ha <- filter(mag_classification, classification == "Human-associated")
pc <- filter(mag_classification, classification == "Probable contaminant")
print(mag_classification %>%
        group_by(classification) %>%
        summarise(med_len = median(genome_length_mbp), med_gc = median(gc_pct),
                  med_is = median(n_IS_elements), .groups = "drop"))
stopifnot(abs(median(ha$genome_length_mbp) - 2.0) < 0.05, abs(median(pc$genome_length_mbp) - 5.0) < 0.05,
          abs(median(ha$gc_pct) - 43.0) < 0.05,           abs(median(pc$gc_pct) - 64.0) < 0.05,
          median(ha$n_IS_elements) == 6,                  median(pc$n_IS_elements) == 12)
rho_ha <- suppressWarnings(cor.test(ha$genome_length_mbp, ha$n_IS_elements, method = "spearman"))$estimate
rho_pc <- suppressWarnings(cor.test(pc$genome_length_mbp, pc$n_IS_elements, method = "spearman"))$estimate
stopifnot(abs(rho_ha - 0.434) < 0.001, abs(rho_pc - 0.415) < 0.001)
cat("Confirmed: genome size 2.0/5.0 Mbp, GC 43.0/64.0%, IS 6/12, Spearman rho 0.434/0.415.\n")

# Optional: compare with the published Tables 1, S1 and S3 (TSV) if placed in published_tables/
ref_file <- function(name) {
  f <- file.path(REFDIR, name)
  if (!file.exists(f)) message("Skipped comparison with ", name, " (not in ", REFDIR, "/)")
  f
}

f <- ref_file("Table_1_contaminant_genera.tsv")
if (file.exists(f)) {
  table1_counts <- read_tsv(f, show_col_types = FALSE, quote = "") %>%
    transmute(genus_key = str_extract(genus, "^[A-Za-z]+(_[A-Z])?"), mags_n)
  pc_counts <- pc %>%
    count(genus_key = str_extract(Genus, "^[A-Za-z]+(_[A-Z])?"), name = "mags_n")
  stopifnot(isTRUE(all.equal(arrange(table1_counts, genus_key), arrange(pc_counts, genus_key),
                             check.attributes = FALSE)))
  cat("Confirmed: matches Table 1.\n")
}

f <- ref_file("Table_S1_hartwig_mag_classification.tsv")
if (file.exists(f)) {
  table_s1 <- read_tsv(f, show_col_types = FALSE)
  stopifnot(isTRUE(all.equal(table_s1[names(mag_classification)], mag_classification,
                             check.attributes = FALSE)))
  cat("Confirmed: matches Table S1.\n")
}

f <- ref_file("Table_S3_screening_allowlists.tsv")
if (file.exists(f)) {
  table_s3 <- read_tsv(f, show_col_types = FALSE)
  stopifnot(isTRUE(all.equal(arrange(table_s3[names(allowlists)], list, rank, taxon),
                             arrange(allowlists, list, rank, taxon), check.attributes = FALSE)))
  cat("Confirmed: matches Table S3.\n")
}

write_tsv(mag_classification, file.path(OUTDIR, "hartwig_mag_classification.tsv"))
cat("Saved: hartwig_mag_classification.tsv (", nrow(mag_classification), "rows)\n")
