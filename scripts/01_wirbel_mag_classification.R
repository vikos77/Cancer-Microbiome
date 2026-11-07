# Wirbel CRC MAG classification, using the same workflow as the Hartwig cohort.
#
# Starts from all 785 Wirbel MAGs and applies the three steps from
# 00_hartwig_mag_classification.R: taxonomic screen, iPHoP host ecology and
# literature review of flagged MAGs.
#
# Input:  screening_allowlists.tsv (run 00_hartwig_mag_classification.R first)
# Output: output/mag_classification/wirbel_mag_classification.tsv

library(tidyverse)

# Run from the repository root: Rscript scripts/01_wirbel_mag_classification.R
stopifnot("Run this script from the repository root" = dir.exists("data"))
WIRBEL_DIR <- file.path("data", "wirbel")
OUTDIR     <- file.path("output", "mag_classification")
REFDIR     <- "published_tables"

allowlists <- read_tsv(file.path(OUTDIR, "screening_allowlists.tsv"), show_col_types = FALSE)
pick <- function(l, r) allowlists %>% filter(list == l, rank == r) %>% pull(taxon)
SCREEN_PHYLA      <- pick("Taxonomic screen", "Phylum")
SCREEN_GENERA     <- pick("Taxonomic screen", "Genus")
gut_phyla_pattern <- paste(pick("iPHoP host classification", "Phylum"), collapse = "|")
GUT_GENERA        <- pick("iPHoP host classification", "Genus")

# -----------------------------------------------------------------------------
# The 785 Wirbel MAGs
# -----------------------------------------------------------------------------

strip_gtdb_prefix <- function(x) str_remove(as.character(x), "^[a-z]__")

canonicalise_phylum <- function(x) {
  x <- str_remove(as.character(x), "^p__")
  case_when(
    is.na(x) | x == ""                                      ~ "Other",
    str_detect(x, "^Bacillota|^Firmicutes")                 ~ "Firmicutes",
    str_detect(x, "^Bacteroidota")                          ~ "Bacteroidota",
    str_detect(x, "^Pseudomonadota|^Proteobacteria")        ~ "Proteobacteria",
    str_detect(x, "^Fusobacteriota")                        ~ "Fusobacteriota",
    str_detect(x, "^Actinomycetota|^Actinobacteriota")      ~ "Actinomycetota",
    str_detect(x, "^Campylobacterota")                      ~ "Campylobacterota",
    str_detect(x, "^Spirochaetota")                         ~ "Spirochaetota",
    str_detect(x, "^Verrucomicrobiota")                     ~ "Verrucomicrobiota",
    str_detect(x, "^Desulfobacterota")                      ~ "Desulfobacterota",
    str_detect(x, "^Methanobacteriota")                     ~ "Methanobacteriota",
    str_detect(x, "^Thermoplasmatota")                      ~ "Thermoplasmatota",
    str_detect(x, "^Cyanobacteriota|^Cyanobacteria")        ~ "Cyanobacteria",
    TRUE                                                     ~ x
  )
}

wirbel_tax <- read_tsv(file.path(WIRBEL_DIR, "gtdbtk_taxonomy_wirbel.tsv"),
                       show_col_types = FALSE) %>%
  transmute(
    mag_id       = gtdb_genome_id,
    phylum_clean = canonicalise_phylum(strip_gtdb_prefix(coalesce(phylum_name, phylum))),
    genus_gtdb   = strip_gtdb_prefix(coalesce(genus_name, genus)),
    order_name,
    unclassified = gtdb_ani_status == "no_result"
  )
stopifnot(nrow(wirbel_tax) == 785, !anyDuplicated(wirbel_tax$mag_id))

# -----------------------------------------------------------------------------
# Step 1: taxonomic screen
# -----------------------------------------------------------------------------

KNOWN_CONTAMINANT_GENERA <- "^Ralstonia|^Bradyrhizobium|^Pelomonas|^Novosphingobium|^Sphingomonas|^Pseudomonas_E|^Curvibacter|^Stagnimonas"

wirbel_tax <- wirbel_tax %>%
  mutate(taxonomic_screen = case_when(
    str_detect(genus_gtdb, regex(KNOWN_CONTAMINANT_GENERA, ignore_case = TRUE)) %in% TRUE
                                                                   ~ "Flagged: known contaminant genus",
    phylum_clean %in% SCREEN_PHYLA | genus_gtdb %in% SCREEN_GENERA ~ "Passed",
    TRUE                                                           ~ "Flagged: not on allowlist"
  ))

# -----------------------------------------------------------------------------
# Step 2: iPHoP host ecology per MAG
# -----------------------------------------------------------------------------

iphop_pct_by_mag <- read_tsv(file.path(WIRBEL_DIR, "iphop_wirbel.tsv"),
                             show_col_types = FALSE) %>%
  transmute(mag_id    = paste(sample_id, bin_id, sep = "_"),
            host_gtdb = `Host genus`,
            host_genus = str_extract(host_gtdb, "g__([^;,]+)") %>% str_remove("g__"),
            is_gut_host = str_detect(host_gtdb, gut_phyla_pattern) | host_genus %in% GUT_GENERA) %>%
  group_by(mag_id) %>%
  summarise(n_iphop_predictions = n(),
            pct_environmental_host    = round(100 * sum(!is_gut_host, na.rm = TRUE) / n(), 1),
            pct_human_associated_host = round(100 * sum(is_gut_host,  na.rm = TRUE) / n(), 1),
            .groups = "drop")
stopifnot(all(iphop_pct_by_mag$mag_id %in% wirbel_tax$mag_id))

# -----------------------------------------------------------------------------
# Step 3: literature review of flagged MAGs (summarised by genus in Table 2)
# -----------------------------------------------------------------------------

# Every flagged genus was confirmed as a gut lineage, so none is environmental
GUT_LINEAGE_GENERA <- c("Scatacola_A", "Scatocola", "51-20", "MGBC133411", "CAZU01",
                        "Duodenibacillus", "Mesosutterella", "Aphodousia", "CAG-196",
                        "Limenecus", "Scatousia", "Stercorousia", "Methanarcanum", "UBA4248")
ENVIRONMENTAL_GENERA <- character(0)

# MAGs with no GTDB-Tk taxonomy cannot be reviewed and are retained.

# -----------------------------------------------------------------------------
# Completeness, genome size, GC and IS counts
# -----------------------------------------------------------------------------

# CheckM2 v1.0.1 (database v2), run on the final 785 MAGs
checkm2 <- read_tsv(file.path(WIRBEL_DIR, "checkm2_wirbel.tsv"), show_col_types = FALSE)
stopifnot(nrow(checkm2) == 785, setequal(checkm2$mag_id, wirbel_tax$mag_id),
          all(checkm2$checkm2_version == "1.0.1"))
checkm2 <- select(checkm2, mag_id, Completeness, Contamination)

# Genome length and GC were computed from the 785 assembled MAG FASTA files
genome_props <- read_tsv(file.path(WIRBEL_DIR, "genome_properties_wirbel.tsv"),
                         show_col_types = FALSE) %>%
  transmute(mag_id = paste(source_sample, source_bin, sep = "_"),
            genome_length_mbp = round(genome_size_mbp, 3), gc_pct)
stopifnot(nrow(genome_props) == 785, all(wirbel_tax$mag_id %in% genome_props$mag_id))

is_counts <- read_tsv(file.path(WIRBEL_DIR, "isescan_wirbel.tsv"),
                      show_col_types = FALSE) %>%
  count(mag_id = paste(sample_id, bin_id, sep = "_"), name = "n_IS_elements")

# -----------------------------------------------------------------------------
# Final classification (same rule as Hartwig)
# -----------------------------------------------------------------------------

wirbel_full <- wirbel_tax %>%
  left_join(iphop_pct_by_mag, by = "mag_id") %>%
  left_join(checkm2, by = "mag_id") %>%
  mutate(n_iphop_predictions = replace_na(n_iphop_predictions, 0L)) %>%
  transmute(
    mag_id,
    phylum_clean = if_else(unclassified, "Unclassified", phylum_clean),
    Order = coalesce(order_name, "Unclassified"),
    Genus = coalesce(genus_gtdb, "Unclassified"),
    Completeness, Contamination,
    taxonomic_screen,
    n_iphop_predictions, pct_environmental_host, pct_human_associated_host,
    iphop_evidence = case_when(
      taxonomic_screen == "Passed" ~ "Not assessed (passed taxonomic screen)",
      n_iphop_predictions == 0     ~ "No iPHoP host predictions",
      pct_environmental_host >= 75 & pct_human_associated_host == 0 ~ "Environmental host profile",
      TRUE                         ~ "Human-associated host signal"
    ),
    classification = if_else(taxonomic_screen != "Passed" & Genus %in% ENVIRONMENTAL_GENERA,
                             "Probable contaminant", "Human-associated")
  ) %>%
  left_join(genome_props, by = "mag_id") %>%
  left_join(is_counts, by = "mag_id") %>%
  mutate(n_IS_elements = replace_na(n_IS_elements, 0L))

# -----------------------------------------------------------------------------
# Checks against the counts reported in the manuscript
# -----------------------------------------------------------------------------

flagged <- filter(wirbel_full, taxonomic_screen != "Passed")
print(count(wirbel_full, taxonomic_screen, iphop_evidence))
stopifnot(nrow(wirbel_full) == 785,
          sum(wirbel_full$taxonomic_screen == "Passed") == 753,
          sum(wirbel_full$taxonomic_screen == "Flagged: not on allowlist") == 32,
          sum(wirbel_full$taxonomic_screen == "Flagged: known contaminant genus") == 0,
          sum(flagged$Genus %in% GUT_LINEAGE_GENERA) == 30,
          sum(flagged$Genus == "Unclassified") == 2,
          all(wirbel_full$classification == "Human-associated"))
stopifnot(sum(wirbel_full$iphop_evidence == "Environmental host profile") == 0,
          sum(wirbel_full$iphop_evidence == "Human-associated host signal") == 29,
          sum(wirbel_full$iphop_evidence == "No iPHoP host predictions") == 3)
cat("Confirmed: 785 MAGs = 753 passed + 32 flagged (30 confirmed gut lineages, 2 unresolved); all human-associated.\n")

stopifnot(abs(median(wirbel_full$genome_length_mbp) - 2.41) < 0.005,
          abs(median(wirbel_full$gc_pct) - 48.7) < 0.05,
          median(wirbel_full$n_IS_elements) == 5)
cat("Confirmed: median genome 2.41 Mbp, GC 48.7%, 5 IS elements per MAG.\n")

# CheckM2 scores are reported, not used to drop MAGs
below <- filter(wirbel_full, Completeness < 90 | Contamination > 5)
cat(sprintf("CheckM2: median completeness %.2f, median contamination %.2f; %d of 785 outside >=90/<=5.\n",
            median(wirbel_full$Completeness), median(wirbel_full$Contamination), nrow(below)))
if (nrow(below) > 0) print(select(below, mag_id, Genus, Completeness, Contamination))
# These three MAGs are named in the Table S2 caption
stopifnot(setequal(below$mag_id, c("CCMD12232071ST-21-0_bin.37", "CCMD79987997ST-21-0_bin.23",
                                   "CCMD87156761ST-21-0_bin.39")))

# Optional: compare with the published Table S2 (TSV) if placed in published_tables/
f <- file.path(REFDIR, "Table_S2_wirbel_785_classification.tsv")
if (file.exists(f)) {
  table_s2 <- read_tsv(f, show_col_types = FALSE)
  stopifnot(isTRUE(all.equal(table_s2[names(wirbel_full)], wirbel_full, check.attributes = FALSE)))
  cat("Confirmed: matches Table S2.\n")
} else {
  message("Skipped comparison with Table S2 (not in ", REFDIR, "/)")
}

write_tsv(wirbel_full, file.path(OUTDIR, "wirbel_mag_classification.tsv"))
cat("Saved: wirbel_mag_classification.tsv (", nrow(wirbel_full), "rows)\n")
