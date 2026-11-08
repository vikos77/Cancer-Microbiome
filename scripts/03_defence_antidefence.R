# Defence and anti-defence systems in the Hartwig and Wirbel MAGs.
#
#   1. Prevalence of every defence and anti-defence subtype, with Fisher's exact tests of
#      Hartwig human-associated MAGs against probable contaminants and against Wirbel CRC
#      (Tables S4 and S5).
#   2. Candidate defence markers: shared Hartwig/Wirbel systems ranked by how rare they are
#      in probable contaminants and in reference genomes (Table 3).
#   3. Figure 3
#      A  Hartwig defence, human-associated vs probable contaminant
#      B  Wirbel CRC, top 20 defence systems
#      C  Top 20 defence systems of Hartwig human-associated or Wirbel CRC, side by side
#      D  Hartwig anti-defence, human-associated vs probable contaminant
#      E  Wirbel CRC, top 10 anti-defence systems
#      F  Anti-defence systems shared between Hartwig human-associated and Wirbel CRC
#
# Input:  MAG classifications from 00_hartwig_mag_classification.R and 01_wirbel_mag_classification.R
# Output: output/tables/ (Tables 3, S4, S5 and Table 3 with full statistics),
#         output/figures/figure_3/, output/figure_data/figure_3/

library(tidyverse)
library(patchwork)
library(scales)

# Run from the repository root: Rscript scripts/03_defence_antidefence.R
stopifnot("Run this script from the repository root" = dir.exists("data"))
CLASS_DIR  <- file.path("output", "mag_classification")
WIRBEL_DIR <- file.path("data", "wirbel")
REF_DIR    <- file.path("data", "reference_genomes")
TAB_DIR    <- file.path("output", "tables")
FIG_DIR    <- file.path("output", "figures", "figure_3")
DATA_DIR   <- file.path("output", "figure_data", "figure_3")
for (d in c(TAB_DIR, FIG_DIR, DATA_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

save_fig <- function(plot, name, w, h) {
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), plot, width = w, height = h)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), plot, width = w, height = h, dpi = 600)
  message("Saved figure: ", file.path(FIG_DIR, name))
}

save_data <- function(df, name) {
  write_tsv(df, file.path(DATA_DIR, paste0(name, ".tsv")))
  message("Saved data:   ", file.path(DATA_DIR, paste0(name, ".tsv")))
}

# -----------------------------------------------------------------------------
# MAG groups and defence / anti-defence calls
# -----------------------------------------------------------------------------

hartwig_class <- read_tsv(file.path(CLASS_DIR, "hartwig_mag_classification.tsv"), show_col_types = FALSE)
contaminant_ids <- hartwig_class$mag_id[hartwig_class$classification == "Probable contaminant"]

N_HARTWIG <- sum(hartwig_class$classification == "Human-associated")
N_CONTAM  <- length(contaminant_ids)
N_WIRBEL  <- nrow(read_tsv(file.path(CLASS_DIR, "wirbel_mag_classification.tsv"), show_col_types = FALSE))
stopifnot(N_HARTWIG == 383, N_CONTAM == 131, N_WIRBEL == 785)

group_of <- function(mag_ids) {
  factor(if_else(mag_ids %in% contaminant_ids, "Probable contaminant", "Human-associated"),
         levels = c("Human-associated", "Probable contaminant"))
}

df_hartwig_full <- read_tsv(file.path("data", "hartwig", "defensefinder_hartwig.tsv"),
                            show_col_types = FALSE) %>%
  filter(activity != "Antidefense") %>%
  mutate(group = group_of(mag_id))

adf_hartwig_full <- read_tsv(file.path("data", "hartwig", "antidefensefinder_hartwig.tsv"),
                             show_col_types = FALSE) %>%
  filter(activity == "Antidefense") %>%
  mutate(group = group_of(mag_id))

# Wirbel defence and anti-defence calls come from one DefenseFinder table
wirbel_df_raw <- read_tsv(file.path(WIRBEL_DIR, "defensefinder_wirbel.tsv"),
                          show_col_types = FALSE) %>%
  mutate(mag_id = paste(sample_id, bin_id, sep = "_"))
wirbel_defence_raw <- wirbel_df_raw %>% filter(activity == "Defense")
wirbel_adf_raw     <- wirbel_df_raw %>% filter(activity == "Antidefense")

# -----------------------------------------------------------------------------
# Statistics
# -----------------------------------------------------------------------------

# Fisher's exact test per subtype, Hartwig human-associated against another group
# (OR > 1 = more common in human-associated). Subtypes that neither group carries are not
# tested, and BH correction covers only the tested subtypes.
fisher_family <- function(n_ha, n_ha_total, n_other, n_other_total) {
  res <- map2_dfr(n_ha, n_other, function(a, b) {
    ft <- fisher.test(matrix(c(a, n_ha_total - a, b, n_other_total - b), nrow = 2))
    tibble(or = unname(ft$estimate), p = ft$p.value)
  })
  testable <- n_ha + n_other > 0
  res$or[!testable] <- NA_real_
  res$p[!testable]  <- NA_real_
  res$p_adj <- NA_real_
  res$p_adj[testable] <- p.adjust(res$p[testable], method = "BH")
  res
}

add_family <- function(tbl, other_col, n_other_total, suffix) {
  res <- fisher_family(tbl$n_human_associated, N_HARTWIG, tbl[[other_col]], n_other_total)
  names(res) <- paste0(c("or_", "p_", "p_adj_"), suffix)
  bind_cols(tbl, res)
}

# AntiDefenseFinder model IDs as protein names. Models that group close homologs keep every
# member, since a hit can be to any of them.
ADF_NAMES <- c(
  "Aca_alone" = "Aca", "acb2" = "Acb2", "acric5" = "AcrIC5", "acriii1" = "AcrIII-1",
  "Adnd_p0020_p0021" = "AdnD", "apyc1" = "ApyC1", "arda_ardu" = "ArdA/ArdU",
  "ardb_klca_acric11" = "ArdB/KlcA/AcrIC11", "ardc" = "ArdC", "ardk" = "ArdK", "dam" = "Dam",
  "gam" = "Gam", "hia5_hin1523_nma1821" = "Hia5/Hin1523/Nma1821", "lar" = "Lar", "mom" = "Mom",
  "NARP1" = "NARP1", "NARP2" = "NARP2", "psiab" = "PsiAB", "racc" = "RacC", "tad1" = "Tad1",
  "tad2_acriia7" = "Tad2/AcrIIA7", "ulx" = "Ulx"
)
fmt_adf <- function(x) {
  out <- if_else(x %in% names(ADF_NAMES), ADF_NAMES[x], str_replace(x, "^acriia(\\d+)$", "AcrIIA\\1"))
  stopifnot(!str_detect(out, "^[a-z]"))
  unname(out)
}

# -----------------------------------------------------------------------------
# Tables S4 and S5: prevalence of every subtype in the three groups
# -----------------------------------------------------------------------------

prevalence_table <- function(hartwig_detail, wirbel_detail) {
  hartwig_counts <- hartwig_detail %>%
    group_by(subtype, group) %>%
    summarise(n_mags = n_distinct(mag_id), .groups = "drop") %>%
    pivot_wider(id_cols = subtype, names_from = group, values_from = n_mags, values_fill = 0) %>%
    rename(n_human_associated = `Human-associated`, n_contaminant = `Probable contaminant`)

  wirbel_counts <- wirbel_detail %>%
    group_by(subtype) %>%
    summarise(n_wirbel = n_distinct(mag_id), .groups = "drop")

  # One DefenseFinder type per subtype (the most common), Hartwig first
  main_type <- function(detail) detail %>%
    count(subtype, type, sort = TRUE) %>%
    distinct(subtype, .keep_all = TRUE) %>%
    select(subtype, type)
  subtype_type <- bind_rows(main_type(hartwig_detail), main_type(wirbel_detail)) %>%
    distinct(subtype, .keep_all = TRUE)

  hartwig_counts %>%
    full_join(wirbel_counts, by = "subtype") %>%
    mutate(across(c(n_human_associated, n_contaminant, n_wirbel), ~ replace_na(., 0))) %>%
    left_join(subtype_type, by = "subtype") %>%
    mutate(
      pct_human_associated = round(n_human_associated / N_HARTWIG * 100, 2),
      pct_contaminant      = round(n_contaminant      / N_CONTAM  * 100, 2),
      pct_wirbel           = round(n_wirbel           / N_WIRBEL  * 100, 2)
    ) %>%
    select(subtype, type, n_human_associated, pct_human_associated,
           n_contaminant, pct_contaminant, n_wirbel, pct_wirbel) %>%
    add_family("n_contaminant", N_CONTAM, "contam_vs_human") %>%
    add_family("n_wirbel",      N_WIRBEL, "wirbel_vs_human") %>%
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.) | is.infinite(.), NA_real_, .))) %>%
    arrange(desc(pct_human_associated))
}

table_s4 <- prevalence_table(df_hartwig_full, wirbel_defence_raw)
table_s5 <- prevalence_table(adf_hartwig_full, wirbel_adf_raw)

stopifnot(nrow(table_s4) == n_distinct(c(df_hartwig_full$subtype, wirbel_defence_raw$subtype)),
          nrow(table_s5) == n_distinct(c(adf_hartwig_full$subtype, wirbel_adf_raw$subtype)))

# Values reported in the Results
pick <- function(tbl, s, col) tbl[[col]][tbl$subtype == s]
stopifnot(abs(pick(table_s4, "RM_Type_II", "pct_human_associated") - 45.43) < 0.05,
          abs(pick(table_s4, "RM_Type_II", "pct_wirbel") - 42.3) < 0.05,
          abs(pick(table_s4, "RM_Type_I",  "pct_contaminant") - 16.03) < 0.1,
          abs(pick(table_s5, "acriia21", "pct_human_associated") - 13.05) < 0.1,
          abs(pick(table_s5, "acriia21", "pct_wirbel") - 19.1) < 0.1,
          abs(pick(table_s5, "ardc", "pct_contaminant") - 32.06) < 0.1,
          pick(table_s5, "acriia21", "type") == "Anti_CRISPR",
          pick(table_s5, "ardc", "type") == "Anti_RM")

write_tsv(table_s4, file.path(TAB_DIR, "Table_S4_defence_subtype_prevalence.tsv"))
write_tsv(mutate(table_s5, name = fmt_adf(subtype), .after = subtype),
          file.path(TAB_DIR, "Table_S5_antidefence_subtype_prevalence.tsv"))
cat(sprintf("Table S4: %d defence subtypes, %d differ between human-associated and probable contaminants (FDR < 0.05)\n",
            nrow(table_s4), sum(table_s4$p_adj_contam_vs_human < 0.05, na.rm = TRUE)))
cat(sprintf("Table S5: %d anti-defence subtypes, %d differ (FDR < 0.05)\n",
            nrow(table_s5), sum(table_s5$p_adj_contam_vs_human < 0.05, na.rm = TRUE)))

# -----------------------------------------------------------------------------
# Table 3: candidate defence markers
# -----------------------------------------------------------------------------

# Candidates are the defence systems in the top 20 of either Hartwig human-associated or
# Wirbel CRC (the systems in Figure 3C)
hartwig_defence_clean_prev <- df_hartwig_full %>%
  filter(group == "Human-associated") %>%
  distinct(mag_id, subtype) %>%
  count(subtype, name = "n_hartwig") %>%
  mutate(pct_hartwig = n_hartwig / N_HARTWIG * 100)

wirbel_defence_prev <- wirbel_defence_raw %>%
  distinct(mag_id, subtype) %>%
  count(subtype, name = "n_wirbel") %>%
  mutate(pct_wirbel = n_wirbel / N_WIRBEL * 100)

top20_hartwig_subtypes <- hartwig_defence_clean_prev %>%
  slice_max(pct_hartwig, n = 20, with_ties = FALSE) %>% pull(subtype)
top20_wirbel_subtypes <- wirbel_defence_prev %>%
  slice_max(pct_wirbel, n = 20, with_ties = FALSE) %>% pull(subtype)
shared_defence_subtypes <- union(top20_hartwig_subtypes, top20_wirbel_subtypes)
stopifnot(length(shared_defence_subtypes) == 26)

# Reference genomes (chromosome only), pooled within each arm: the human-associated genera
# shown in Figure 6 and the probable-contaminant genera with at least 10 MAGs
HUMAN_ARM_GENERA  <- c("Prevotella", "Fusobacterium", "Bacteroides", "Streptococcus")
CONTAM_ARM_GENERA <- c("Ralstonia", "Bradyrhizobium", "Pelomonas")

ref_summary <- read_tsv(file.path(REF_DIR, "defensefinder_summary_refgenome.tsv"),
                        show_col_types = FALSE) %>%
  filter(chromosome_length_bp > 0, source_category %in% c("human", "environmental"))

ref_detail <- read_tsv(file.path(REF_DIR, "defensefinder_detail_refgenome.tsv"),
                       show_col_types = FALSE) %>%
  filter(activity != "Antidefense") %>%
  inner_join(select(ref_summary, genus, accession, source_category), by = c("genus", "accession")) %>%
  distinct(genus, accession, source_category, subtype)

ref_groups <- tribble(
  ~suffix,         ~genera,           ~source,
  "HAarm_human",    HUMAN_ARM_GENERA,  "human",
  "HAarm_nonhuman", HUMAN_ARM_GENERA,  "environmental",
  "COarm_nonhuman", CONTAM_ARM_GENERA, "environmental",
  "COarm_human",    CONTAM_ARM_GENERA, "human"
)

candidates <- df_hartwig_full %>%
  distinct(mag_id, subtype, group) %>%
  count(subtype, group) %>%
  pivot_wider(names_from = group, values_from = n, values_fill = 0L) %>%
  rename(n_human_associated = `Human-associated`, n_contaminant = `Probable contaminant`) %>%
  full_join(wirbel_defence_raw %>% distinct(mag_id, subtype) %>% count(subtype, name = "n_wirbel"),
            by = "subtype") %>%
  mutate(across(c(n_human_associated, n_contaminant, n_wirbel), ~ replace_na(., 0L))) %>%
  filter(subtype %in% shared_defence_subtypes) %>%
  mutate(pct_human_associated = round(100 * n_human_associated / N_HARTWIG, 2),
         pct_contaminant      = round(100 * n_contaminant      / N_CONTAM,  2),
         pct_wirbel           = round(100 * n_wirbel           / N_WIRBEL,  2)) %>%
  add_family("n_contaminant", N_CONTAM, "contam_vs_human") %>%
  add_family("n_wirbel",      N_WIRBEL, "wirbel_vs_human")

for (i in seq_len(nrow(ref_groups))) {
  g <- ref_groups[i, ]
  n_col <- paste0("n_ref_", g$suffix)
  n_total <- sum(ref_summary$genus %in% g$genera[[1]] & ref_summary$source_category == g$source)
  counts <- ref_detail %>%
    filter(genus %in% g$genera[[1]], source_category == g$source) %>%
    distinct(accession, subtype) %>%
    count(subtype, name = n_col)
  candidates <- candidates %>%
    left_join(counts, by = "subtype") %>%
    mutate(!!n_col := replace_na(.data[[n_col]], 0L),
           !!paste0("pct_ref_", g$suffix) := round(100 * .data[[n_col]] / n_total, 2))
}
for (i in seq_len(nrow(ref_groups))) {
  g <- ref_groups[i, ]
  n_total <- sum(ref_summary$genus %in% g$genera[[1]] & ref_summary$source_category == g$source)
  candidates <- add_family(candidates, paste0("n_ref_", g$suffix), n_total, paste0(g$suffix, "_vs_human"))
}

# Rank by lowest probable-contaminant prevalence, then lowest reference prevalence, then the
# better-supported of Hartwig/Wirbel. CAS_Cluster is DefenseFinder's unclassified CRISPR-Cas
# category, not a characterised system.
markers_full <- candidates %>%
  mutate(max_ref_pct = pmax(pct_ref_HAarm_human, pct_ref_HAarm_nonhuman,
                            pct_ref_COarm_nonhuman, pct_ref_COarm_human)) %>%
  arrange(pct_contaminant, max_ref_pct, desc(pmin(pct_human_associated, pct_wirbel))) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.) | is.infinite(.), NA_real_, .))) %>%
  filter(subtype != "CAS_Cluster") %>%
  slice_head(n = 10) %>%
  select(system = subtype, n_human_associated, n_contaminant, n_wirbel,
         pct_human_associated, pct_contaminant, pct_wirbel,
         ends_with("contam_vs_human"), ends_with("wirbel_vs_human") & !starts_with("n_"),
         starts_with("n_ref_"), starts_with("pct_ref_"), matches("arm_.*_vs_human$"))

stopifnot(nrow(markers_full) == 10,
          abs(markers_full$pct_human_associated[markers_full$system == "Wadjet_III"] - 2.36) < 0.1,
          abs(markers_full$pct_wirbel[markers_full$system == "Wadjet_III"] - 6.24) < 0.1,
          abs(markers_full$pct_ref_COarm_nonhuman[markers_full$system == "Wadjet_III"] - 0.51) < 0.1,
          abs(markers_full$pct_ref_HAarm_human[markers_full$system == "AbiQ"] - 3.22) < 0.1)

write_tsv(markers_full, file.path(TAB_DIR, "Table_3_defence_candidate_markers_full_statistics.tsv"))

# Reference-group sizes for the column headers
ref_n <- setNames(
  map_int(seq_len(nrow(ref_groups)), function(i) {
    sum(ref_summary$genus %in% ref_groups$genera[[i]] & ref_summary$source_category == ref_groups$source[i])
  }),
  ref_groups$suffix)
stopifnot(ref_n[["HAarm_human"]] == 1553, ref_n[["HAarm_nonhuman"]] == 476,
          ref_n[["COarm_nonhuman"]] == 584, ref_n[["COarm_human"]] == 8)

hdr <- function(label, n) sprintf("%s (n=%d, %%)", label, n)

table_3 <- markers_full %>%
  transmute(
    System = str_replace(system, "^CAS_Class[12]-Subtype-", "CRISPR ") %>% str_replace_all("_", "-"),
    !!hdr("Hartwig human-associated", N_HARTWIG)                              := pct_human_associated,
    !!hdr("Hartwig probable-contaminant", N_CONTAM)                           := pct_contaminant,
    !!hdr("Wirbel CRC", N_WIRBEL)                                             := pct_wirbel,
    !!hdr("Ref. genomes, human-assoc. genera, human-source", ref_n[["HAarm_human"]])  := pct_ref_HAarm_human,
    !!hdr("Ref. genomes, human-assoc. genera, non-human", ref_n[["HAarm_nonhuman"]])  := pct_ref_HAarm_nonhuman,
    !!hdr("Ref. genomes, contaminant genera, non-human", ref_n[["COarm_nonhuman"]])   := pct_ref_COarm_nonhuman,
    !!hdr("Ref. genomes, contaminant genera, human-source", ref_n[["COarm_human"]])   := pct_ref_COarm_human
  )
write_tsv(table_3, file.path(TAB_DIR, "Table_3_defence_candidate_markers.tsv"))
cat("Table 3:", paste(table_3$System, collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# Figure 3
# -----------------------------------------------------------------------------

MIRROR_COLOURS <- c(
  "Probable-contaminant-exclusive" = "#B36F00",
  "Probable-contaminant-enriched"  = "#E69F00",
  "Human-associated-enriched"      = "#56B4E9",
  "Human-associated-exclusive"     = "#0072B2",
  "Not significant"                = "#BDC3C7"
)

HW_WB_COLOURS <- c(
  "Hartwig human-associated" = "#0072B2",
  "Wirbel CRC"                = "#4D4D4D"
)

WIRBEL_BAR_COLOUR <- "#4D4D4D"

MIRROR_LABELS <- c(
  "Probable-contaminant-exclusive" = "Probable-contaminant-exclusive  (absent from all human-associated MAGs)",
  "Probable-contaminant-enriched"  = "Probable-contaminant-enriched   (present in both, FDR < 0.05, higher in probable contaminants)",
  "Human-associated-enriched"      = "Human-associated-enriched   (present in both, FDR < 0.05, higher in human-associated)",
  "Human-associated-exclusive"     = "Human-associated-exclusive  (absent from all probable-contaminant MAGs)",
  "Not significant"                = "Not significant   (present in both, FDR >= 0.05)"
)

fmt_subtype <- function(x) {
  x %>%
    str_replace("CAS_Class1-Subtype-", "CRISPR ") %>%
    str_replace("CAS_Class2-Subtype-", "CRISPR ") %>%
    str_replace("CAS_Class1-Type-",    "CRISPR (general) ") %>%
    str_replace("CAS_Class2-Type-",    "CRISPR (general) ") %>%
    str_replace("CAS_Cluster",         "CRISPR (cluster)") %>%
    str_replace("^RM_Type_",           "RM Type ")
}

SHARED_THEME <- theme_bw(base_size = 11) +
  theme(
    legend.text        = element_text(size = 8.5),
    axis.text.y        = element_text(size = 9, face = "bold"),
    axis.title         = element_text(size = 11),
    panel.grid.major.y = element_blank(),
    plot.title         = element_text(size = 12, face = "bold")
  )

# Human-associated vs probable-contaminant prevalence per subtype, with the tests from
# Tables S4/S5. Exclusive (absent from the other group) is descriptive and not gated on
# significance; enrichment is called only at FDR < 0.05.
hartwig_comparison <- function(detail, stats_table) {
  detail %>%
    group_by(group, subtype) %>%
    summarise(n_mags_with = n_distinct(mag_id), .groups = "drop") %>%
    left_join(tibble(group = c("Human-associated", "Probable contaminant"),
                     n_total = c(N_HARTWIG, N_CONTAM)), by = "group") %>%
    mutate(prevalence = n_mags_with / n_total * 100) %>%
    pivot_wider(id_cols = subtype, names_from = group,
                values_from = c(n_mags_with, prevalence), values_fill = 0) %>%
    rename(n_human_associated       = `n_mags_with_Human-associated`,
           n_probable_contaminant   = `n_mags_with_Probable contaminant`,
           pct_human_associated     = `prevalence_Human-associated`,
           pct_probable_contaminant = `prevalence_Probable contaminant`) %>%
    left_join(select(stats_table, subtype, fisher_p = p_contam_vs_human, fdr = p_adj_contam_vs_human),
              by = "subtype") %>%
    mutate(
      pct_difference = pct_human_associated - pct_probable_contaminant,
      status = case_when(
        pct_human_associated     == 0 ~ "Probable-contaminant-exclusive",
        pct_probable_contaminant == 0 ~ "Human-associated-exclusive",
        fdr >= 0.05                   ~ "Not significant",
        pct_difference > 0            ~ "Human-associated-enriched",
        TRUE                          ~ "Probable-contaminant-enriched"
      )
    )
}

mirror_long <- function(plot_data, label_fun) {
  plot_data %>%
    mutate(display_name = factor(label_fun(subtype), levels = rev(label_fun(subtype)))) %>%
    pivot_longer(cols = c(pct_human_associated, pct_probable_contaminant),
                 names_to = "grp", values_to = "prevalence") %>%
    mutate(
      grp      = if_else(grp == "pct_human_associated", "Human-associated", "Probable contaminant"),
      x_val    = if_else(grp == "Probable contaminant", -prevalence, prevalence),
      fill_cat = case_when(
        status != "Human-associated-enriched" & status != "Probable-contaminant-enriched" ~ status,
        grp == "Probable contaminant" ~ "Probable-contaminant-enriched",
        TRUE                          ~ "Human-associated-enriched"
      ),
      fill_cat = factor(fill_cat, levels = names(MIRROR_COLOURS))
    )
}

# Panel A: top 20 per group, plus any exclusive system at >= 3% in its group
sys_wide <- hartwig_comparison(df_hartwig_full, table_s4)

top_ha      <- sys_wide %>% slice_max(pct_human_associated,     n = 20) %>% pull(subtype)
top_contam  <- sys_wide %>% slice_max(pct_probable_contaminant, n = 20) %>% pull(subtype)
excl_ha     <- sys_wide %>% filter(pct_probable_contaminant == 0, pct_human_associated >= 3) %>% pull(subtype)
excl_contam <- sys_wide %>% filter(pct_human_associated == 0, pct_probable_contaminant >= 3) %>% pull(subtype)

sys_plot <- sys_wide %>%
  filter(subtype %in% unique(c(top_ha, top_contam, excl_ha, excl_contam))) %>%
  arrange(pct_difference)

n_sys_A <- nrow(sys_plot)
max_x_A <- max(sys_plot$pct_human_associated, sys_plot$pct_probable_contaminant)
sys_mirror <- mirror_long(sys_plot, fmt_subtype)

panel_A <- ggplot(sys_mirror, aes(x = x_val, y = display_name, fill = fill_cat)) +
  geom_col(width = 0.70, colour = NA) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.55) +
  geom_text(
    data = sys_mirror %>% filter(abs(x_val) > 3),
    aes(label = sprintf("%.0f%%", abs(x_val)), hjust = if_else(x_val < 0, 1.18, -0.18)),
    size = 2.65, colour = "grey20"
  ) +
  annotate("text", x = -(max_x_A / 2), y = n_sys_A + 1.1,
           label = sprintf("<- Probable contaminants  (n = %d)", N_CONTAM),
           hjust = 0.5, size = 3.0, colour = "#B36F00", fontface = "bold") +
  annotate("text", x =  (max_x_A / 2), y = n_sys_A + 1.1,
           label = sprintf("Human-associated  (n = %d) ->", N_HARTWIG),
           hjust = 0.5, size = 3.0, colour = "#0072B2", fontface = "bold") +
  scale_fill_manual(values = MIRROR_COLOURS, name = NULL, labels = MIRROR_LABELS) +
  scale_x_continuous(breaks = seq(-50, 50, 10), labels = function(x) paste0(abs(x), "%"),
                     limits = c(-(max_x_A + 10), max_x_A + 10), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.5, n_sys_A + 1.8), clip = "off") +
  labs(title = "Hartwig defence system composition: human-associated vs probable-contaminant MAGs",
       x = "% of MAGs carrying system", y = NULL) +
  SHARED_THEME +
  theme(legend.position = "bottom", legend.direction = "vertical",
        legend.key.size = unit(0.38, "cm"),
        plot.margin = margin(t = 25, r = 10, b = 5, l = 5))

save_data(sys_wide %>% arrange(desc(pct_human_associated)), "Figure_3A_hartwig_defence_composition_data")

# Panel B
top20_wirbel_defence <- wirbel_defence_prev %>%
  slice_max(n_wirbel, n = 20, with_ties = FALSE) %>%
  arrange(n_wirbel) %>%
  mutate(subtype_short = factor(fmt_subtype(subtype), levels = fmt_subtype(subtype)))

panel_B <- ggplot(top20_wirbel_defence, aes(x = pct_wirbel, y = subtype_short)) +
  geom_col(width = 0.7, fill = WIRBEL_BAR_COLOUR) +
  geom_text(aes(label = sprintf("%.1f%%  (n=%d)", pct_wirbel, n_wirbel)),
            hjust = -0.1, size = 2.8, colour = "grey20") +
  scale_x_continuous(limits = c(0, max(top20_wirbel_defence$pct_wirbel) * 1.3),
                     expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  labs(title = "Wirbel CRC: top 20 defence systems", x = "% of MAGs with this system", y = NULL) +
  SHARED_THEME

save_data(wirbel_defence_prev %>% arrange(desc(n_wirbel)), "Figure_3B_wirbel_defence_composition_data")

# Panel C
shared_defence_df <- full_join(hartwig_defence_clean_prev, wirbel_defence_prev, by = "subtype") %>%
  filter(subtype %in% shared_defence_subtypes) %>%
  mutate(across(c(n_hartwig, pct_hartwig, n_wirbel, pct_wirbel), ~replace_na(.x, 0)),
         subtype_short = fmt_subtype(subtype), combined = n_hartwig + n_wirbel) %>%
  arrange(desc(combined))

N_SHARED_C <- nrow(shared_defence_df)

shared_mirror_C <- shared_defence_df %>%
  mutate(display_name = factor(subtype_short, levels = rev(subtype_short))) %>%
  pivot_longer(cols = c(pct_hartwig, pct_wirbel), names_to = "grp", values_to = "prevalence") %>%
  mutate(
    grp = case_when(grp == "pct_hartwig" ~ "Hartwig human-associated", TRUE ~ "Wirbel CRC"),
    x_val = if_else(grp == "Wirbel CRC", prevalence, -prevalence),
    n_label = if_else(grp == "Hartwig human-associated", n_hartwig, n_wirbel)
  )

max_x_C  <- max(abs(shared_mirror_C$x_val))
n_rows_C <- n_distinct(shared_mirror_C$display_name)

panel_C <- ggplot(shared_mirror_C, aes(x = x_val, y = display_name, fill = grp)) +
  geom_col(width = 0.72, colour = NA) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.55) +
  geom_text(aes(label = sprintf("%.1f%%  (n=%d)", abs(x_val), n_label),
                hjust = if_else(x_val < 0, 1.05, -0.05)), size = 2.5, colour = "grey20") +
  annotate("text", x = -(max_x_C + 18) / 2, y = n_rows_C + 2.2,
           label = sprintf("<- Hartwig human-associated  (n = %d)", N_HARTWIG),
           hjust = 0.5, size = 3.0, colour = "#0072B2", fontface = "bold") +
  annotate("text", x =  (max_x_C + 18) / 2, y = n_rows_C + 2.2,
           label = sprintf("Wirbel CRC  (n = %d) ->", N_WIRBEL),
           hjust = 0.5, size = 3.0, colour = "#4D4D4D", fontface = "bold") +
  scale_fill_manual(values = HW_WB_COLOURS, name = NULL) +
  scale_x_continuous(labels = function(x) paste0(abs(x), "%"),
                     limits = c(-(max_x_C + 18), max_x_C + 18), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.5, n_rows_C + 3), clip = "off") +
  labs(title = sprintf("Defence systems: Hartwig human-associated vs Wirbel CRC (n = %d, top-20-per-dataset union)", N_SHARED_C),
       x = "% of MAGs with this system", y = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 8, face = "bold"),
        panel.grid.major.y = element_blank(), plot.title = element_text(size = 12, face = "bold"))

save_data(shared_defence_df, "Figure_3C_shared_defence_composition_data")

# Panel D: top 10 per group, plus any exclusive system at >= 3% in its group.
# Zero-prevalence subtypes are dropped before ranking, as most of the small anti-defence
# catalogue ties at 0% on one side.
adf_wide <- hartwig_comparison(adf_hartwig_full, table_s5)

top_ha_D      <- adf_wide %>% filter(pct_human_associated > 0) %>%
  slice_max(pct_human_associated, n = 10, with_ties = FALSE) %>% pull(subtype)
top_contam_D  <- adf_wide %>% filter(pct_probable_contaminant > 0) %>%
  slice_max(pct_probable_contaminant, n = 10, with_ties = FALSE) %>% pull(subtype)
excl_ha_D     <- adf_wide %>% filter(pct_probable_contaminant == 0, pct_human_associated >= 3) %>% pull(subtype)
excl_contam_D <- adf_wide %>% filter(pct_human_associated == 0, pct_probable_contaminant >= 3) %>% pull(subtype)

adf_plot <- adf_wide %>%
  filter(subtype %in% unique(c(top_ha_D, top_contam_D, excl_ha_D, excl_contam_D))) %>%
  arrange(pct_difference)

n_sys_D <- nrow(adf_plot)
max_x_D <- max(adf_plot$pct_human_associated, adf_plot$pct_probable_contaminant, na.rm = TRUE)
adf_mirror_D <- mirror_long(adf_plot, fmt_adf)

panel_D <- ggplot(adf_mirror_D, aes(x = x_val, y = display_name, fill = fill_cat)) +
  geom_col(width = 0.70, colour = NA) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.55) +
  geom_text(data = adf_mirror_D %>% filter(x_val != 0),
            aes(label = sprintf("%.0f%%", abs(x_val)), hjust = if_else(x_val < 0, 1.18, -0.18)),
            size = 2.65, colour = "grey20") +
  annotate("text", x = -(max_x_D + 8), y = n_sys_D + 1.1,
           label = sprintf("<- Probable contaminants  (n = %d)", N_CONTAM),
           hjust = 0, size = 2.7, colour = "#B36F00", fontface = "bold") +
  annotate("text", x =  (max_x_D + 8), y = n_sys_D + 1.1,
           label = sprintf("Human-associated  (n = %d) ->", N_HARTWIG),
           hjust = 1, size = 2.7, colour = "#0072B2", fontface = "bold") +
  scale_fill_manual(values = MIRROR_COLOURS, name = NULL, labels = MIRROR_LABELS) +
  scale_x_continuous(breaks = seq(-50, 50, 10), labels = function(x) paste0(abs(x), "%"),
                     limits = c(-(max_x_D + 10), max_x_D + 10), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.5, n_sys_D + 1.8), clip = "off") +
  labs(title = "Hartwig anti-defence composition:\nhuman-associated vs probable contaminant",
       x = "% of MAGs carrying system", y = NULL) +
  SHARED_THEME +
  theme(legend.position = "bottom", legend.direction = "vertical",
        legend.key.size = unit(0.38, "cm"), axis.text.y = element_text(size = 9.5, face = "bold"),
        plot.margin = margin(t = 25, r = 10, b = 5, l = 5))

save_data(adf_wide %>% arrange(desc(pct_human_associated)), "Figure_3D_hartwig_adf_composition_data")

# Panel E
wirbel_adf_prev <- wirbel_adf_raw %>%
  distinct(mag_id, subtype) %>%
  count(subtype, name = "n_wirbel") %>%
  mutate(pct_wirbel = n_wirbel / N_WIRBEL * 100)

N_TOPN_E <- min(10, nrow(wirbel_adf_prev))

top_wirbel_adf <- wirbel_adf_prev %>%
  slice_max(n_wirbel, n = N_TOPN_E, with_ties = FALSE) %>%
  arrange(n_wirbel) %>%
  mutate(display_name = factor(fmt_adf(subtype), levels = fmt_adf(subtype)))

panel_E <- ggplot(top_wirbel_adf, aes(x = pct_wirbel, y = display_name)) +
  geom_col(width = 0.7, fill = WIRBEL_BAR_COLOUR) +
  geom_text(aes(label = sprintf("%.1f%%  (n=%d)", pct_wirbel, n_wirbel)),
            hjust = -0.1, size = 2.8, colour = "grey20") +
  scale_x_continuous(limits = c(0, max(top_wirbel_adf$pct_wirbel) * 1.3),
                     expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  labs(title = sprintf("Wirbel CRC: top %d anti-defence systems", N_TOPN_E),
       x = "% of MAGs with this system", y = NULL) +
  SHARED_THEME

save_data(wirbel_adf_prev %>% arrange(desc(n_wirbel)), "Figure_3E_wirbel_adf_composition_data")

# Panel F: every anti-defence subtype found in both datasets
adf_hartwig_clean_prev <- adf_hartwig_full %>%
  filter(group == "Human-associated") %>%
  distinct(mag_id, subtype) %>%
  count(subtype, name = "n_hartwig") %>%
  mutate(pct_hartwig = n_hartwig / N_HARTWIG * 100)

shared_adf_df <- inner_join(adf_hartwig_clean_prev, wirbel_adf_prev, by = "subtype") %>%
  mutate(combined = n_hartwig + n_wirbel) %>%
  arrange(desc(combined))

N_SHARED_F <- nrow(shared_adf_df)

adf_mirror_F <- shared_adf_df %>%
  mutate(display_name = factor(fmt_adf(subtype), levels = rev(fmt_adf(subtype)))) %>%
  pivot_longer(cols = c(pct_hartwig, pct_wirbel), names_to = "grp", values_to = "prevalence") %>%
  mutate(
    grp = case_when(grp == "pct_hartwig" ~ "Hartwig human-associated", TRUE ~ "Wirbel CRC"),
    x_val = if_else(grp == "Wirbel CRC", prevalence, -prevalence),
    n_label = if_else(grp == "Hartwig human-associated", n_hartwig, n_wirbel)
  )

max_x_F  <- max(abs(adf_mirror_F$x_val))
n_rows_F <- n_distinct(adf_mirror_F$display_name)

panel_F <- ggplot(adf_mirror_F, aes(x = x_val, y = display_name, fill = grp)) +
  geom_col(width = 0.72, colour = NA) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.55) +
  geom_text(aes(label = sprintf("%.1f%%  (n=%d)", abs(x_val), n_label),
                hjust = if_else(x_val < 0, 1.05, -0.05)), size = 2.5, colour = "grey20") +
  annotate("text", x = -(max_x_F + 18) / 2, y = n_rows_F + 1.4,
           label = sprintf("<- Hartwig human-associated  (n = %d)", N_HARTWIG),
           hjust = 0.5, size = 3.0, colour = "#0072B2", fontface = "bold") +
  annotate("text", x =  (max_x_F + 18) / 2, y = n_rows_F + 1.4,
           label = sprintf("Wirbel CRC  (n = %d) ->", N_WIRBEL),
           hjust = 0.5, size = 3.0, colour = "#4D4D4D", fontface = "bold") +
  scale_fill_manual(values = HW_WB_COLOURS, name = NULL) +
  scale_x_continuous(labels = function(x) paste0(abs(x), "%"),
                     limits = c(-(max_x_F + 18), max_x_F + 18), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.5, n_rows_F + 2), clip = "off") +
  labs(title = sprintf("Shared anti-defence: Hartwig vs Wirbel CRC (n = %d)", N_SHARED_F),
       x = "% of MAGs with this system", y = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 9, face = "bold"),
        panel.grid.major.y = element_blank(), plot.title = element_text(size = 12, face = "bold"))

save_data(shared_adf_df, "Figure_3F_shared_adf_composition_data")

# Assemble and save
save_fig(panel_A, "Figure_3A_hartwig_defence_composition", w = 10, h = 11)
save_fig(panel_B, "Figure_3B_wirbel_defence_composition",  w = 8,  h = 8.5)
save_fig(panel_C, "Figure_3C_shared_defence_composition",  w = 9,  h = 10)
save_fig(panel_D, "Figure_3D_hartwig_adf_composition",     w = 10, h = 8)
save_fig(panel_E, "Figure_3E_wirbel_adf_composition",      w = 8,  h = 8.5)
save_fig(panel_F, "Figure_3F_shared_adf_composition",      w = 9,  h = 6.5)

figure_3_top    <- (panel_A | panel_B) + plot_layout(widths = c(1, 0.8))
figure_3_bottom <- (panel_D | panel_E | panel_F) + plot_layout(widths = c(0.9, 0.85, 1))

# Tag styling must be applied with & theme(); plot_annotation(theme = ) ignores plot.tag.
figure_3 <- (figure_3_top / panel_C / figure_3_bottom) +
  plot_layout(heights = c(1.35, 1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))

save_fig(figure_3, "Figure_3", w = 20, h = 26)

cat(sprintf("Panel A: %d defence subtypes shown\n", n_sys_A))
cat(sprintf("Panel B: top 20 of %d Wirbel defence subtypes\n", nrow(wirbel_defence_prev)))
cat(sprintf("Panel C: %d shared defence subtypes\n", N_SHARED_C))
cat(sprintf("Panel D: %d anti-defence subtypes shown\n", n_sys_D))
cat(sprintf("Panel E: top %d of %d Wirbel anti-defence subtypes\n", N_TOPN_E, nrow(wirbel_adf_prev)))
cat(sprintf("Panel F: %d shared anti-defence subtypes\n", N_SHARED_F))
