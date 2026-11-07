# Insertion sequence (IS) families in the Hartwig and Wirbel MAGs.
#
#   1. Prevalence of every IS family, with Fisher's exact tests of Hartwig human-associated
#      MAGs against probable contaminants and against Wirbel CRC (Table S6).
#   2. Figure 4
#      A  Hartwig IS families, human-associated vs probable contaminant
#      B  Wirbel CRC, top 20 IS families
#      C  IS families shared between Hartwig human-associated and Wirbel CRC
#
# Input:  MAG classifications from 00_hartwig_mag_classification.R and 01_wirbel_mag_classification.R
# Output: output/tables/ (Table S6),
#         output/figures/figure_4/, output/figure_data/figure_4/

library(tidyverse)
library(patchwork)
library(scales)

# Run from the repository root: Rscript scripts/04_insertion_sequences.R
stopifnot("Run this script from the repository root" = dir.exists("data"))
CLASS_DIR <- file.path("output", "mag_classification")
TAB_DIR   <- file.path("output", "tables")
FIG_DIR   <- file.path("output", "figures", "figure_4")
DATA_DIR  <- file.path("output", "figure_data", "figure_4")
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
# MAG groups and IS calls
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

is_hartwig_full <- read_tsv(file.path("data", "hartwig", "isescan_hartwig.tsv"),
                            show_col_types = FALSE) %>%
  mutate(group = group_of(mag_id))

is_wirbel_raw <- read_tsv(file.path("data", "wirbel", "isescan_wirbel.tsv"),
                          show_col_types = FALSE) %>%
  mutate(mag_id = paste(sample_id, bin_id, sep = "_"))

# -----------------------------------------------------------------------------
# Table S6: prevalence of every IS family in the three groups
# -----------------------------------------------------------------------------

# Fisher's exact test per family, Hartwig human-associated against another group
# (OR > 1 = more common in human-associated). Families that neither group carries are not
# tested, and BH correction covers only the tested families.
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

table_s6 <- is_hartwig_full %>%
  distinct(mag_id, family, group) %>%
  count(family, group, name = "n_mags") %>%
  pivot_wider(id_cols = family, names_from = group, values_from = n_mags, values_fill = 0) %>%
  rename(n_human_associated = `Human-associated`, n_contaminant = `Probable contaminant`) %>%
  full_join(is_wirbel_raw %>% distinct(mag_id, family) %>% count(family, name = "n_wirbel"),
            by = "family") %>%
  mutate(across(c(n_human_associated, n_contaminant, n_wirbel), ~ replace_na(., 0L))) %>%
  mutate(
    pct_human_associated = round(n_human_associated / N_HARTWIG * 100, 2),
    pct_contaminant      = round(n_contaminant      / N_CONTAM  * 100, 2),
    pct_wirbel           = round(n_wirbel           / N_WIRBEL  * 100, 2)
  ) %>%
  select(family, n_human_associated, pct_human_associated,
         n_contaminant, pct_contaminant, n_wirbel, pct_wirbel) %>%
  add_family("n_contaminant", N_CONTAM, "contam_vs_human") %>%
  add_family("n_wirbel",      N_WIRBEL, "wirbel_vs_human") %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.) | is.infinite(.), NA_real_, .))) %>%
  arrange(desc(pct_human_associated))

# 19 of the 25 families differ between human-associated and probable contaminants (Results)
stopifnot(nrow(table_s6) == n_distinct(c(is_hartwig_full$family, is_wirbel_raw$family)),
          nrow(table_s6) == 25,
          sum(table_s6$p_adj_contam_vs_human < 0.05) == 19)

write_tsv(table_s6, file.path(TAB_DIR, "Table_S6_IS_family_prevalence.tsv"))
cat(sprintf("Table S6: %d IS families, %d differ between human-associated and probable contaminants (FDR < 0.05)\n",
            nrow(table_s6), sum(table_s6$p_adj_contam_vs_human < 0.05)))

# -----------------------------------------------------------------------------
# Figure 4
# -----------------------------------------------------------------------------

IS_MIRROR_COLOURS_5 <- c(
  "Probable-contaminant-exclusive" = "#B36F00",
  "Probable-contaminant-enriched"  = "#E69F00",
  "Not significant (FDR >= 0.05)"  = "#BDC3C7",
  "Human-associated-enriched"      = "#56B4E9",
  "Human-associated-exclusive"     = "#0072B2"
)

HW_WB_COLOURS <- c(
  "Hartwig human-associated" = "#0072B2",
  "Wirbel CRC"                = "#4D4D4D"
)

WIRBEL_BAR_COLOUR <- "#4D4D4D"

SHARED_THEME <- theme_bw(base_size = 11) +
  theme(
    legend.text        = element_text(size = 8.5),
    axis.text.y        = element_text(size = 9, face = "bold"),
    axis.title         = element_text(size = 11),
    panel.grid.major.y = element_blank(),
    plot.title         = element_text(size = 12, face = "bold")
  )

# Panel A: every IS family, with the tests from Table S6
is_family_diff <- is_hartwig_full %>%
  group_by(group, family) %>%
  summarise(n_mags_with = n_distinct(mag_id), .groups = "drop") %>%
  pivot_wider(id_cols = family, names_from = group, values_from = n_mags_with, values_fill = 0) %>%
  rename(n_human_associated = `Human-associated`, n_probable_contaminant = `Probable contaminant`) %>%
  mutate(
    pct_human_associated     = 100 * n_human_associated     / N_HARTWIG,
    pct_probable_contaminant = 100 * n_probable_contaminant / N_CONTAM,
    pct_difference           = pct_human_associated - pct_probable_contaminant
  ) %>%
  left_join(select(table_s6, family, fisher_p = p_contam_vs_human, fdr = p_adj_contam_vs_human),
            by = "family") %>%
  arrange(pct_difference) %>%
  mutate(family = factor(family, levels = family))

n_is_A <- nrow(is_family_diff)

# Every colour category, exclusive included, requires FDR < 0.05
is_mirror_long <- is_family_diff %>%
  mutate(
    status = case_when(
      fdr >= 0.05                   ~ "Not significant (FDR >= 0.05)",
      pct_probable_contaminant == 0 ~ "Human-associated-exclusive",
      pct_human_associated == 0     ~ "Probable-contaminant-exclusive",
      TRUE                          ~ "enriched"
    ),
    display_name = factor(family, levels = rev(family))
  ) %>%
  pivot_longer(cols = c(pct_human_associated, pct_probable_contaminant),
               names_to = "grp", values_to = "prevalence") %>%
  mutate(
    grp      = if_else(grp == "pct_human_associated", "Human-associated", "Probable contaminant"),
    x_val    = if_else(grp == "Probable contaminant", -prevalence, prevalence),
    fill_cat = case_when(
      status != "enriched"          ~ status,
      grp == "Probable contaminant" ~ "Probable-contaminant-enriched",
      TRUE                          ~ "Human-associated-enriched"
    ),
    fill_cat = factor(fill_cat, levels = names(IS_MIRROR_COLOURS_5))
  )

max_x_A <- max(is_family_diff$pct_human_associated, is_family_diff$pct_probable_contaminant)

panel_A <- ggplot(is_mirror_long, aes(x = x_val, y = display_name, fill = fill_cat)) +
  geom_col(width = 0.72, colour = NA) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.55) +
  geom_text(
    data = is_mirror_long %>% filter(x_val != 0),
    aes(label = sprintf("%.0f%%", abs(x_val)), hjust = if_else(x_val < 0, 1.18, -0.18)),
    size = 2.65, colour = "grey20"
  ) +
  annotate("text", x = -(max_x_A / 2), y = n_is_A + 1.1,
           label = sprintf("<- Probable contaminants  (n = %d)", N_CONTAM),
           hjust = 0.5, size = 3.0, colour = "#B36F00", fontface = "bold") +
  annotate("text", x =  (max_x_A / 2), y = n_is_A + 1.1,
           label = sprintf("Human-associated  (n = %d) ->", N_HARTWIG),
           hjust = 0.5, size = 3.0, colour = "#0072B2", fontface = "bold") +
  scale_fill_manual(
    values = IS_MIRROR_COLOURS_5, name = NULL,
    labels = c(
      "Probable-contaminant-exclusive" = "Probable-contaminant-exclusive  (0% in human-associated MAGs)",
      "Probable-contaminant-enriched"  = "Probable-contaminant-enriched   (FDR < 0.05)",
      "Not significant (FDR >= 0.05)"  = "Not significant  (FDR >= 0.05)",
      "Human-associated-enriched"      = "Human-associated-enriched   (FDR < 0.05)",
      "Human-associated-exclusive"     = "Human-associated-exclusive  (0% in probable-contaminant MAGs)"
    )
  ) +
  scale_x_continuous(breaks = seq(-80, 80, 20), labels = function(x) paste0(abs(x), "%"),
                     limits = c(-(max_x_A + 16), max_x_A + 16), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.5, n_is_A + 1.8), clip = "off") +
  labs(title = "Hartwig IS family composition: human-associated vs probable-contaminant MAGs",
       x = "% of MAGs carrying IS family", y = NULL) +
  SHARED_THEME +
  theme(legend.position = "bottom", legend.direction = "vertical",
        legend.key.size = unit(0.38, "cm"),
        plot.margin = margin(t = 25, r = 10, b = 5, l = 5))

save_data(is_family_diff %>% arrange(desc(pct_human_associated)), "Figure_4A_hartwig_IS_composition_data")

# Panel B
wirbel_is_prev <- is_wirbel_raw %>%
  distinct(mag_id, family) %>%
  count(family, name = "n_wirbel") %>%
  mutate(pct_wirbel = n_wirbel / N_WIRBEL * 100)

top20_wirbel_is <- wirbel_is_prev %>%
  slice_max(n_wirbel, n = 20, with_ties = FALSE) %>%
  arrange(n_wirbel) %>%
  mutate(family = factor(family, levels = family))

panel_B <- ggplot(top20_wirbel_is, aes(x = pct_wirbel, y = family)) +
  geom_col(width = 0.7, fill = WIRBEL_BAR_COLOUR) +
  geom_text(aes(label = sprintf("%.1f%%  (n=%d)", pct_wirbel, n_wirbel)),
            hjust = -0.1, size = 2.8, colour = "grey20") +
  scale_x_continuous(limits = c(0, max(top20_wirbel_is$pct_wirbel) * 1.3),
                     expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  labs(title = "Wirbel CRC: top 20 IS families", x = "% of MAGs with this IS family", y = NULL) +
  SHARED_THEME

save_data(wirbel_is_prev %>% arrange(desc(n_wirbel)), "Figure_4B_wirbel_IS_composition_data")

# Panel C: every IS family found in both datasets
hartwig_is_ha_prev <- is_hartwig_full %>%
  filter(group == "Human-associated") %>%
  distinct(mag_id, family) %>%
  count(family, name = "n_hartwig") %>%
  mutate(pct_hartwig = n_hartwig / N_HARTWIG * 100)

shared_is_df <- inner_join(hartwig_is_ha_prev, wirbel_is_prev, by = "family") %>%
  mutate(combined = n_hartwig + n_wirbel) %>%
  arrange(desc(combined))

N_SHARED_C <- nrow(shared_is_df)

shared_mirror_C <- shared_is_df %>%
  mutate(display_name = factor(family, levels = rev(family))) %>%
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
  labs(title = sprintf("Shared IS families: Hartwig human-associated vs Wirbel CRC (n = %d)", N_SHARED_C),
       x = "% of MAGs with this IS family", y = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 8.5, face = "bold"),
        panel.grid.major.y = element_blank(), plot.title = element_text(size = 12, face = "bold"))

save_data(shared_is_df, "Figure_4C_shared_IS_composition_data")

# Assemble and save
save_fig(panel_A, "Figure_4A_hartwig_IS_composition", w = 10, h = 11)
save_fig(panel_B, "Figure_4B_wirbel_IS_composition",  w = 8,  h = 8.5)
save_fig(panel_C, "Figure_4C_shared_IS_composition",  w = 9,  h = 9)

figure_4_top <- (panel_A | panel_B) + plot_layout(widths = c(1, 0.8))

# Tag styling must be applied with & theme(); plot_annotation(theme = ) ignores plot.tag.
figure_4 <- (figure_4_top / panel_C) +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))

save_fig(figure_4, "Figure_4", w = 18, h = 18)

cat(sprintf("Panel A: %d IS families shown\n", n_is_A))
cat(sprintf("Panel B: top 20 of %d Wirbel IS families\n", nrow(wirbel_is_prev)))
cat(sprintf("Panel C: %d shared IS families\n", N_SHARED_C))
