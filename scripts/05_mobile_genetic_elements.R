# Mobile genetic elements in Hartwig human-associated, probable-contaminant and Wirbel CRC MAGs
# (Figure 5).
#   A  Provirus burden per MAG (CheckV)
#   B  Plasmid burden per MAG (geNomad)
#   C  Phage lifecycle of proviral contigs (Pharokka)
#   D  CDS per proviral contig (Pharokka)
# Each panel compares the three groups pairwise (Mann-Whitney, or Fisher's exact for C), with
# BH correction across the three comparisons in that panel.
#
# Input:  MAG classifications from 00_hartwig_mag_classification.R and 01_wirbel_mag_classification.R
# Output: output/figures/figure_5/, output/figure_data/figure_5/

library(tidyverse)
library(patchwork)
library(scales)

# Run from the repository root: Rscript scripts/05_mobile_genetic_elements.R
stopifnot("Run this script from the repository root" = dir.exists("data"))
CLASS_DIR <- file.path("output", "mag_classification")
HBASE     <- file.path("data", "hartwig")
WBASE     <- file.path("data", "wirbel")
FIG_DIR   <- file.path("output", "figures", "figure_5")
DATA_DIR  <- file.path("output", "figure_data", "figure_5")
dir.create(FIG_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

save_fig <- function(plot, name, w, h) {
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), plot, width = w, height = h)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), plot, width = w, height = h, dpi = 600)
  message("Saved figure: ", file.path(FIG_DIR, name))
}

save_data <- function(df, name) {
  write_tsv(df, file.path(DATA_DIR, paste0(name, ".tsv")))
  message("Saved data:   ", file.path(DATA_DIR, paste0(name, ".tsv")))
}

GROUP3_LEVELS <- c("Human-associated", "Probable contaminant", "Wirbel CRC")

GROUP3_COLOURS <- c(
  "Human-associated"     = "#0072B2",
  "Probable contaminant" = "#E69F00",
  "Wirbel CRC"           = "#4D4D4D"
)
LIFESTYLE_COLOURS <- c("Temperate (lysogenic)" = "#6A3D9A", "Lytic" = "#E69F00")

SHARED_THEME <- theme_bw(base_size = 11) +
  theme(
    legend.text        = element_text(size = 8.5),
    axis.text.y        = element_text(size = 9,  face = "bold"),
    axis.text.x        = element_text(size = 9,  face = "bold"),
    axis.title         = element_text(size = 11),
    plot.title         = element_text(size = 12, face = "bold")
  )

sig_label <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")

# All three pairwise brackets on every panel: the two adjacent pairs level at y1, the
# full span (human-associated vs Wirbel) above at y2
sig_brackets_3 <- function(p12, p23, p13, y1, y2, label_size = 3.4) {
  list(
    annotate("segment", x = 1,    xend = 1.95, y = y1, yend = y1, colour = "grey40", linewidth = 0.45),
    annotate("text",    x = 1.5,  y = y1, label = sig_label(p12), vjust = -0.35, size = label_size, colour = "grey20", fontface = "bold"),
    annotate("segment", x = 2.05, xend = 3,    y = y1, yend = y1, colour = "grey40", linewidth = 0.45),
    annotate("text",    x = 2.5,  y = y1, label = sig_label(p23), vjust = -0.35, size = label_size, colour = "grey20", fontface = "bold"),
    annotate("segment", x = 1,    xend = 3,    y = y2, yend = y2, colour = "grey40", linewidth = 0.45),
    annotate("text",    x = 2,    y = y2, label = sig_label(p13), vjust = -0.35, size = label_size, colour = "grey20", fontface = "bold")
  )
}

# Three pairwise tests for one panel, BH-corrected within the panel
pairwise_3 <- function(test_fun, panel) {
  pairs <- list(c("Human-associated", "Probable contaminant"),
                c("Human-associated", "Wirbel CRC"),
                c("Probable contaminant", "Wirbel CRC"))
  res <- map_dfr(pairs, function(p) {
    t <- test_fun(p[1], p[2])
    tibble(panel = panel, group_1 = p[1], group_2 = p[2],
           odds_ratio = if (is.null(t$estimate)) NA_real_ else unname(t$estimate),
           p = t$p.value)
  })
  res$p_adj <- p.adjust(res$p, method = "BH")
  res
}

mann_whitney <- function(df, value) {
  function(g1, g2) wilcox.test(as.formula(paste(value, "~ group")),
                               data = df %>% filter(group %in% c(g1, g2)) %>% droplevels(),
                               exact = FALSE)
}

p_of <- function(tests, g1, g2) tests$p_adj[tests$group_1 == g1 & tests$group_2 == g2]

# -----------------------------------------------------------------------------
# MAG groups
# -----------------------------------------------------------------------------

hartwig_class <- read_tsv(file.path(CLASS_DIR, "hartwig_mag_classification.tsv"), show_col_types = FALSE)
wirbel_ids    <- read_tsv(file.path(CLASS_DIR, "wirbel_mag_classification.tsv"), show_col_types = FALSE)$mag_id
contaminant_ids <- hartwig_class$mag_id[hartwig_class$classification == "Probable contaminant"]

N_HARTWIG <- sum(hartwig_class$classification == "Human-associated")
N_CONTAM  <- length(contaminant_ids)
N_WIRBEL  <- length(wirbel_ids)
stopifnot(N_HARTWIG == 383, N_CONTAM == 131, N_WIRBEL == 785)

group_of <- function(mag_ids) {
  factor(if_else(mag_ids %in% contaminant_ids, "Probable contaminant", "Human-associated"),
         levels = GROUP3_LEVELS)
}

x_labels3 <- c(
  "Human-associated"     = sprintf("Human-\nassociated\n(n = %d)", N_HARTWIG),
  "Probable contaminant" = sprintf("Probable\ncontaminant\n(n = %d)", N_CONTAM),
  "Wirbel CRC"           = sprintf("Wirbel CRC\n(n = %d)", N_WIRBEL)
)

# -----------------------------------------------------------------------------
# Panel A: provirus burden
# -----------------------------------------------------------------------------

hartwig_provirus <- read_tsv(file.path(HBASE, "checkv_summary_hartwig.tsv"),
                             show_col_types = FALSE) %>%
  right_join(select(hartwig_class, mag_id), by = "mag_id") %>%
  mutate(n_proviruses = replace_na(n_proviruses, 0), group = group_of(mag_id)) %>%
  select(mag_id, n_proviruses, group)

wirbel_checkv_raw <- read_tsv(file.path(WBASE, "checkv_contigs_wirbel.tsv"), show_col_types = FALSE) %>%
  mutate(mag_id = paste(sample_id, bin_id, sep = "_"))

wirbel_provirus <- tibble(mag_id = wirbel_ids) %>%
  left_join(wirbel_checkv_raw %>% filter(provirus == "Yes") %>% count(mag_id, name = "n_proviruses"),
            by = "mag_id") %>%
  replace_na(list(n_proviruses = 0)) %>%
  mutate(group = factor("Wirbel CRC", levels = GROUP3_LEVELS))

provirus_df <- bind_rows(hartwig_provirus, wirbel_provirus)
stopifnot(sum(hartwig_provirus$n_proviruses) == 212)

tests_A <- pairwise_3(mann_whitney(provirus_df, "n_proviruses"), "A")

pv_max <- quantile(provirus_df$n_proviruses, 0.98)
pv_y1  <- pv_max * 1.12
pv_y2  <- pv_max * 1.40

set.seed(42)
panel_A <- ggplot(provirus_df, aes(x = group, y = n_proviruses, colour = group, fill = group)) +
  geom_boxplot(alpha = 0.20, outlier.shape = NA, linewidth = 0.6, width = 0.45) +
  geom_jitter(alpha = 0.18, size = 0.8, width = 0.14, show.legend = FALSE) +
  sig_brackets_3(p_of(tests_A, "Human-associated", "Probable contaminant"),
                 p_of(tests_A, "Probable contaminant", "Wirbel CRC"),
                 p_of(tests_A, "Human-associated", "Wirbel CRC"), pv_y1, pv_y2) +
  scale_fill_manual(values = GROUP3_COLOURS,   guide = "none") +
  scale_colour_manual(values = GROUP3_COLOURS, guide = "none") +
  scale_x_discrete(labels = x_labels3) +
  coord_cartesian(ylim = c(0, pv_max * 1.55)) +
  labs(title = "Provirus burden (CheckV)", x = NULL, y = "CheckV-confirmed proviruses / MAG") +
  SHARED_THEME

save_data(provirus_df, "Figure_5A_provirus_burden_data")

# -----------------------------------------------------------------------------
# Panel B: plasmid burden
# -----------------------------------------------------------------------------

hartwig_plasmid <- select(hartwig_class, mag_id) %>%
  left_join(read_tsv(file.path(HBASE, "genomad_plasmid_hartwig.tsv"),
                     show_col_types = FALSE, col_types = cols(.default = "c")) %>%
              count(mag_id, name = "n_plasmids"),
            by = "mag_id") %>%
  mutate(n_plasmids = replace_na(n_plasmids, 0), group = group_of(mag_id))

wirbel_plasmid <- tibble(mag_id = wirbel_ids) %>%
  left_join(read_tsv(file.path(WBASE, "genomad_plasmid_wirbel.tsv"), show_col_types = FALSE) %>%
              count(mag_id = paste(sample_id, bin_id, sep = "_"), name = "n_plasmids"),
            by = "mag_id") %>%
  replace_na(list(n_plasmids = 0)) %>%
  mutate(group = factor("Wirbel CRC", levels = GROUP3_LEVELS))

plasmid_df <- bind_rows(hartwig_plasmid, wirbel_plasmid)

tests_B <- pairwise_3(mann_whitney(plasmid_df, "n_plasmids"), "B")

pl_max_val <- max(plasmid_df$n_plasmids)
pl_y1 <- pl_max_val * 1.35
pl_y2 <- pl_max_val * 1.75

set.seed(42)
panel_B <- ggplot(plasmid_df, aes(x = group, y = n_plasmids, colour = group, fill = group)) +
  geom_boxplot(alpha = 0.20, outlier.shape = NA, linewidth = 0.6, width = 0.45) +
  geom_jitter(alpha = 0.18, size = 0.8, width = 0.14, show.legend = FALSE) +
  sig_brackets_3(p_of(tests_B, "Human-associated", "Probable contaminant"),
                 p_of(tests_B, "Probable contaminant", "Wirbel CRC"),
                 p_of(tests_B, "Human-associated", "Wirbel CRC"), pl_y1, pl_y2) +
  scale_fill_manual(values = GROUP3_COLOURS,   guide = "none") +
  scale_colour_manual(values = GROUP3_COLOURS, guide = "none") +
  scale_x_discrete(labels = x_labels3) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0, 1, 5, 10, 25, 50, 100, 250, 500),
                     labels = c("0", "1", "5", "10", "25", "50", "100", "250", "500"),
                     limits = c(0, pl_max_val * 2.2)) +
  labs(title = "Plasmid burden (geNomad, log1p scale)", x = NULL,
       y = "geNomad plasmid sequences / MAG  [log(n+1)]") +
  SHARED_THEME

save_data(plasmid_df, "Figure_5B_plasmid_burden_data")

# -----------------------------------------------------------------------------
# Panels C and D: proviral contigs annotated by Pharokka
# -----------------------------------------------------------------------------

# Per-contig CheckV calls for Hartwig were regenerated with CheckV v1.0 (DIAMOND 2.0.4,
# database v1.0), as only the per-MAG summary of the original run was kept. They reproduce the
# original provirus counts in every MAG.
hartwig_checkv_pv <- read_tsv(file.path(HBASE, "checkv_contigs_hartwig.tsv"),
                              show_col_types = FALSE) %>%
  filter(provirus == "Yes") %>%
  select(mag_id, contig_id)
stopifnot(hartwig_provirus %>%
            left_join(count(hartwig_checkv_pv, mag_id, name = "n_contigs"), by = "mag_id") %>%
            with(all(n_proviruses == replace_na(n_contigs, 0L))))

contig_to_mag <- read_tsv(file.path(HBASE, "genomad_virus_hartwig.tsv"),
                          show_col_types = FALSE, col_types = cols(.default = "c")) %>%
  distinct(seq_name, mag_id)

hartwig_mapped <- read_tsv(file.path(HBASE, "pharokka_cds_hartwig.tsv"),
                           show_col_types = FALSE) %>%
  left_join(contig_to_mag, by = c("contig" = "seq_name")) %>%
  filter(!is.na(mag_id)) %>%
  inner_join(hartwig_checkv_pv, by = c("mag_id", "contig" = "contig_id")) %>%
  mutate(group = group_of(mag_id),
         category_clean = if_else(is.na(category) | category == "", "unknown function", category))

wirbel_mapped <- read_tsv(file.path(WBASE, "pharokka_cds_wirbel.tsv"), show_col_types = FALSE,
                          col_select = c(sample_id, bin_id, contig, category)) %>%
  inner_join(wirbel_checkv_raw %>% filter(provirus == "Yes") %>% select(sample_id, bin_id, contig_id),
             by = c("sample_id", "bin_id", "contig" = "contig_id")) %>%
  mutate(group = factor("Wirbel CRC", levels = GROUP3_LEVELS),
         category_clean = if_else(is.na(category) | category == "", "unknown function", category))
# Contigs are grouped by name below, which is safe only while proviral contig names are unique
stopifnot(n_distinct(wirbel_mapped$contig) == n_distinct(paste(wirbel_mapped$sample_id, wirbel_mapped$contig)))

# A contig is temperate if Pharokka annotates any integration/excision gene on it
per_contig <- function(mapped, ...) mapped %>% group_by(contig, group) %>% summarise(..., .groups = "drop")
lifestyle_all <- bind_rows(
  per_contig(hartwig_mapped, is_temperate = any(str_detect(category_clean, "integration and excision"))),
  per_contig(wirbel_mapped,  is_temperate = any(str_detect(category_clean, "integration and excision")))
)

lifestyle_pct <- lifestyle_all %>%
  count(group, is_temperate) %>%
  group_by(group) %>%
  mutate(n_total = sum(n), pct = n / n_total * 100,
         lifestyle = factor(if_else(is_temperate, "Temperate (lysogenic)", "Lytic"),
                            levels = c("Temperate (lysogenic)", "Lytic"))) %>%
  ungroup()

# Values reported in the Figure 5 caption: 42 of 44 human-associated proviral contigs temperate
ha_lifestyle <- filter(lifestyle_pct, group == "Human-associated", is_temperate)
stopifnot(ha_lifestyle$n == 42, ha_lifestyle$n_total == 44)

lifecycle_fisher <- function(g1, g2) {
  tally <- function(g) lifestyle_all %>% filter(group == g) %>%
    summarise(T = sum(is_temperate), L = sum(!is_temperate))
  t1 <- tally(g1); t2 <- tally(g2)
  fisher.test(matrix(c(t1$T, t1$L, t2$T, t2$L), nrow = 2))
}
tests_C <- pairwise_3(lifecycle_fisher, "C")

bar_labels <- lifestyle_pct %>% filter(is_temperate) %>%
  mutate(label = sprintf("%.1f%%\n(%d/%d)", pct, n, n_total))

contig_n <- lifestyle_pct %>% distinct(group, n_total) %>% deframe()
x_labels_contigs <- function(n) c(
  "Human-associated"     = sprintf("Human-\nassociated\n(%d contigs)", n[["Human-associated"]]),
  "Probable contaminant" = sprintf("Probable\ncontaminant\n(%d contigs)", n[["Probable contaminant"]]),
  "Wirbel CRC"           = sprintf("Wirbel CRC\n(%d contigs)", n[["Wirbel CRC"]])
)

panel_C <- ggplot(lifestyle_pct, aes(x = group, y = pct, fill = lifestyle)) +
  geom_col(width = 0.55, colour = "white", linewidth = 0.3) +
  geom_text(data = bar_labels, aes(x = group, y = 100 - pct / 2, label = label),
            colour = "white", fontface = "bold", size = 3.3, inherit.aes = FALSE) +
  sig_brackets_3(p_of(tests_C, "Human-associated", "Probable contaminant"),
                 p_of(tests_C, "Probable contaminant", "Wirbel CRC"),
                 p_of(tests_C, "Human-associated", "Wirbel CRC"), 112, 128) +
  scale_fill_manual(values = LIFESTYLE_COLOURS, name = NULL) +
  scale_x_discrete(labels = x_labels_contigs(contig_n)) +
  scale_y_continuous(labels = label_percent(scale = 1), limits = c(0, 140), expand = c(0, 0),
                     breaks = seq(0, 100, 25)) +
  coord_cartesian(ylim = c(0, 140)) +
  labs(title = "Phage lifecycle (Pharokka, CheckV proviruses)", x = NULL,
       y = "Percentage of proviral contigs") +
  SHARED_THEME +
  theme(legend.position = "bottom", legend.key.size = unit(0.4, "cm"),
        panel.grid.major.x = element_blank())

save_data(lifestyle_pct, "Figure_5C_phage_lifestyle_data")

# Panel D
cds_df <- bind_rows(per_contig(hartwig_mapped, n_cds = n()), per_contig(wirbel_mapped, n_cds = n()))

tests_D <- pairwise_3(mann_whitney(cds_df, "n_cds"), "D")

cds_max <- quantile(cds_df$n_cds, 0.97)
cds_y1  <- cds_max * 1.12
cds_y2  <- cds_max * 1.40

set.seed(42)
panel_D <- ggplot(cds_df, aes(x = group, y = n_cds, colour = group, fill = group)) +
  geom_boxplot(alpha = 0.20, outlier.shape = NA, linewidth = 0.6, width = 0.45) +
  geom_jitter(alpha = 0.18, size = 0.8, width = 0.14, show.legend = FALSE) +
  sig_brackets_3(p_of(tests_D, "Human-associated", "Probable contaminant"),
                 p_of(tests_D, "Probable contaminant", "Wirbel CRC"),
                 p_of(tests_D, "Human-associated", "Wirbel CRC"), cds_y1, cds_y2) +
  scale_fill_manual(values = GROUP3_COLOURS,   guide = "none") +
  scale_colour_manual(values = GROUP3_COLOURS, guide = "none") +
  scale_x_discrete(labels = x_labels_contigs(table(cds_df$group))) +
  coord_cartesian(ylim = c(0, cds_max * 1.55)) +
  labs(title = "Phage gene content (Pharokka)", x = NULL, y = "CDS per proviral contig") +
  SHARED_THEME

save_data(cds_df, "Figure_5D_cds_per_phage_data")

# -----------------------------------------------------------------------------
# Pairwise tests and assembly
# -----------------------------------------------------------------------------

pairwise_tests <- bind_rows(tests_A, tests_B, tests_C, tests_D) %>%
  mutate(test = if_else(panel == "C", "Fisher's exact", "Mann-Whitney"), .after = panel)
save_data(pairwise_tests, "Figure_5_pairwise_tests")
print(pairwise_tests %>% mutate(across(c(p, p_adj), ~ signif(., 3))), n = Inf)

save_fig(panel_A, "Figure_5A_provirus_burden", w = 8, h = 8)
save_fig(panel_B, "Figure_5B_plasmid_burden",  w = 8, h = 8)
save_fig(panel_C, "Figure_5C_phage_lifestyle", w = 8, h = 8.5)
save_fig(panel_D, "Figure_5D_cds_per_phage",   w = 8, h = 8)

# Tag styling must be applied with & theme(); plot_annotation(theme = ) ignores plot.tag.
figure_5 <- (panel_A | panel_B) / (panel_C | panel_D) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))

save_fig(figure_5, "Figure_5", w = 16, h = 16)

cat(sprintf("Panels C, D: %d / %d / %d proviral contigs\n",
            contig_n[["Human-associated"]], contig_n[["Probable contaminant"]], contig_n[["Wirbel CRC"]]))
