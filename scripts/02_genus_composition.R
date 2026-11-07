# Genus composition of the Hartwig and Wirbel MAGs (Figure 2).
#   A  Hartwig human-associated vs probable-contaminant MAGs
#   B  Wirbel CRC, top 20 genera
#   C  Genera shared between Hartwig human-associated and Wirbel CRC
# Input:  MAG classifications from 00_hartwig_mag_classification.R and 01_wirbel_mag_classification.R
# Output: output/figures/figure_2/, output/figure_data/figure_2/

library(tidyverse)
library(patchwork)
library(scales)

# Run from the repository root: Rscript scripts/02_genus_composition.R
stopifnot("Run this script from the repository root" = dir.exists("data"))
FIG_DIR    <- file.path("output", "figures", "figure_2")
DATA_DIR   <- file.path("output", "figure_data", "figure_2")
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

phylum_colours <- c(
  "Firmicutes"        = "#4477AA",
  "Bacteroidota"      = "#66CCEE",
  "Proteobacteria"    = "#228833",
  "Fusobacteriota"    = "#CCBB44",
  "Actinomycetota"    = "#EE6677",
  "Campylobacterota"  = "#AA3377",
  "Spirochaetota"     = "#332288",
  "Synergistota"      = "#000000",
  "Verrucomicrobiota" = "#004D25",
  "Other"             = "#BBBBBB"
)

phylum_clean <- function(phylum_raw) {
  case_when(
    str_detect(phylum_raw, "^Bacillota|^Firmicutes")          ~ "Firmicutes",
    str_detect(phylum_raw, "^Bacteroidota")                   ~ "Bacteroidota",
    str_detect(phylum_raw, "^Pseudomonadota|^Proteobacteria") ~ "Proteobacteria",
    str_detect(phylum_raw, "^Fusobacteriota")                 ~ "Fusobacteriota",
    str_detect(phylum_raw, "^Actinomycetota")                 ~ "Actinomycetota",
    str_detect(phylum_raw, "^Campylobacterota")               ~ "Campylobacterota",
    str_detect(phylum_raw, "^Spirochaetota")                  ~ "Spirochaetota",
    str_detect(phylum_raw, "^Synergistota")                   ~ "Synergistota",
    str_detect(phylum_raw, "^Verrucomicrobiota")              ~ "Verrucomicrobiota",
    TRUE                                                       ~ "Other"
  )
}

MIRROR_COLOURS_A <- c(
  "Probable-contaminant-exclusive"       = "#B36F00",
  "Probable-contaminant-enriched"        = "#E69F00",
  "Human-associated-enriched"  = "#56B4E9",
  "Human-associated-exclusive" = "#0072B2"
)

HW_WB_COLOURS <- c(
  "Hartwig human-associated" = "#0072B2",
  "Wirbel CRC"                = "#4D4D4D"
)

# Latin genus names are italicised. GTDB placeholder names (e.g. UBA11524, CAG-177) stay roman.
is_latin_genus <- function(x) grepl("^[A-Z][a-z]+(_[A-Z0-9]+)?$", x)

italic_genus_labels <- function(genus_levels) {
  lapply(genus_levels, function(g) {
    if (is_latin_genus(g)) bquote(bolditalic(.(g))) else bquote(bold(.(g)))
  })
}

SHARED_THEME <- theme_bw(base_size = 11) +
  theme(
    legend.text        = element_text(size = 8.5),
    axis.text.y        = element_text(size = 9, face = "bold"),
    axis.title         = element_text(size = 11),
    panel.grid.major.y = element_blank(),
    plot.title         = element_text(size = 12, face = "bold"),
    plot.subtitle      = element_text(size = 9, colour = "grey40")
  )

CLASS_DIR <- file.path("output", "mag_classification")
hartwig_class <- read_tsv(file.path(CLASS_DIR, "hartwig_mag_classification.tsv"),
                          show_col_types = FALSE) %>%
  select(mag_id, mag_group = classification)
wirbel_class <- read_tsv(file.path(CLASS_DIR, "wirbel_mag_classification.tsv"),
                         show_col_types = FALSE) %>%
  select(mag_id, mag_group = classification)

hartwig_tax <- read_csv(file.path("data", "hartwig", "gtdbtk_checkm2_hartwig.csv"),
                        show_col_types = FALSE) %>%
  inner_join(hartwig_class, by = c("bin_id" = "mag_id")) %>%
  mutate(group    = mag_group,
         phylum_c = phylum_clean(Phylum))

N_HARTWIG <- sum(hartwig_tax$group == "Human-associated")
N_CONTAM  <- sum(hartwig_tax$group == "Probable contaminant")
stopifnot(N_HARTWIG == 383, N_CONTAM == 131)

wirbel_tax <- read_tsv(
  file.path("data", "wirbel", "gtdbtk_taxonomy_wirbel.tsv"),
  show_col_types = FALSE
) %>%
  inner_join(wirbel_class, by = c("gtdb_genome_id" = "mag_id")) %>%
  filter(mag_group == "Human-associated") %>%
  mutate(phylum_c = phylum_clean(phylum_name))

N_WIRBEL <- nrow(wirbel_tax)
stopifnot(N_WIRBEL == 785)

cat(sprintf("Hartwig: %d human-associated + %d probable-contaminant MAGs\n", N_HARTWIG, N_CONTAM))
cat(sprintf("Wirbel:  %d MAGs\n", N_WIRBEL))

genus_wide <- hartwig_tax %>%
  filter(Genus != "", !is.na(Genus)) %>%
  count(group, Genus) %>%
  complete(group, Genus, fill = list(n = 0)) %>%
  left_join(tibble(group = c("Human-associated", "Probable contaminant"),
                    n_total = c(N_HARTWIG, N_CONTAM)),
            by = "group") %>%
  mutate(prevalence = n / n_total * 100) %>%
  pivot_wider(id_cols = Genus, names_from = group, values_from = c(n, prevalence),
              values_fill = 0) %>%
  rename(n_human_associated       = `n_Human-associated`,
         n_probable_contaminant   = `n_Probable contaminant`,
         pct_human_associated     = `prevalence_Human-associated`,
         pct_probable_contaminant = `prevalence_Probable contaminant`) %>%
  mutate(
    pct_difference = pct_human_associated - pct_probable_contaminant,
    status = case_when(
      pct_human_associated     == 0 ~ "Probable-contaminant-exclusive",
      pct_probable_contaminant == 0 ~ "Human-associated-exclusive",
      pct_difference > 0            ~ "Human-associated-enriched",
      TRUE                          ~ "Probable-contaminant-enriched"
    )
  )

N_DISTINCT_GENERA <- nrow(genus_wide)

# Panel A: every genus found among probable contaminants, topped up to 20 with the
# most frequent human-associated-only genera. Remaining genera are not shown.
TOP_N_GENUS  <- 20
contam_genera <- genus_wide %>% filter(n_probable_contaminant > 0) %>% pull(Genus)
top_ha_only_genera <- genus_wide %>%
  filter(!Genus %in% contam_genera) %>%
  arrange(desc(n_human_associated)) %>%
  slice_head(n = max(TOP_N_GENUS - length(contam_genera), 0)) %>%
  pull(Genus)
top_genera <- c(contam_genera, top_ha_only_genera)

panel_A_data <- genus_wide %>%
  filter(Genus %in% top_genera) %>%
  arrange(pct_difference) %>%
  mutate(display_name = factor(Genus, levels = Genus))

max_x_A  <- max(panel_A_data$pct_human_associated, panel_A_data$pct_probable_contaminant)
n_rows_A <- nrow(panel_A_data)

panel_A_mirror <- panel_A_data %>%
  mutate(display_name = factor(Genus, levels = rev(Genus))) %>%
  pivot_longer(cols = c(pct_human_associated, pct_probable_contaminant),
               names_to = "grp", values_to = "prevalence") %>%
  mutate(
    grp      = if_else(grp == "pct_human_associated", "Human-associated", "Probable contaminant"),
    x_val    = if_else(grp == "Probable contaminant", -prevalence, prevalence),
    n_label  = if_else(grp == "Human-associated", n_human_associated, n_probable_contaminant),
    fill_cat = case_when(
      str_ends(status, "exclusive")  ~ status,
      grp == "Probable contaminant"  ~ "Probable-contaminant-enriched",
      TRUE                           ~ "Human-associated-enriched"
    ),
    fill_cat = factor(fill_cat, levels = names(MIRROR_COLOURS_A))
  )

panel_A <- ggplot(panel_A_mirror, aes(x = x_val, y = display_name, fill = fill_cat)) +
  geom_col(width = 0.70, colour = NA) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.55) +
  geom_text(
    data = panel_A_mirror %>% filter(x_val != 0),
    aes(label = sprintf("%.0f%%  (n=%d)", abs(x_val), n_label),
        hjust = if_else(x_val < 0, 1.05, -0.05)),
    size = 2.65, colour = "grey20"
  ) +
  annotate("text", x = -(max_x_A / 2), y = n_rows_A + 1.1,
           label = sprintf("<- Probable contaminants  (n = %d)", N_CONTAM),
           hjust = 0.5, size = 3.0, colour = "#B36F00", fontface = "bold") +
  annotate("text", x =  (max_x_A / 2), y = n_rows_A + 1.1,
           label = sprintf("Human-associated  (n = %d) ->", N_HARTWIG),
           hjust = 0.5, size = 3.0, colour = "#0072B2", fontface = "bold") +
  scale_fill_manual(
    values = MIRROR_COLOURS_A, name = NULL,
    labels = c(
      "Probable-contaminant-exclusive"       = "Probable-contaminant-exclusive  (absent from all human-associated MAGs)",
      "Probable-contaminant-enriched"        = "Probable-contaminant-enriched   (present in both, higher in probable contaminants)",
      "Human-associated-enriched"  = "Human-associated-enriched   (present in both, higher in human-associated)",
      "Human-associated-exclusive" = "Human-associated-exclusive  (absent from all probable-contaminant MAGs)"
    )
  ) +
  scale_x_continuous(
    breaks = seq(-100, 100, 10),
    labels = function(x) paste0(abs(x), "%"),
    limits = c(-(max_x_A + 18), max_x_A + 18),
    expand = c(0, 0)
  ) +
  scale_y_discrete(labels = italic_genus_labels(levels(panel_A_mirror$display_name))) +
  coord_cartesian(ylim = c(0.5, n_rows_A + 1.8), clip = "off") +
  labs(
    title = "Hartwig genus composition: human-associated vs probable-contaminant MAGs (top 20)",
    x = "% of MAGs in this genus", y = NULL
  ) +
  SHARED_THEME +
  theme(legend.position = "bottom", legend.direction = "vertical",
        legend.key.size = unit(0.38, "cm"),
        plot.margin = margin(t = 25, r = 10, b = 5, l = 5))

save_data(genus_wide %>% arrange(desc(n_human_associated + n_probable_contaminant)),
          "Figure_2A_hartwig_genus_composition_data")

# Panel B
wirbel_genus_counts <- wirbel_tax %>%
  filter(genus_name != "", !is.na(genus_name)) %>%
  count(genus_name, phylum_c, name = "n_wirbel") %>%
  rename(Genus = genus_name) %>%
  mutate(pct_wirbel = n_wirbel / N_WIRBEL * 100)

top20_wirbel <- wirbel_genus_counts %>%
  slice_max(n_wirbel, n = 20, with_ties = FALSE) %>%
  arrange(n_wirbel) %>%
  mutate(Genus = factor(Genus, levels = Genus))

panel_B <- ggplot(top20_wirbel, aes(x = pct_wirbel, y = Genus, fill = phylum_c)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%  (n=%d)", pct_wirbel, n_wirbel)),
            hjust = -0.1, size = 2.8, colour = "grey20") +
  scale_x_continuous(limits = c(0, max(top20_wirbel$pct_wirbel) * 1.3),
                     expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  scale_y_discrete(labels = italic_genus_labels(levels(top20_wirbel$Genus))) +
  scale_fill_manual(values = phylum_colours, name = "Phylum") +
  labs(title = "Wirbel CRC: top 20 genera by MAG count",
       x = "% of MAGs in this genus", y = NULL) +
  SHARED_THEME +
  theme(legend.position = "bottom")

save_data(wirbel_genus_counts %>% arrange(desc(n_wirbel)), "Figure_2B_wirbel_genus_composition_data")

# Panel C: all genera present in both datasets
hartwig_genus_counts <- hartwig_tax %>%
  filter(group == "Human-associated", Genus != "", !is.na(Genus)) %>%
  count(Genus, phylum_c, name = "n_hartwig") %>%
  mutate(pct_hartwig = n_hartwig / N_HARTWIG * 100)

shared_genus_df <- inner_join(hartwig_genus_counts, wirbel_genus_counts,
                              by = "Genus", suffix = c("_h", "_w")) %>%
  mutate(phylum_c = coalesce(phylum_c_h, phylum_c_w),
         combined = n_hartwig + n_wirbel) %>%
  arrange(desc(combined))

N_SHARED <- nrow(shared_genus_df)

shared_mirror <- shared_genus_df %>%
  mutate(display_name = factor(Genus, levels = rev(Genus))) %>%
  pivot_longer(cols = c(pct_hartwig, pct_wirbel), names_to = "grp", values_to = "prevalence") %>%
  mutate(
    grp     = case_when(grp == "pct_hartwig" ~ "Hartwig human-associated", TRUE ~ "Wirbel CRC"),
    x_val   = if_else(grp == "Wirbel CRC", prevalence, -prevalence),
    n_label = if_else(grp == "Hartwig human-associated", n_hartwig, n_wirbel)
  )

max_x_C  <- max(abs(shared_mirror$x_val))
n_rows_C <- n_distinct(shared_mirror$display_name)

panel_C <- ggplot(shared_mirror, aes(x = x_val, y = display_name, fill = grp)) +
  geom_col(width = 0.72, colour = NA) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.55) +
  geom_text(
    aes(label = sprintf("%.1f%%  (n=%d)", abs(x_val), n_label),
        hjust = if_else(x_val < 0, 1.05, -0.05)),
    size = 2.5, colour = "grey20"
  ) +
  annotate("text", x = -(max_x_C / 2), y = n_rows_C + 1.4,
           label = sprintf("<- Hartwig human-associated  (n = %d)", N_HARTWIG),
           hjust = 0.5, size = 3.0, colour = "#0072B2", fontface = "bold") +
  annotate("text", x =  (max_x_C / 2), y = n_rows_C + 1.4,
           label = sprintf("Wirbel CRC  (n = %d) ->", N_WIRBEL),
           hjust = 0.5, size = 3.0, colour = "#4D4D4D", fontface = "bold") +
  scale_fill_manual(values = HW_WB_COLOURS, name = NULL) +
  scale_x_continuous(
    labels = function(x) paste0(abs(x), "%"),
    limits = c(-(max_x_C + 4), max_x_C + 4),
    expand = c(0, 0)
  ) +
  scale_y_discrete(labels = italic_genus_labels(levels(shared_mirror$display_name))) +
  coord_cartesian(ylim = c(0.5, n_rows_C + 2), clip = "off") +
  labs(
    title = sprintf("Shared genera: Hartwig human-associated vs Wirbel CRC (n = %d genera)", N_SHARED),
    x = "% of MAGs in this genus", y = NULL
  ) +
  SHARED_THEME +
  theme(legend.position = "bottom")

save_data(shared_genus_df, "Figure_2C_shared_genera_composition_data")

# Assemble and save
save_fig(panel_A, "Figure_2A_hartwig_genus_composition",   w = 10, h = 10)
save_fig(panel_B, "Figure_2B_wirbel_genus_composition",    w = 10, h = 8.5)
save_fig(panel_C, "Figure_2C_shared_genera_composition",   w = 9,  h = 15.5)

figure_2_top <- (panel_A | panel_B) + plot_layout(widths = c(1, 0.8))

figure_2 <- (figure_2_top / panel_C) +
  plot_layout(heights = c(1, 1.5)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))
# Tag styling must be applied with & theme(); plot_annotation(theme = ) ignores plot.tag.

save_fig(figure_2, "Figure_2", w = 18, h = 20.5)

cat(sprintf("Panel A: %d of %d Hartwig genera shown\n", n_rows_A, N_DISTINCT_GENERA))
cat(sprintf("Panel B: top 20 of %d Wirbel genera\n", nrow(wirbel_genus_counts)))
cat(sprintf("Panel C: %d shared genera\n", N_SHARED))
