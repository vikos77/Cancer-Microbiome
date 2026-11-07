# Predicted phage host ecology of Hartwig human-associated and probable-contaminant MAGs
# (Figure S1).
#
# Each iPHoP-predicted phage host is called human-associated or not with the same host lists as
# the classification (screening_allowlists.tsv). Each MAG is then placed in one category: at
# least one human-associated host predicted, only non-human-associated hosts predicted, or no
# prediction. The two MAG groups are compared with Fisher's exact test (human-associated host
# predicted vs the rest).
#
# Input:  MAG classification and host lists from 00_hartwig_mag_classification.R
# Output: output/figures/figure_S1/, output/figure_data/figure_S1/

library(tidyverse)
library(scales)

# Run from the repository root: Rscript scripts/07_phage_host_prediction.R
stopifnot("Run this script from the repository root" = dir.exists("data"))
CLASS_DIR <- file.path("output", "mag_classification")
FIG_DIR   <- file.path("output", "figures", "figure_S1")
DATA_DIR  <- file.path("output", "figure_data", "figure_S1")
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

# -----------------------------------------------------------------------------
# MAG groups and host lists
# -----------------------------------------------------------------------------

hartwig_class <- read_tsv(file.path(CLASS_DIR, "hartwig_mag_classification.tsv"), show_col_types = FALSE) %>%
  transmute(mag_id, group = factor(classification, levels = c("Human-associated", "Probable contaminant")))
N_HUMAN  <- sum(hartwig_class$group == "Human-associated")
N_CONTAM <- sum(hartwig_class$group == "Probable contaminant")
stopifnot(N_HUMAN == 383, N_CONTAM == 131)

allowlists   <- read_tsv(file.path(CLASS_DIR, "screening_allowlists.tsv"), show_col_types = FALSE) %>%
  filter(list == "iPHoP host classification")
IPHOP_PHYLA  <- allowlists$taxon[allowlists$rank == "Phylum"]
IPHOP_GENERA <- allowlists$taxon[allowlists$rank == "Genus"]

# -----------------------------------------------------------------------------
# Host ecology per MAG
# -----------------------------------------------------------------------------

iphop <- read_csv(file.path("data", "hartwig", "iphop_hartwig.csv"), show_col_types = FALSE) %>%
  transmute(mag_id,
            host_genus          = str_extract(`Host genus`, "g__([^;,]+)") %>% str_remove("g__"),
            is_human_assoc_host = str_detect(`Host genus`, paste(IPHOP_PHYLA, collapse = "|")) |
                                  host_genus %in% IPHOP_GENERA)

CAT_LEVELS <- c("Human-associated-host predicted", "Non-human-associated-host predicted", "No host prediction")

iphop_mag_summary <- hartwig_class %>%
  left_join(iphop %>%
              group_by(mag_id) %>%
              summarise(human_host    = any(is_human_assoc_host,  na.rm = TRUE),
                        nonhuman_host = any(!is_human_assoc_host, na.rm = TRUE),
                        .groups = "drop"),
            by = "mag_id") %>%
  replace_na(list(human_host = FALSE, nonhuman_host = FALSE)) %>%
  mutate(category = factor(case_when(human_host    ~ "Human-associated-host predicted",
                                     nonhuman_host ~ "Non-human-associated-host predicted",
                                     TRUE          ~ "No host prediction"),
                           levels = CAT_LEVELS))
save_data(iphop_mag_summary, "Figure_S1_per_mag_iphop_host_ecology")

# TRUE column first so the odds ratio reads human-associated vs probable contaminant
ft <- fisher.test(table(iphop_mag_summary$group,
                        factor(iphop_mag_summary$category == "Human-associated-host predicted",
                               levels = c(TRUE, FALSE))))
save_data(tibble(comparison = "Human-associated vs Probable contaminant (human-associated-host predicted vs rest)",
                 odds_ratio = unname(ft$estimate), conf_low = ft$conf.int[1], conf_high = ft$conf.int[2],
                 p_value = ft$p.value),
          "Figure_S1_fisher_test_human_associated_host_predicted")

print(count(iphop_mag_summary, group, category) %>% group_by(group) %>% mutate(pct = round(100 * n / sum(n), 1)))
cat(sprintf("Fisher's exact: OR = %.3g, p = %.3g\n", unname(ft$estimate), ft$p.value))

# Values reported in the Figure S1 caption
pct_human_host <- iphop_mag_summary %>%
  group_by(group) %>%
  summarise(pct = 100 * mean(category == "Human-associated-host predicted"), .groups = "drop")
stopifnot(round(pct_human_host$pct[pct_human_host$group == "Human-associated"]) == 77,
          pct_human_host$pct[pct_human_host$group == "Probable contaminant"] < 1,
          round(unname(ft$estimate)) == 430, signif(ft$p.value, 2) == 5.7e-61)

# -----------------------------------------------------------------------------
# Figure S1
# -----------------------------------------------------------------------------

CAT_COLOURS <- c(
  "Human-associated-host predicted"     = "#0072B2",
  "Non-human-associated-host predicted" = "#E69F00",
  "No host prediction"                  = "#BBBBBB"
)

# Label positions computed from the stacking order (first category on top)
plot_data <- iphop_mag_summary %>%
  count(group, category, .drop = FALSE) %>%
  group_by(group) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  arrange(group, category) %>%
  mutate(ymax = rev(cumsum(rev(pct))), ymin = ymax - pct, ymid = (ymax + ymin) / 2) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, "\n(n = ", if_else(group == "Human-associated", N_HUMAN, N_CONTAM), ")"))

# Segments under 4% are too thin for an inside label, so their label sits above the bar
LABEL_PCT_THRESHOLD <- 4

figure_S1 <- ggplot(plot_data, aes(x = group_label, y = pct, fill = category)) +
  geom_col(width = 0.55, colour = "white", linewidth = 0.4) +
  geom_text(data = plot_data %>% filter(n > 0, pct >= LABEL_PCT_THRESHOLD),
            aes(y = ymid, label = paste0(round(pct, 0), "%\n(n=", n, ")")),
            size = 3.3, colour = "white", fontface = "bold") +
  geom_text(data = plot_data %>% filter(n > 0, pct < LABEL_PCT_THRESHOLD),
            aes(y = ymax + 2.6,
                label = paste0(if_else(pct < 1, "<1", as.character(round(pct, 0))), "%  (n=", n, ")")),
            size = 3.0, colour = "grey20", fontface = "bold") +
  scale_fill_manual(values = CAT_COLOURS, name = "iPHoP host-ecology prediction") +
  scale_y_continuous(breaks = seq(0, 100, 25), labels = label_percent(scale = 1),
                     expand = expansion(mult = c(0, 0.10))) +
  labs(x = NULL, y = "% of MAGs") +
  theme_bw(base_size = 12) +
  theme(legend.text      = element_text(size = 9.5),
        legend.title     = element_text(size = 10, face = "bold"),
        axis.text        = element_text(size = 10),
        axis.text.x      = element_text(face = "bold"),
        axis.title       = element_text(size = 11.5),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(size = 13, face = "bold"))

save_fig(figure_S1, "Figure_S1", w = 8, h = 6.5)
