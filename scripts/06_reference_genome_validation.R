# Reference-genome validation of the MAG defence and IS results.
#
# MAGs are compared with complete RefSeq genomes of the same genus (chromosome only), split by
# NCBI BioSample source into human and non-human (environmental) reference genomes. A genus or
# reference group is shown only if it has at least 10 genomes.
#   1. Burden: defence systems and IS elements per Mbp, MAG vs reference (Figure S2)
#   2. Composition: the ten most prevalent defence subtypes and IS families per genus in MAGs,
#      and how often they occur in reference genomes
#        Figure 6   human-associated genera and Wirbel Alistipes
#        Figure S3  probable-contaminant genera (built for reference, not part of the paper)
#
# Input:  MAG classifications from 00_hartwig_mag_classification.R;
#         reference-genome DefenseFinder/ISEScan tables in data/reference_genomes/
# Output: output/figures/ and output/figure_data/ (figure_6, figure_S2, figure_S3)

library(tidyverse)
library(patchwork)

# Run from the repository root: Rscript scripts/06_reference_genome_validation.R
stopifnot("Run this script from the repository root" = dir.exists("data"))
HARTWIG <- file.path("data", "hartwig")
REFDIR  <- file.path("data", "reference_genomes")
WIRBEL  <- file.path("data", "wirbel")
OUT <- list(
  fig6  = list(fig = file.path("output", "figures", "figure_6"),
               data = file.path("output", "figure_data", "figure_6")),
  figS2 = list(fig = file.path("output", "figures", "figure_S2"),
               data = file.path("output", "figure_data", "figure_S2")),
  figS3 = list(fig = file.path("output", "figures", "figure_S3"),
               data = file.path("output", "figure_data", "figure_S3"))
)
walk(unlist(OUT), dir.create, showWarnings = FALSE, recursive = TRUE)

save_fig <- function(plot, out, name, w, h) {
  ggsave(file.path(out$fig, paste0(name, ".pdf")), plot, width = w, height = h, limitsize = FALSE)
  ggsave(file.path(out$fig, paste0(name, ".png")), plot, width = w, height = h, dpi = 600, limitsize = FALSE)
  message("Saved figure: ", file.path(out$fig, name))
}

save_data <- function(df, out, name) {
  write_tsv(df, file.path(out$data, paste0(name, ".tsv")))
  message("Saved data:   ", file.path(out$data, paste0(name, ".tsv")))
}

GENUS_MIN_N <- 10
TOP_N       <- 10

# -----------------------------------------------------------------------------
# Genus scope
# -----------------------------------------------------------------------------

hartwig_class <- read_tsv(file.path("output", "mag_classification", "hartwig_mag_classification.tsv"),
                          show_col_types = FALSE)
contaminant_ids <- hartwig_class$mag_id[hartwig_class$classification == "Probable contaminant"]

gtdb <- read_csv(file.path(HARTWIG, "gtdbtk_checkm2_hartwig.csv"), show_col_types = FALSE) %>%
  rename(mag_id = bin_id)
gtdb_human  <- gtdb %>% filter(!mag_id %in% contaminant_ids)
gtdb_contam <- gtdb %>% filter(mag_id %in% contaminant_ids)
stopifnot(nrow(gtdb_human) == 383, nrow(gtdb_contam) == 131)

human_genus_n  <- gtdb_human  %>% count(Genus, name = "n_mag") %>% filter(n_mag >= GENUS_MIN_N) %>% arrange(desc(n_mag))
contam_genus_n <- gtdb_contam %>% count(Genus, name = "n_mag") %>% filter(n_mag >= GENUS_MIN_N) %>% arrange(desc(n_mag))
human_genera  <- human_genus_n$Genus
contam_genera <- contam_genus_n$Genus
stopifnot(identical(human_genera,  c("Prevotella", "Fusobacterium", "Bacteroides", "Streptococcus")),
          identical(contam_genera, c("Ralstonia", "Bradyrhizobium", "Pelomonas")))

# Wirbel Alistipes: 33 MAGs, 26 human reference genomes, no non-human reference genomes on NCBI
alistipes <- read_tsv(file.path(WIRBEL, "gtdbtk_taxonomy_wirbel.tsv"), show_col_types = FALSE) %>%
  filter(genus_name == "Alistipes") %>%
  select(source_sample, source_bin) %>%
  mutate(mag_id = paste(source_sample, source_bin, sep = "_"))
stopifnot(nrow(alistipes) == 33)

# MAG-side calls
df_hartwig <- read_tsv(file.path(HARTWIG, "defensefinder_hartwig.tsv"), show_col_types = FALSE) %>%
  filter(activity != "Antidefense")
is_hartwig <- read_tsv(file.path(HARTWIG, "isescan_hartwig.tsv"), show_col_types = FALSE)
df_wirbel  <- read_tsv(file.path(WIRBEL, "defensefinder_wirbel.tsv"), show_col_types = FALSE) %>%
  filter(activity == "Defense") %>%
  mutate(mag_id = paste(sample_id, bin_id, sep = "_"))
is_wirbel  <- read_tsv(file.path(WIRBEL, "isescan_wirbel.tsv"), show_col_types = FALSE) %>%
  mutate(mag_id = paste(sample_id, bin_id, sep = "_"))

# Reference-side calls (chromosome only)
ref_df_summary <- read_tsv(file.path(REFDIR, "defensefinder_summary_refgenome.tsv"), show_col_types = FALSE)
ref_is_summary <- read_tsv(file.path(REFDIR, "isescan_summary_refgenome.tsv"),        show_col_types = FALSE)
ref_df_detail  <- read_tsv(file.path(REFDIR, "defensefinder_detail_refgenome.tsv"),  show_col_types = FALSE) %>%
  filter(activity != "Antidefense")
ref_is_detail  <- read_tsv(file.path(REFDIR, "isescan_detail_refgenome.tsv"),        show_col_types = FALSE)

# -----------------------------------------------------------------------------
# Shared plotting style
# -----------------------------------------------------------------------------

dataset_colours <- c("MAG" = "#0072B2", "Human reference" = "#EBC9DC", "Non-human reference" = "#E49659")

FS_TITLE      <- 11
FS_AXIS_TITLE <- 10
FS_AXIS_TEXT  <- 9
FS_TAG        <- 13
# Tag position must be set on each panel; plot_annotation(theme = ) ignores plot.tag.location
tag_theme <- theme(plot.tag = element_text(face = "bold", size = FS_TAG),
                   plot.tag.position = "topleft", plot.tag.location = "plot")

LEGEND_STRIP_IN <- 1.5 / 2.54

build_legend_strip <- function(colours) {
  legend_src <- ggplot(tibble(dataset = factor(names(colours), levels = names(colours))),
                       aes(x = 1, y = 1, fill = dataset)) +
    geom_point(shape = 22, size = 5, colour = NA) +
    scale_fill_manual(values = colours, name = NULL, drop = FALSE) +
    guides(fill = guide_legend(nrow = 1, override.aes = list(shape = 22, size = 5))) +
    theme_void() +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.text = element_text(size = 10), legend.key.size = unit(0.45, "cm"))
  gt <- ggplotGrob(legend_src)
  idx <- which(grepl("guide-box", vapply(gt$grobs, function(g) g$name, character(1))))
  idx <- idx[vapply(gt$grobs[idx], function(g) !inherits(g, "zeroGrob"), logical(1))]
  stopifnot(length(idx) == 1)
  gt$grobs[[idx[1]]]
}

sys_short <- function(x) {
  mapping <- c(
    "CAS_Class1-Subtype-I-B"   = "CRISPR I-B",   "CAS_Class1-Subtype-I-C"   = "CRISPR I-C",
    "CAS_Class1-Subtype-I-E"   = "CRISPR I-E",   "CAS_Class1-Subtype-I-F"   = "CRISPR I-F",
    "CAS_Class2-Subtype-II-A"  = "CRISPR II-A",  "CAS_Class2-Subtype-II-C"  = "CRISPR II-C",
    "CAS_Class2-Subtype-VI-A"  = "CRISPR VI-A",
    "CAS_Class1-Subtype-III-A" = "CRISPR III-A", "CAS_Class1-Subtype-III-B" = "CRISPR III-B",
    "CAS_Class1-Subtype-III-D" = "CRISPR III-D",
    "RM_Type_I" = "RM Type I", "RM_Type_II" = "RM Type II", "RM_Type_III" = "RM Type III",
    "RM_Type_IIG" = "RM Type IIG", "RM_Type_IIG_2" = "RM Type IIG-2", "RM_Type_IV" = "RM Type IV",
    "RM_Type_IV_1" = "RM Type IV-1",
    "BREX_I" = "BREX-I", "PARIS_II" = "PARIS-II", "PARIS_I_merge" = "PARIS-I",
    "MazEF" = "MazEF (TA)", "Wadjet_I" = "Wadjet-I", "Wadjet_III" = "Wadjet-III",
    "Thoeris_I" = "Thoeris-I", "DarTG" = "DarTG",
    "Druantia_TypeI" = "Druantia-I", "Druantia_TypeIII" = "Druantia-III",
    "Zorya_TypeII" = "Zorya-II", "Zorya_TypeIII" = "Zorya-III", "SoFic" = "SoFic",
    "Hachiman" = "Hachiman", "Hachiman_II" = "Hachiman-II", "Prometheus" = "Prometheus",
    "Gabija" = "Gabija", "Lamassu_Fam1" = "Lamassu", "Lamassu-Mrr" = "Lamassu-Mrr"
  )
  ifelse(x %in% names(mapping), mapping[x], x)
}

# =============================================================================
# 1. Burden: systems per Mbp, MAG vs reference (Figure S2)
# =============================================================================

ref_pool <- ref_df_summary %>%
  select(genus, accession, source_category, chromosome_length_bp, n_defense) %>%
  left_join(ref_is_summary %>% select(genus, accession, n_IS_elements), by = c("genus", "accession")) %>%
  filter(chromosome_length_bp > 0) %>%
  mutate(defense_per_Mbp = n_defense / (chromosome_length_bp / 1e6),
         IS_per_Mbp      = n_IS_elements / (chromosome_length_bp / 1e6))

# Human-associated genera are compared with human reference genomes, and with non-human ones
# where at least 10 exist (not Fusobacterium). Contaminant genera are compared with their
# non-human reference genomes.
nonhuman_n <- ref_pool %>% filter(source_category == "environmental") %>% count(genus, name = "n_nonhuman")
human_genera_with_nonhuman <- nonhuman_n %>%
  filter(genus %in% human_genera, n_nonhuman >= GENUS_MIN_N) %>% pull(genus)
stopifnot(setequal(human_genera_with_nonhuman, c("Prevotella", "Bacteroides", "Streptococcus")),
          all(nonhuman_n$n_nonhuman[nonhuman_n$genus %in% contam_genera] >= GENUS_MIN_N))

ref_hartwig <- ref_pool %>%
  filter(genus %in% c(human_genera, contam_genera)) %>%
  mutate(genus_group = if_else(genus %in% contam_genera, "Probable-contaminant genera", "Human-associated genera")) %>%
  filter((genus_group == "Human-associated genera" & source_category == "human") |
         (genus_group == "Human-associated genera" & source_category == "environmental" &
            genus %in% human_genera_with_nonhuman) |
         (genus_group == "Probable-contaminant genera" & source_category == "environmental")) %>%
  mutate(dataset = if_else(source_category == "human", "Human reference", "Non-human reference"))

mag_hartwig <- gtdb %>%
  filter(Genus %in% c(human_genera, contam_genera)) %>%
  select(mag_id, genus = Genus, Length) %>%
  left_join(count(df_hartwig, mag_id, name = "n_defense"), by = "mag_id") %>%
  left_join(count(is_hartwig, mag_id, name = "n_IS_elements"), by = "mag_id") %>%
  mutate(n_defense       = replace_na(n_defense, 0),
         n_IS_elements   = replace_na(n_IS_elements, 0),
         defense_per_Mbp = n_defense / (Length / 1e6),
         IS_per_Mbp      = n_IS_elements / (Length / 1e6),
         dataset         = "MAG",
         genus_group     = if_else(genus %in% contam_genera, "Probable-contaminant genera", "Human-associated genera"))

mag_wirbel <- alistipes %>%
  left_join(read_tsv(file.path(WIRBEL, "genome_properties_wirbel.tsv"), show_col_types = FALSE) %>%
              select(source_sample, source_bin, genome_size_mbp),
            by = c("source_sample", "source_bin")) %>%
  left_join(count(df_wirbel, sample_id, bin_id, name = "n_defense"),
            by = c("source_sample" = "sample_id", "source_bin" = "bin_id")) %>%
  left_join(count(is_wirbel, sample_id, bin_id, name = "n_IS_elements"),
            by = c("source_sample" = "sample_id", "source_bin" = "bin_id")) %>%
  mutate(n_defense       = replace_na(n_defense, 0),
         n_IS_elements   = replace_na(n_IS_elements, 0),
         defense_per_Mbp = n_defense / genome_size_mbp,
         IS_per_Mbp      = n_IS_elements / genome_size_mbp,
         dataset         = "MAG",
         genus           = "Alistipes",
         genus_group     = "Wirbel genus (Alistipes)")
stopifnot(!anyNA(mag_wirbel$genome_size_mbp))

ref_wirbel <- ref_pool %>%
  filter(genus == "Alistipes", source_category == "human") %>%
  mutate(dataset = "Human reference", genus_group = "Wirbel genus (Alistipes)")
stopifnot(nrow(ref_wirbel) == 26)

comparison_long <- bind_rows(ref_hartwig, mag_hartwig, ref_wirbel, mag_wirbel) %>%
  select(genus, genus_group, dataset, defense_per_Mbp, IS_per_Mbp) %>%
  mutate(
    genus       = factor(genus, levels = c(human_genera, contam_genera, "Alistipes")),
    genus_group = factor(genus_group, levels = c("Human-associated genera", "Probable-contaminant genera",
                                                 "Wirbel genus (Alistipes)")),
    dataset     = factor(dataset, levels = c("MAG", "Human reference", "Non-human reference"))
  )

genus_summary <- comparison_long %>%
  group_by(genus_group, genus, dataset) %>%
  summarise(n = n(),
            median_defense_per_Mbp = median(defense_per_Mbp, na.rm = TRUE),
            median_IS_per_Mbp      = median(IS_per_Mbp, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(genus_group, genus, dataset)
save_data(genus_summary, OUT$figS2, "Figure_S2_genus_level_summary")

# Wilcoxon rank-sum test per genus, MAG vs that arm's reference group, BH within each arm
wilcox_p <- function(x, y) tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_)

compute_wilcox_for_arm <- function(df, ref_label) {
  df %>%
    group_by(genus_group, genus) %>%
    summarise(n_mag = sum(dataset == "MAG"), n_ref = sum(dataset == ref_label),
              p_defense = wilcox_p(defense_per_Mbp[dataset == "MAG"], defense_per_Mbp[dataset == ref_label]),
              p_IS      = wilcox_p(IS_per_Mbp[dataset == "MAG"],      IS_per_Mbp[dataset == ref_label]),
              .groups = "drop")
}

wilcox_results <- bind_rows(
  compute_wilcox_for_arm(comparison_long %>% filter(genus_group == "Human-associated genera"),     "Human reference"),
  compute_wilcox_for_arm(comparison_long %>% filter(genus_group == "Probable-contaminant genera"), "Non-human reference"),
  compute_wilcox_for_arm(comparison_long %>% filter(genus_group == "Wirbel genus (Alistipes)"),    "Human reference")
) %>%
  group_by(genus_group) %>%
  mutate(p_defense_BH = p.adjust(p_defense, method = "BH"),
         p_IS_BH      = p.adjust(p_IS,      method = "BH")) %>%
  ungroup()
save_data(wilcox_results, OUT$figS2, "Figure_S2_wilcoxon_tests")
print(wilcox_results, n = Inf)

# Human-associated arm: all three pairwise comparisons, each BH-corrected across the four genera
compute_wilcox_pair <- function(df, label_a, label_b) {
  df %>%
    group_by(genus) %>%
    summarise(p_defense = wilcox_p(defense_per_Mbp[dataset == label_a], defense_per_Mbp[dataset == label_b]),
              p_IS      = wilcox_p(IS_per_Mbp[dataset == label_a],      IS_per_Mbp[dataset == label_b]),
              .groups = "drop") %>%
    mutate(p_defense_BH = p.adjust(p_defense, method = "BH"),
           p_IS_BH      = p.adjust(p_IS,      method = "BH"))
}

human_arm_long      <- comparison_long %>% filter(genus_group == "Human-associated genera") %>% droplevels()
row1_mag_human      <- compute_wilcox_pair(human_arm_long, "MAG", "Human reference")
row1_mag_nonhuman   <- compute_wilcox_pair(human_arm_long, "MAG", "Non-human reference")
row1_human_nonhuman <- compute_wilcox_pair(human_arm_long, "Human reference", "Non-human reference")

save_data(bind_rows(
  row1_mag_human      %>% mutate(comparison = "MAG vs Human reference"),
  row1_mag_nonhuman   %>% mutate(comparison = "MAG vs Non-human reference"),
  row1_human_nonhuman %>% mutate(comparison = "Human reference vs Non-human reference")
), OUT$figS2, "Figure_S2_row1_pairwise_wilcoxon_tests")

sig_stars <- function(p) case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")

# All genera here are Latin names, so the whole axis is italic
GENUS_AXIS_THEME <- theme(axis.text.x = element_text(angle = 40, hjust = 1, face = "bold.italic", size = FS_AXIS_TEXT))

sig_brackets_3_at <- function(x_center, p_mag_human, p_human_nonhuman, p_mag_nonhuman,
                              y1, y2, dodge = 0.25, gap = 0.03, label_size = 2.5) {
  list(
    annotate("segment", x = x_center - dodge, xend = x_center - gap, y = y1, yend = y1,
             colour = "grey40", linewidth = 0.4),
    annotate("text", x = x_center - dodge / 2, y = y1, label = sig_stars(p_mag_human),
             vjust = -0.35, size = label_size, fontface = "bold", colour = "grey20"),
    annotate("segment", x = x_center + gap, xend = x_center + dodge, y = y1, yend = y1,
             colour = "grey40", linewidth = 0.4),
    annotate("text", x = x_center + dodge / 2, y = y1, label = sig_stars(p_human_nonhuman),
             vjust = -0.35, size = label_size, fontface = "bold", colour = "grey20"),
    annotate("segment", x = x_center - dodge, xend = x_center + dodge, y = y2, yend = y2,
             colour = "grey40", linewidth = 0.4),
    annotate("text", x = x_center, y = y2, label = sig_stars(p_mag_nonhuman),
             vjust = -0.35, size = label_size, fontface = "bold", colour = "grey20")
  )
}

# Two-box genera (no non-human reference group): ggplot dodges the two boxes to +/- 0.1875
sig_bracket_2_at <- function(x_center, p_mag_human, y1, dodge = 0.1875, label_size = 2.5) {
  list(
    annotate("segment", x = x_center - dodge, xend = x_center + dodge, y = y1, yend = y1,
             colour = "grey40", linewidth = 0.4),
    annotate("text", x = x_center, y = y1, label = sig_stars(p_mag_human),
             vjust = -0.35, size = label_size, fontface = "bold", colour = "grey20")
  )
}

burden_theme <- function() {
  theme_bw(base_size = 11) +
    GENUS_AXIS_THEME +
    theme(plot.title  = element_text(face = "bold", size = FS_TITLE),
          axis.title  = element_text(size = FS_AXIS_TITLE),
          axis.text.y = element_text(size = FS_AXIS_TEXT))
}

# Human-associated arm: three boxes per genus (two for Fusobacterium), all pairwise brackets
make_panel_3way <- function(data_long, group_name, metric_col, y_lab, title,
                            y_trans = "identity", y_breaks = waiver(), y_expand_top = 0.5) {
  d <- data_long %>% filter(genus_group == group_name) %>% droplevels()
  genus_levels <- levels(d$genus)
  genus_max <- d %>% group_by(genus) %>% summarise(genus_max = max(.data[[metric_col]], na.rm = TRUE), .groups = "drop")
  n_labels <- d %>% count(genus, dataset) %>% left_join(genus_max, by = "genus") %>% mutate(y = genus_max * 1.05)
  p_col <- if (metric_col == "defense_per_Mbp") "p_defense_BH" else "p_IS_BH"

  tier_y <- function(gm) {
    if (identical(y_trans, "identity")) return(c(y1 = gm * 1.18, y2 = gm * 1.36))
    t <- y_trans$transform(gm)
    c(y1 = y_trans$inverse(t + 0.22), y2 = y_trans$inverse(t + 0.42))
  }
  p_at <- function(tbl, g) tbl[[p_col]][tbl$genus == g]

  brackets <- map(seq_along(genus_levels), function(i) {
    g  <- genus_levels[i]
    ys <- tier_y(genus_max$genus_max[genus_max$genus == g])
    if ("Non-human reference" %in% d$dataset[d$genus == g]) {
      sig_brackets_3_at(i, p_at(row1_mag_human, g), p_at(row1_human_nonhuman, g), p_at(row1_mag_nonhuman, g),
                        y1 = ys[["y1"]], y2 = ys[["y2"]])
    } else {
      sig_bracket_2_at(i, p_at(row1_mag_human, g), y1 = ys[["y1"]])
    }
  })

  ggplot(d, aes(x = genus, y = .data[[metric_col]], fill = dataset)) +
    geom_boxplot(outlier.size = 0.6, position = position_dodge(width = 0.75), width = 0.6) +
    geom_text(data = n_labels, aes(y = y, label = paste0("n=", n), group = dataset),
              position = position_dodge(width = 0.75), size = 2.6, vjust = 0) +
    brackets +
    scale_y_continuous(trans = y_trans, breaks = y_breaks, expand = expansion(mult = c(0.02, y_expand_top))) +
    scale_fill_manual(values = dataset_colours, drop = FALSE) +
    labs(x = NULL, y = y_lab, fill = NULL, title = title) +
    burden_theme()
}

# Contaminant and Alistipes arms: two boxes per genus, one star
make_panel <- function(data_long, group_name, metric_col, y_lab, title,
                       y_trans = "identity", y_breaks = waiver(), y_expand_top = 0.32) {
  d <- data_long %>% filter(genus_group == group_name) %>% droplevels()
  genus_max <- d %>% group_by(genus) %>% summarise(genus_max = max(.data[[metric_col]], na.rm = TRUE), .groups = "drop")
  n_labels <- d %>% count(genus, dataset) %>% left_join(genus_max, by = "genus") %>% mutate(y = genus_max * 1.08)
  p_col <- if (metric_col == "defense_per_Mbp") "p_defense_BH" else "p_IS_BH"
  sig_labels <- wilcox_results %>% filter(genus_group == group_name) %>%
    mutate(label = sig_stars(.data[[p_col]])) %>%
    left_join(genus_max, by = "genus") %>% mutate(y = genus_max * 1.22)

  ggplot(d, aes(x = genus, y = .data[[metric_col]], fill = dataset)) +
    geom_boxplot(outlier.size = 0.6, position = position_dodge(width = 0.75), width = 0.6) +
    geom_text(data = n_labels, aes(y = y, label = paste0("n=", n), group = dataset),
              position = position_dodge(width = 0.75), size = 2.6, vjust = 0) +
    geom_text(data = sig_labels, aes(x = genus, y = y, label = label), inherit.aes = FALSE, size = 4, fontface = "bold") +
    scale_y_continuous(trans = y_trans, breaks = y_breaks, expand = expansion(mult = c(0.02, y_expand_top))) +
    scale_fill_manual(values = dataset_colours, drop = FALSE) +
    labs(x = NULL, y = y_lab, fill = NULL, title = title) +
    burden_theme()
}

# IS panels use a log-like scale that is defined at zero
IS_Y_TRANS  <- scales::pseudo_log_trans(sigma = 1, base = 10)
IS_Y_LAB    <- "IS elements / Mbp (log-like scale)"
DEF_Y_LAB   <- "Defence systems / Mbp"

burden_panels <- list(
  A = make_panel_3way(comparison_long, "Human-associated genera", "defense_per_Mbp", DEF_Y_LAB,
                      "Hartwig human-associated genera - defence"),
  B = make_panel_3way(comparison_long, "Human-associated genera", "IS_per_Mbp", IS_Y_LAB,
                      "Hartwig human-associated genera - IS elements", y_trans = IS_Y_TRANS,
                      y_breaks = c(0, 1, 2, 5, 10, 20, 50, 100), y_expand_top = 0.28),
  C = make_panel(comparison_long, "Probable-contaminant genera", "defense_per_Mbp", DEF_Y_LAB,
                 "Hartwig probable-contaminant genera - defence"),
  D = make_panel(comparison_long, "Probable-contaminant genera", "IS_per_Mbp", IS_Y_LAB,
                 "Hartwig probable-contaminant genera - IS elements", y_trans = IS_Y_TRANS,
                 y_breaks = c(0, 1, 2, 5, 10, 20, 50, 100), y_expand_top = 0.10),
  E = make_panel(comparison_long, "Wirbel genus (Alistipes)", "defense_per_Mbp", DEF_Y_LAB,
                 "Wirbel CRC, Alistipes - defence"),
  F = make_panel(comparison_long, "Wirbel genus (Alistipes)", "IS_per_Mbp", IS_Y_LAB,
                 "Wirbel CRC, Alistipes - IS elements", y_trans = IS_Y_TRANS,
                 y_breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200, 300))
)

burden_files <- c(A = "human_associated_defence", B = "human_associated_IS",
                  C = "contaminant_defence",      D = "contaminant_IS",
                  E = "wirbel_alistipes_defence", F = "wirbel_alistipes_IS")
for (tag in names(burden_panels)) {
  save_fig(burden_panels[[tag]], OUT$figS2, paste0("Figure_S2", tag, "_", burden_files[[tag]]),
           w = if (tag %in% c("A", "B")) 11 else 8.5, h = 5.5)
}

ROW_HEIGHT_S2 <- 5.5
figure_S2 <- (wrap_plots(imap(burden_panels, ~ .x + labs(tag = .y) + tag_theme + theme(legend.position = "none")),
                         ncol = 2) /
                wrap_elements(full = build_legend_strip(dataset_colours))) +
  plot_layout(heights = c(3 * ROW_HEIGHT_S2, LEGEND_STRIP_IN))
save_fig(figure_S2, OUT$figS2, "Figure_S2", w = 20, h = 3 * ROW_HEIGHT_S2 + LEGEND_STRIP_IN)

# =============================================================================
# 2. Composition: top ten systems per genus, MAG vs reference (Figures 6 and S3)
# =============================================================================

# ---- Figure 6: human-associated genera and Wirbel Alistipes ----------------------------------

# MAG on the left; human reference and (where at least 10 exist) non-human reference dodged
# within the same row on the right. Rows are the genus's top ten by MAG prevalence only.
fig6_genera <- c(human_genera, "Alistipes")

ref_meta6 <- ref_df_summary %>%
  select(genus, accession, source_category) %>%
  filter(genus %in% fig6_genera, source_category %in% c("human", "environmental")) %>%
  mutate(ref_group = if_else(source_category == "human", "Human", "NonHuman"))
ref_df6 <- ref_df_detail %>% filter(accession %in% ref_meta6$accession) %>% left_join(ref_meta6, by = c("genus", "accession"))
ref_is6 <- ref_is_detail %>% filter(accession %in% ref_meta6$accession) %>% left_join(ref_meta6, by = c("genus", "accession"))

mag_ids6 <- c(map(human_genera, ~ gtdb %>% filter(Genus == .x) %>% pull(mag_id)) %>% set_names(human_genera),
              list(Alistipes = alistipes$mag_id))

denom6 <- tibble(genus = fig6_genera, n_mag = map_int(mag_ids6[fig6_genera], length)) %>%
  left_join(ref_meta6 %>% count(genus, ref_group) %>% pivot_wider(names_from = ref_group, values_from = n, values_fill = 0),
            by = "genus") %>%
  mutate(NonHuman = replace_na(NonHuman, 0), Human = replace_na(Human, 0),
         NonHuman_shown = if_else(NonHuman >= GENUS_MIN_N, NonHuman, 0L))
print(denom6)

prevalence <- function(detail, id_col, type_col, n) {
  detail %>% distinct(.data[[id_col]], .data[[type_col]]) %>%
    count(.data[[type_col]], name = "n_with") %>%
    mutate(prev = n_with / n * 100) %>% select(-n_with)
}

build_prevalence6 <- function(genus, type_col) {
  dn <- denom6 %>% filter(genus == !!genus)
  mag_detail <- if (genus == "Alistipes") {
    if (type_col == "subtype") df_wirbel else is_wirbel
  } else {
    if (type_col == "subtype") df_hartwig else is_hartwig
  }
  ref_detail <- if (type_col == "subtype") ref_df6 else ref_is6

  wide <- full_join(
    prevalence(filter(mag_detail, mag_id %in% mag_ids6[[genus]]), "mag_id", type_col, dn$n_mag) %>%
      rename(prev_MAG = prev),
    prevalence(filter(ref_detail, genus == !!genus, ref_group == "Human"), "accession", type_col, dn$Human) %>%
      rename(prev_Human = prev),
    by = type_col) %>%
    mutate(across(c(prev_MAG, prev_Human), ~ replace_na(.x, 0)))

  has_nh <- dn$NonHuman_shown > 0
  if (has_nh) {
    wide <- full_join(wide,
                      prevalence(filter(ref_detail, genus == !!genus, ref_group == "NonHuman"), "accession",
                                 type_col, dn$NonHuman_shown) %>% rename(prev_NonHuman = prev),
                      by = type_col) %>%
      mutate(prev_NonHuman = replace_na(prev_NonHuman, 0))
  } else {
    wide$prev_NonHuman <- NA_real_
  }

  top_of <- function(col) wide %>% slice_max(.data[[col]], n = TOP_N, with_ties = FALSE) %>% pull(.data[[type_col]])
  top_mag <- top_of("prev_MAG")

  list(data = wide %>% filter(.data[[type_col]] %in% top_mag) %>% rename(type = all_of(type_col)) %>%
         arrange(desc(prev_MAG)),
       n_mag = dn$n_mag, n_human = dn$Human, n_nonhuman = dn$NonHuman_shown, genus = genus,
       overlap_mag_human    = length(intersect(top_mag, top_of("prev_Human"))),
       overlap_mag_nonhuman = if (has_nh) length(intersect(top_mag, top_of("prev_NonHuman"))) else NA_integer_)
}

colours6      <- c("MAG" = "#0072B2", "Human" = "#EBC9DC", "NonHuman" = "#E49659")
colours6_text <- c("MAG" = "#0072B2", "Human" = "#A8447F", "NonHuman" = "#B54900")

make_panel6 <- function(res, label_fn, y_lab, tag) {
  d <- res$data
  n_types <- nrow(d)
  d <- d %>% mutate(row = row_number(), display_name = factor(label_fn(type), levels = rev(label_fn(type))))
  has_nh <- res$n_nonhuman > 0
  offset <- if (has_nh) 0.19 else 0

  long <- bind_rows(
    d %>% transmute(row, x = -prev_MAG,  grp = "MAG",   y_off = 0,      w = 0.72),
    d %>% transmute(row, x =  prev_Human, grp = "Human", y_off = offset, w = if (has_nh) 0.34 else 0.72),
    if (has_nh) d %>% transmute(row, x = prev_NonHuman, grp = "NonHuman", y_off = -offset, w = 0.34) else NULL
  ) %>%
    mutate(y = (n_types - row + 1) + y_off, grp = factor(grp, levels = c("MAG", "Human", "NonHuman")))
  max_x <- max(abs(long$x), na.rm = TRUE)

  ref_labels <- if (has_nh) list(
    annotate("text", x = max_x / 2, y = n_types + 1.15, hjust = 0.5, size = 2.7, fontface = "bold",
             colour = colours6_text[["Human"]], label = sprintf("Human (n=%d)", res$n_human)),
    annotate("text", x = max_x / 2, y = n_types + 0.65, hjust = 0.5, size = 2.7, fontface = "bold",
             colour = colours6_text[["NonHuman"]], label = sprintf("Non-human (n=%d) ->", res$n_nonhuman))
  ) else
    annotate("text", x = max_x / 2, y = n_types + 0.9, hjust = 0.5, size = 2.7, fontface = "bold",
             colour = colours6_text[["Human"]], label = sprintf("Human\n(n=%d) ->", res$n_human))

  ggplot(long, aes(x = x, y = y, fill = grp, width = w)) +
    geom_col(orientation = "y", colour = NA) +
    geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.5) +
    scale_y_continuous(breaks = n_types:1, labels = d$display_name[order(d$row)],
                       limits = c(0.3, n_types + 1.3), expand = c(0, 0)) +
    scale_x_continuous(breaks = c(-100, -50, 0, 50, 100), labels = function(x) paste0(abs(x), "%"),
                       limits = c(-105, 105), expand = c(0, 0)) +
    scale_fill_manual(values = colours6, name = NULL, drop = FALSE) +
    annotate("text", x = -(max_x / 2), y = n_types + 0.9, hjust = 0.5, size = 2.7, fontface = "bold",
             colour = colours6_text[["MAG"]], label = sprintf("<- MAG\n(n=%d)", res$n_mag)) +
    ref_labels +
    labs(title = res$genus, x = "Prevalence among genomes (%)", y = y_lab, tag = tag) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          axis.text.y = element_text(size = FS_AXIS_TEXT),
          axis.text.x = element_text(size = FS_AXIS_TEXT),
          axis.title  = element_text(size = FS_AXIS_TITLE),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(size = FS_TITLE, face = "bold.italic"),
          plot.margin = margin(t = 22, r = 6, b = 6, l = 6)) +
    tag_theme
}

res6_defence <- map(fig6_genera, ~ build_prevalence6(.x, "subtype")) %>% set_names(fig6_genera)
res6_IS      <- map(fig6_genera, ~ build_prevalence6(.x, "family"))  %>% set_names(fig6_genera)

overlap_summary6 <- tibble(
  genus      = fig6_genera,
  n_mag      = map_dbl(res6_defence, "n_mag"),
  n_human    = map_dbl(res6_defence, "n_human"),
  n_nonhuman = map_dbl(res6_defence, "n_nonhuman"),
  defence_human_overlap    = sprintf("%d/%d", map_dbl(res6_defence, "overlap_mag_human"), TOP_N),
  defence_nonhuman_overlap = ifelse(n_nonhuman > 0, sprintf("%d/%d", map_dbl(res6_defence, "overlap_mag_nonhuman"), TOP_N), "n/a"),
  IS_human_overlap         = sprintf("%d/%d", map_dbl(res6_IS, "overlap_mag_human"), TOP_N),
  IS_nonhuman_overlap      = ifelse(n_nonhuman > 0, sprintf("%d/%d", map_dbl(res6_IS, "overlap_mag_nonhuman"), TOP_N), "n/a")
)
save_data(overlap_summary6, OUT$fig6, "Figure_6_overlap_summary")
print(overlap_summary6)

# Row 1 (A-E) defence, row 2 (F-J) IS, one column per genus
tags6_defence <- LETTERS[seq_along(fig6_genera)]
tags6_IS      <- LETTERS[length(fig6_genera) + seq_along(fig6_genera)]
panels6_defence <- map2(res6_defence, tags6_defence, ~ make_panel6(.x, sys_short, "Defence subtype", .y))
panels6_IS      <- map2(res6_IS,      tags6_IS,      ~ make_panel6(.x, identity,  "IS family",       .y))

for (i in seq_along(fig6_genera)) {
  g <- fig6_genera[i]
  save_fig(panels6_defence[[g]], OUT$fig6, sprintf("Figure_6%s_%s_defence", tags6_defence[i], g), w = 4.6, h = 4.3)
  save_fig(panels6_IS[[g]],      OUT$fig6, sprintf("Figure_6%s_%s_IS",      tags6_IS[i],      g), w = 4.6, h = 4.3)
  save_data(res6_defence[[g]]$data, OUT$fig6, sprintf("Figure_6%s_%s_defence_data", tags6_defence[i], g))
  save_data(res6_IS[[g]]$data,      OUT$fig6, sprintf("Figure_6%s_%s_IS_data",      tags6_IS[i],      g))
}

ROW_HEIGHT_6 <- 6
figure_6 <- (wrap_plots(c(unname(panels6_defence), unname(panels6_IS)), nrow = 2) /
               wrap_elements(full = build_legend_strip(dataset_colours))) +
  plot_layout(heights = c(2 * ROW_HEIGHT_6, LEGEND_STRIP_IN))
save_fig(figure_6, OUT$fig6, "Figure_6", w = 20.5, h = 2 * ROW_HEIGHT_6 + LEGEND_STRIP_IN)

# ---- Figure S3: probable-contaminant genera ----------------------------------------------------

# Contaminant genera are compared with their non-human reference genomes only, so each panel
# is a plain two-sided butterfly. Rows are the genus's top ten by MAG prevalence; a genus with
# fewer than ten IS families (Ralstonia, 8) shows fewer rows.
colours_S3 <- c("MAG" = "#0072B2", "Non-human reference" = "#E49659")

mag_lookup_S3 <- gtdb_contam %>% filter(Genus %in% contam_genera) %>% select(mag_id, genus = Genus)
n_mag_S3 <- count(mag_lookup_S3, genus, name = "n_mag")
stopifnot(sum(n_mag_S3$n_mag) == sum(contam_genus_n$n_mag))

acc_env_S3 <- ref_df_summary %>% select(genus, accession, source_category) %>%
  filter(genus %in% contam_genera, source_category == "environmental")
n_ref_S3 <- count(acc_env_S3, genus, name = "n_ref")

df_mag_S3 <- df_hartwig %>% filter(mag_id %in% mag_lookup_S3$mag_id) %>% left_join(mag_lookup_S3, by = "mag_id")
is_mag_S3 <- is_hartwig %>% filter(mag_id %in% mag_lookup_S3$mag_id) %>% left_join(mag_lookup_S3, by = "mag_id") %>%
  rename(subtype = family)
df_ref_S3 <- ref_df_detail %>% filter(accession %in% acc_env_S3$accession)
is_ref_S3 <- ref_is_detail %>% filter(accession %in% acc_env_S3$accession) %>% rename(subtype = family)

compare_S3 <- function(mag_detail, ref_detail, genus_name, n_mag, n_ref) {
  mag_prev <- mag_detail %>% filter(genus == genus_name) %>%
    distinct(mag_id, subtype) %>% count(subtype, name = "n_mag_with") %>%
    mutate(mag_pct = round(n_mag_with / n_mag * 100, 1))
  ref_prev <- ref_detail %>% filter(genus == genus_name) %>%
    distinct(accession, subtype) %>% count(subtype, name = "n_ref_with") %>%
    mutate(ref_pct = round(n_ref_with / n_ref * 100, 1))

  top_mag <- mag_prev %>% slice_max(mag_pct, n = TOP_N, with_ties = FALSE) %>% pull(subtype)
  top_ref <- ref_prev %>% slice_max(ref_pct, n = TOP_N, with_ties = FALSE) %>% pull(subtype)

  data <- tibble(subtype = top_mag) %>%
    left_join(mag_prev, by = "subtype") %>%
    left_join(ref_prev, by = "subtype") %>%
    mutate(across(c(n_mag_with, mag_pct, n_ref_with, ref_pct), ~ replace_na(.x, 0)),
           prev_mag = mag_pct, prev_ref = ref_pct, genus = genus_name) %>%
    arrange(desc(prev_mag))
  list(data = data, overlap_n = length(intersect(top_mag, top_ref)))
}

make_panel_S3 <- function(plot_df, label_fn, n_mag, n_ref, title, x_lab, y_lab, tag) {
  mirror <- plot_df %>%
    mutate(display_name = factor(label_fn(subtype), levels = label_fn(subtype))) %>%
    pivot_longer(cols = c(prev_mag, prev_ref), names_to = "grp", values_to = "prevalence") %>%
    mutate(grp   = factor(if_else(grp == "prev_mag", "MAG", "Non-human reference"), levels = names(colours_S3)),
           x_val = if_else(grp == "MAG", -prevalence, prevalence))
  n_types <- nrow(plot_df)
  max_x   <- max(mirror$prevalence)

  ggplot(mirror, aes(x = x_val, y = display_name, fill = grp)) +
    geom_col(width = 0.70, colour = NA, show.legend = FALSE) +
    geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.55) +
    geom_text(data = mirror %>% filter(abs(x_val) > 3),
              aes(label = sprintf("%.0f%%", abs(x_val)), hjust = if_else(x_val < 0, 1.18, -0.18)),
              size = 2.6, colour = "grey20") +
    annotate("text", x = -(max_x / 2), y = n_types + 1.1, label = sprintf("<- MAG  (n = %d)", n_mag),
             hjust = 0.5, size = 3.0, colour = "#0072B2", fontface = "bold") +
    annotate("text", x = max_x / 2, y = n_types + 1.1, label = sprintf("Non-human reference  (n = %d) ->", n_ref),
             hjust = 0.5, size = 3.0, colour = "#B54900", fontface = "bold") +
    scale_fill_manual(values = colours_S3, name = NULL, drop = FALSE) +
    scale_x_continuous(breaks = c(-100, -50, 0, 50, 100), labels = function(x) paste0(abs(x), "%"),
                       limits = c(-105, 105), expand = c(0, 0)) +
    coord_cartesian(ylim = c(0.5, n_types + 1.8), clip = "off") +
    labs(title = title, x = x_lab, y = y_lab, tag = tag) +
    theme_bw(base_size = 11) +
    theme(axis.text.y        = element_text(size = FS_AXIS_TEXT),
          axis.title         = element_text(size = FS_AXIS_TITLE),
          panel.grid.major.y = element_blank(),
          plot.title         = element_text(size = FS_TITLE, face = "bold.italic"),
          plot.margin        = margin(t = 25, r = 10, b = 5, l = 5)) +
    tag_theme
}

panels_S3   <- list()
row_heights <- numeric(0)
summary_S3  <- list()
for (i in seq_along(contam_genera)) {
  g <- contam_genera[i]
  n_mag_g <- n_mag_S3$n_mag[n_mag_S3$genus == g]
  n_ref_g <- n_ref_S3$n_ref[n_ref_S3$genus == g]
  tag_def <- LETTERS[2 * i - 1]; tag_is <- LETTERS[2 * i]

  def <- compare_S3(df_mag_S3, df_ref_S3, g, n_mag_g, n_ref_g)
  isf <- compare_S3(is_mag_S3, is_ref_S3, g, n_mag_g, n_ref_g)
  p_def <- make_panel_S3(def$data, sys_short, n_mag_g, n_ref_g, g, "% of MAGs / genomes carrying subtype",
                         "Defence subtype", tag_def)
  p_is  <- make_panel_S3(isf$data, identity,  n_mag_g, n_ref_g, g, "% of MAGs / genomes carrying family",
                         "IS family", tag_is)

  save_fig(p_def, OUT$figS3, sprintf("Figure_S3%s_%s_defence", tag_def, g), w = 8.5, h = 7 + nrow(def$data) * 0.14)
  save_fig(p_is,  OUT$figS3, sprintf("Figure_S3%s_%s_IS",      tag_is,  g), w = 8.5, h = 7 + nrow(isf$data) * 0.14)
  save_data(def$data, OUT$figS3, sprintf("Figure_S3%s_%s_defence_data", tag_def, g))
  save_data(isf$data, OUT$figS3, sprintf("Figure_S3%s_%s_IS_data",      tag_is,  g))

  panels_S3 <- c(panels_S3, list(p_def, p_is))
  row_heights <- c(row_heights, 4 + max(nrow(def$data), nrow(isf$data)) * 0.11)
  summary_S3[[i]] <- tibble(genus = g, n_mag = n_mag_g, n_ref = n_ref_g,
                            defence_overlap = sprintf("%d/%d", def$overlap_n, TOP_N),
                            IS_overlap      = sprintf("%d/%d", isf$overlap_n, TOP_N))
}

figure_S3 <- (wrap_plots(panels_S3, ncol = 2, heights = row_heights) /
                wrap_elements(full = build_legend_strip(colours_S3))) +
  plot_layout(heights = c(sum(row_heights), LEGEND_STRIP_IN))
save_fig(figure_S3, OUT$figS3, "Figure_S3", w = 18, h = sum(row_heights) + LEGEND_STRIP_IN)

summary_S3 <- bind_rows(summary_S3)
save_data(summary_S3, OUT$figS3, "Figure_S3_overlap_summary")
print(summary_S3)
