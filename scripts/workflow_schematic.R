# Figure 1: graphical abstract and study workflow.
# Builds the figure as SVG and converts it to PDF/PNG with rsvg-convert (brew install librsvg).
# Input:  MAG classification from 00_hartwig_mag_classification.R (for the retained/excluded counts)
# Output: output/figures/figure_1/, output/figure_data/figure_1/

library(tidyverse)
library(systemfonts)

# Run from the repository root: Rscript scripts/workflow_schematic.R
stopifnot("Run this script from the repository root" = dir.exists("data"))
FIG_DIR  <- file.path("output", "figures",     "figure_1")
DATA_DIR <- file.path("output", "figure_data", "figure_1")
dir.create(FIG_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

stopifnot("rsvg-convert not found on PATH; install it with `brew install librsvg`" =
            nchar(Sys.which("rsvg-convert")) > 0)

mag_class  <- read_tsv(file.path("output", "mag_classification", "hartwig_mag_classification.tsv"),
                       show_col_types = FALSE)
n_total    <- nrow(mag_class)
n_retained <- sum(mag_class$classification == "Human-associated")
n_excluded <- sum(mag_class$classification == "Probable contaminant")
stopifnot(n_total == 514, n_retained == 383, n_excluded == 131)

# Icons from Lucide (https://lucide.dev), ISC licence.
LUCIDE <- list(
  dna = '<path d="m10 16 1.5 1.5"/><path d="m14 8-1.5-1.5"/><path d="M15 2c-1.798 1.998-2.518 3.995-2.807 5.993"/><path d="m16.5 10.5 1 1"/><path d="m17 6-2.891-2.891"/><path d="M2 15c6.667-6 13.333 0 20-6"/><path d="m20 9 .891.891"/><path d="M3.109 14.109 4 15"/><path d="m6.5 12.5 1 1"/><path d="m7 18 2.891 2.891"/><path d="M9 22c1.798-1.998 2.518-3.995 2.807-5.993"/>',
  network = '<rect x="16" y="16" width="6" height="6" rx="1"/><rect x="2" y="16" width="6" height="6" rx="1"/><rect x="9" y="2" width="6" height="6" rx="1"/><path d="M5 16v-3a1 1 0 0 1 1-1h12a1 1 0 0 1 1 1v3"/><path d="M12 12V8"/>',
  virus = '<path d="M12 14.991h.01"/><path d="M12 22v-3"/><path d="M12 2v3"/><path d="M13 22h-2"/><path d="M13 2h-2"/><path d="M13.99 10H14"/><path d="m16.5 19.794-1-1.733"/><path d="m16.5 4.205-1 1.732"/><path d="m19.794 16.5-1.732-1"/><path d="m19.794 7.5-1.732 1"/><path d="M2 12h3"/><path d="M2 13v-2"/><path d="M22 12h-3"/><path d="M22 13v-2"/><path d="m4.206 16.5 1.732-1"/><path d="m4.206 7.5 1.732 1"/><path d="m7.5 19.794 1-1.733"/><path d="m7.5 4.205 1 1.732"/><path d="M9 12h.01"/><circle cx="12" cy="12" r="7"/>',
  "scan-search" = '<path d="M3 7V5a2 2 0 0 1 2-2h2"/><path d="M17 3h2a2 2 0 0 1 2 2v2"/><path d="M21 17v2a2 2 0 0 1-2 2h-2"/><path d="M7 21H5a2 2 0 0 1-2-2v-2"/><circle cx="12" cy="12" r="3"/><path d="m16 16-1.9-1.9"/>',
  shield = '<path d="M20 13c0 5-3.5 7.5-7.66 8.95a1 1 0 0 1-.67-.01C7.5 20.5 4 18 4 13V6a1 1 0 0 1 1-1c2 0 4.5-1.2 6.24-2.72a1.17 1.17 0 0 1 1.52 0C14.51 3.81 17 5 19 5a1 1 0 0 1 1 1z"/>',
  puzzle = '<path d="M15.39 4.39a1 1 0 0 0 1.68-.474 2.5 2.5 0 1 1 3.014 3.015 1 1 0 0 0-.474 1.68l1.683 1.682a2.414 2.414 0 0 1 0 3.414L19.61 15.39a1 1 0 0 1-1.68-.474 2.5 2.5 0 1 0-3.014 3.015 1 1 0 0 1 .474 1.68l-1.683 1.682a2.414 2.414 0 0 1-3.414 0L8.61 19.61a1 1 0 0 0-1.68.474 2.5 2.5 0 1 1-3.014-3.015 1 1 0 0 0 .474-1.68l-1.683-1.682a2.414 2.414 0 0 1 0-3.414L4.39 8.61a1 1 0 0 1 1.68.474 2.5 2.5 0 1 0 3.014-3.015 1 1 0 0 1-.474-1.68l1.683-1.682a2.414 2.414 0 0 1 3.414 0z"/>',
  refresh_cw = '<path d="M3 12a9 9 0 0 1 9-9 9.75 9.75 0 0 1 6.74 2.74L21 8"/><path d="M21 3v5h-5"/><path d="M21 12a9 9 0 0 1-9 9 9.75 9.75 0 0 1-6.74-2.74L3 16"/><path d="M8 16H3v5"/>',
  scale = '<path d="M12 3v18"/><path d="m19 8 3 8a5 5 0 0 1-6 0zV7"/><path d="M3 7h1a17 17 0 0 0 8-2 17 17 0 0 0 8 2h1"/><path d="m5 8 3 8a5 5 0 0 1-6 0zV7"/><path d="M7 21h10"/>',
  shield_check = '<path d="M20 13c0 5-3.5 7.5-7.66 8.95a1 1 0 0 1-.67-.01C7.5 20.5 4 18 4 13V6a1 1 0 0 1 1-1c2 0 4.5-1.2 6.24-2.72a1.17 1.17 0 0 1 1.52 0C14.51 3.81 17 5 19 5a1 1 0 0 1 1 1z"/><path d="m9 12 2 2 4-4"/>',
  leaf = '<path d="M11 20a10 10 0 0010-10 25.9 25.9 0 00-1.04-7.281 1 1 0 00-1.755-.325C15.833 5.5 13 5.5 9.8 6.1A7 7 0 0011 20"/><path d="M2 21a5 5 0 012.911-4.544C7.613 15.212 8.351 15.24 11 13"/>',
  users = '<path d="M16 21v-2a4 4 0 0 0-4-4H6a4 4 0 0 0-4 4v2"/><path d="M16 3.128a4 4 0 0 1 0 7.744"/><path d="M22 21v-2a4 4 0 0 0-3-3.87"/><circle cx="9" cy="7" r="4"/>',
  merge = '<path d="m8 6 4-4 4 4"/><path d="M12 2v10.3a4 4 0 0 1-1.172 2.872L4 22"/><path d="m20 22-5-5"/>',
  split = '<path d="M16 3h5v5"/><path d="M8 3H3v5"/><path d="M12 22v-8.3a4 4 0 0 0-1.172-2.872L3 3"/><path d="m15 9 6-6"/>'
)

lucide_icon <- function(name, cx, cy, size = 22, colour = "#0072B2", stroke_w = 1.8, rotate180 = FALSE) {
  scale_f <- size / 24
  tx <- cx - size / 2
  ty <- cy - size / 2
  rot <- if (rotate180) ' rotate(180,12,12)' else ''
  sprintf('<g transform="translate(%.2f,%.2f) scale(%.4f)%s" fill="none" stroke="%s" stroke-width="%.2f" stroke-linecap="round" stroke-linejoin="round">%s</g>',
          tx, ty, scale_f, rot, colour, stroke_w / scale_f, LUCIDE[[name]])
}

BADGE_R <- 14
badge <- function(x, y, n, fill = "#0B2545", text_col = "#ffffff") {
  cx <- x + 18; cy <- y
  sprintf('<circle cx="%.1f" cy="%.1f" r="%d" fill="%s"/><text x="%.1f" y="%.1f" fill="%s" font-size="14" text-anchor="middle" font-weight="bold">%s</text>',
          cx, cy, BADGE_R, fill, cx, cy + 5, text_col, n)
}

icon_swatch <- function(x, y, icon_name, swatch_fill = "#EAF2FA", icon_colour = "#0072B2", shape = "square", size = 34, rotate180 = FALSE) {
  cx <- x + size / 2; cy <- y + size / 2
  bg <- if (shape == "circle")
    sprintf('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s"/>', cx, cy, size / 2, swatch_fill)
  else
    sprintf('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="7" fill="%s"/>', x, y, size, size, swatch_fill)
  paste0(bg, lucide_icon(icon_name, cx, cy, size = size * 0.62, colour = icon_colour, rotate180 = rotate180))
}

FONT_FAMILY <- "Helvetica"
stopifnot("Helvetica not resolvable by systemfonts on this machine" =
            length(systemfonts::match_fonts(FONT_FAMILY)$path) > 0)
# Measured at 600 dpi: at 72 dpi systemfonts fails with a FreeType error below ~8 pt.
MEASURE_RES  <- 600
UNITS_PER_IN <- 100
PT_PER_UNIT  <- 72 / UNITS_PER_IN

text_width_units <- function(txt, size_units, weight = "normal", italic = FALSE) {
  pt <- size_units * PT_PER_UNIT
  px <- systemfonts::string_width(txt, family = FONT_FAMILY, size = pt,
                                   res = MEASURE_RES, weight = weight, italic = italic)
  (px / MEASURE_RES) * UNITS_PER_IN
}

ital <- function(x) { attr(x, "italic") <- TRUE; x }
is_ital <- function(x) isTRUE(attr(x, "italic"))

LEFT_PAD    <- 15
ICON_SIZE   <- 34
ICON_GAP    <- 12
RIGHT_PAD   <- 20
ICON_TOP_PAD <- 18
SUB_GAP1    <- 15
LINE_H      <- 15.5
BOTTOM_PAD  <- 16

stopifnot("badge would overlap icon swatch" = ICON_TOP_PAD > BADGE_R)

box_size <- function(title, sublines, title_size = 15, sub_size = 11.5) {
  title_w <- text_width_units(title, title_size, weight = "bold")
  sub_ws  <- vapply(sublines, function(s) text_width_units(s, sub_size, italic = is_ital(s)), numeric(1))
  content_w <- max(LEFT_PAD + ICON_SIZE + ICON_GAP + title_w,
                    LEFT_PAD + if (length(sub_ws)) max(sub_ws) else 0)
  w <- content_w + RIGHT_PAD
  h <- ICON_TOP_PAD + ICON_SIZE + SUB_GAP1 + max(0, length(sublines) - 1) * LINE_H + BOTTOM_PAD
  list(w = w, h = h)
}

content <- list(
  input    = list(title = "Input metagenomes", icon = "dna",
                   sub = list("Battaglia et al. 2024 (Cell) — 4,160 tumour biopsies", "→ 514 high- and medium-quality MAGs")),
  taxonomy = list(title = "Taxonomic classification", icon = "network",
                   sub = list("GTDB-Tk v2.3.0 (R214)", ital("primary classification trigger"))),
  viral    = list(title = "Viral sequence identification", icon = "virus",
                   sub = list("geNomad v1.9.7", "prophage + extrachromosomal, all 514 MAGs")),
  iphop    = list(title = "Phage-host prediction", icon = "scan-search",
                   sub = list("iPHoP v1.4.2", ital("corroborating evidence"))),
  checkv   = list(title = "CheckV v1.0.1", icon = "shield_check", title_size = 13, sub_size = 9.5,
                   sub = list(ital("≥50% completeness, low contamination"), ital("(not a classification input)"))),
  classify = list(title = "Classification decision", icon = "merge", icon_rot = TRUE,
                   sub = list("GTDB-Tk screen + iPHoP evidence", "final call by literature review")),
  pharokka = list(title = "Pharokka v1.9.1", icon = "refresh_cw", title_size = 13, sub_size = 10,
                   sub = list(ital("lifecycle: temperate vs lytic"), ital("read at inference stage"))),
  outcome  = list(title = "MAG grouping outcome", icon = "split", icon_rot = TRUE,
                   sub = list(sprintf("%d MAGs (human-associated) retained", n_retained),
                              sprintf("%d MAGs (probable contaminant) excluded", n_excluded))),
  br_contam = list(title = "MAG (probable contaminant)", icon = "leaf", title_size = 14,
                    sub = list(sprintf("n = %d", n_excluded))),
  br_human  = list(title = "MAG (human-associated)", icon = "users", title_size = 14,
                    sub = list(sprintf("n = %d", n_retained))),
  val_def  = list(title = "(i) Defence arms race", icon = "shield", title_size = 13, sub_size = 10.5,
                   sub = list("Bakta v1.12.0 → DefenseFinder v2.0.1", "+ AntiDefenseFinder (DefenseFinder v3.0.0)", ital("compared between groups"))),
  val_is   = list(title = "(ii) IS composition", icon = "puzzle", title_size = 13, sub_size = 10.5,
                   sub = list("ISEScan v1.7.2.3", ital("compared between groups"))),
  val_phage = list(title = "(iii) Phage lifecycle", icon = "refresh_cw", title_size = 13, sub_size = 10.5,
                    sub = list("temperate vs lytic ecology,", "Pharokka calls read here", ital("compared between groups"))),
  inference = list(title = "Comparative genomic inference", icon = "scale",
                    sub = list("Compare defence, IS-element, and phage-lifecycle signatures between",
                               "the two groups — confirmatory biological validation, NOT an input",
                               "to the classification decision above"))
)

SZ <- lapply(content, function(b) {
  box_size(b$title, b$sub,
           title_size = if (!is.null(b$title_size)) b$title_size else 15,
           sub_size   = if (!is.null(b$sub_size))   b$sub_size   else 11.5)
})

MARGIN   <- 30
GAP_LANE <- 30
GAP_ROW  <- 45

lane_w <- c(
  L = max(SZ$taxonomy$w, SZ$val_def$w,  SZ$br_contam$w),
  C = max(SZ$viral$w,    SZ$iphop$w,    SZ$val_is$w, SZ$br_human$w),
  R = max(SZ$checkv$w,   SZ$pharokka$w, SZ$val_phage$w)
)
lane_x0 <- c(L = MARGIN, C = MARGIN + lane_w[["L"]] + GAP_LANE, R = NA)
lane_x0[["R"]] <- lane_x0[["C"]] + lane_w[["C"]] + GAP_LANE

W_LC   <- lane_w[["L"]] + GAP_LANE + lane_w[["C"]]
W_full <- W_LC + GAP_LANE + lane_w[["R"]]

grow_needed_LC   <- max(0, SZ$input$w - W_LC, SZ$classify$w - W_LC, SZ$outcome$w - W_LC)
if (grow_needed_LC > 0) {
  lane_w[["C"]] <- lane_w[["C"]] + grow_needed_LC
  lane_x0[["R"]] <- lane_x0[["C"]] + lane_w[["C"]] + GAP_LANE
  W_LC   <- lane_w[["L"]] + GAP_LANE + lane_w[["C"]]
  W_full <- W_LC + GAP_LANE + lane_w[["R"]]
}
grow_needed_full <- max(0, SZ$inference$w - W_full)
if (grow_needed_full > 0) {
  lane_w[["R"]] <- lane_w[["R"]] + grow_needed_full
  W_full <- W_LC + GAP_LANE + lane_w[["R"]]
}

W <- lane_x0[["R"]] + lane_w[["R"]] + MARGIN

place_in_lane <- function(lane, box_w) {
  x0 <- lane_x0[[lane]]; lw <- lane_w[[lane]]
  x0 + (lw - box_w) / 2
}
place_spanning <- function(span_w, box_w) {
  MARGIN + (span_w - box_w) / 2
}

Y <- list()
Y$input    <- 20
Y$taxonomy <- Y$input + SZ$input$h + GAP_ROW
Y$viral    <- Y$taxonomy
Y$iphop    <- Y$taxonomy + max(SZ$taxonomy$h, SZ$viral$h) + GAP_ROW
Y$checkv   <- Y$iphop
Y$classify <- Y$iphop + max(SZ$iphop$h, SZ$checkv$h) + GAP_ROW
Y$pharokka <- Y$classify
Y$outcome  <- Y$classify + max(SZ$classify$h, SZ$pharokka$h) + GAP_ROW
Yb         <- Y$outcome + SZ$outcome$h + 45
Y$br_contam <- Yb
Y$br_human  <- Yb
BR_H <- max(SZ$br_contam$h, SZ$br_human$h)
Y$val_def   <- Yb + BR_H + 55
Y$val_is    <- Y$val_def
Y$val_phage <- Y$val_def
VAL_H <- max(SZ$val_def$h, SZ$val_is$h, SZ$val_phage$h)
Y$inference <- Y$val_def + VAL_H + 55
H <- Y$inference + SZ$inference$h + 110

R <- list()
R$input    <- list(x = place_spanning(W_LC, SZ$input$w),    y = Y$input,    w = SZ$input$w,    h = SZ$input$h)
R$taxonomy <- list(x = place_in_lane("L", SZ$taxonomy$w),   y = Y$taxonomy, w = SZ$taxonomy$w, h = SZ$taxonomy$h)
R$viral    <- list(x = place_in_lane("C", SZ$viral$w),      y = Y$viral,    w = SZ$viral$w,    h = SZ$viral$h)
R$iphop    <- list(x = place_in_lane("C", SZ$iphop$w),      y = Y$iphop,    w = SZ$iphop$w,    h = SZ$iphop$h)
R$checkv   <- list(x = place_in_lane("R", SZ$checkv$w),     y = Y$checkv,   w = SZ$checkv$w,   h = SZ$checkv$h)
R$classify <- list(x = place_spanning(W_LC, SZ$classify$w), y = Y$classify, w = SZ$classify$w, h = SZ$classify$h)
R$pharokka <- list(x = place_in_lane("R", SZ$pharokka$w),   y = Y$pharokka, w = SZ$pharokka$w, h = SZ$pharokka$h)
R$outcome  <- list(x = place_spanning(W_LC, SZ$outcome$w),  y = Y$outcome,  w = SZ$outcome$w,  h = SZ$outcome$h)
R$br_contam <- list(x = place_in_lane("L", SZ$br_contam$w), y = Y$br_contam, w = SZ$br_contam$w, h = SZ$br_contam$h)
R$br_human  <- list(x = place_in_lane("C", SZ$br_human$w),  y = Y$br_human,  w = SZ$br_human$w,  h = SZ$br_human$h)
R$val_def   <- list(x = place_in_lane("L", SZ$val_def$w),   y = Y$val_def,   w = SZ$val_def$w,   h = SZ$val_def$h)
R$val_is    <- list(x = place_in_lane("C", SZ$val_is$w),    y = Y$val_is,    w = SZ$val_is$w,    h = SZ$val_is$h)
R$val_phage <- list(x = place_in_lane("R", SZ$val_phage$w), y = Y$val_phage, w = SZ$val_phage$w, h = SZ$val_phage$h)
R$inference <- list(x = place_spanning(W_full, SZ$inference$w), y = Y$inference, w = SZ$inference$w, h = SZ$inference$h)

cx <- function(id) R[[id]]$x + R[[id]]$w / 2
cy <- function(id) R[[id]]$y + R[[id]]$h / 2
top    <- function(id) R[[id]]$y
bottom <- function(id) R[[id]]$y + R[[id]]$h
left   <- function(id) R[[id]]$x
right  <- function(id) R[[id]]$x + R[[id]]$w

aabb_overlap <- function(a, b) (a$x < b$x + b$w) && (b$x < a$x + a$w) && (a$y < b$y + b$h) && (b$y < a$y + a$h)

all_rects <- R
for (id in names(R)) {
  bx <- R[[id]]$x + 18 - BADGE_R; by <- R[[id]]$y - BADGE_R
  all_rects[[paste0(id, "_badge")]] <- list(x = bx, y = by, w = 2 * BADGE_R, h = 2 * BADGE_R)
}
ids <- names(all_rects)
for (i in seq_along(ids)) for (j in seq_along(ids)) if (i < j) {
  a_id <- ids[i]; b_id <- ids[j]
  if (identical(a_id, paste0(b_id, "_badge")) || identical(b_id, paste0(a_id, "_badge"))) next
  if (aabb_overlap(all_rects[[a_id]], all_rects[[b_id]])) {
    stop(sprintf("Layout bug: '%s' overlaps '%s'.", a_id, b_id))
  }
}

smallest_pt_native <- 9.5 * PT_PER_UNIT
for (mm in c(180, 85)) {
  scale <- (mm / 25.4) / (W / UNITS_PER_IN)
  message(sprintf("Smallest text at %g mm width: %.2f pt", mm, smallest_pt_native * scale))
}

card <- function(id, n = NULL, border = "#0072B2", border_w = 1.8, fill = "#ffffff",
                  swatch_fill = "#EAF2FA", icon_colour = "#0072B2",
                  title_col = "#12233b", sub_col = "#4a5a6a",
                  dashed = FALSE, shadow = TRUE, icon_shape = "square") {
  b <- content[[id]]; r <- R[[id]]
  title_size <- if (!is.null(b$title_size)) b$title_size else 15
  sub_size   <- if (!is.null(b$sub_size))   b$sub_size   else 11.5
  dash <- if (dashed) ' stroke-dasharray="5 4"' else ''
  shadow_rect <- if (shadow)
    sprintf('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="12" fill="#000000" fill-opacity="0.13"/>',
            r$x + 2.2, r$y + 3, r$w, r$h)
  else ""
  rect <- paste0(shadow_rect,
           sprintf('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="12" fill="%s" stroke="%s" stroke-width="%.1f"%s/>',
                   r$x, r$y, r$w, r$h, fill, border, border_w, dash))
  bdg  <- if (!is.null(n)) badge(r$x, r$y, n) else ""
  ix <- r$x + LEFT_PAD; iy <- r$y + ICON_TOP_PAD
  swatch <- icon_swatch(ix, iy, b$icon, swatch_fill, icon_colour, shape = icon_shape, size = ICON_SIZE,
                         rotate180 = isTRUE(b$icon_rot))
  tx <- ix + ICON_SIZE + ICON_GAP
  ty <- iy + ICON_SIZE / 2 - 2
  ttl <- sprintf('<text x="%.1f" y="%.1f" text-anchor="start" font-weight="bold" font-size="%.1f" fill="%s">%s</text>',
                  tx, ty, title_size, title_col, b$title)
  sub_y0 <- iy + ICON_SIZE + SUB_GAP1
  subs <- imap_chr(b$sub, function(line, i) {
    sprintf('<text x="%.1f" y="%.1f" text-anchor="start" font-size="%.1f" font-style="%s" fill="%s">%s</text>',
            ix, sub_y0 + (i - 1) * LINE_H, sub_size, if (is_ital(line)) "italic" else "normal", sub_col, line)
  })
  paste(c(rect, bdg, swatch, ttl, subs), collapse = "\n")
}

orthconn <- function(pts, colour = "#0072B2", w = 1.6, dashed = FALSE, arrow = TRUE) {
  d <- paste(sprintf("%s%.1f,%.1f", c("M", rep("L", length(pts) / 2 - 1)),
                      pts[seq(1, length(pts), 2)], pts[seq(2, length(pts), 2)]),
             collapse = " ")
  dash <- if (dashed) ' stroke-dasharray="6 4"' else ''
  marker_id <- if (colour == "#6b7480") "arrowd" else if (colour == "#009E73") "arrowg" else "arrow"
  amark <- if (arrow) sprintf(' marker-end="url(#%s)"', marker_id) else ''
  sprintf('<path d="%s" fill="none" stroke="%s" stroke-width="%.2f"%s%s/>', d, colour, w, dash, amark)
}

assert_lands <- function(pts, from = NULL, to = NULL, tol = 0.6) {
  n <- length(pts)
  if (!is.null(to)) {
    lx <- pts[n - 1]; ly <- pts[n]
    if (abs(ly - top(to)) > tol || lx < left(to) - tol || lx > right(to) + tol) {
      stop(sprintf("Connector does not land on '%s': point (%.1f,%.1f), box x=[%.1f,%.1f] top=%.1f",
                    to, lx, ly, left(to), right(to), top(to)))
    }
  }
  if (!is.null(from)) {
    fx <- pts[1]; fy <- pts[2]
    if (abs(fy - bottom(from)) > tol || fx < left(from) - tol || fx > right(from) + tol) {
      stop(sprintf("Connector does not start on '%s': point (%.1f,%.1f), box x=[%.1f,%.1f] bottom=%.1f",
                    from, fx, fy, left(from), right(from), bottom(from)))
    }
  }
  invisible(TRUE)
}

BLUE <- "#0072B2"; GREEN <- "#009E73"; NAVY <- "#0B2545"; GREY <- "#8a94a0"; GREYD <- "#6b7480"

svg_parts <- c(
sprintf('<svg viewBox="0 0 %.1f %.1f" width="%.2fin" height="%.2fin" xmlns="http://www.w3.org/2000/svg" font-family="Helvetica, Arial, sans-serif">', W, H, W/100, H/100),
'<defs>
  <marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6.4" markerHeight="6.4" orient="auto-start-reverse">
    <path d="M0,0 L10,5 L0,10 z" fill="#0072B2"/>
  </marker>
  <marker id="arrowd" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6.4" markerHeight="6.4" orient="auto-start-reverse">
    <path d="M0,0 L10,5 L0,10 z" fill="#6b7480"/>
  </marker>
  <marker id="arrowg" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6.4" markerHeight="6.4" orient="auto-start-reverse">
    <path d="M0,0 L10,5 L0,10 z" fill="#009E73"/>
  </marker>
</defs>',
sprintf('<rect x="0" y="0" width="%.1f" height="%.1f" fill="#ffffff"/>', W, H),

card("input", 1),

{p1 <- c(cx("input"),bottom("input"), cx("input"),bottom("input")+25, cx("taxonomy"),bottom("input")+25, cx("taxonomy"),top("taxonomy"))
 assert_lands(p1, from = "input", to = "taxonomy"); orthconn(p1)},
{p2 <- c(cx("input"),bottom("input"), cx("input"),bottom("input")+25, cx("viral"),bottom("input")+25, cx("viral"),top("viral"))
 assert_lands(p2, from = "input", to = "viral"); orthconn(p2)},

card("taxonomy", 2),
card("viral", 3),

{LAND_PRIMARY <- left("classify") + R$classify$w * 0.25
 LAND_IPHOP   <- left("classify") + R$classify$w * 0.75
 jog_y <- top("classify") - 20
 p3 <- c(cx("taxonomy"),bottom("taxonomy"), cx("taxonomy"),jog_y, LAND_PRIMARY,jog_y, LAND_PRIMARY,top("classify"))
 assert_lands(p3, from = "taxonomy", to = "classify")
 orthconn(p3, w = 3.4)},
sprintf('<text x="%.0f" y="%.0f" font-size="11" fill="#0072B2" font-weight="bold">primary</text>', cx("taxonomy")-60, (bottom("taxonomy")+top("classify"))/2),
{p4 <- c(cx("viral"),bottom("viral"), cx("viral"),bottom("viral")+25, cx("iphop"),bottom("viral")+25, cx("iphop"),top("iphop"))
 assert_lands(p4, from = "viral", to = "iphop"); orthconn(p4)},
{p5 <- c(cx("viral"),bottom("viral"), cx("viral"),bottom("viral")+25, cx("checkv"),bottom("viral")+25, cx("checkv"),top("checkv"))
 assert_lands(p5, from = "viral", to = "checkv"); orthconn(p5, colour = GREYD)},

card("iphop", 4),
card("checkv", NULL, border = GREY, dashed = TRUE, shadow = FALSE, swatch_fill = "#F1F2F4", icon_colour = GREYD),

{p6 <- c(cx("iphop"),bottom("iphop"), cx("iphop"),bottom("iphop")+35, LAND_IPHOP,bottom("iphop")+35, LAND_IPHOP,top("classify"))
 assert_lands(p6, from = "iphop", to = "classify"); orthconn(p6, colour = GREYD, dashed = TRUE)},
{p7 <- c(cx("checkv"),bottom("checkv"), cx("checkv"),top("pharokka"))
 assert_lands(p7, from = "checkv", to = "pharokka"); orthconn(p7, colour = GREYD)},

card("classify", 5, border_w = 2.2),
card("pharokka", NULL, border = GREY, dashed = TRUE, shadow = FALSE, swatch_fill = "#F1F2F4", icon_colour = GREYD),

{p8 <- c(cx("classify"),bottom("classify"), cx("classify"),top("outcome"))
 assert_lands(p8, from = "classify", to = "outcome"); orthconn(p8, w = 3.6)},
sprintf('<text x="%.0f" y="%.0f" font-size="11" fill="#0072B2" font-weight="bold">the split</text>', cx("classify")+15, (bottom("classify")+top("outcome"))/2),
{p9 <- c(cx("pharokka"),bottom("pharokka"), cx("pharokka"),top("val_phage"))
 assert_lands(p9, from = "pharokka", to = "val_phage"); orthconn(p9, colour = GREYD)},

card("outcome", 6, border_w = 2.2),

{po1 <- c(cx("outcome")-80,bottom("outcome"), cx("outcome")-80,bottom("outcome")+22, cx("br_contam"),bottom("outcome")+22, cx("br_contam"),top("br_contam"))
 assert_lands(po1, from = "outcome", to = "br_contam"); orthconn(po1, colour = GREEN)},
{po2 <- c(cx("outcome")+80,bottom("outcome"), cx("outcome")+80,bottom("outcome")+22, cx("br_human"),bottom("outcome")+22, cx("br_human"),top("br_human"))
 assert_lands(po2, from = "outcome", to = "br_human"); orthconn(po2, colour = BLUE)},

card("br_contam", NULL, border = GREEN, border_w = 2.4, swatch_fill = "#E4F5F0", icon_colour = GREEN,
     title_col = "#00543d", sub_col = "#00543d", icon_shape = "circle"),
card("br_human", NULL, border = BLUE, border_w = 2.4, swatch_fill = "#E7F0FA", icon_colour = BLUE,
     title_col = "#00304d", sub_col = "#00304d", icon_shape = "circle"),

{
  bus_y <- bottom("br_contam") + 22
  val_phage_entry_x <- cx("val_phage") - 45
  stopifnot("val_phage bus entry point falls outside the box" =
              val_phage_entry_x >= left("val_phage") && val_phage_entry_x <= right("val_phage"))
  d1 <- c(cx("val_def"),bus_y, cx("val_def"),top("val_def"))
  d2 <- c(cx("val_is"),bus_y, cx("val_is"),top("val_is"))
  d3 <- c(val_phage_entry_x,bus_y, val_phage_entry_x,top("val_phage"))
  assert_lands(d1, to = "val_def"); assert_lands(d2, to = "val_is"); assert_lands(d3, to = "val_phage")
  paste(
    orthconn(c(cx("br_contam"),bottom("br_contam"), cx("br_contam"),bus_y), colour = GREYD, dashed = TRUE, arrow = FALSE),
    orthconn(c(cx("br_human"),bottom("br_human"), cx("br_human"),bus_y), colour = GREYD, dashed = TRUE, arrow = FALSE),
    orthconn(c(min(cx("val_def"), cx("br_contam")), bus_y, val_phage_entry_x, bus_y), colour = GREYD, dashed = TRUE, arrow = FALSE),
    orthconn(d1, colour = GREYD, dashed = TRUE),
    orthconn(d2, colour = GREYD, dashed = TRUE),
    orthconn(d3, colour = GREYD, dashed = TRUE),
    sep = "\n"
  )
},

card("val_def", 7, border = GREY, border_w = 1.5, swatch_fill = "#F1F2F4", icon_colour = "#4a5a6a"),
card("val_is", 7, border = GREY, border_w = 1.5, swatch_fill = "#F1F2F4", icon_colour = "#4a5a6a"),
card("val_phage", 7, border = GREY, border_w = 1.5, swatch_fill = "#F1F2F4", icon_colour = "#4a5a6a"),

{
  bus2_y <- max(bottom("val_def"), bottom("val_is"), bottom("val_phage")) + 25
  entry_x <- cx("inference")
  assert_lands(c(entry_x,bus2_y, entry_x,top("inference")), to = "inference")
  paste(
    orthconn(c(cx("val_def"),bottom("val_def"), cx("val_def"),bus2_y), colour = GREYD, dashed = TRUE, arrow = FALSE),
    orthconn(c(cx("val_is"),bottom("val_is"), cx("val_is"),bus2_y), colour = GREYD, dashed = TRUE, arrow = FALSE),
    orthconn(c(cx("val_phage"),bottom("val_phage"), cx("val_phage"),bus2_y), colour = GREYD, dashed = TRUE, arrow = FALSE),
    orthconn(c(cx("val_def"),bus2_y, cx("val_phage"),bus2_y), colour = GREYD, dashed = TRUE, arrow = FALSE),
    orthconn(c(entry_x,bus2_y, entry_x,top("inference")), colour = GREYD, dashed = TRUE),
    sep = "\n"
  )
},

card("inference", 8, border_w = 2.6, icon_shape = "circle"),

{ly <- bottom("inference")+40; sprintf('<text x="30" y="%.0f" font-size="11" font-weight="bold" fill="#1a1a1a">Legend</text>', ly)},
{ly <- bottom("inference")+55; sprintf('<line x1="30" y1="%.0f" x2="70" y2="%.0f" stroke="%s" stroke-width="3.4"/><text x="80" y="%.0f" font-size="10" fill="#333333">heavy = determines the split</text>', ly, ly, BLUE, ly+4)},
{ly <- bottom("inference")+72; sprintf('<line x1="30" y1="%.0f" x2="70" y2="%.0f" stroke="%s" stroke-width="1.6"/><text x="80" y="%.0f" font-size="10" fill="#333333">plain solid = ordinary data flow</text>', ly, ly, BLUE, ly+4)},
{ly <- bottom("inference")+89; sprintf('<line x1="30" y1="%.0f" x2="70" y2="%.0f" stroke="%s" stroke-width="1.5" stroke-dasharray="6 4"/><text x="80" y="%.0f" font-size="10" fill="#333333">dashed = corroborating / confirmatory, not determinative</text>', ly, ly, GREYD, ly+4)},

'</svg>'
)

svg_content <- paste(svg_parts, collapse = "\n")

svg_path <- file.path(FIG_DIR, "Figure_1.svg")
writeLines(svg_content, svg_path)
message("Saved SVG: ", svg_path)

pdf_path <- file.path(FIG_DIR, "Figure_1.pdf")
png_path <- file.path(FIG_DIR, "Figure_1.png")

system2("rsvg-convert", args = c("-f", "pdf", "-o", shQuote(pdf_path), shQuote(svg_path)))
system2("rsvg-convert", args = c("-f", "png", "--dpi-x", "600", "--dpi-y", "600",
                                  "-o", shQuote(png_path), shQuote(svg_path)))

# rsvg-convert does not record the resolution in the PNG; sips (macOS only) sets it to 600 dpi
if (nzchar(Sys.which("sips"))) {
  system2("sips", args = c("-s", "dpiWidth", "600", "-s", "dpiHeight", "600", shQuote(png_path)),
          stdout = FALSE, stderr = FALSE)
}

stopifnot(file.exists(pdf_path), file.exists(png_path))
message("Saved figure: ", pdf_path)
message("Saved figure: ", png_path)

nodes <- tribble(
  ~step, ~label,                        ~tool_version,
  "1",   "Input metagenomes",           "Battaglia et al. 2024 (Cell)",
  "2",   "Taxonomic classification",    "GTDB-Tk v2.3.0 (R214)",
  "3",   "Viral sequence identification","geNomad v1.9.7",
  "4",   "Phage-host prediction",       "iPHoP v1.4.2",
  "-",   "Viral quality filter",        "CheckV v1.0.1",
  "5",   "Classification decision",     "GTDB-Tk screen + iPHoP evidence, final call by literature review",
  "-",   "Phage lifecycle annotation",  "Pharokka v1.9.1",
  "6",   "MAG grouping outcome",        sprintf("%d retained / %d excluded", n_retained, n_excluded),
  "7i",  "Defence arms race (validation)", "Bakta v1.12.0 -> DefenseFinder v2.0.1 + AntiDefenseFinder (DefenseFinder v3.0.0)",
  "7ii", "IS composition (validation)", "ISEScan v1.7.2.3",
  "7iii","Phage lifecycle (validation)","Pharokka v1.9.1 (read, not recomputed)",
  "8",   "Comparative genomic inference", "confirmatory, not a classification input"
)
write_tsv(nodes, file.path(DATA_DIR, "Figure_1_node_data.tsv"))

edges <- tribble(
  ~from,       ~to,         ~style,
  "input",     "taxonomy",  "solid",
  "input",     "viral",     "solid",
  "taxonomy",  "classify",  "primary",
  "viral",     "iphop",     "solid",
  "viral",     "checkv",    "solid",
  "iphop",     "classify",  "dashed",
  "checkv",    "pharokka",  "solid",
  "classify",  "outcome",   "primary",
  "pharokka",  "val_phage", "solid",
  "outcome",   "br_contam", "solid",
  "outcome",   "br_human",  "solid",
  "br_contam", "val_def",   "dashed",
  "br_contam", "val_is",    "dashed",
  "br_contam", "val_phage", "dashed",
  "br_human",  "val_def",   "dashed",
  "br_human",  "val_is",    "dashed",
  "br_human",  "val_phage", "dashed",
  "val_def",   "inference", "dashed",
  "val_is",    "inference", "dashed",
  "val_phage", "inference", "dashed"
)
write_tsv(edges, file.path(DATA_DIR, "Figure_1_edge_data.tsv"))
message("Saved data: ", file.path(DATA_DIR, "Figure_1_node_data.tsv"))
message("Saved data: ", file.path(DATA_DIR, "Figure_1_edge_data.tsv"))

