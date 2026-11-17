# Figure 4

# Order of the script:
#   Setup    paths, print sizes, helper functions
#   Part 1   mate-pair synchrony: the label-permutation test
#   Part 2   panel B (mate-pair synchrony), drawn from the test results
#   then     panel C (scatter), panel D (bubble plot), the snRNA-seq row
#            (E, F, G), and the assembly
#   Panel letters as in Figures/Figure4_11122026.ai / .png

# Inputs: Data/phosphotrap/ (edgeR tables, TMM expression matrix, sample
#         sheets)
# Output: Figures/Figure4_11112026.{pdf,png}
#         Supplementary/derived_tables/PairSynchrony_results.csv, Supplementary/derived_tables/PairSynchrony_perPair.csv
#         Figures/Figure4_panels/: the synchrony panels on their own (all
#         three regions, and the POA panel)

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr)
  library(ggplot2); library(patchwork); library(cowplot)
  library(scales); library(tibble)
})

repo_dir <- here::here()
fig_dir        <- file.path(repo_dir, "Figures")
ptrap_data     <- file.path(repo_dir, "Data", "phosphotrap")
panel_dir      <- file.path(fig_dir, "Figure4_panels")     # stand-alone copies of the synchrony panels
out_dir        <- file.path(repo_dir, "Supplementary", "derived_tables")           # the permutation-test results
dir.create(panel_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Panel A is a drawing: Figures/artwork/Figure4_artwork.png, exported from
# Figure4_artwork.ai and cropped to the drawing.
IMG_PANEL_A <- file.path(repo_dir, "Figures", "artwork", "Figure4_artwork.png")

# ---------------------------------------------------------------------------
# Print size and fonts. The figure is composed at the width it is printed
# (double column), so every font size below is its printed size. Every ggplot
# panel goes through print_ready(), which sets theme text and in-plot text
# layers to these sizes.
# ---------------------------------------------------------------------------
FIG_W_IN <- 7.5

GEOM_SCALE_DEFAULTS <- 0.6   # theme-default line widths / point sizes at print size

# Nothing in the plots is smaller than 6 pt at the printed size.
# (Panel A's drawing has its own, smaller text and is not matched.)
FONT_PT  <- 7      # axis titles
SMALL_PT <- 6      # tick labels, in-plot labels, legends, subtitles

font_theme <- theme(
  line              = element_line(linewidth = 0.3),
  # legend keys/spacing in pt, not "lines" (which scale with each source
  # plot's own base font size and left big gaps between legend rows)
  legend.key.size   = grid::unit(FONT_PT + 1, "pt"),
  legend.key.spacing = grid::unit(2, "pt"),
  legend.spacing    = grid::unit(2, "pt"),
  legend.margin     = margin(1, 1, 1, 1),
  axis.text         = element_text(size = SMALL_PT),
  axis.text.x       = element_text(size = SMALL_PT),
  axis.text.y       = element_text(size = SMALL_PT),
  axis.title        = element_text(size = FONT_PT),
  axis.title.x      = element_text(size = FONT_PT),
  axis.title.y      = element_text(size = FONT_PT),
  plot.title        = element_text(size = FONT_PT + 1),
  plot.subtitle     = element_text(size = SMALL_PT),
  strip.text        = element_text(size = SMALL_PT),
  strip.text.y.left = element_text(size = SMALL_PT),
  legend.text       = element_text(size = SMALL_PT),
  legend.title      = element_text(size = SMALL_PT)
)

# ggplot2 >= 4.0 takes default line widths / point sizes from the theme
# (geoms whose size wasn't set explicitly, e.g. Panel C's dashed threshold
# lines, B's bar outlines). Scale those defaults down too.
if ("element_geom" %in% getNamespaceExports("ggplot2")) {
  font_theme <- font_theme + theme(geom = ggplot2::element_geom(
    linewidth   = 0.5 * GEOM_SCALE_DEFAULTS,
    borderwidth = 0.5 * GEOM_SCALE_DEFAULTS,
    pointsize   = 1.5 * GEOM_SCALE_DEFAULTS))
}

harmonise_text <- function(p) {
  # In-plot text layers (geom_text, annotate("text")) take their size in mm,
  # not pt. Layers are ggproto objects, so this edits them in place.
  for (lyr in p$layers) {
    if (inherits(lyr$geom, c("GeomText", "GeomLabel")) &&
        !is.null(lyr$aes_params$size)) {
      lyr$aes_params$size <- (SMALL_PT) / ggplot2::.pt
    }
  }
  p + font_theme
}

GEOM_SCALE <- 0.45
scale_geoms <- function(p, f = GEOM_SCALE) {
  plot_mapped <- names(p$mapping)
  for (lyr in p$layers) {
    if (inherits(lyr$geom, c("GeomText", "GeomLabel"))) next
    mapped <- c(names(lyr$mapping), if (isTRUE(lyr$inherit.aes)) plot_mapped)
    for (a in c("size", "linewidth", "stroke")) {
      if (a %in% mapped) next
      cur <- lyr$aes_params[[a]]
      if (is.null(cur)) cur <- lyr$geom$default_aes[[a]]
      if (is.numeric(cur)) lyr$aes_params[[a]] <- cur * f
    }
  }
  p
}

# Everything that goes into the figure gets both, and loses its plot title
# (panels are identified by letter). Subtitles, e.g. D's stats line, are kept.
print_ready <- function(p) scale_geoms(harmonise_text(p)) + labs(title = NULL)

# Drawn in place of an image file that is missing
placeholder_panel <- function(note) {
  ggplot() + theme_void() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
             fill = NA, color = "grey70", linetype = "dashed", linewidth = 0.3) +
    annotate("text", x = 0.5, y = 0.5, color = "grey55", fontface = "italic",
             size = SMALL_PT / ggplot2::.pt, label = note) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE)
}

# Panel A - drawing
has_A <- file.exists(IMG_PANEL_A)
if (!has_A) message("NOTE, Figure 4: panel A image not found, drawing a placeholder: ", IMG_PANEL_A)
img_A <- if (has_A) png::readPNG(IMG_PANEL_A) else NULL
ar_A  <- if (has_A) ncol(img_A) / nrow(img_A) else 0.70   # 0.70 = shape of the drawing, used for the placeholder

#panel b Mate-pair synchrony: the label-permutation test
#
# Are mates' baseline (Input) transcriptomes more alike
#   1. Pearson r of log2(TMM+1) Input expression for every pair of same-region
#      samples.
#   2. Mean r across the N true mate pairs (one female and one male of the
#      same breeding pair, both with usable Input data).
#   3. Null: male identity shuffled across the same N pairs, 10,000 times.
#   4. p = share of shuffles with mean mate r at or above the observed one.
# A Mann-Whitney test of mate r against all non-mate r would treat the
# non-mate values as independent; they are not (each sample is in about N-1
# of them). Shuffling labels keeps the real dependence structure.
#
# Run for the POA (the panel), pallium and subpallium.
# Writes Supplementary/derived_tables/PairSynchrony_results.csv and PairSynchrony_perPair.csv, and
# a stand-alone view of all three regions in Figures/Figure4_panels/.
# ============================================================================
set.seed(123)

# ---------- Paths ----------
fp_expr <- file.path(ptrap_data, "output.gene.TMM.EXPR.matrix")
fp_list <- file.path(ptrap_data, "pTRAP List.xlsx")
fp_meta <- file.path(ptrap_data, "phosphoTRAP data.xlsx")
stopifnot(file.exists(fp_expr), file.exists(fp_list), file.exists(fp_meta))

# Colors of the "real mate pairs" bar and of its dots (dark outline)
SYNC_GREEN      <- "#5aa04e"
SYNC_GREEN_DARK <- "#2c5a25"

# ---------- Theme ----------
theme_nature <- theme_classic(base_size = 9, base_family = "sans") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8))
theme_set(theme_nature)

# ---------- Build sample annotation: sid -> pair / region / sex ----------
ptrap_list <- read_excel(fp_list, sheet = "collectionsheet") |>
  rename(sid = Sample_id, region = Brain_region, sex = Sex,
         beg = BegNoBeg, frog_id = FrogID, testid = TestID) |>
  mutate(sid = as.integer(sid))

ptrap_meta <- read_excel(fp_meta, sheet = "Sheet1") |>
  rename(testid = `...1`, pair = pair) |>
  select(testid, pair) |> distinct()

samp <- ptrap_list |>
  left_join(ptrap_meta, by = "testid") |>
  filter(!is.na(pair))

cat("Loaded", nrow(samp), "annotated samples across",
    length(unique(samp$pair)), "breeding pairs\n")

# ---------- Load TMM matrix ----------
expr <- read.delim(fp_expr, header = TRUE, row.names = 1, check.names = FALSE,
                   stringsAsFactors = FALSE)
expr <- as.matrix(expr)

# Parse column names like "13i_R1_001" -> sid=13, fraction='i'
parse_colname <- function(cn) {
  m <- regmatches(cn, regexec("^(\\d+)([ip])_", cn))[[1]]
  if (length(m) < 3) return(c(NA, NA))
  c(sid = as.integer(m[2]), fraction = m[3])
}
parsed <- t(sapply(colnames(expr), parse_colname))
parsed_df <- data.frame(col = colnames(expr),
                        sid = as.integer(parsed[, 1]),
                        fraction = parsed[, 2],
                        stringsAsFactors = FALSE) |>
  filter(!is.na(sid))

# ---------- Per-region permutation test ----------
run_region <- function(region_name, n_perm = 10000) {
  region_sids <- samp |> filter(region == region_name) |> pull(sid)
  inputs <- parsed_df |>
    filter(sid %in% region_sids, fraction == "i")
  if (nrow(inputs) == 0)
    stop("No Input columns found for region: ", region_name)

  M <- log2(expr[, inputs$col, drop = FALSE] + 1)
  variable <- apply(M, 1, var, na.rm = TRUE) > 0
  M <- M[variable, , drop = FALSE]
  colnames(M) <- inputs$sid

  R <- cor(M, method = "pearson")

  region_samp <- samp |> filter(region == region_name, sid %in% inputs$sid)
  pair_table <- region_samp |>
    group_by(pair) |>
    filter(n_distinct(sex) == 2) |>
    summarise(F_sid = sid[sex == "Female"][1],
              M_sid = sid[sex == "Male"][1],
              .groups = "drop") |>
    filter(!is.na(F_sid), !is.na(M_sid))
  N <- nrow(pair_table)

  mate_rs <- mapply(function(f, m) R[as.character(f), as.character(m)],
                    pair_table$F_sid, pair_table$M_sid)
  obs_mean <- mean(mate_rs)

  females <- pair_table$F_sid
  males   <- pair_table$M_sid
  null_means <- numeric(n_perm)
  for (k in seq_len(n_perm)) {
    perm_males <- sample(males)
    null_means[k] <- mean(mapply(function(f, m) R[as.character(f), as.character(m)],
                                  females, perm_males))
  }

  list(region = region_name,
       N_pairs = N,
       observed = obs_mean,
       null_mean = mean(null_means),
       null_sd = sd(null_means),
       null_means = null_means,
       p = mean(null_means >= obs_mean),
       pair_table = pair_table |> mutate(r = mate_rs))
}

cat("\nRunning label-permutation tests (10,000 perms each) ...\n")
res <- lapply(c("POA", "Pallium", "Subpallium"), run_region)
names(res) <- c("POA", "Pallium", "Subpallium")

summary_tbl <- do.call(rbind, lapply(res, function(r) {
  data.frame(region = r$region,
             N_mate_pairs = r$N_pairs,
             observed_r = round(r$observed, 4),
             null_mean = round(r$null_mean, 4),
             null_sd = round(r$null_sd, 4),
             z = round((r$observed - r$null_mean) / r$null_sd, 2),
             p_label_permutation = round(r$p, 4))
}))
cat("\n"); print(summary_tbl, row.names = FALSE)

#region     N_mate_pairs   observed_r null_mean null_sd zscore p_label_permutation
#POA              9        0.5219    0.5030  0.0099 1.91              0.0245
#Pallium            8       0.7312    0.7281  0.0035 0.89              0.2033
#Subpallium            9     0.6767    0.6720  0.0048 0.99              0.1800

make_bar_panel <- function(r, n_breaks = 5) {
  bar_df <- data.frame(
    group = factor(c("Real\nmate pairs", "Null\n(shuffled)"),
                   levels = c("Real\nmate pairs", "Null\n(shuffled)")),
    mean  = c(r$observed, r$null_mean),
    err   = c(NA, 1.96 * r$null_sd)
  )
  dot_df <- data.frame(group = "Real\nmate pairs", r_val = r$pair_table$r)

  y_lo <- min(c(r$pair_table$r, r$null_mean - 2.5 * r$null_sd)) - 0.005
  y_hi <- max(c(r$pair_table$r, r$null_mean + 2.5 * r$null_sd)) + 0.005

  ggplot(bar_df, aes(group, mean, fill = group)) +
    geom_col(color = "black", linewidth = 0.5, width = 0.55, alpha = 0.7) +
    geom_errorbar(aes(ymin = mean - err, ymax = mean + err),
                  width = 0.18, linewidth = 0.6, na.rm = TRUE) +
    geom_jitter(data = dot_df, aes(group, r_val),
                shape = 21, fill = SYNC_GREEN, color = SYNC_GREEN_DARK, size = 1.8, stroke = 0.5,
                width = 0.10, height = 0, inherit.aes = FALSE) +
    scale_fill_manual(values = c("Real\nmate pairs" = SYNC_GREEN,
                                 "Null\n(shuffled)" = "#888888"),
                      guide = "none") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = n_breaks)) +
    coord_cartesian(ylim = c(y_lo, y_hi)) +
    labs(x = NULL, y = "Mean mate-pair r",
         title = sprintf("%s  (n = %d pairs)", r$region, r$N_pairs),
         subtitle = sprintf("obs = %.3f, null = %.3f, p = %.3f",
                            r$observed, r$null_mean, r$p))
}

fig <- (make_bar_panel(res$POA, n_breaks = 6) |
        make_bar_panel(res$Pallium) |
        make_bar_panel(res$Subpallium)) +
  plot_annotation(
    title = "Input transcriptional synchrony between mate pairs",
    subtitle = "Bars: observed (red) vs label-permutation null (gray, ± 1.96·SD); dots = individual mate pairs",
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 9, color = "grey40")))


ggsave(file.path(panel_dir, "Pair_synchrony_INPUT_R.png"), fig,
       width = 6, height = 4, dpi = 300, bg = "white")
ggsave(file.path(panel_dir, "Pair_synchrony_INPUT_R.pdf"), fig,
       width = 6, height = 4)

write.csv(summary_tbl,
          file.path(out_dir, "PairSynchrony_results.csv"),
          row.names = FALSE)

per_pair_tbl <- do.call(rbind, lapply(res, function(r) {
  data.frame(region = r$region, r$pair_table)
}))
write.csv(per_pair_tbl,
          file.path(out_dir, "PairSynchrony_perPair.csv"),
          row.names = FALSE)

sync_results <- read.csv(file.path(out_dir, "PairSynchrony_results.csv"))
sync_perPair <- read.csv(file.path(out_dir, "PairSynchrony_perPair.csv"))

poa  <- sync_results |> filter(region == "POA")
dots <- sync_perPair |> filter(region == "POA")

bar <- data.frame(
  group = factor(c("Real\nmate pairs", "Null\n(shuffled)"),
                 levels = c("Real\nmate pairs", "Null\n(shuffled)")),
  mean  = c(poa$observed_r, poa$null_mean),
  err   = c(NA, 1.96 * poa$null_sd))

p_txt <- if (poa$p_label_permutation < 0.001) "p < 0.001" else
         sprintf("p = %.3f", poa$p_label_permutation)
y_lo <- min(c(dots$r, poa$null_mean - 2.5 * poa$null_sd)) - 0.005
y_hi <- max(c(dots$r, poa$null_mean + 2.5 * poa$null_sd)) + 0.005

poa_panel <- ggplot(bar, aes(group, mean, fill = group)) +
  geom_col(color = "black", linewidth = 0.5, width = 0.55, alpha = 0.7) +
  geom_errorbar(aes(ymin = mean - err, ymax = mean + err),
                width = 0.18, linewidth = 0.6, na.rm = TRUE) +
  # one dot per mate pair: the bar's green, outlined in a darker shade of it
  # (the dot style of the Figure 2 and Figure 3 box panels)
  geom_jitter(data = dots, aes(x = "Real\nmate pairs", y = r), inherit.aes = FALSE,
              shape = 21, fill = SYNC_GREEN, color = SYNC_GREEN_DARK,
              size = 2.4, stroke = 0.9, width = 0.10, height = 0) +
  scale_fill_manual(values = c("Real\nmate pairs" = SYNC_GREEN,
                               "Null\n(shuffled)" = "#888888"), guide = "none") +
  coord_cartesian(ylim = c(y_lo, y_hi)) +
  labs(x = NULL, y = "Mean mate-pair r", title = "POA synchrony",
       subtitle = sprintf("obs %.3f vs null %.3f\n%s  (n = %d pairs)",
                          poa$observed_r, poa$null_mean, p_txt, poa$N_mate_pairs)) +
  theme_classic(base_size = 9) +   # ~matches Figure3 top-row fonts
  theme(panel.grid = element_blank(),
        axis.line  = element_line(color = "black", linewidth = 0.4),
        axis.ticks = element_line(color = "black", linewidth = 0.4),
        legend.position = "none",
        plot.title    = element_text(face = "bold", size = 9.5, hjust = 0.5),
        plot.subtitle = element_text(size = 7, hjust = 0.5))

poa_panel

ggsave(file.path(panel_dir, "POA_synchrony_panel.pdf"), poa_panel, width = 3.3, height = 4.4)
ggsave(file.path(panel_dir, "POA_synchrony_panel.png"), poa_panel,
       width = 3.3, height = 4.4, dpi = 300, bg = "white")
cat(sprintf("Saved POA_synchrony_panel.{pdf,png}  (POA obs r = %.3f, p = %.4f, n = %d)\n",
            poa$observed_r, poa$p_label_permutation, poa$N_mate_pairs))
panelB <- print_ready(poa_panel)

d_sub <- panelB$labels$subtitle
if (!is.null(d_sub)) {
  m <- regmatches(d_sub, gregexpr("[0-9]*\\.?[0-9]+", d_sub))[[1]]
  # m = obs, null, p, n (numbers in the order the subtitle prints them)
  if (length(m) >= 4) d_sub <- sprintf("%s vs %s\np = %s, n = %s", m[1], m[2], m[3], m[4])   # two lines: B is narrow
}

#Panel C -- POA, sexes combined: begging against no-begging log2 fold change,
# from the two edgeR tables (IP vs Input)
# ============================================================================
read_edgeR <- function(filename) read.delim(file.path(ptrap_data, filename), header = TRUE)

nobeg_poa <- read_edgeR("IP_Input_nobeg_bothsex_poa_edgeR_Anno.txt")
beg_poa   <- read_edgeR("IP_Input_beg_bothsex_poa_edgeR_Anno.txt")

poa_compare <- nobeg_poa[, c("Gene", "Gene.symbol", "logFC", "padj")]
colnames(poa_compare)[colnames(poa_compare) == "logFC"] <- "logFC_nobeg"
colnames(poa_compare)[colnames(poa_compare) == "padj"]  <- "padj_nobeg"

beg_tmp <- beg_poa[, c("Gene", "logFC", "padj")]  # Gene.symbol dropped here --
                                                    # keeping it on both sides made
                                                    # merge() rename to Gene.symbol.x/.y,
                                                    # so poa_compare$Gene.symbol later
                                                    # resolved to NULL (0-length) and
                                                    # broke the $highlight assignment.
colnames(beg_tmp)[colnames(beg_tmp) == "logFC"] <- "logFC_beg"
colnames(beg_tmp)[colnames(beg_tmp) == "padj"]  <- "padj_beg"

poa_compare <- merge(poa_compare, beg_tmp, by = "Gene")

poa_compare <- poa_compare %>%
  mutate(
    sig_nobeg = padj_nobeg < 0.05 & abs(logFC_nobeg) > 0.6,
    sig_beg   = padj_beg   < 0.05 & abs(logFC_beg)   > 0.6,
    DE_class = case_when(
      sig_nobeg & sig_beg  ~ "DT in Both",
      sig_nobeg & !sig_beg ~ "DT in No-Begging",
      !sig_nobeg & sig_beg ~ "DT in Begging only",
      TRUE                 ~ "Not DT"
    ),
    DE_class = factor(DE_class, levels = c("DT in Begging only", "DT in No-Begging",
                                            "Not DT", "DT in Both"))
  )

cat("Panel C | POA genes by class (padj < 0.05 and |log2FC| > 0.6):\n")
print(table(poa_compare$DE_class))
#DT in Begging only   DT in No-Begging             Not DT         DT in Both 
#289               1093               6150                141

# Genes labelled in the panel, by gene symbol.
genes_to_label <- c("NEXMIF", "FAAH1", "CRF", "OPRM", "PDYN", "CNR1", "VIP", "PENKA",
                     "PENK", "GAL", "OXT", "OT", "AVP", "CRFR1", "COMTA", "BOMB", "PNOC",
                     "GPR83")
poa_compare$highlight <- poa_compare$Gene.symbol %in% genes_to_label
# CRFR1 and GPR83 each match two duplicate-annotated transcripts in this assembly (see the note above panel D's gene table on this same issue). For each, only the stronger/more significant transcript gets labeled.
#   CRFR1:  rb_8081_pilon (|logFC| ~3-4, padj_nobeg = 0.038) kept over
#           rb_598_pilon  (|logFC| ~1-2, padj_nobeg = 0.224)
#   GPR83:  rb_14905_pilon (padj_beg = 1.1e-4, the begging-only hit) kept over
#           rb_20664_pilon (not significant in either condition)
poa_compare$highlight[poa_compare$Gene.symbol == "CRFR1" & poa_compare$Gene != "rb_8081_pilon"] <- FALSE
poa_compare$highlight[poa_compare$Gene.symbol == "GPR83" & poa_compare$Gene != "rb_14905_pilon"] <- FALSE

# Both axes show -logFC, so that positive = enriched in the IP fraction.

panelC_x    <- -poa_compare$logFC_nobeg
panelC_xrng <- range(panelC_x, na.rm = TRUE)
panelC_xpad <- diff(panelC_xrng)

poa_compare$label_color_cat <- as.character(poa_compare$DE_class)
poa_compare$label_color_cat[poa_compare$DE_class == "Not DT"] <- "Not DT (label)"

LABEL_SPREAD  <- c(0.04, 0.96)     # bottom and top of the label column, as fractions of the y range
panelC_yrng   <- range(-poa_compare$logFC_beg, na.rm = TRUE)
panelC_labels <- subset(poa_compare, highlight)
panelC_labels <- panelC_labels[order(panelC_labels$logFC_beg), ]          # highest point first
panelC_labels$x_lab <- panelC_xrng[2] + 0.38 * panelC_xpad
panelC_labels$y_lab <- seq(panelC_yrng[1] + LABEL_SPREAD[2] * diff(panelC_yrng),
                           panelC_yrng[1] + LABEL_SPREAD[1] * diff(panelC_yrng),
                           length.out = nrow(panelC_labels))

lines_cross <- function(i, j, d) {
  px <- -d$logFC_nobeg; py <- -d$logFC_beg
  side <- function(ax, ay, bx, by, cx, cy) sign((bx - ax) * (cy - ay) - (by - ay) * (cx - ax))
  side(px[i], py[i], d$x_lab[i], d$y_lab[i], px[j], py[j]) != side(px[i], py[i], d$x_lab[i], d$y_lab[i], d$x_lab[j], d$y_lab[j]) &&
    side(px[j], py[j], d$x_lab[j], d$y_lab[j], px[i], py[i]) != side(px[j], py[j], d$x_lab[j], d$y_lab[j], d$x_lab[i], d$y_lab[i])
}
for (pass in seq_len(nrow(panelC_labels))) {
  swapped <- FALSE
  for (i in seq_len(nrow(panelC_labels) - 1)) if (lines_cross(i, i + 1, panelC_labels)) {
    y_two <- panelC_labels$y_lab[c(i, i + 1)]
    panelC_labels[c(i, i + 1), ] <- panelC_labels[c(i + 1, i), ]
    panelC_labels$y_lab[c(i, i + 1)] <- y_two
    swapped <- TRUE
  }
  if (!swapped) break
}

panelC_core <- ggplot(poa_compare, aes(x = -logFC_nobeg, y = -logFC_beg, color = DE_class)) +
  geom_point(size = 1.0, alpha = 0.7) +   # x GEOM_SCALE at print size (~0.45)
  geom_segment(data = panelC_labels,
               aes(xend = x_lab - 0.01 * panelC_xpad, yend = y_lab, color = label_color_cat),
               linewidth = 0.3 / GEOM_SCALE, show.legend = FALSE) +   # print_ready() scales it back to 0.3
  geom_text(data = panelC_labels,
            aes(x = x_lab, y = y_lab, label = Gene.symbol, color = label_color_cat),
            size = 3, fontface = "italic", hjust = 0, vjust = 0.5, show.legend = FALSE) +
  geom_hline(yintercept = c(-0.6, 0.6), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-0.6, 0.6), linetype = "dashed", color = "grey40") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40") +
  scale_color_manual(values = c(
    "DT in Begging only" = "#D32F2F",
    "DT in No-Begging"   = "#1565C0",
    "Not DT"             = "grey80",
    "DT in Both"         = "grey30",
    "Not DT (label)"     = "black"),
    breaks = c("DT in Begging only", "DT in No-Begging", "Not DT", "DT in Both"),
    name = NULL) +
  coord_cartesian(xlim = c(panelC_xrng[1] - 0.05 * panelC_xpad,
                           panelC_xrng[2] + 0.70 * panelC_xpad),
                  clip = "off") +
  labs(x = expression(Log[2] * " Fold Change (No begging)"),
       y = expression(Log[2] * " Fold Change (Begging)"),
       title = "POA (combined sexes)") +
  theme_classic(base_size = 11)

panelC_core

panelC_core <- print_ready(panelC_core) +
  guides(color = guide_legend(nrow = 1, override.aes = list(size = 1.8, alpha = 1)))

panelC   <- panelC_core + theme(legend.position = "none")
legend_C <- cowplot::get_legend(
  panelC_core + theme(legend.position = "bottom", legend.title = element_blank())
)

# panel D
# Input: Data/phosphotrap/IP_Input_{nobeg,beg}_{male,female}_poa_edgeR*.txt,
# edgeR results (IP vs Input) for the POA, one file per sex and stimulus. The
# no-begging female file has a different column layout from the other three;
# both are read. logFC is negated so that positive = enriched in the IP
# fraction, as in panel C. Several gene symbols match more than one transcript
# in the annotation; the transcript used for each symbol is fixed in the
# table below.
# ---- canonical symbol -> transcript ID, and system category ----------------
canon <- tribble(
  ~gene.symbol, ~Gene,                          ~category,
  "OPRM",       "rb_16192_pilon",               "Opioid system",
  "PDYN",       "rb_18826_pilon",               "Opioid system",
  "PENKA",      "rb_7418_pilon",                "Opioid system",
  "PNOC",       "rb_20033_pilon",               "Opioid system",
  "FAAH1",      "rb_23602_pilon",               "Endocannabinoid",
  "AVP",        "rb_31702_pilon_billieadded",   "Neuropeptides",
  "OXT",        "rb_31695_pilon_billieadded",   "Neuropeptides",
  "CRF",        "rb_13255_pilon",               "Neuropeptides",
  "CRFR1",      "rb_8081_pilon",                "Neuropeptides",   # the same CRFR1 transcript that panel C labels
  "UCN3",       "rb_31696_pilon_billieadded",   "Neuropeptides",
  "TRH",        "rb_9464_pilon",                "Neuropeptides",
  "TRHDE",      "rb_15346_pilon",               "Neuropeptides",
  "GNRH",       "rb_31703_pilon_billieadded",   "Neuropeptides",
  "GAL",        "rb_31694_pilon_billieadded",   "Neuropeptides",
  "NPY",        "rb_19415_pilon",               "Neuropeptides",
  "VIP",        "rb_3280_pilon",                "Neuropeptides",
  "PRRP",       "rb_8360_pilon",                "Neuropeptides",
  "CART",       "rb_13238_pilon",               "Neuropeptides",
)
canon$category <- factor(canon$category,
                          levels = c("Opioid system","Endocannabinoid","Neuropeptides"))

# ---- read one per-sex POA edgeR file (either column layout) --------------
# Both layouts share Gene/logFC/padj; we join on the pinned transcript ID
# (canon$Gene) rather than any symbol column, so the symbol-naming
# inconsistency between layouts never matters here.
read_poa_sexed <- function(filename) {
  df <- read.delim(file.path(ptrap_data, filename), header = TRUE)
  stopifnot(all(c("Gene", "logFC", "padj") %in% names(df)))
  df[, c("Gene", "logFC", "padj")]
}

conditions <- tribble(
  ~treatment, ~sex,     ~file,
  "nobeg",    "male",   "IP_Input_nobeg_male_poa_edgeR.GeneName_Anno.txt",
  "beg",      "male",   "IP_Input_beg_male_poa_edgeR.GeneName_Anno.txt",
  "beg",      "female", "IP_Input_beg_female_poa_edgeR.GeneName_Anno.txt",
  "nobeg",    "female", "IP_Input_nobeg_female_poa_edgeR_Anno.txt",
)

poa_only_df <- dplyr::bind_rows(lapply(seq_len(nrow(conditions)), function(i) {
  cond <- conditions[i, ]
  raw  <- read_poa_sexed(cond$file)
  merge(canon, raw, by = "Gene") |>
    distinct(gene.symbol, .keep_all = TRUE) |>   # collapse duplicate synonym rows per transcript
    mutate(treatment = cond$treatment, sex = cond$sex, brainregion = "poa")
}))

# Negate logFC to keep the positive = IP-enriched convention
poa_only_df <- poa_only_df |>
  mutate(
    logFC        = -logFC,
    neglog10padj = -log10(pmax(padj, 1e-10)),
    # Plain "M"/"F" rather than Mars/Venus glyphs (<U+2642>/<U+2640>): those
    # Unicode symbols don't render in every graphics device/font -- they came
    # out as invisible empty glyphs here, which silently made the male and
    # female condition columns look identical ("POA NoBEG" x2, "POA BEG" x2)
    # instead of 4 distinct columns. Found while widening this panel.
    sex_lab      = factor(sex, levels = c("male","female"), labels = c("M","F")),
    tx_lab       = factor(treatment, levels = c("nobeg","beg"), labels = c("NoBEG","BEG")),
    condition    = paste("POA", sex_lab, tx_lab),
    category     = factor(category,
                           levels = c("Opioid system","Endocannabinoid","Neuropeptides"),
                           labels = c("Endogenous opioids","Endocannabinoid","Neuropeptides")),
    gene.symbol  = factor(gene.symbol, levels = canon$gene.symbol)
  )

cond_levels <- c()
for (sex in c("M","F"))
  for (tx in c("NoBEG","BEG"))
    cond_levels <- c(cond_levels, paste("POA", sex, tx))
poa_only_df$condition <- factor(poa_only_df$condition, levels = cond_levels)

bubble_poa <- ggplot(poa_only_df, aes(x = condition, y = gene.symbol)) +
  geom_point(aes(size = neglog10padj, fill = logFC),
             shape = 21, color = "black", stroke = 0.3, alpha = 0.95) +
  scale_fill_gradient2(low = "#542788", mid = "#f7f7f7", high = "#b35806",
                       midpoint = 0, limits = c(-8, 8), oob = scales::squish,
                       name = "logFC\n(orange = IP-enriched\npurple = IP-depleted)") +
  scale_size_continuous(range = c(1, 6),
                        breaks = -log10(c(0.5, 0.05, 0.01, 0.001)),
                        labels = c("0.5","0.05","0.01","0.001"),
                        name = "padj") +
  # No facet_grid by category: faceting put FAAH1 (the lone Endocannabinoid
  # entry) in its own single-row strip between Opioid system and
  # Neuropeptides, which read as a visual "split". gene.symbol's factor
  # levels (set above from canon$gene.symbol, in category order) already
  # keep the three groups in the right order on one continuous axis.
  scale_y_discrete(limits = rev) +
  theme_classic(base_size = 9) +
  theme(axis.text.x       = element_text(angle = 0, hjust = 0.5, size = 8),
        axis.text.y       = element_text(size = 8),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
        plot.title         = element_text(face = "bold", size = 11),
        plot.subtitle      = element_text(size = 9, color = "grey30")) +
  labs(x = NULL, y = NULL,
       title    = "POA — Endogenous opioids, Endocannabinoid, Neuropeptides",
       subtitle = "Color = logFC (orange = IP-enriched, purple = IP-depleted); size = −log10(padj)\nRebuilt directly from per-sex edgeR output (see header note)")

panelD_core <- bubble_poa +
  # "POA " prefix dropped from the column labels -- the panel title already
  # says POA, and at print size the full labels collide.
  scale_x_discrete(labels = function(x) sub("^POA ", "", x)) +
  # No facet_grid here either (see the note in the plot above) -- one
  # continuous axis, so the single-gene Endocannabinoid group (FAAH1) no
  # longer reads as split off from Opioid system/Neuropeptides above/below
  # it. Order is preserved via gene.symbol's factor levels.
  # Bubble sizes: same breaks/labels as the standalone script, range scaled
  # down to print size like the other geometry.
  scale_size_continuous(range = c(1, 6) * GEOM_SCALE,
                        breaks = -log10(c(0.5, 0.05, 0.01, 0.001)),
                        labels = c("0.5", "0.05", "0.01", "0.001"),
                        name = "padj") +
  labs(subtitle = NULL) +
    guides(fill = guide_colorbar(title = "logFC", direction = "horizontal", title.position = "left",
                                title.vjust = 0.9,
                                barwidth = grid::unit(0.55, "in"), barheight = grid::unit(0.06, "in"),
                                order = 1),
         size = guide_legend(direction = "horizontal", title.position = "left", order = 2)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "bottom", legend.box = "horizontal",
        legend.box.spacing = grid::unit(0, "pt"),
        legend.spacing.y = grid::unit(1, "pt"),
        plot.title.position = "plot")
panelD_core <- print_ready(panelD_core)
panelD   <- panelD_core + theme(legend.position = "none")
legend_D <- cowplot::get_legend(panelD_core)
panelD

#snRNAseq
#These panels are NOT drawn here. They were made in Python by Prateek Kalakuntla(the Jupyter
# notebooks in Data/Analysis_ForBillie_Figure4/, which need scanpy, SAMap
# and the .h5ad expression files). This script only places the exported
# images, left to right:
#   t-SNE of the nuclei colored by cluster | Penk expression on the t-SNE |
#   cumulative distribution of the clusters' enrichment ranks for begging
#   marker genes, cluster 106 (PVpo-VMPO-MPN Hmx2 Gaba) marked
#   (snRNAseq_enrichment_rank_ECDF.png; no notebook or data for this one is
#   in the folder)
# Each image keeps its own shape; the columns are sized so all three have the
# same height. If an image is missing, a labelled placeholder fills the row.
# ============================================================================
snrna_dir    <- file.path(repo_dir, "Data", "Analysis_ForBillie_Figure4")
SNRNA_IMAGES <- file.path(snrna_dir, c("RI_sss_subclass_v5_10012026.png",
                                       "Imitator_Frog_Penk_10012026.png",
                                       "snRNAseq_enrichment_rank_ECDF.png"))
panel_snRNAseq <- if (all(file.exists(SNRNA_IMAGES))) {
  imgs   <- lapply(SNRNA_IMAGES, png::readPNG)
  aspect <- vapply(imgs, function(m) ncol(m) / nrow(m), numeric(1))     # width / height
  wrap_plots(lapply(imgs, function(m)
    ggplot() + theme_void() + annotation_custom(grid::rasterGrob(m, interpolate = TRUE))),
    nrow = 1, widths = aspect)
} else {
  message("NOTE, Figure 4: snRNA-seq image(s) not found, drawing a placeholder: ",
          paste(basename(SNRNA_IMAGES[!file.exists(SNRNA_IMAGES)]), collapse = ", "))
  placeholder_panel("[ single-nucleus RNA-seq panels: images not found ]")
}

panel_snRNAseq

# Assemble
#
# Panel letters as in Figures/Figure4_11122026.ai / .png:
#   A  schematic (drawing)           B  mate-pair synchrony, under A
#   C  scatter, sexes combined       D  bubble plot, sexes separated
#   E, F, G  the three snRNA-seq images (second row)
#
#   Row 1:  [A, with B underneath it] | C | D
#   Row 2:  E | F | G; the legends of C and D sit in a strip at the top of
#           this row, each under its own panel

# C and D are built separately and have different margins around their plot
# areas; align_plots() lines up the top and bottom edges of the plot areas so
# the two x-axes sit level. B is in the call because its margins take part:
# without it C and D's plot areas come out a different height.
aligned_row1   <- cowplot::align_plots(panelC, panelD, panelB, align = "hv", axis = "tb")
panelC_aligned <- wrap_elements(full = aligned_row1[[1]])
panelD_aligned <- wrap_elements(full = aligned_row1[[2]])

# ---- sizes (inches) -----------------------------------------------------------
ROW1_H_IN <- 4.3          # A (with B under it) | C | D
ROW2_H_IN <- 2.85         # snRNA-seq images E to G (about 2.2 in tall) plus the
                          # legend strip above them and the margin below
FIG_H_IN  <- ROW1_H_IN + ROW2_H_IN
A_W_IN    <- 1.6          # printed width of A; the left column is exactly this wide
B_W_IN    <- 1.3          # width of B (synchrony) under A
W_C       <- 3.1          # C and D share the width right of the A/B column
W_D       <- 1.4          #   in this proportion

# ---- left column: A on top at its true shape, B in the space underneath ------
panelA  <- if (has_A) {
  cowplot::ggdraw() +
    cowplot::draw_grob(grid::rasterGrob(img_A, y = 1, vjust = 1, width = grid::unit(1, "npc")))
} else placeholder_panel("[ panel A:\nimage not found ]")
panelAB <- (panelA / ((panelB | plot_spacer()) +
                        plot_layout(widths = grid::unit.c(grid::unit(B_W_IN, "in"), grid::unit(1, "null"))))) +
  plot_layout(heights = grid::unit.c(grid::unit(A_W_IN / ar_A, "in"), grid::unit(1, "null")))

row1 <- (panelAB | panelC_aligned | panelD_aligned) +
  plot_layout(widths = grid::unit.c(grid::unit(A_W_IN, "in"),
                                    grid::unit(W_C, "null"), grid::unit(W_D, "null")))

# Row 2 (E, F, G): the top margin is the strip that holds the legends of C and D
row2 <- wrap_elements(full = panel_snRNAseq) +
  theme(plot.margin = margin(t = 26, r = 5.5, b = 22, l = 5.5))

final_layout <- (row1 / row2) + plot_layout(heights = c(ROW1_H_IN, ROW2_H_IN))

# ---- legends of C and D, drawn in the strip at the top of row 2 --------------
# D's legend (color bar and two circle sizes) needs a fixed width whatever
# the width of panel D; C's legend gets the rest of the width right of A.
LEGEND_OFFSET_IN <- 0.10   # legend centre, this far above the top of row 2
LEGEND_D_W_IN    <- 2.0
x_legC <- A_W_IN / FIG_W_IN                 # left edge of C's legend
x_legD <- 1 - LEGEND_D_W_IN / FIG_W_IN      # left edge of D's legend
leg_y  <- (ROW2_H_IN + LEGEND_OFFSET_IN) / FIG_H_IN
leg_h  <- 0.2 / FIG_H_IN
final_fig <- cowplot::ggdraw(final_layout) +
  cowplot::draw_grob(legend_C, x = x_legC, y = leg_y - leg_h / 2, width = x_legD - x_legC, height = leg_h) +
  cowplot::draw_grob(legend_D, x = x_legD, y = leg_y - leg_h / 2, width = 1 - x_legD, height = leg_h)

ggsave(file.path(fig_dir, "Figure4_11112026.pdf"), final_fig, width = FIG_W_IN, height = FIG_H_IN)
ggsave(file.path(fig_dir, "Figure4_11112026.png"), final_fig, width = FIG_W_IN, height = FIG_H_IN,
       dpi = 300, bg = "white")
cat("Saved Figure4_11112026.{pdf,png} to", fig_dir,
    sprintf("| %.1f x %.2f in\n", FIG_W_IN, FIG_H_IN))
