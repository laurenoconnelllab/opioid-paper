# FIGURE 2 -- coordination as a pair-level care phenotype
# Input : Supplementary/derived_tables/VisitDurations_byType_derived_v4.csv, Data/recording_metadata.csv
#         (and Figures/artwork/Figure2_firstpanel_v2.png if SHOW_PHOTOS is TRUE)
# Output: Figures/Figure2_11102026.pdf / .png
#         Supplementary/derived_tables/CarePCA_sixMetrics_loadings.csv, CarePCA_sixMetrics_scores.csv
#         Supplementary/derived_tables/CarePC1_classification.csv, CarePC1_rule.csv
#         Supplementary/derived_tables/CareSummary_bySex.csv, AccumulatingVisits_bySex.csv

suppressPackageStartupMessages({
  library(dplyr); library(tidyr)
  library(ggplot2); library(patchwork); library(lme4)
})

# ---------- Panel letters (change here to re-letter the whole figure) ----------
LETTER <- c(pca = "A", index = "B", growth = "C", photos = "D",
            father_visits  = "E", mother_visits  = "F", coord_visits  = "G",
            father_minutes = "H", mother_minutes = "I", coord_minutes = "J")

# ---------- Theme + palettes ----------
theme_nature <- theme_classic(base_size = 10, base_family = "sans") +
  theme(panel.grid = element_blank(),
        axis.line  = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.ticks.length   = unit(2, "pt"),
        plot.title          = element_text(face = "bold", size = 10, hjust = 0, margin = margin(b = 2)),
        plot.title.position = "plot",
        legend.background   = element_blank(),
        legend.key.size     = unit(0.75, "lines"),
        legend.text         = element_text(size = 9),
        axis.title          = element_text(size = 10),
        axis.text           = element_text(size = 9))
theme_set(theme_nature)

CARE_PALETTE <- c("Low" = "#B2E4DB", "High" = "#005F52")
SEX_CARE_PALETTE <- c(
  "Male.High" = "#0077BB",  "Male.Low" = "#7FBBDD",
  "Female.High" = "#EE7733","Female.Low" = "#F6BB99",
  "Family.High" = "#005F52","Family.Low" = "#B2E4DB")
sex_care_pal <- function(s) c("High" = SEX_CARE_PALETTE[[paste0(s,".High")]],
                               "Low"  = SEX_CARE_PALETTE[[paste0(s,".Low")]])

darken_hex <- function(hex, amount = 0.35) {
  m  <- col2rgb(hex)                    # rgb2hsv() expects 0-255, not 0-1 (fixed in v20)
  h  <- rgb2hsv(m[1, ], m[2, ], m[3, ])
  v2 <- pmax(h["v", ] * (1 - amount), 0)
  out <- hsv(h["h", ], h["s", ], v2)
  names(out) <- names(hex)
  out
}
CARE_PALETTE_DARK     <- darken_hex(CARE_PALETTE)
SEX_CARE_PALETTE_DARK <- darken_hex(SEX_CARE_PALETTE)
sex_care_pal_dark <- function(s) c("High" = SEX_CARE_PALETTE_DARK[[paste0(s,".High")]],
                                    "Low"  = SEX_CARE_PALETTE_DARK[[paste0(s,".Low")]])

# ---------- Where things live ----------

repo_dir <- here::here()
derived   <- file.path(repo_dir, "Supplementary", "derived_tables")
fig_dir   <- file.path(repo_dir, "Figures")
meta_file <- file.path(repo_dir, "Data", "recording_metadata.csv")
dir.create(derived, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# visits analysis
# written by Figure1_ethogram_panels_11162026.R
visit_file <- file.path(derived, "VisitDurations_byType_derived_v4.csv")
if (!file.exists(visit_file)) stop("Run Figure1_ethogram_panels_11162026.R first: it writes ", visit_file)

LONG_SOLO_VISIT_MIN <- 360                      # 6 h cutoff for individual visits. this is because there's a very strong possibility that parents may fall asleep when the lights turn off 
visitors <- c("Male", "Female", "Family")

visits <- read.csv(visit_file, stringsAsFactors = FALSE)
meta   <- read.csv(meta_file,  stringsAsFactors = FALSE)
recordings <- meta$ObservationID
full_days  <- setNames(meta$NumberofFullDays, recordings)
stopifnot(all(visits$tank %in% recordings))

# 6-variables
long_solo <- visits$visitor %in% c("Male", "Female") & visits$duration_min > LONG_SOLO_VISIT_MIN
kept      <- visits[!long_solo, ]
cat(sprintf("Left out of the per-day rates: %d of %d single-parent visits longer than %d min\n",
            sum(long_solo), nrow(visits), LONG_SOLO_VISIT_MIN))

#you should see that 28 uniparental visits are kicked off
count_by <- function(d) table(factor(d$tank, recordings), factor(d$visitor, visitors))
sum_by   <- function(d) {
  s <- tapply(d$duration_min, list(factor(d$tank, recordings), factor(d$visitor, visitors)), sum)
  s[is.na(s)] <- 0
  s
}
n_visits <- count_by(kept)
minutes  <- sum_by(kept)

care <- data.frame(ObservationID = recordings, Tank = meta$Pair,
                   NumberofFullDays = meta$NumberofFullDays, stringsAsFactors = FALSE)
for (v in visitors) care[[paste0("Visits_", v, "_perDay")]]  <- as.numeric(n_visits[, v]) / care$NumberofFullDays
for (v in visitors) care[[paste0("Minutes_", v, "_perDay")]] <- as.numeric(minutes[, v])  / care$NumberofFullDays

# 6-metric PCA
pop_sd      <- function(x) sqrt(mean((x - mean(x))^2))        # divide by n, not n - 1
z_score     <- function(X) sweep(sweep(X, 2, colMeans(X), "-"), 2, apply(X, 2, pop_sd), "/")
six_metrics <- c("Visits_Male_perDay", "Visits_Female_perDay", "Visits_Family_perDay",
                 "Minutes_Male_perDay", "Minutes_Female_perDay", "Minutes_Family_perDay")
pca6      <- prcomp(z_score(as.matrix(care[, six_metrics])), center = FALSE, scale. = FALSE)
pca6 #you should see that male visits epr day and family visits per day are around 0.52+
loadings6 <- pca6$rotation
scores6   <- pca6$x
if (loadings6["Visits_Family_perDay", 1] < 0) {               # PC1 points towards more coordinated visits
  loadings6[, 1] <- -loadings6[, 1]; scores6[, 1] <- -scores6[, 1]
}
var6 <- pca6$sdev^2 / sum(pca6$sdev^2)

cat(sprintf("Six-metric PCA: PC1 explains %.1f%% variance, PC2 %.1f%%. Loadings on PC1, largest first", #then run next line
            100 * var6[1], 100 * var6[2]))
print(round(sort(loadings6[, 1], decreasing = TRUE), 3))

#Visits_Family_perDay    Visits_Male_perDay   Minutes_Male_perDay Minutes_Female_perDay 
#0.535                 0.517                 0.406                 0.367 
#Visits_Female_perDay Minutes_Family_perDay 
#0.312                 0.223 

pca_loadings <- data.frame(Variable = six_metrics, PC1 = loadings6[, 1], PC2 = loadings6[, 2],
                            VarExplained_PC1 = var6[1], VarExplained_PC2 = var6[2], row.names = NULL)
pca_scores   <- data.frame(ObservationID = care$ObservationID, PC1 = scores6[, 1], PC2 = scores6[, 2])

# To keep each parent's solo care free to be tested against coordination, the CAREPC1 index uses the two coordinated rates only.
family_vars <- c("Visits_Family_perDay", "Minutes_Family_perDay")
X       <- as.matrix(care[, family_vars])
centre  <- colMeans(X)                                         # kept: the other cohorts are scored with these
spread  <- apply(X, 2, pop_sd)
Z       <- z_score(X)
pca     <- prcomp(Z, center = FALSE, scale. = FALSE)
loading <- pca$rotation[, 1]
if (loading[1] < 0) loading <- -loading                        # higher = more coordination
care$CarePC1 <- as.numeric(Z %*% loading)
var_explained <- pca$sdev^2 / sum(pca$sdev^2)
var_explained

#establishing the care rule from the carepc1 index
set.seed(42) #42 is canon in ML. there is no reason for it. comes from hitchhiker's guide to galaxy
km   <- kmeans(care$CarePC1, centers = 2, nstart = 20)
high <- which.max(tapply(care$CarePC1, km$cluster, mean))
care$CareLevel <- ifelse(km$cluster == high, "High", "Low")
cutoff <- (max(care$CarePC1[care$CareLevel == "Low"]) + min(care$CarePC1[care$CareLevel == "High"])) / 2
# cutoff is 0.1558056
care$bodylength    <- meta$BodyLength_mm
care$tadpoleweight <- meta$TadpoleWeight_g
care$gosnerstage   <- meta$GosnerStage
pca_scores <- pca_scores[order(care$Tank, care$ObservationID), ]
care <- care[order(care$Tank, care$ObservationID), ]
#keep in mind since you only uploaded 2 variables, the loadings are balanced
cat(sprintf("CarePC1 (coordinated visits and minutes only): explains %.1f%% of their variance; loadings %+.3f (visits), %+.3f (minutes)\n",
            100 * var_explained[1], loading[1], loading[2]))
cat(sprintf("Cutoff = %+.3f; High = %d, Low = %d\n", cutoff,
            sum(care$CareLevel == "High"), sum(care$CareLevel == "Low")))
#cutoff is 0.156, high = 11, low = 15.
rule <- data.frame(Variable = family_vars, Mean = centre, SD = spread, Loading = loading,
                   Threshold = c(cutoff, NA), VarExplained_PC1 = c(var_explained[1], NA),
                   row.names = NULL)

all_n   <- count_by(visits)
all_min <- sum_by(visits)
sm_long <- do.call(rbind, lapply(recordings, function(id) data.frame(
  ObservationID = id, Sex = visitors,
  TotalVisits = as.numeric(all_n[id, visitors]), TotalMinutes = as.numeric(all_min[id, visitors]),
  NumberofFullDays = full_days[[id]], stringsAsFactors = FALSE)))
sm_long$VisitsPerDay  <- sm_long$TotalVisits  / sm_long$NumberofFullDays
sm_long$MinutesPerDay <- sm_long$TotalMinutes / sm_long$NumberofFullDays

gap_days <- lapply(setNames(meta$GapDays, recordings), function(s) {
  if (is.na(s) || s == "") integer(0) else as.integer(strsplit(as.character(s), ";")[[1]])
})
last_day <- sapply(recordings, function(id)
  max(visits$day[visits$tank == id], full_days[[id]], na.rm = TRUE))
max_day  <- max(last_day)

acc_rows <- list()
for (id in recordings) for (v in visitors) {
  days    <- visits$day[visits$tank == id & visits$visitor == v]
  per_day <- tabulate(days[!is.na(days) & days >= 1], nbins = last_day[[id]])
  cum     <- as.numeric(cumsum(per_day))
  cum[seq_along(cum) %in% gap_days[[id]]] <- NA                # gap days are unknown, not zero
  cum     <- c(cum, rep(NA, max_day - length(cum)))
  acc_rows[[length(acc_rows) + 1]] <- data.frame(
    ObservationID = id, Sex = v, NumberofFullDays = full_days[[id]], TotalVisits = sum(per_day),
    setNames(as.list(cum), paste0("Day", seq_len(max_day))), stringsAsFactors = FALSE)
}
acc_wide <- do.call(rbind, acc_rows)

level <- care[, c("ObservationID", "CareLevel", "CarePC1")]
sm_long  <- merge(sm_long,  level, by = "ObservationID", sort = FALSE)
acc_wide <- merge(acc_wide, level, by = "ObservationID", sort = FALSE)

write.csv(pca_loadings, file.path(derived, "CarePCA_sixMetrics_loadings.csv"), row.names = FALSE)
write.csv(pca_scores,   file.path(derived, "CarePCA_sixMetrics_scores.csv"),   row.names = FALSE)
write.csv(care,         file.path(derived, "CarePC1_classification.csv"),      row.names = FALSE)
write.csv(rule,         file.path(derived, "CarePC1_rule.csv"),                row.names = FALSE)
write.csv(sm_long,      file.path(derived, "CareSummary_bySex.csv"),           row.names = FALSE)
write.csv(acc_wide,     file.path(derived, "AccumulatingVisits_bySex.csv"),    row.names = FALSE)

#figure 2
care$CareLevel     <- factor(care$CareLevel,     levels = c("Low","High"))
sm_long$CareLevel  <- factor(sm_long$CareLevel,  levels = c("Low","High"))
acc_wide$CareLevel <- factor(acc_wide$CareLevel, levels = c("Low","High"))
KMEANS_THRESHOLD   <- cutoff

#helper scripts
#fmt_p generates the p-values the figures will print. 
fmt_p   <- function(p) {
  if (is.na(p)) return("n.a.")
  if (p < 0.001) return(sprintf("%.2e", p))
  sprintf("%.3f", p)
}
#sig_code is giving you the significance asterisks
sig_code <- function(p) {
  if (is.na(p)) return("n.a.")
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "n.s."
}
# Minimal bottom-right Spearman-only callout used by the body-length panel.
# Built as a plotmath expression (parse=TRUE at the call site) so the
# rho glyph renders via the symbol font instead of a raw UTF-8 character,
# which this graphics device silently fails to convert (mbcsToSbcs).
spear_label_expr <- function(x, y, n) {
  rho <- cor(x, y, method="spearman")
  p   <- suppressWarnings(cor.test(x, y, method="spearman", exact=FALSE)$p.value)
  # Each line wrapped in displaystyle() so nesting atop() inside atop()
  # doesn't auto-shrink the inner two lines relative to the first --
  # otherwise "p = ..." renders visibly larger than "rho = ..."/"n = ..."
  # even though all three use the same "==" spacing in the source.
  sprintf("atop(displaystyle(italic(p) == %s), atop(displaystyle(rho == %.2f), displaystyle(italic(n) == %d)))",
          fmt_p(p), rho, n)
}

#panel A
load   <- pca_loadings; rownames(load) <- load$Variable
ev_pct <- 100 * c(load$VarExplained_PC1[1], load$VarExplained_PC2[1])
ev_pct
d_pca  <- merge(pca_scores, care[, c("ObservationID", "CareLevel")], by = "ObservationID", sort = FALSE)
d_pca
mw  <- suppressWarnings(wilcox.test(PC1 ~ CareLevel, data=d_pca))
#Wilcoxon rank sum exact test
#data:  PC1 by CareLevel
#W = 14, p-value = 0.0001297
#alternative hypothesis: true location shift is not equal to 0

.auc <- function(a,b) { p <- expand.grid(a=a,b=b); mean(p$a>p$b)+0.5*mean(p$a==p$b) }
auc_val <- .auc(d_pca$PC1[d_pca$CareLevel=="High"], d_pca$PC1[d_pca$CareLevel=="Low"])
#auc_val is 0.915

# Proper "1 x 10^-4"-style plotmath label instead of "p=1e-04"
pval <- mw$p.value
pca_stat_label <- if (pval < 1e-3) {
  pv_exp  <- floor(log10(pval))
  pv_mant <- pval / 10^pv_exp
  if (abs(pv_mant - round(pv_mant)) < 0.05) {
    sprintf("italic(p) == %.0f %%*%% 10^%d", round(pv_mant), pv_exp)
  } else {
    sprintf("italic(p) == %.1f %%*%% 10^%d", pv_mant, pv_exp)
  }
} else {
  sprintf("italic(p) == %.3f", pval)
}
pca_stat_label <- c(sprintf("AUC == %.2f", max(auc_val, 1 - auc_val)), pca_stat_label)   # two lines: AUC, then p

cat(sprintf("[Panel %s | Six-metric PCA: PC1 %.1f%%, PC2 %.1f%% of variance; loadings on PC1: %s]\n",
    LETTER[["pca"]], ev_pct[1], ev_pct[2],
    paste(sprintf("%s %.2f", load$Variable, load$PC1)[order(-load$PC1)], collapse = ", ")))
#Six-metric PCA: PC1 37.7%, PC2 24.6% of variance; loadings on PC1: Visits_Family_perDay 0.53, Visits_Male_perDay 0.52, Minutes_Male_perDay 0.41, Minutes_Female_perDay 0.37, Visits_Female_perDay 0.31, Minutes_Family_perDay 0.22]
cat(sprintf("  High vs Low care (CarePC1) on PC1 (n=%d): AUC=%.3f, Wilcoxon p=%.3g\n\n",
    nrow(d_pca), max(auc_val, 1-auc_val), mw$p.value))

SHOW_VARS  <- c("Visits_Male_perDay","Visits_Female_perDay","Visits_Family_perDay")
SHOW_LBLS  <- c("Visits Male","Visits Female","Visits Family")
HIGHLIGHT  <- c("Visits Male","Visits Family")
ARROW_SCALE <- 2.4
loadings_show <- data.frame(var=SHOW_LBLS, L1=load[SHOW_VARS,"PC1"], L2=load[SHOW_VARS,"PC2"]) |>
  mutate(highlight=var %in% HIGHLIGHT,
         x_end=L1*ARROW_SCALE, y_end=L2*ARROW_SCALE,
         x_lab=x_end*1.18,     y_lab=y_end*1.18,
         # the father and coordinated arrows nearly coincide: move their labels apart
         y_lab=y_lab + c("Visits Male"=-0.75, "Visits Female"=0.35, "Visits Family"=0.45)[var])

panel_pca <- ggplot(d_pca, aes(PC1, PC2)) +
  stat_ellipse(aes(fill=CareLevel, color=CareLevel),
               geom="polygon", level=0.95, alpha=0.18, linewidth=0.5) +
  geom_hline(yintercept=0, color="grey80", linewidth=0.2) +
  geom_vline(xintercept=0, color="grey80", linewidth=0.2) +
  geom_point(aes(fill=CareLevel, color=CareLevel), shape=21, size=2.4, stroke=0.4) +
  geom_segment(data=loadings_show,
               aes(x=0, y=0, xend=x_end, yend=y_end,
                   color=ifelse(highlight,"highlight","other"),
                   linewidth=ifelse(highlight,1.0,0.5)),
               arrow=arrow(length=unit(0.13,"cm"), type="closed"),
               inherit.aes=FALSE, show.legend=FALSE) +
  geom_label(data=loadings_show,
             aes(x=x_lab, y=y_lab, label=var,
                 fontface=ifelse(highlight,"bold","plain")),
             color="black", fill="white", size=2.5,
             label.size=0.2, label.padding=unit(0.10,"lines"),
             inherit.aes=FALSE, show.legend=FALSE) +
  scale_fill_manual(values=CARE_PALETTE, guide="none") +
  scale_color_manual(values=c(highlight="#CC3311", other="black",
                                Low=CARE_PALETTE_DARK[["Low"]], High=CARE_PALETTE_DARK[["High"]]),
                      guide="none") +
  scale_linewidth_identity() +
  scale_y_continuous(expand=expansion(mult=c(0.05, 0.45))) +          # head-room for the AUC and p labels
  annotate("text", x=-Inf, y=Inf, hjust=-0.08, vjust=c(1.4, 2.7),
           label=pca_stat_label, parse=TRUE, size=3.0, color="black") +
  labs(x=sprintf("PC1: Coordination\nintensity (%.1f%%)",ev_pct[1]),          # plain ":" (a long dash does not print in the PDF)
       y=sprintf("PC2: Maternal- vs.\npaternal-weighted (%.1f%%)",ev_pct[2]),
       title=paste0(LETTER[["pca"]], ".")) +
  theme(legend.position="none",
        axis.title.x=element_text(size=7.3, lineheight=0.95),
        axis.title.y=element_text(size=7.3, lineheight=0.95))
#panel b
d_index <- care |> select(ObservationID, CarePC1, CareLevel) |> arrange(CarePC1) |>
  mutate(ObservationID = factor(ObservationID, levels=ObservationID))
index_range <- range(d_index$CarePC1)
panel_index <- ggplot(d_index, aes(ObservationID, CarePC1, fill=CareLevel, color=CareLevel)) +
  geom_col(linewidth=0.3) +
  geom_hline(yintercept=KMEANS_THRESHOLD, color="red", linetype="dashed", linewidth=0.5) +
  annotate("text", x=0.6, y=KMEANS_THRESHOLD+0.12,
           label=sprintf("cutoff\n%+.2f", KMEANS_THRESHOLD), color="red", size=2.8,
           hjust=0, vjust=0, lineheight=0.9) +
  annotate("text", x=nlevels(d_index$ObservationID)*0.18, y=index_range[2]*0.82,
           label="High", color=CARE_PALETTE_DARK[["High"]],
           fontface="bold", size=2.9, hjust=0) +
  annotate("segment", x=nlevels(d_index$ObservationID)*0.10, xend=nlevels(d_index$ObservationID)*0.10,
           y=index_range[2]*0.45, yend=index_range[2]*0.78,
           color=CARE_PALETTE_DARK[["High"]], linewidth=0.6,
           arrow=arrow(length=unit(0.10,"cm"), type="closed")) +
  annotate("text", x=nlevels(d_index$ObservationID)*0.70, y=index_range[1]*0.80,
           label="Low", color=CARE_PALETTE_DARK[["Low"]],
           fontface="bold", size=2.9, hjust=0) +
  annotate("segment", x=nlevels(d_index$ObservationID)*0.64, xend=nlevels(d_index$ObservationID)*0.64,
           y=index_range[1]*0.40, yend=index_range[1]*0.75,
           color=CARE_PALETTE_DARK[["Low"]], linewidth=0.6,
           arrow=arrow(length=unit(0.10,"cm"), type="closed")) +
  scale_fill_manual(values=CARE_PALETTE, guide="none") +
  scale_color_manual(values=CARE_PALETTE_DARK, guide="none") +
  labs(x="Recordings", y="CarePC1", title=paste0(LETTER[["index"]], ".")) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

# panel c
d_growth <- care |> select(Minutes_Family_perDay, bodylength, CareLevel) |> drop_na()
panel_growth <- ggplot(d_growth, aes(Minutes_Family_perDay, bodylength)) +
  geom_smooth(method="lm", se=FALSE, color="black", linetype="dashed", linewidth=0.5) +
  geom_point(aes(fill=CareLevel, color=CareLevel), shape=21, size=2.6, stroke=0.4) +
  annotate("text", x=Inf, y=-Inf, hjust=1.05, vjust=-0.15,
           label=spear_label_expr(d_growth$Minutes_Family_perDay, d_growth$bodylength, nrow(d_growth)),
           parse=TRUE, size=2.9, lineheight=0.9) +
  scale_fill_manual(values=CARE_PALETTE, guide="none") +
  scale_color_manual(values=CARE_PALETTE_DARK, guide="none") +
  labs(x="Coordinated visit\nduration (min/day)", y="Body length (mm)", title=paste0(LETTER[["growth"]], ".")) +
  theme(legend.position="none")

cat(sprintf(paste0("[Panel ", LETTER[["growth"]], " | Body length x coordinated minutes per day  (n=%d)]\n",
    "  Pearson:  r=%+.3f, p=%.3g\n",
    "  Spearman: rho=%+.3f, p=%.3g\n\n"),
  nrow(d_growth),
  cor(d_growth$Minutes_Family_perDay, d_growth$bodylength),
  cor.test(d_growth$Minutes_Family_perDay, d_growth$bodylength)$p.value,
  cor(d_growth$Minutes_Family_perDay, d_growth$bodylength, method="spearman"),
  suppressWarnings(cor.test(d_growth$Minutes_Family_perDay, d_growth$bodylength,
                            method="spearman", exact=FALSE)$p.value)))

#[Panel C | Body length x coordinated minutes per day  (n=23)]
#Pearson:  r=+0.423, p=0.0442
#Spearman: rho=+0.579, p=0.00381 


# panel D 
photo_file <- file.path(repo_dir, "Figures", "artwork", "Figure2_firstpanel_v2.png")
photo_panel <- function(piece, letter = NULL) {
  p <- ggplot() + theme_void()
  if (!is.null(piece))
    p <- p + annotation_custom(grid::rasterGrob(piece, height = grid::unit(1, "npc"), interpolate = TRUE))
  if (!is.null(letter))
    p <- p + labs(title = paste0(letter, ".")) +
      theme(plot.title = element_text(face = "bold", size = 10, hjust = 0, margin = margin(b = 2)))
  p
}
SHOW_PHOTOS <- FALSE   # FALSE: leave a placeholder
                       # TRUE : draw the photos from the artwork file
if (SHOW_PHOTOS && file.exists(photo_file)) {
  strip  <- png::readPNG(photo_file)
  inked  <- apply(strip[, , 1:3], 2, function(column) any(column < 0.98))   # columns that are not blank white
  runs   <- rle(inked)
  ends   <- cumsum(runs$lengths); starts <- ends - runs$lengths + 1
  keep   <- which(runs$values)
  stopifnot(length(keep) == 3)                                              # father, mother, coordinating
  pieces <- lapply(keep, function(i) strip[, starts[i]:ends[i], , drop = FALSE])
} else {
  pieces <- list(NULL, NULL, NULL)                                          # blank row, lettered D
}
panel_photo_father <- photo_panel(pieces[[1]], LETTER[["photos"]])
panel_photo_mother <- photo_panel(pieces[[2]])
panel_photo_coord  <- photo_panel(pieces[[3]])

# panels e-j
DAY_LIM  <- 30
day_cols <- grep("^Day[0-9]+$", names(acc_wide), value=TRUE)
acc_long <- acc_wide |>
  select(ObservationID, Sex, CareLevel, all_of(day_cols)) |>
  pivot_longer(all_of(day_cols), names_to="Day", values_to="Cum") |>
  mutate(Day=as.integer(sub("Day","",Day)),
         Sex=factor(Sex, levels=c("Male","Female","Family"))) |>
  filter(Day <= DAY_LIM)

acc_summary <- acc_long |>
  group_by(Sex, CareLevel, Day) |>
  summarise(med=median(Cum,na.rm=TRUE),
            q1=quantile(Cum,0.25,na.rm=TRUE),
            q3=quantile(Cum,0.75,na.rm=TRUE), .groups="drop")

acc_mm_dat <- acc_long |>
  filter(!is.na(Cum)) |>
  arrange(ObservationID, Sex, Day) |>
  group_by(ObservationID, Sex) |>
  mutate(VisitsPerDay = pmax(Cum - dplyr::lag(Cum, default=0), 0)) |>
  ungroup() |>
  mutate(DayC=Day-1, CareLevel=factor(CareLevel, levels=c("Low","High")))

# Visits per day: negative-binomial mixed model per visitor.
#   main model:        visits ~ care level + day + (1 + day | recording)
#   interaction model: visits ~ care level * day + (1 + day | recording)
cat("[Per-day visits | NB-GLMM -- CarePC1]\n")
glmm_stats <- data.frame(Sex=c("Male","Female","Family"),
                          p_care=NA_real_, RR_care=NA_real_,
                          p_day=NA_real_,  p_inter=NA_real_,
                          stringsAsFactors=FALSE)

for (sx in c("Male","Female","Family")) {
  d_sx <- acc_mm_dat |> filter(Sex == sx)

  fit_main <- suppressMessages(suppressWarnings(
    glmer.nb(VisitsPerDay ~ CareLevel + DayC + (1 + DayC | ObservationID),
             data=d_sx, control=glmerControl(optimizer="bobyqa",
                                              optCtrl=list(maxfun=200000)))))
  cf <- summary(fit_main)$coefficients
  i_care <- grep("^CareLevelHigh$", rownames(cf))[1]
  i_day  <- grep("^DayC$",          rownames(cf))[1]

  fit_inter <- tryCatch(suppressMessages(suppressWarnings(
    glmer.nb(VisitsPerDay ~ CareLevel * DayC + (1 + DayC | ObservationID),
             data=d_sx, control=glmerControl(optimizer="bobyqa",
                                              optCtrl=list(maxfun=200000))))),
    error=function(e) NULL)
  p_inter <- if (!is.null(fit_inter)) {
    cf_i <- summary(fit_inter)$coefficients
    cf_i[grep(":", rownames(cf_i))[1], "Pr(>|z|)"]
  } else NA_real_

  glmm_stats[glmm_stats$Sex==sx,] <- list(sx,
    cf[i_care,"Pr(>|z|)"], exp(cf[i_care,"Estimate"]),
    cf[i_day, "Pr(>|z|)"], p_inter)

  cat(sprintf("  %-6s  Care: RR=%.2fx, p=%.3g | Day: p=%.3g | Day x Care: p=%.3g\n",
    sx, exp(cf[i_care,"Estimate"]), cf[i_care,"Pr(>|z|)"],
    cf[i_day,"Pr(>|z|)"], p_inter))
}

get_p_care  <- function(s) glmm_stats$p_care[ glmm_stats$Sex==s]
get_p_day   <- function(s) glmm_stats$p_day[  glmm_stats$Sex==s]
get_p_inter <- function(s) glmm_stats$p_inter[ glmm_stats$Sex==s]

# Cumulative visits over 30 days (median and quartiles by care level)
panel_traj <- function(sex_lbl, letter, pretty=sex_lbl) {
  d   <- acc_summary |> filter(Sex==sex_lbl)
  pal <- sex_care_pal(sex_lbl)
  pmm <- get_p_care(sex_lbl)
  dp  <- get_p_day(sex_lbl);  ip <- get_p_inter(sex_lbl)
  stats_lab <- paste(
    sprintf("Care %s", sig_code(pmm)),
    sprintf("Day %s", sig_code(dp)),
    sprintf("Day \u00d7 Care %s", sig_code(ip)),                               # \u00d7 = multiplication sign, written this way so it prints in the PDF on any computer
    sep="\n")
  ggplot(d, aes(Day, med, color=CareLevel, fill=CareLevel, linetype=CareLevel)) +
    geom_ribbon(aes(ymin=q1, ymax=q3), color=NA, alpha=0.25) +
    geom_line(linewidth=0.9) +
    annotate("text", x=1, y=Inf, hjust=0, vjust=1.3,
             label=stats_lab, size=2.9, fontface="bold", lineheight=1.25) +
    scale_color_manual(values=pal, guide="none") +
    scale_fill_manual(values=pal, guide="none") +
    scale_linetype_manual(values=c(High="solid",Low="dashed"), guide="none") +
    scale_x_continuous(breaks=c(1,10,20,30), limits=c(1,DAY_LIM)) +
    scale_y_continuous(expand=expansion(mult=c(0.05, 0.50))) +             # head-room for the three test results
    labs(x="Day", y="Cumulative visits", title=sprintf("%s.", letter))
}

# Minutes per day by care level

panel_min_box <- function(sex_lbl, letter, pretty=sex_lbl) {
  d  <- sm_long |> filter(Sex==sex_lbl) |> drop_na(MinutesPerDay)
  hi <- d$MinutesPerDay[d$CareLevel=="High"]; lo <- d$MinutesPerDay[d$CareLevel=="Low"]
  pw <- suppressWarnings(wilcox.test(d$MinutesPerDay ~ d$CareLevel)$p.value)
  sig  <- sig_code(pw)
  ymax <- max(d$MinutesPerDay, na.rm=TRUE)
  pal <- sex_care_pal(sex_lbl); pal_dark <- sex_care_pal_dark(sex_lbl)
  ticks <- c(Low=sprintf("Low\n(n=%d)",length(lo)), High=sprintf("High\n(n=%d)",length(hi)))
  ggplot(d, aes(CareLevel, MinutesPerDay, fill=CareLevel, color=CareLevel)) +
    geom_boxplot(width=0.55, alpha=0.55, outlier.shape=NA, linewidth=0.6, fatten=1.4) +
    geom_point(shape=21, stroke=0.4, size=1.6,
               position=position_jitter(width=0.07, height=0, seed=1)) +   # seed: same sideways jitter every run
    annotate("text", x=1.5, y=ymax*1.08, label=sig, size=4.5, fontface="bold") +
    scale_fill_manual(values=pal, guide="none") +
    scale_color_manual(values=pal_dark, guide="none") +
    scale_y_continuous(expand=expansion(mult=c(0.05,0.14))) +
    scale_x_discrete(labels=ticks) +
    labs(x=NULL, y=sprintf("%s minutes / day", pretty), title=sprintf("%s.", letter))
}

panel_father_visits  <- panel_traj("Male",      LETTER[["father_visits"]],  "Father")
panel_father_minutes <- panel_min_box("Male",   LETTER[["father_minutes"]], "Father")
panel_mother_visits  <- panel_traj("Female",    LETTER[["mother_visits"]],  "Mother")
panel_mother_minutes <- panel_min_box("Female", LETTER[["mother_minutes"]], "Mother")
panel_coord_visits   <- panel_traj("Family",    LETTER[["coord_visits"]],   "Family")
panel_coord_minutes  <- panel_min_box("Family", LETTER[["coord_minutes"]],  "Family")

cat("[Minutes/day | Wilcoxon rank-sum -- CarePC1]\n")
for (box_row in list(c("Male",LETTER[["father_minutes"]],"Father"), c("Female",LETTER[["mother_minutes"]],"Mother"), c("Family",LETTER[["coord_minutes"]],"Family"))) {
  d_box  <- sm_long |> filter(Sex == box_row[1]) |> drop_na(MinutesPerDay)
  hi_box <- d_box$MinutesPerDay[d_box$CareLevel == "High"]
  lo_box <- d_box$MinutesPerDay[d_box$CareLevel == "Low"]
  pw_box <- suppressWarnings(wilcox.test(d_box$MinutesPerDay ~ d_box$CareLevel)$p.value)
  d_eff  <- (mean(hi_box) - mean(lo_box)) /
    sqrt(((length(hi_box)-1)*var(hi_box) + (length(lo_box)-1)*var(lo_box)) /
         (length(hi_box) + length(lo_box) - 2))
  cat(sprintf("  %s (Panel %s): p=%.3g, d=%+.2f  %s\n",
      box_row[3], box_row[2], pw_box, d_eff, sig_code(pw_box)))
}
cat("\n")

row1 <- (panel_pca | panel_index | panel_growth) + plot_layout(widths=c(1.8, 1.2, 1.5))

# Rows 2 to 4: one column per visit type -- photo, cumulative visits, minutes per day. Built as one grid so each photo is centered over its own plots
columns <- (panel_photo_father  | panel_photo_mother  | panel_photo_coord) /
           (panel_father_visits | panel_mother_visits | panel_coord_visits) /
           (panel_father_minutes | panel_mother_minutes | panel_coord_minutes) +
  plot_layout(heights=c(1.6, 2.2, 2.0))

# Sizes in inches
FIG_W <- 5.2043       # inches (was 8.5)
ROW_H <- c(index = 2.5, columns = 1.6 + 2.2 + 2.0)
fig2 <- (wrap_elements(full = row1) / wrap_elements(full = columns)) +
  plot_layout(heights = ROW_H) & theme(plot.margin=margin(3,5,3,5))

fig2

ggsave(file.path(fig_dir,"Figure2_11102026.png"), fig2, width=FIG_W, height=sum(ROW_H), dpi=300, bg="white")
ggsave(file.path(fig_dir,"Figure2_11102026.pdf"), fig2, width=FIG_W, height=sum(ROW_H))
cat("Saved Figure2_11102026 to", fig_dir, "\n")
