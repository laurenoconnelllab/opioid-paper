#Figure5 
# The micrographs are not drawn here: the top row places the finished image Figures/artwork/Figure5_micrographs.png.
# Bar annotation convention:
#   bold star  = BH-corrected q (q_MW) within that metric's regional screen
#   grey (star) in parens = uncorrected p (p_MW)
#   S = ANOVA Sex main effect p<0.05   I = ANOVA Sex x Condition p<0.05
# Dot convention (every per-animal/per-pair panel):
#   filled, colored by Sex  = Coord
#   hollow, outline colored by Sex = Control
#   Sex: Male = blue, Female = orange, Mixed (averaged breeding pair) = grey
# Pairs: mates are not independent, so the CarePC1 panels use one point per
# breeding pair (the mean of its two animals, or the one animal that was imaged). Which two samples form a pair comes from Supplementary/derived_tables/RNAscope_CarePC1_perPair.csv.
# Input:
#   Data/rnascope/Cleaned_Data_ForBillie_9-22-2026.csv   one row per image:
#       cell counts per marker and marker combination, and the area counted
#   Data/rnascope/rnascope_pairs.csv, Data/rnascope/event_logs/*.xlsx
#       the breeding pairs and their nest-visit event logs
#   Data/pharmacology/opioid-behavior-spreadsheet.xlsx  sheet "rawdata"
#       naltrexone experiment, one row per nest visit
#   Supplementary/derived_tables/CarePC1_rule.csv   the CarePC1 rule (Figure2_11162026.R; run that first)
#   Figures/artwork/Figure5_naltrexone_timeline_stacked.png   the injection timeline drawing
#   Figures/artwork/Figure5_micrographs.png   the micrograph panel
#   (if one of the two images is missing, a labelled placeholder is drawn in
#   its place and the rest of the figure is still made)
# Output:
#   Figures/Figure5_11112026.pdf / .png
#   Supplementary/derived_tables/RNAscope_CarePC1_perPair.csv          CarePC1 per breeding pair
#   Supplementary/derived_tables/Pharmacology_perTank.csv              naltrexone experiment, per tank
#   Supplementary/derived_tables/Figure5_region_stats_allmetrics.csv   per-region statistics, all 4 metrics
#   Supplementary/derived_tables/Figure5_CarePC1_region_stats.csv      pair-level Spearman, every region
#
# Open the .Rproj file (or set the working directory to the project folder) before running.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(patchwork); library(lmerTest)
})
JITTER_SEED <- 1   # fixes the sideways jitter of the points: every plot, the PDF and the PNG get the same positions on every run

repo_dir <- here::here()
fig_dir  <- file.path(repo_dir, "Figures")
out_dir  <- file.path(repo_dir, "Supplementary", "derived_tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

CSV   <- file.path(repo_dir, "Data", "rnascope", "Cleaned_Data_ForBillie_9-22-2026.csv")
PAIRS <- file.path(out_dir, "RNAscope_CarePC1_perPair.csv")
stopifnot(file.exists(CSV))

# =============================================================================
# 0. CarePC1 for the two cohorts of this figure. Both are scored with the rule
#    fitted on the 26 reference recordings (Supplementary/derived_tables/CarePC1_rule.csv), so High
#    and Low mean the same thing as in Figure 2.
#    Writes Supplementary/derived_tables/RNAscope_CarePC1_perPair.csv and Pharmacology_perTank.csv.
# =============================================================================
local({
  raw_dir <- file.path(repo_dir, "Data")
  source(file.path(repo_dir, "Functions", "shared_functions.R"), local = TRUE)
  rule <- load_care_rule(file.path(out_dir, "CarePC1_rule.csv"))
  LONG_SOLO_VISIT_MIN <- 360                     # 6 h, as for the reference recordings

  # ---- RNAscope breeding pairs ---------------------------------------------------
  # rnascope_pairs.csv lists each breeding pair, its two brain samples, its event
  # log and the number of days in its recording window. The event logs are read
  # with exactly the same rules as the reference recordings
  # (shared_functions.R). Coordinated visits and minutes are divided by the
  # window length.
  pairs <- read.csv(file.path(raw_dir, "rnascope", "rnascope_pairs.csv"), stringsAsFactors = FALSE)
  log_dir <- file.path(raw_dir, "rnascope", "event_logs")

  family_totals <- function(event_log) {
    v <- read_visits(event_log, log_dir)
    family <- v[v$visitor == "Family", ]
    c(Visits_Family = nrow(family), Minutes_Family = sum(family$duration_min))
  }
  totals <- t(vapply(pairs$EventLog, family_totals, numeric(2)))
  pairs$Visits_Family  <- totals[, "Visits_Family"]
  pairs$Minutes_Family <- totals[, "Minutes_Family"]
  pairs$Visits_Family_perDay  <- pairs$Visits_Family  / pairs$WindowDays
  pairs$Minutes_Family_perDay <- pairs$Minutes_Family / pairs$WindowDays
  pairs$CarePC1   <- rule$care_pc1(pairs$Visits_Family_perDay, pairs$Minutes_Family_perDay)
  pairs$CareLevel <- rule$care_level(pairs$CarePC1)
  write.csv(pairs[, c("Pair", "Tank", "Cohort", "Sample1", "Sample2", "WindowDays",
                      "Visits_Family", "Minutes_Family", "Visits_Family_perDay",
                      "Minutes_Family_perDay", "CarePC1", "CareLevel")],
            PAIRS, row.names = FALSE)
  cat(sprintf("RNAscope: %d breeding pairs scored\n", nrow(pairs)))

  # ---- Naltrexone experiment: seven days per treatment --------------------------
  # One row per visit. Treatment F = saline, E = naltrexone; each tank received
  # both. Minutes per day = total minutes / 7 days. As for the reference
  # recordings, single-parent visits longer than 6 h are left out. Each tank's
  # CarePC1 and care level come from its saline week (its baseline). Tank B16 has
  # no parental visits and is left out.
  N_DAYS <- 7
  pharm <- suppressMessages(as.data.frame(
    read_excel(file.path(raw_dir, "pharmacology", "opioid-behavior-spreadsheet.xlsx"), sheet = "rawdata")))
  pharm <- data.frame(Tank = pharm$`Tank ID`,
                      Treatment = ifelse(pharm$Treatment == "F", "Saline", "Naltrexone"),
                      Visitor = pharm$Sex, minutes = as.numeric(pharm$TimeDifference_inminutes),
                      stringsAsFactors = FALSE)
  pharm <- pharm[pharm$Tank != "B16", ]
  pharm <- pharm[!(pharm$Visitor %in% c("Male", "Female") & pharm$minutes > LONG_SOLO_VISIT_MIN), ]

  ph <- expand.grid(Tank = sort(unique(pharm$Tank)), Treatment = c("Saline", "Naltrexone"),
                    Visitor = c("Male", "Female", "Family"), stringsAsFactors = FALSE)
  key <- function(d) paste(d$Tank, d$Treatment, d$Visitor)
  ph$Visits  <- as.numeric(table(factor(key(pharm), key(ph))))
  ph$Minutes <- as.numeric(tapply(pharm$minutes, factor(key(pharm), key(ph)), sum))
  ph$Minutes[is.na(ph$Minutes)] <- 0
  ph$VisitsPerDay  <- ph$Visits  / N_DAYS
  ph$MinutesPerDay <- ph$Minutes / N_DAYS

  baseline <- ph[ph$Treatment == "Saline" & ph$Visitor == "Family", ]
  baseline$CarePC1   <- rule$care_pc1(baseline$VisitsPerDay, baseline$MinutesPerDay)
  baseline$CareLevel <- rule$care_level(baseline$CarePC1)
  ph <- merge(ph, baseline[, c("Tank", "CarePC1", "CareLevel")], by = "Tank")
  ph <- ph[order(ph$Tank, ph$Treatment, ph$Visitor), ]
  write.csv(ph, file.path(out_dir, "Pharmacology_perTank.csv"), row.names = FALSE)
  cat(sprintf("Naltrexone experiment: %d tanks; saline-week care level High = %d, Low = %d\n",
              nrow(baseline), sum(baseline$CareLevel == "High"), sum(baseline$CareLevel == "Low")))
})

# ---------- sizes ----------
# Final printed sizes (the figure is saved FIG_W = 14 in wide). Text in points,
# line widths in ggplot units. Change them here to restyle the whole figure.
PT_BASE   <- 13.8     # theme base size, plot titles, legend keys
PT_TITLE  <- 12.08    # axis titles, facet titles
PT_TEXT   <- 10.35    # tick labels, legend text, text inside the naltrexone panels
PT_TEXT_X <- 11.21    # x tick labels of the box plots in B and F
LW_AXIS   <- 0.49     # axis lines, ticks, bar outlines
LW_THIN   <- 0.37     # box outlines and brackets
MARGIN_BOX <- margin(t = 1.23, r = 12.32, b = 4.93, l = 4.93)   # around the box plots in B and F

# ---------- palettes & theme ----------
SEX_PALETTE <- c(Male = "#0077BB", Female = "#EE7733", Mixed = "grey45")
DIR_PAL     <- c("Coord higher" = "#1C2B39", "Control higher" = "white")  # ink navy

theme_nature <- theme_classic(base_size = PT_BASE, base_family = "sans") +
  theme(
    panel.grid        = element_blank(),
    axis.line         = element_line(linewidth = LW_AXIS),
    axis.ticks        = element_line(linewidth = LW_AXIS),
    axis.ticks.length = unit(2.46, "pt"),
    axis.title        = element_text(size = PT_TITLE),
    axis.text         = element_text(size = PT_TEXT),
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold", size = PT_TITLE),
    legend.background = element_blank(),
    legend.key.size   = unit(PT_BASE, "pt"),
    legend.text       = element_text(size = PT_TEXT),
    legend.title      = element_text(size = PT_TEXT),
    plot.title        = element_text(size = PT_BASE, face = "plain"))
theme_set(theme_nature)

# =============================================================================
# 1. Load & tidy — one row per SampleID x Sex x Condition x BrainRegion,
#    all 4 metrics. Region-code duplicates collapsed; Fos density excludes
#    images with a blank area field (see comment below).
# =============================================================================
collapse_map <- c(LP2 = "LP", LS2 = "LS", MP2 = "MP", MPOA2 = "MPOA",
                  NAC = "NACC", DPVPV = "DPVP",
                  STR2 = "STR", BST2 = "BST",
                  APOA2 = "APOA", DP2 = "DP",
                  DPVP2 = "DPVP", NACC2 = "NACC",
                  VH. = "VH", MERGED_VH = "VH")

dat <- read_csv(CSV, show_col_types = FALSE, name_repair = "unique") %>%
  rename(SampleID    = `Sample ID...1`,
         BrainRegion = brainregion,
         Fos_Penk    = `Fos+Penk`,
         Fos_Gal     = `Fos+Gal`,
         Penk_Gal    = `Penk + Gal`,
         All3        = `All.3`,
         Area        = area)

# Counts are nested: a cell positive for all three markers is also counted as Fos+Penk, Fos+Gal and Penk+Gal. 
#The syn37 APOA row is the sum of six images of File8; in one of them (Image33) the counting workbook has 18 triple-positive cells but 0 Fos+Penk cells, so the Fos+Penk type was not saved for that image. Its Fos+Penk count is set to the minimum it can be, 18 (the triple-positive cells), which makes the row's Fos+Penk 29 instead of 11.
fix_row <- which(dat$SampleID == "syn37" & dat$BrainRegion == "APOA" &
                 dat$Fos_Penk == 11 & dat$All3 == 20)
if (length(fix_row) == 1) {
  dat$Fos_Penk[fix_row] <- 29
  cat("Correction: syn37 APOA Fos+Penk 11 -> 29 (File8 Image33: Fos+Penk set to its minimum, the 18 triple-positive cells)\n")
} else {
  cat("Correction for syn37 APOA not applied: the row is not as expected (already corrected in the data file?)\n")
}

dat <- dat %>%
  mutate(Condition   = recode(str_to_title(str_trim(Condition)), Treated = "Coord"),
         Sex         = str_to_title(str_trim(Sex)),
         BrainRegion = recode(BrainRegion, !!!collapse_map),
         across(c(Fos, Penk, Gal, Fos_Penk, Fos_Gal, Penk_Gal, All3, Area), as.numeric)) %>%
  group_by(SampleID, Sex, Condition, BrainRegion) %>%
  # Fos density uses only images that HAVE an area
  # The ratio metrics don't need area, so they keep every image.
  summarise(Fos_wArea = sum(Fos[!is.na(Area)], na.rm = TRUE),
            across(c(Fos, Penk, Gal, Fos_Penk, Fos_Gal, Penk_Gal, All3, Area),
                   ~ sum(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(Fos_density   = ifelse(Area > 0,     Fos_wArea / Area * 1e6, NA_real_),
         FosPenk_ratio = ifelse(Penk > 0,     Fos_Penk / Penk,       NA_real_),
         FosGal_ratio  = ifelse(Gal > 0,      Fos_Gal  / Gal,        NA_real_),
         All3_PenkGal  = ifelse(Penk_Gal > 0, All3     / Penk_Gal,   NA_real_))

cat(sprintf("Loaded %d animal x region rows, %d animals, %d regions\n",
            nrow(dat), n_distinct(dat$SampleID), n_distinct(dat$BrainRegion)))

#loaded 289 animal x region rows, 33 animals, 16 regions
# =============================================================================
# 2. Per-region nonparametric stats, all 4 metrics:
#    Mann-Whitney (Coord vs Control) + Hodges-Lehmann shift + CI,
#    BH correction across regions WITHIN each metric's screen,
#    plus a 2-way ANOVA (Sex x Condition) per region for the S / I tags.
#    Regions are only included where every Sex x Condition cell has n >= 3.
# =============================================================================
per_region_HL <- function(d, value_col, min_n = 2) {
  d %>%
    filter(!is.na(.data[[value_col]])) %>%
    group_by(BrainRegion) %>%
    group_modify(~ {
      dd <- .x
      x <- dd[[value_col]][dd$Condition == "Coord"]
      y <- dd[[value_col]][dd$Condition == "Control"]
      n_Coord <- length(x); n_Control <- length(y)
      if (n_Coord < min_n || n_Control < min_n) return(tibble::tibble())
      wt <- suppressWarnings(wilcox.test(x, y, conf.int = TRUE, exact = FALSE))
      tibble::tibble(n_Coord = n_Coord, n_Control = n_Control,
                     mean_Coord = mean(x), mean_Control = mean(y),
                     HL = unname(wt$estimate), HL_lo = wt$conf.int[1], HL_hi = wt$conf.int[2],
                     p_MW = wt$p.value)
    }) %>% ungroup() %>%
    mutate(metric = value_col, q_MW = p.adjust(p_MW, method = "BH"),
           star = case_when(p_MW < 0.001 ~ "***", p_MW < 0.01 ~ "**",
                            p_MW < 0.05 ~ "*", p_MW < 0.10 ~ ".", TRUE ~ ""))
}
anova_terms <- function(data, value_col, min_cell_n = 2) {
  data %>%
    filter(!is.na(.data[[value_col]])) %>%
    group_by(BrainRegion) %>%
    group_modify(~ {
      d <- .x
      cell_n <- d %>% count(Sex, Condition)
      ok <- nrow(cell_n) == 4 && all(cell_n$n >= min_cell_n)
      if (!ok) return(tibble::tibble(p_sex = NA_real_, p_cond = NA_real_, p_int = NA_real_))
      fit <- tryCatch(lm(d[[value_col]] ~ Sex * Condition, data = d), error = function(e) NULL)
      if (is.null(fit)) return(tibble::tibble(p_sex = NA_real_, p_cond = NA_real_, p_int = NA_real_))
      a <- tryCatch(stats::anova(fit), error = function(e) NULL)
      if (is.null(a)) return(tibble::tibble(p_sex = NA_real_, p_cond = NA_real_, p_int = NA_real_))
      tibble::tibble(p_sex = a["Sex", "Pr(>F)"], p_cond = a["Condition", "Pr(>F)"], p_int = a["Sex:Condition", "Pr(>F)"])
    }) %>% ungroup()
}
regions_min3 <- function(data, value_col, min_cell_n = 3) {
  data %>%
    filter(!is.na(.data[[value_col]])) %>%
    count(BrainRegion, Sex, Condition) %>%
    group_by(BrainRegion) %>%
    summarise(n_cells = n(), min_n = min(n), .groups = "drop") %>%
    filter(n_cells == 4, min_n >= min_cell_n) %>%
    pull(BrainRegion)
}
cliffs_d_fn <- function(a, b) {
  a <- a[!is.na(a)]; b <- b[!is.na(b)]
  if (length(a) == 0 || length(b) == 0) return(NA_real_)
  mean(outer(a, b, ">")) - mean(outer(a, b, "<"))
}

build_panel_data <- function(value_col) {
  st    <- per_region_HL(dat, value_col)
  terms <- anova_terms(dat, value_col)
  keep  <- regions_min3(dat, value_col)
  cd    <- dat %>% filter(!is.na(.data[[value_col]]), BrainRegion %in% keep) %>%
    group_by(BrainRegion) %>%
    summarise(cliffs_delta = cliffs_d_fn(.data[[value_col]][Condition == "Coord"],
                                         .data[[value_col]][Condition == "Control"]),
              .groups = "drop")
  st %>% filter(BrainRegion %in% keep) %>%
    left_join(terms, by = "BrainRegion") %>%
    left_join(cd, by = "BrainRegion") %>%
    mutate(direction = factor(ifelse(HL >= 0, "Coord higher", "Control higher"),
                              levels = c("Control higher", "Coord higher")),
           sex_sig = !is.na(p_sex) & p_sex < 0.05,
           int_sig = !is.na(p_int) & p_int < 0.05,
           tag = paste0(ifelse(sex_sig, " S", ""), ifelse(int_sig, " I", "")))
}
sd_ <- build_panel_data("Fos_density")
sr_ <- build_panel_data("FosPenk_ratio")
sg_ <- build_panel_data("FosGal_ratio")
sa_ <- build_panel_data("All3_PenkGal")

# ---- supplementary table: all 4 metrics, same schema ----
supp_all <- bind_rows(
  sd_ %>% mutate(Metric = "Fos density"),
  sr_ %>% mutate(Metric = "(Fos+Penk)/Penk"),
  sg_ %>% mutate(Metric = "(Fos+Gal)/Gal"),
  sa_ %>% mutate(Metric = "All3/(Penk+Gal)")
) %>%
  transmute(Metric, BrainRegion, n_Coord, n_Control,
            mean_Coord = round(mean_Coord, 3), mean_Control = round(mean_Control, 3),
            HL = round(HL, 3), HL_lo = round(HL_lo, 3), HL_hi = round(HL_hi, 3),
            cliffs_delta = round(cliffs_delta, 3),
            p_condition_MW = signif(p_MW, 3), q_condition_BH = signif(q_MW, 3),
            p_sex_ANOVA = signif(p_sex, 3), p_interaction_ANOVA = signif(p_int, 3),
            star_condition = star,
            sex_sig = ifelse(sex_sig, "Y", ""), interaction_sig = ifelse(int_sig, "Y", "")) %>%
  arrange(Metric, desc(abs(HL)))
write.csv(supp_all, file.path(out_dir, "Figure5_region_stats_allmetrics.csv"), row.names = FALSE)
cat("\n=== Per-region stats, all 4 metrics (BH-corrected within metric) ===\n")
print(supp_all %>% filter(q_condition_BH < 0.10 | sex_sig == "Y" | interaction_sig == "Y"), n = 50)

# Panel builders

# ---- 3a. Single region, sexes combined on x-axis, filled/hollow + sex color ----
cond_sex_panel <- function(region, metric, title, ylab) {
  d <- dat %>% filter(BrainRegion == region, !is.na(.data[[metric]])) %>%
    mutate(Condition = factor(Condition, levels = c("Control", "Coord")),
           point_fill = ifelse(Condition == "Coord", SEX_PALETTE[Sex], "white"),
           point_col  = SEX_PALETTE[Sex])
  wt <- suppressWarnings(wilcox.test(d[[metric]] ~ d$Condition, conf.int = TRUE, exact = FALSE))
  ymax <- max(d[[metric]], na.rm = TRUE)
  # wilcox.test gives Control - Coord (the order of the factor levels); the sign is
  # flipped so HL is Coord - Control, as in the bar charts
  lab <- sprintf("p=%.3f\nHL=%+.3g", wt$p.value, -unname(wt$estimate))
  ggplot(d, aes(x = Condition, y = .data[[metric]])) +
    geom_boxplot(aes(group = Condition), outlier.shape = NA, width = 0.5, alpha = 0.0,
                 color = "grey30", linewidth = LW_THIN) +
    geom_jitter(aes(fill = point_fill, color = point_col),
                position = position_jitter(width = 0.08, height = 0, seed = JITTER_SEED), shape = 21, size = 2.46,
                stroke = 0.62, alpha = 0.95) +
    scale_fill_identity() + scale_color_identity() +
    geom_segment(aes(x = 1, xend = 2, y = ymax * 1.08, yend = ymax * 1.08),
                 inherit.aes = FALSE, linewidth = LW_THIN) +
    annotate("text", x = 1.5, y = ymax * 1.26, label = lab, size = 3.62, lineheight = 0.9) +
    annotate("text", x = 1.5, y = ymax * 1.47, label = "Male",   color = SEX_PALETTE["Male"],   size = 3.62, hjust = 1.1) +
    annotate("text", x = 1.5, y = ymax * 1.47, label = "Female", color = SEX_PALETTE["Female"], size = 3.62, hjust = -0.1) +
    scale_y_continuous(limits = c(0, ymax * 1.55), expand = expansion(mult = c(0.02, 0.04))) +   # room above the points for the test result and the sex key
    scale_x_discrete(labels = c(Control = "Ctl", Coord = "Coord"), expand = expansion(add = 0.7)) +
    coord_cartesian(clip = "off") +          # the sex key may run past the narrow panel's edge
    labs(title = NULL, x = NULL, y = ylab) +
    theme_nature +
    theme(axis.text.x = element_text(size = PT_TEXT_X),
          plot.margin = MARGIN_BOX)
}
p_apoa_all3 <- cond_sex_panel("APOA", "All3_PenkGal", "APOA All3/(Penk+Gal)", "All3 / (Penk+Gal)")
p_apoa_all3
p_vh_all3   <- cond_sex_panel("VH",   "All3_PenkGal", "VH All3/(Penk+Gal)",   "All3 / (Penk+Gal)")
p_vh_all3

# BST, sex-faceted (the one Sex x Condition interaction in the screen) ----
d_bst <- dat %>% filter(BrainRegion == "BST", !is.na(Fos_density)) %>%
  mutate(Condition = factor(Condition, levels = c("Control", "Coord")),
         point_fill = ifelse(Condition == "Coord", SEX_PALETTE[Sex], "white"),
         point_col  = SEX_PALETTE[Sex])
p_male_bst   <- suppressWarnings(wilcox.test(Fos_density ~ Condition, data = d_bst,
                                              subset = Sex == "Male", exact = FALSE))$p.value
p_female_bst <- suppressWarnings(wilcox.test(Fos_density ~ Condition, data = d_bst,
                                              subset = Sex == "Female", exact = FALSE))$p.value
ymax_bst <- max(d_bst$Fos_density, na.rm = TRUE)
ann_bst <- data.frame(Sex = c("Male", "Female"),
                       label = c(sprintf("p=%.3f", p_male_bst), sprintf("p=%.3f", p_female_bst)))

p_bst_d <- ggplot(d_bst, aes(x = Condition, y = Fos_density)) +
  geom_boxplot(aes(group = Condition), outlier.shape = NA, width = 0.5, alpha = 0.0,
               color = "grey30", linewidth = LW_THIN) +
  geom_jitter(aes(fill = point_fill, color = point_col),
              position = position_jitter(width = 0.08, height = 0, seed = JITTER_SEED), shape = 21, size = 2.46, stroke = 0.62, alpha = 0.95) +
  scale_fill_identity() + scale_color_identity() +
  geom_segment(aes(x = 1, xend = 2, y = ymax_bst * 1.08, yend = ymax_bst * 1.08),
               inherit.aes = FALSE, linewidth = LW_THIN) +
  geom_text(data = ann_bst, aes(x = 1.5, y = ymax_bst * 1.18, label = label),
            inherit.aes = FALSE, size = 3.62) +
  facet_wrap(~ Sex) +
  scale_x_discrete(labels = c(Control = "Ctl", Coord = "Coord")) +
  scale_y_continuous(limits = c(0, ymax_bst * 1.3), expand = expansion(mult = c(0.02, 0.05))) +
  labs(title = NULL, x = NULL,
       y = expression("Fos density (cells/mm"^2*")")) +
  theme_nature +
  theme(axis.text.x = element_text(size = PT_TEXT_X), plot.title = element_text(size = PT_TITLE),
        strip.text = element_text(size = PT_TEXT_X),
        plot.margin = MARGIN_BOX)

p_bst_d
# Regional screen diverging bars

add_tiers <- function(d) d %>% mutate(
  star_bh  = case_when(q_MW < 0.001 ~ "***", q_MW < 0.01 ~ "**", q_MW < 0.05 ~ "*", TRUE ~ ""),
  star_unc = case_when(p_MW < 0.001 ~ "***", p_MW < 0.01 ~ "**", p_MW < 0.05 ~ "*", p_MW < 0.10 ~ ".", TRUE ~ ""))
sd2 <- add_tiers(sd_); sr2 <- add_tiers(sr_); sg2 <- add_tiers(sg_)

min_bar <- function(d) {
  rng <- range(d$HL, na.rm = TRUE)
  eps <- diff(rng) * 0.015
  d %>% mutate(HL_draw = ifelse(abs(HL) < eps, eps * sign(ifelse(HL == 0, 1, HL)), HL))
}
sd2 <- min_bar(sd2); sr2 <- min_bar(sr2); sg2 <- min_bar(sg2)

make_bar4 <- function(d, title, xlab, show_legend = FALSE) {
  d <- d %>% mutate(BrainRegion = fct_reorder(BrainRegion, HL))
  rng <- range(d$HL, na.rm = TRUE)
  pad <- diff(rng) * 0.3
  d <- d %>% mutate(lbl_bh = star_bh,
                     lbl_unc = ifelse(star_unc == "", "", paste0("(", star_unc, ")")),
                     lbl_tag = tag)
  ggplot(d, aes(y = BrainRegion, x = HL_draw, fill = direction)) +
    geom_col(width = 0.65, color = DIR_PAL["Coord higher"], linewidth = LW_AXIS) +
    geom_text(aes(label = lbl_bh,  x = HL_draw + sign(ifelse(HL_draw==0,1,HL_draw)) * (pad*0.25 + 0.01*max(abs(rng)))),
              size = 3.80, color = "black", hjust = 0, fontface = "bold") +
    geom_text(aes(label = lbl_unc, x = HL_draw + sign(ifelse(HL_draw==0,1,HL_draw)) * (pad*1.00 + 0.01*max(abs(rng)))),
              size = 3.45, color = "grey50", hjust = 0, fontface = "plain") +
    geom_text(aes(label = lbl_tag, x = HL_draw + sign(ifelse(HL_draw==0,1,HL_draw)) * (pad*2.35 + 0.01*max(abs(rng)))),
              size = 3.80, color = "black", hjust = 0, fontface = "bold") +
    geom_vline(xintercept = 0, color = "grey40", linewidth = LW_AXIS) +
    scale_fill_manual(values = DIR_PAL, name = NULL, limits = names(DIR_PAL), drop = FALSE) +
    scale_x_continuous(limits = c(rng[1] - pad, rng[2] + pad*3.0)) +
    labs(title = NULL, x = xlab, y = NULL) +
    theme_nature +
    theme(legend.position = if (show_legend) "inside" else "none",
          legend.position.inside = c(0.78, 0.22),
          legend.justification.inside = c(0.5, 0.5),
          legend.background = element_rect(fill = alpha("white", 0.75), color = NA),
          legend.margin = margin(t = 1, r = 2, b = 1, l = 2),
          axis.text.y = element_text(size = PT_TEXT))
}
# HL shift = Hodges-Lehmann shift, Coord - Control (the axis titles are short: five panel groups share the row)
p_bar_d <- make_bar4(sd2, "Fos density", expression("HL shift (cells/mm"^2*")"), show_legend = TRUE)
p_bar_d
p_bar_r <- make_bar4(sr2, "(Fos+Penk) / Penk", "HL shift (ratio)")
p_bar_r
p_bar_g <- make_bar4(sg2, "(Fos+Gal) / Gal", "HL shift (ratio)")
p_bar_g
# CarePC1 vs Fos activity, pair-averaged (mates are not independent, so
#    each breeding pair contributes one point per region
CARE_DENSITY_REGIONS <- c("STR", "NACC")
CARE_RECRUIT_REGIONS <- c("BST", "DPVP")

# one row per sample: the breeding pair it belongs to and that pair's CarePC1
care <- read.csv(PAIRS, stringsAsFactors = FALSE) %>%
  pivot_longer(c(Sample1, Sample2), values_to = "SampleID_match") %>%
  filter(SampleID_match != "") %>%
  transmute(SampleID_match, pair = Pair, CarePC1)

mer <- dat %>%
  mutate(SampleID_match = gsub("^syn(\\d+)$", "syn_\\1", SampleID)) %>%
  inner_join(care, by = "SampleID_match")

pair_metric <- function(col) {
  mer %>%
    filter(!is.na(.data[[col]])) %>%
    group_by(BrainRegion, pair) %>%
    summarise(value   = mean(.data[[col]]),
              CarePC1 = mean(CarePC1), .groups = "drop")
}
pair_dens    <- pair_metric("Fos_density")
pair_recruit <- pair_metric("FosPenk_ratio")

pair_meta <- mer %>%
  group_by(pair) %>%
  summarise(Condition = first(Condition),
            Sex_label = ifelse(n_distinct(Sex) == 2, "Mixed", first(Sex)),
            .groups = "drop")
pair_dens    <- pair_dens    %>% left_join(pair_meta, by = "pair")
pair_recruit <- pair_recruit %>% left_join(pair_meta, by = "pair")

spearman_safe <- function(x, y) {
  ok <- complete.cases(x, y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 5 || sd(x) == 0 || sd(y) == 0)
    return(list(n = length(x), rho = NA_real_, p = NA_real_))
  ct <- tryCatch(suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)),
                 error = function(e) NULL)
  if (is.null(ct)) return(list(n = length(x), rho = NA_real_, p = NA_real_))
  list(n = length(x), rho = unname(ct$estimate), p = ct$p.value)
}
region_stats <- function(pd) {
  pd %>% group_by(BrainRegion) %>%
    group_modify(~ { s <- spearman_safe(.x$CarePC1, .x$value)
                     tibble::tibble(n = s$n, rho = s$rho, p = s$p) }) %>%
    ungroup() %>% filter(!is.na(rho)) %>% arrange(p)
}
care_stats_dens    <- region_stats(pair_dens)
care_stats_recruit <- region_stats(pair_recruit)
care_stats <- bind_rows(mutate(care_stats_dens,    metric = "Fos_density"),
                        mutate(care_stats_recruit, metric = "FosPenk_recruitment"))
write.csv(care_stats, file.path(out_dir, "Figure5_CarePC1_region_stats.csv"), row.names = FALSE)
cat("\n=== CarePC1 vs Fos density / (Fos+Penk)/Penk, pair-level Spearman ===\n")
print(care_stats, n = 100)

care_scatter2 <- function(pd, st, regions, ylab) {
  d  <- pd %>% filter(BrainRegion %in% regions) %>%
    mutate(BrainRegion = factor(BrainRegion, levels = regions),
           Condition = factor(Condition, levels = c("Control", "Coord")),
           point_fill = ifelse(Condition == "Coord", SEX_PALETTE[Sex_label], "white"),
           point_col  = SEX_PALETTE[Sex_label])
  sl <- st %>% filter(BrainRegion %in% regions) %>%
    mutate(BrainRegion = factor(BrainRegion, levels = regions),
           lab = sprintf("rho = %+.2f\np %s\n%d pairs", rho,
                         ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p)), n))
  ggplot(d, aes(x = CarePC1, y = value)) +
    geom_smooth(method = "lm", se = FALSE, color = "grey40",
                linetype = "dashed", linewidth = LW_AXIS) +
    geom_point(aes(fill = point_fill, color = point_col),
               shape = 21, size = 2.22, stroke = 0.55, alpha = 0.95) +
    scale_fill_identity() + scale_color_identity() +
    facet_wrap(~ BrainRegion, scales = "free_y", nrow = 1,
               labeller = as_labeller(c(STR = "STR", NACC = "NAcc", BST = "BST",
                                        DPVP = "DPVP"), default = identity)) +
    geom_text(data = sl, aes(x = -Inf, y = Inf, label = lab),
              hjust = -0.08, vjust = 1.15, size = 3.45, lineheight = 0.9,
              inherit.aes = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +
    labs(x = "CarePC1", y = ylab) +
    theme_nature +
    theme(legend.position = "none", strip.text = element_text(size = PT_TITLE))
}
pcap_dens    <- care_scatter2(pair_dens,    care_stats_dens,    CARE_DENSITY_REGIONS,
                              expression("Fos density (cells/mm"^2*")"))
pcap_recruit <- care_scatter2(pair_recruit, care_stats_recruit, CARE_RECRUIT_REGIONS,
                              "(Fos+Penk) / Penk")

ROW1_H <- 3.278       # plot-area height of the row with B to F (the bar charts need the height for their region labels)
ROW2_H <- 1.60        # plot-area height of the row with G and H
FIG_W   <- 14.0       # width the figure is saved at: every row is squeezed into the width of the micrograph strip

p_care <- (pcap_dens | pcap_recruit)   # native vector patchwork composition, no raster round-trip


# "A" is part of the micrograph image. The rest are added here, one per panel
# group, at the top left of the group's first plot.
LETTER <- c(triple = "B", fos_density = "C", penk = "D", gal = "E", bst = "F",
            care_density = "G", care_penk = "H", timeline = "I", ntx_boxes = "J", ntx_care = "K")
lettered <- function(p, key) p + labs(tag = LETTER[[key]]) +
  theme(plot.tag = element_text(face = "bold", size = 20, color = "black"),
        plot.tag.position = "topleft")
p_apoa_all3  <- lettered(p_apoa_all3,  "triple")
p_bar_d      <- lettered(p_bar_d,      "fos_density")
p_bar_r      <- lettered(p_bar_r,      "penk")
p_bar_g      <- lettered(p_bar_g,      "gal")
p_bst_d      <- lettered(p_bst_d,      "bst")
pcap_dens    <- lettered(pcap_dens,    "care_density")
pcap_recruit <- lettered(pcap_recruit, "care_penk")
p_care       <- (pcap_dens | pcap_recruit)

row1 <- (p_apoa_all3 | p_vh_all3 | p_bar_d | p_bar_r | p_bar_g | p_bst_d) +
  plot_layout(widths = c(1.7, 1.7, 2.95, 2.95, 2.95, 3.3))
row1
row2 <- p_care

# naltrexone and parental care

# Input: Supplementary/derived_tables/Pharmacology_perTank.csv (one row per tank, treatment and
#        visitor, with the tank's saline-week CarePC1 and care level) and
#        Supplementary/derived_tables/CarePC1_rule.csv (the High/Low cutoff). The first is written
#        in section 0 of this script, the second by Figure2_11162026.R.
# The six panels are also saved on their own in Figures/Figure5_panels/.
panel_dir <- file.path(fig_dir, "Figure5_panels")
dir.create(panel_dir, showWarnings = FALSE, recursive = TRUE)

cutoff  <- read.csv(file.path(repo_dir, "Supplementary", "derived_tables", "CarePC1_rule.csv"))$Threshold[1]
dat_fam <- read.csv(file.path(repo_dir, "Supplementary", "derived_tables", "Pharmacology_perTank.csv"), stringsAsFactors = FALSE) |>
  rename(Sex = Visitor, Treatment_label = Treatment) |>
  mutate(Treatment_label = factor(Treatment_label, levels = c("Saline", "Naltrexone")),
         CareLevel       = factor(CareLevel,       levels = c("High", "Low")))

# ---- Figure-1 palettes (SEX_CARE_PALETTE / CARE_PALETTE) for the family row ----
SEX_CARE  <- c("Male.High" = "#0077BB", "Male.Low" = "#7FBBDD",
               "Family.High" = "#005F52", "Family.Low" = "#B2E4DB")
sex_dark  <- c(Male = "#0077BB", Family = "#005F52")   # sex hue (naltrexone / outline)
THEME_ROW <- theme_classic(base_size = 9) +
  theme(panel.border = element_blank(), panel.grid = element_blank(),
        axis.line  = element_line(color = "black", linewidth = LW_AXIS),
        axis.ticks = element_line(color = "black", linewidth = LW_AXIS),
        axis.title = element_text(size = PT_TITLE),
        axis.text  = element_text(size = PT_TEXT),
        axis.text.x = element_text(size = PT_TEXT),
        legend.position = "none")
NTX_TEXT <- PT_TEXT / ggplot2::.pt             # in-plot text of the naltrexone panels (mm)


#stats

cat("Saline-week CarePC1 and care level per tank:\n")

print(unique(dat_fam[, c("Tank", "CarePC1", "CareLevel")]), row.names = FALSE)
#Tank    CarePC1 CareLevel
#A9 -0.5056289       Low
#B15  2.0100082      High
#B7 -1.5363606       Low
#B9 -0.2893995       Low
#D10 -1.5363606       Low
#D15 -0.3139595       Low
#D16 -0.2784840       Low
#D8  2.2921238      High
#E5 -1.3873126       Low
#E6 -1.5363606       Low
#E8 -1.5363606       Low
#F2  0.3921736      High
#F7 -0.8702586       Low


for (sex_value in c("Male", "Family")) for (care_value in c("High", "Low")) {
  wide <- dat_fam |> filter(Sex == sex_value, CareLevel == care_value) |>
    pivot_wider(id_cols = Tank, names_from = Treatment_label, values_from = MinutesPerDay)
  tt <- t.test(wide$Naltrexone, wide$Saline, paired = TRUE)
  cat(sprintf("%-6s %-4s care: n=%d tanks, saline %.1f vs naltrexone %.1f min/day, paired t=%+.2f, p=%.4f\n",
              ifelse(sex_value == "Male", "Father", "Coord."), care_value, nrow(wide),
              mean(wide$Saline), mean(wide$Naltrexone), tt$statistic, tt$p.value))
}

#Father High care: n=3 tanks, saline 34.0 vs naltrexone 13.3 min/day, paired t=-5.10, p=0.0364
#Father Low  care: n=10 tanks, saline 14.9 vs naltrexone 15.7 min/day, paired t=+0.18, p=0.8598
#Coord. High care: n=3 tanks, saline 114.6 vs naltrexone 32.7 min/day, paired t=-3.40, p=0.0767
#Coord. Low  care: n=10 tanks, saline 6.7 vs naltrexone 10.3 min/day, paired t=+1.59, p=0.1463

for (sex_value in c("Male", "Family")) {
  co <- summary(lmerTest::lmer(MinutesPerDay ~ CarePC1 * Treatment_label + (1 | Tank),
                               data = filter(dat_fam, Sex == sex_value)))$coefficients
  i  <- grep("CarePC1:Treatment", rownames(co))
  cat(sprintf("%-6s CarePC1 x treatment: slope under saline %+.2f, change under naltrexone %+.2f, t(%.0f)=%+.2f, p=%.4f\n",
              ifelse(sex_value == "Male", "Father", "Coord."), co["CarePC1", "Estimate"],
              co[i, "Estimate"], co[i, "df"], co[i, "t value"], co[i, "Pr(>|t|)"]))
}

#Father CarePC1 x treatment: slope under saline +9.25, change under naltrexone -9.14, t(11)=-3.68, p=0.0036
#Coord. CarePC1 x treatment: slope under saline +38.06, change under naltrexone -25.51, t(11)=-4.33, p=0.0012

# Boxplot panel: Naltrexone solid / Saline hollow. show_y: draw the y-axis title
# (only the first of the four panels has it).
make_box <- function(sex_value, care_value, sexkey, ymax, show_y = FALSE) {
  d <- dat_fam |> filter(Sex == sex_value, CareLevel == care_value)
  idcol <- SEX_CARE[[paste0(sexkey, ".", care_value)]]
  outline <- sex_dark[[sexkey]]
  wide <- d |> pivot_wider(id_cols = Tank, names_from = Treatment_label,
                           values_from = MinutesPerDay) |> drop_na(Saline, Naltrexone)
  p_label <- if (nrow(wide) >= 2 && sd(wide$Naltrexone - wide$Saline) > 0) {
    pv <- t.test(wide$Naltrexone, wide$Saline, paired = TRUE)$p.value
    if (pv < 0.01) "p < 0.01" else if (pv < 0.05) sprintf("p = %.3f", pv) else "ns"
  } else "n.s."
  by_ <- ymax * 0.85; bt_ <- ymax * 0.02; ty_ <- ymax * 0.93
  ggplot(d, aes(Treatment_label, MinutesPerDay, fill = Treatment_label)) +
    geom_boxplot(outlier.shape = NA, color = outline, linewidth = 0.7, width = 0.6, alpha = 0.95) +
    geom_line(aes(group = Tank), color = "#9E9E9E", alpha = 0.6, linewidth = 0.4) +
    geom_jitter(position = position_jitter(width = 0.07, height = 0, seed = JITTER_SEED), size = 1.9, shape = 21,
                fill = idcol, color = outline, stroke = 0.4, alpha = 0.85) +
    stat_summary(fun = mean, geom = "point", shape = 23, fill = "white",
                 color = "black", size = 3, stroke = 0.9) +
    scale_fill_manual(values = c("Saline" = "white", "Naltrexone" = idcol)) +
    annotate("segment", x = 1, xend = 2, y = by_, yend = by_, linewidth = 0.5) +
    annotate("segment", x = 1, xend = 1, y = by_ - bt_, yend = by_, linewidth = 0.5) +
    annotate("segment", x = 2, xend = 2, y = by_ - bt_, yend = by_, linewidth = 0.5) +
    annotate("text", x = 1.5, y = ty_, label = p_label, size = NTX_TEXT, vjust = 0.1) +   # vjust lifts it clear of the bracket
    scale_x_discrete(labels = c(Saline = "Veh", Naltrexone = "NTX")) +
    ylim(0, ymax) +
    labs(x = NULL, y = if (show_y) "Minutes / day" else NULL) +
    THEME_ROW
}

# lmerTest CarePC1 x Treatment interaction (Satterthwaite) -> slope-difference p
interaction_p <- function(sex_value) {
  d <- dat_fam |> filter(Sex == sex_value)
  co <- summary(lmerTest::lmer(MinutesPerDay ~ CarePC1 * Treatment_label + (1 | Tank),
                               data = d))$coefficients
  co[grep("CarePC1:Treatment", rownames(co)), "Pr(>|t|)"]
}

interaction_p

# Scatter panel
make_scatter <- function(sex_value, sexkey, ymax) {
  d <- dat_fam |> filter(Sex == sex_value)
  pv <- interaction_p(sex_value)
  p_txt <- if (pv < 0.01) "interaction p < 0.01" else sprintf("interaction p = %.3f", pv)
  col <- sex_dark[[sexkey]]
  sal <- dplyr::filter(d, Treatment_label == "Saline")
  nal <- dplyr::filter(d, Treatment_label == "Naltrexone")
  lx0 <- min(d$CarePC1); lx1 <- lx0 + 0.45; ltx <- lx0 + 0.6
  # The key sits in head-room above the data: the p value, then both treatments side by side
  y_p   <- ymax * 1.36                             # p value
  y_key <- ymax * 1.08                             # Naltrexone and Saline keys
  dx    <- 1.75                                    # Saline key, this far right of the Naltrexone key
  y_top <- ymax * 1.38
  ggplot(d, aes(CarePC1, MinutesPerDay)) +
    geom_smooth(data = sal, method = "lm", se = TRUE,
                color = col, fill = col, alpha = 0.10, linewidth = 0.7) +
    geom_smooth(data = nal, method = "lm", se = TRUE,
                color = col, fill = col, alpha = 0.14, linewidth = 1.6) +
    geom_point(aes(fill = Treatment_label), shape = 21, color = col,
               size = 1.9, stroke = 0.5, alpha = 0.9) +
    geom_vline(xintercept = cutoff, linetype = "dotted", color = "grey60") +
    scale_fill_manual(values = c("Saline" = "white", "Naltrexone" = col)) +
    annotate("text", x = -Inf, y = y_p, hjust = 0, vjust = 1,
             label = p_txt, size = NTX_TEXT) +
    annotate("segment", x = lx0, xend = lx1, y = y_key, yend = y_key,
             color = col, linewidth = 1.4) +
    annotate("point", x = (lx0 + lx1) / 2, y = y_key, shape = 21,
             fill = col, color = col, size = 1.7) +
    annotate("text", x = ltx, y = y_key, hjust = 0, vjust = 0.5,
             label = "Naltrexone", color = col, size = NTX_TEXT) +
    annotate("segment", x = lx0 + dx, xend = lx1 + dx, y = y_key, yend = y_key,
             color = col, linewidth = 0.7) +
    annotate("point", x = (lx0 + lx1) / 2 + dx, y = y_key, shape = 21,
             fill = "white", color = col, size = 1.7) +
    annotate("text", x = ltx + dx, y = y_key, hjust = 0, vjust = 0.5,
             label = "Saline", color = col, size = NTX_TEXT) +
    coord_cartesian(ylim = c(0, y_top)) +   # crop view without deleting ribbon
    labs(x = "Baseline CarePC1", y = NULL) +
    THEME_ROW
}

ymax_father <- max(dat_fam$MinutesPerDay[dat_fam$Sex == "Male"],   na.rm = TRUE) * 1.15
ymax_family <- max(dat_fam$MinutesPerDay[dat_fam$Sex == "Family"], na.rm = TRUE) * 1.15

# The six panels, in the order father High, coordinated High, father Low,
# coordinated Low, father vs CarePC1, coordinated vs CarePC1
pharm <- list(make_box("Male",   "High", "Male",   ymax_father, show_y = TRUE),
              make_box("Family", "High", "Family", ymax_family),
              make_box("Male",   "Low",  "Male",   ymax_father),
              make_box("Family", "Low",  "Family", ymax_family),
              make_scatter("Male",   "Male",   ymax_father),
              make_scatter("Family", "Family", ymax_family))

# A copy of the six panels on their own
fig_fam <- wrap_plots(pharm, nrow = 1, widths = c(1, 1, 1, 1, 1.25, 1.25))
ggsave(file.path(panel_dir, "Pharmacology_panels.pdf"), fig_fam, width = 11.55, height = 3.2)
ggsave(file.path(panel_dir, "Pharmacology_panels.png"), fig_fam, width = 11.55, height = 3.2,
       dpi = 300, bg = "white")
cat("Saved Pharmacology_panels.{pdf,png}\n")

TIMELINE_PNG   <- file.path(repo_dir, "Figures", "artwork", "Figure5_naltrexone_timeline_stacked.png")   # panel I of Figure5_11122026.ai, without its letter
TIMELINE_W     <- 4.3                                  # width of the timeline column, inches
TIMELINE_IMG_W <- 2.981                                # printed width of the drawing, inches (2.565 in the .ai, enlarged
                                                       #   with J and K so its foot stays level with the foot of K)
TIMELINE_X     <- 0.223                                # its left edge, inches into its column (0.49 in from the edge of the figure)
TIMELINE_PAD   <- c(top = 0.115, bottom = 0.195)       # space above and below the drawing, inches; the bottom value
                                                       #   puts its foot level with the foot of panel K
# Drawn in place of an image file that is missing
placeholder_panel <- function(note) {
  ggplot() + theme_void() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
             fill = NA, color = "grey70", linetype = "dashed", linewidth = LW_AXIS) +
    annotate("text", x = 0.5, y = 0.5, color = "grey55", fontface = "italic",
             size = NTX_TEXT, label = note) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE)
}
if (file.exists(TIMELINE_PNG)) {
  img         <- png::readPNG(TIMELINE_PNG)
  timeline_ar <- ncol(img) / nrow(img)                  # width / height
  p_timeline  <- ggplot() + theme_void() +
    annotation_custom(grid::rasterGrob(img, interpolate = TRUE,
                                       x = grid::unit(TIMELINE_X, "in"), hjust = 0,
                                       y = grid::unit(TIMELINE_PAD[["bottom"]], "in"), vjust = 0,
                                       width = grid::unit(TIMELINE_IMG_W, "in")))
} else {
  message("NOTE, Figure 5: timeline image not found, drawing a placeholder: ", TIMELINE_PNG)
  timeline_ar <- 0.694                                  # shape of the drawing, used for the placeholder
  p_timeline  <- placeholder_panel("[ naltrexone timeline:\nimage not found ]")
}
ROW3_H        <- TIMELINE_PAD[["top"]] + TIMELINE_IMG_W / timeline_ar + TIMELINE_PAD[["bottom"]]   # whole naltrexone block
NTX_BOX_H     <- 2.8                                   # relative heights: row of four box plots
NTX_SCATTER_H <- 3.0                                   #   and row of two CarePC1 panels
# Box plots ordered so the two father panels (blue) sit together, then the two
# coordinated-visit panels (green): father High, father Low, coordinated High,
# coordinated Low. (pharm is built in the order father High, coordinated High,
# father Low, coordinated Low.)
p_timeline <- lettered(p_timeline, "timeline")
pharm[[1]] <- lettered(pharm[[1]], "ntx_boxes")
pharm[[5]] <- lettered(pharm[[5]], "ntx_care")
ntx_panels <- wrap_plots(pharm[c(1, 3, 2, 4)], nrow = 1, widths = c(1.2, 1, 1, 1)) /
              wrap_plots(pharm[5:6], nrow = 1) +
              plot_layout(heights = c(NTX_BOX_H, NTX_SCATTER_H))
row3   <- (wrap_elements(full = p_timeline) | wrap_elements(full = ntx_panels)) +
  plot_layout(widths = c(TIMELINE_W, FIG_W - TIMELINE_W))

# ---- Top row: micrographs (a finished image, placed at its own aspect ratio) ---
# Figures/artwork/Figure5_micrographs.png is exported from the Illustrator file
# Figure5_11122026.ai (its panel A: coordinating and non-coordinating
# RNAscope images, each with an enlarged inset, and the color key).
MICRO_PNG <- file.path(repo_dir, "Figures", "artwork", "Figure5_micrographs.png")
if (file.exists(MICRO_PNG)) {
  micro_img <- png::readPNG(MICRO_PNG)
  micro_ar  <- ncol(micro_img) / nrow(micro_img)       # width / height
  p_micro   <- ggplot() + theme_void() + annotation_custom(grid::rasterGrob(micro_img, interpolate = TRUE))
} else {
  message("NOTE, Figure 5: micrograph image not found, drawing a placeholder: ", MICRO_PNG)
  micro_ar <- 3.72                                      # shape of the micrograph panel, used for the placeholder
  p_micro  <- placeholder_panel("[ micrographs: image not found ]")
}
# The micrograph strip spans the whole figure width, edge to edge (it is not
# lined up with the plot areas below, which start inside their axis labels).
row0   <- wrap_elements(full = p_micro)
ROW0_H <- FIG_W / micro_ar

combined <- row0 / row1 / row2 / wrap_elements(full = row3) +
  plot_layout(heights = grid::unit(c(ROW0_H, ROW1_H, ROW2_H, ROW3_H), "in"))

# Figure height: the four rows plus the axis text, titles and margins between them.
# Measured on a blank PDF device, so the height does not depend on whatever
# graphics device happens to be open (RStudio, Quarto, none).
grDevices::pdf(NULL, width = FIG_W, height = 30) 
FIG_H <- sum(grid::convertHeight(patchwork::patchworkGrob(combined)$heights, "in", valueOnly = TRUE))
invisible(grDevices::dev.off())

ggsave(file.path(fig_dir, "Figure5_11112026.pdf"), combined, width = FIG_W, height = FIG_H, units = "in")
ggsave(file.path(fig_dir, "Figure5_11112026.png"), combined, width = FIG_W, height = FIG_H, units = "in",
       dpi = 300, bg = "white")

cat(sprintf("\nFinal composite: %.2f in wide x %.2f in tall (%.0f x %.0f px @ 300 dpi)\n",
            FIG_W, FIG_H, FIG_W * 300, FIG_H * 300))
cat("Saved Figure5_11112026.{pdf,png} to", fig_dir, "and the two statistics tables to", out_dir, "\n")
