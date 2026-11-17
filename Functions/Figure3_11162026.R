# Figure3
# Inputs
#   Data/begging_frequency/combinedbegginpeakfreq_v2.xlsx  sheet "refined"
#       one row per behaviour bout: its peak frequency
#   Data/begging_frequency/averagedbeggingpeaks.xlsx       sheet "Sheet1"
#       one row per individual and behaviour: its mean peak frequency
#   Data/tadbot/TadBot_posvnegctrl_v2.xlsx                 sheet "raw_data"
#       one row per nest visit: tank, treatment, trial day, visitor, duration
#   Supplementary/derived_tables/CarePC1_rule.csv   the CarePC1 rule
# Output: Figures/Figure3_11112026.png + .pdf
#   Supplementary/derived_tables/TadBot_dailyVisits.csv  visits per tank, treatment and trial day (panels D, E)
#   Supplementary/derived_tables/TadBot_first3days.csv   per tank and treatment, first three days (printed tests)

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr)
  library(ggplot2); library(patchwork)   # panels D and E also need lme4 (called as lme4::glmer.nb)
})

# ---------- Theme + palettes ----------
theme_nature <- theme_classic(base_size = 12, base_family = "sans") +
  theme(panel.grid        = element_blank(),
        axis.line         = element_line(linewidth = 0.4),
        axis.ticks        = element_line(linewidth = 0.4),
        axis.ticks.length = unit(2, "pt"),
        plot.title        = element_blank(),
        plot.title.position = "plot",
        legend.background = element_blank(),
        legend.key.size   = unit(0.9, "lines"),
        legend.text       = element_text(size = 10),
        axis.title        = element_text(size = 12),
        axis.text         = element_text(size = 10))
theme_set(theme_nature)

FREQ_PALETTE <- c("Begging" = "#882255", "Wiggling" = "#DDAACC", "Robotad" = "#CC3311")

darken_hex <- function(hex, amount = 0.35) {
  m  <- col2rgb(hex)
  h  <- rgb2hsv(m[1, ], m[2, ], m[3, ])
  v2 <- pmax(h["v", ] * (1 - amount), 0)
  out <- hsv(h["h", ], h["s", ], v2)
  names(out) <- names(hex)
  out
}

C_PALETTE      <- c("Begging" = "#882255", "Wiggling" = "#DDAACC", "TadBot" = "#CC3311")
C_PALETTE_DARK <- darken_hex(C_PALETTE)

repo_dir <- here::here()
fig_dir  <- file.path(repo_dir, "Figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

fp_freq   <- file.path(repo_dir, "Data", "begging_frequency", "combinedbegginpeakfreq_v2.xlsx")
fp_avg    <- file.path(repo_dir, "Data", "begging_frequency", "averagedbeggingpeaks.xlsx")
fp_tadbot <- file.path(repo_dir, "Supplementary", "derived_tables", "TadBot_first3days.csv")
fp_daily  <- file.path(repo_dir, "Supplementary", "derived_tables", "TadBot_dailyVisits.csv")

# TadBot tables
# One row per visit. Each tank was recorded with the TadBot moving at 15 Hz ("TadBot") and not moving ("NegCtrl"). A recorded day with no visits has one row with no visitor; a day with no row at all was not recorded.
local({
  raw_dir <- file.path(repo_dir, "Data")
  source(file.path(repo_dir, "Functions", "shared_functions.R"), local = TRUE)
  rule <- load_care_rule(file.path(repo_dir, "Supplementary", "derived_tables", "CarePC1_rule.csv"))

  # min in duration cell (read by readxl as a clock time)
  clock_minutes <- function(x) as.numeric(format(x, "%H", tz = "UTC")) * 60 +
                               as.numeric(format(x, "%M", tz = "UTC")) +
                               as.numeric(format(x, "%S", tz = "UTC")) / 60

  tadbot <- suppressWarnings(as.data.frame(
    read_excel(file.path(raw_dir, "tadbot", "TadBot_posvnegctrl_v2.xlsx"), sheet = "raw_data")))
  paired_tanks <- names(which(table(unique(tadbot[, c("Tank", "Treatment")])$Tank) == 2))   # recorded under both treatments

  # 1a. Visits on each recorded trial day, days 1 to 14 (Figure 3D and E)
  tb_days <- unique(tadbot[tadbot$Tank %in% paired_tanks & tadbot$day_trial %in% 1:14, c("Tank", "Treatment", "day_trial")])
  count_visits <- function(who) mapply(function(tank, treatment, day)
    sum(tadbot$Tank == tank & tadbot$Treatment == treatment & tadbot$day_trial == day & tadbot$sex %in% who),
    tb_days$Tank, tb_days$Treatment, tb_days$day_trial)
  tb_daily <- data.frame(Tank = tb_days$Tank, Treatment = tb_days$Treatment, TrialDay = tb_days$day_trial,
                         AllVisits     = count_visits(c("Male", "Female", "Family")),
                         Visits_Male   = count_visits("Male"),
                         Visits_Female = count_visits("Female"),
                         Visits_Family = count_visits("Family"))
  tb_daily <- tb_daily[order(tb_daily$Tank, tb_daily$Treatment, tb_daily$TrialDay), ]
  write.csv(tb_daily, fp_daily, row.names = FALSE)

  # 1b. First three trial days. For each tank and treatment, visits on trial days 1 to 3 are counted and divided by the number of those days that were recorded (C7 under TadBot has only 2 days).
  tadbot <- tadbot[tadbot$day_trial %in% 1:3, ]
  tadbot$minutes <- clock_minutes(tadbot$timedifference)
  groups <- unique(tadbot[, c("Tank", "Treatment")])
  tb_rows <- lapply(seq_len(nrow(groups)), function(i) {
    d      <- tadbot[tadbot$Tank == groups$Tank[i] & tadbot$Treatment == groups$Treatment[i], ]
    n_days <- length(unique(d$day_trial))        # a recorded day with no visits has a row with no visitor
    d      <- d[!is.na(d$sex), ]
    family <- d[d$sex == "Family", ]
    data.frame(Tank = groups$Tank[i], Treatment = groups$Treatment[i], DaysRecorded = n_days,
               AllVisits_perDay      = nrow(d) / n_days,
               Visits_Family_perDay  = nrow(family) / n_days,
               Minutes_Family_perDay = sum(family$minutes) / n_days)
  })
  tb <- do.call(rbind, tb_rows)
  tb$CarePC1 <- rule$care_pc1(tb$Visits_Family_perDay, tb$Minutes_Family_perDay)
  tb <- tb[order(tb$Tank, tb$Treatment), ]
  both <- names(which(table(tb$Tank) == 2))      # tanks recorded under both treatments
  tb <- tb[tb$Tank %in% both, ]
  write.csv(tb, fp_tadbot, row.names = FALSE)
  cat(sprintf("TadBot: %d tanks recorded under both treatments\n", length(both)))
})

#TadBot: 6 tanks recorded under both treatments                                                       

#figure 3
stopifnot(file.exists(fp_freq), file.exists(fp_avg), file.exists(fp_tadbot), file.exists(fp_daily))

tadbot   <- read.csv(fp_tadbot)
freq_dat <- read_excel(fp_freq,  sheet = "refined") |>
  filter(!is.na(Behavior), !is.na(X)) |>
  mutate(Behavior = factor(Behavior, levels = c("Begging","Wiggling","Robotad")))
freq_avg <- read_excel(fp_avg,   sheet = "Sheet1") |>
  filter(!is.na(Behavior), !is.na(X)) |>
  mutate(Behavior = factor(Behavior, levels = c("Begging","Wiggling","Robotad")))

tb_pvn <- tadbot |>
  select(Tank, Treatment, AllVisits_perDay) |>
  pivot_wider(names_from = Treatment, values_from = AllVisits_perDay) |>
  drop_na(NegCtrl, TadBot)


cat("Figure 3 -- behaviors: ", nrow(freq_dat), " (",
    paste(names(table(freq_dat$Behavior)), table(freq_dat$Behavior), collapse = ", "), ")\n", sep = "")
#Figure 3 -- behaviors: 236 (Begging 136, Wiggling 74, Robotad 26)
# helperscripts
make_blank <- function(letter, sub = "Illustrator overlay", show_title = TRUE) {
  p <- ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
             fill = NA, color = "#cccccc", linetype = "dashed", linewidth = 0.5) +
    annotate("text", x = 0.5, y = 0.55, label = sprintf("[ %s ]", letter),
             color = "#aaaaaa", fontface = "italic", size = 3.0) +
    annotate("text", x = 0.5, y = 0.45, label = sub,
             color = "#aaaaaa", size = 2.4) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme(axis.line = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), axis.title = element_blank(),
          panel.grid = element_blank(), panel.border = element_blank())
  if (show_title) p <- p + labs(title = letter)
  p
}

sigstr <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "n.s."
}

panelA <- make_blank("A. TadBot apparatus")

freq_n      <- freq_dat |> count(Behavior) |>                       # the data call the TadBot "Robotad"
  mutate(label = sprintf("%s (n=%d)", ifelse(Behavior == "Robotad", "TadBot", as.character(Behavior)), n))
freq_legend <- setNames(freq_n$label, as.character(freq_n$Behavior))

panelB <- ggplot(freq_dat, aes(x = X, fill = Behavior, color = Behavior)) +
  geom_density(alpha = 0.55, linewidth = 0.5, adjust = 0.9) +
  scale_fill_manual(values = FREQ_PALETTE, labels = freq_legend, name = NULL) +
  scale_color_manual(values = FREQ_PALETTE, guide = "none") +
  scale_x_continuous(limits = c(0, max(freq_dat$X, na.rm = TRUE) * 1.05), expand = c(0, 0)) +
  labs(x = "Peak frequency (Hz)", y = "Density",
       title = "B. Peak frequency: Begging, Wiggling, TadBot") +
  theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1))

all_freq_avg <- freq_avg |>
  filter(Behavior %in% c("Begging","Wiggling","Robotad")) |>
  mutate(Behavior = factor(dplyr::recode(Behavior, Robotad = "TadBot"),
                           levels = c("Begging","Wiggling","TadBot")))

# Paired Begging-Wiggling (connecting lines + paired t-test)
beg_wig <- freq_avg |> filter(Behavior %in% c("Begging","Wiggling")) |>
  select(ID, Behavior, X) |>
  pivot_wider(names_from = Behavior, values_from = X) |> drop_na(Begging, Wiggling)
C_bw_p <- t.test(beg_wig$Wiggling, beg_wig$Begging, paired = TRUE)$p.value
# C_bw_p is 0.01803853
beg_wig_long <- beg_wig |>
  pivot_longer(c(Begging, Wiggling), names_to = "Behavior", values_to = "X") |>
  mutate(Behavior = factor(Behavior, levels = c("Begging","Wiggling","TadBot")))


set.seed(42)
bw_ids <- unique(all_freq_avg$ID[all_freq_avg$Behavior != "TadBot"])
jit_bw <- setNames(runif(length(bw_ids), -0.08, 0.08), bw_ids)
all_freq_avg <- all_freq_avg |>
  mutate(x_jit = as.numeric(Behavior) +
           ifelse(Behavior == "TadBot", 0, jit_bw[as.character(ID)]))
beg_wig_long <- beg_wig_long |>
  mutate(x_jit = as.numeric(Behavior) + jit_bw[as.character(ID)])

# Begging vs TadBot (Mann-Whitney)
beg_rob  <- freq_avg |> filter(Behavior %in% c("Begging","Robotad")) |> select(ID, Behavior, X)
C_br_p   <- suppressWarnings(wilcox.test(X ~ Behavior, data = beg_rob))$p.value
ymax_C   <- max(all_freq_avg$X, na.rm = TRUE)
#ymax_C should be 18.45706
panelC <- ggplot(all_freq_avg, aes(Behavior, X, fill = Behavior, color = Behavior)) +
  geom_boxplot(width = 0.55, alpha = 0.55, outlier.shape = NA,
               linewidth = 0.5, fatten = 1.3) +
  scale_color_manual(values = C_PALETTE_DARK, guide = "none") +
  geom_line(data = beg_wig_long, aes(x = x_jit, y = X, group = ID),
            color = "#888888", alpha = 0.5, linewidth = 0.5) +
  geom_point(data = all_freq_avg |> filter(Behavior != "TadBot"),
             aes(x = x_jit, y = X),
             shape = 21, size = 2.0, stroke = 0.4) +          # outline = darker shade of the fill
  geom_jitter(data = all_freq_avg |> filter(Behavior == "TadBot"),
              shape = 21, size = 1.6, stroke = 0.4,
              width = 0.07, height = 0) +
  scale_fill_manual(values = c("Begging" = "#882255", "Wiggling" = "#DDAACC",
                                "TadBot"  = "#CC3311"), guide = "none") +
  # Beg-Wig bracket (lower)
  annotate("segment", x = 1, xend = 2, y = ymax_C*1.06, yend = ymax_C*1.06, linewidth = 0.5) +
  annotate("text",    x = 1.5, y = ymax_C*1.11, label = sigstr(C_bw_p), size = 5.3, fontface = "bold") +
  annotate("text",    x = 1.5, y = ymax_C*1.17,
           label = sprintf("paired t, p=%.3f", C_bw_p), size = 3.0, family = "mono") +
  # Beg-TadBot bracket (higher)
  annotate("segment", x = 1, xend = 3, y = ymax_C*1.27, yend = ymax_C*1.27, linewidth = 0.5) +
  annotate("text",    x = 2, y = ymax_C*1.32, label = sigstr(C_br_p), size = 5.3, fontface = "bold") +
  annotate("text",    x = 2, y = ymax_C*1.38,
           label = sprintf("M-W p=%.3f", C_br_p), size = 3.0, family = "mono") +
  coord_cartesian(ylim = c(NA, ymax_C * 1.45)) +
  labs(x = NULL, y = "Peak frequency (Hz)",
       title = sprintf("C. Begging vs Wiggling vs TadBot (n=%d)", nrow(beg_wig)))

# D and E  Cumulative visits over the 14 trial days, TadBot off vs on
# D  all parental visits
# E  coordinated visits only

# Test: negative-binomial mixed model on the daily counts,
#   visits ~ treatment * day + (1 | tank),  day counted from trial day 1.
#   "TadBot"       = on vs off on trial day 1
#   "Day"          = change per day with the TadBot off
#   "Day x TadBot" = whether that change differs with the TadBot on
# Dashed line in D: the trial day on which the model's visits per day are equal
# in the two conditions, i.e. where the early advantage of the on condition ends.
daily <- read.csv(fp_daily, stringsAsFactors = FALSE) |>
  mutate(Treatment = factor(Treatment, levels = c("NegCtrl", "TadBot")), d0 = TrialDay - 1)
tanks_DE <- unique(daily$Tank)
TREAT_COLORS <- c("NegCtrl" = "#888888", "TadBot" = "#CC3311")

# cumulative curve of one visit column, with its jackknife SE
cumulative_curve <- function(column) {
  curve <- function(d) d |>
    group_by(Treatment, TrialDay) |> summarise(mean_visits = mean(.data[[column]]), .groups = "drop") |>
    group_by(Treatment) |> arrange(TrialDay) |> mutate(Cumulative = cumsum(mean_visits)) |> ungroup()
  se <- bind_rows(lapply(tanks_DE, function(tank) curve(filter(daily, Tank != tank)))) |>
    group_by(Treatment, TrialDay) |>
    summarise(SE = sqrt((n() - 1) / n() * sum((Cumulative - mean(Cumulative))^2)), .groups = "drop")
  bind_rows(data.frame(Treatment = factor(c("NegCtrl", "TadBot"), levels = c("NegCtrl", "TadBot")),
                       TrialDay = 0, mean_visits = 0, Cumulative = 0, SE = 0),
            left_join(curve(daily), se, by = c("Treatment", "TrialDay")))
}
#jacknife SE is the sqrt((n - 1) / n * sum((value - mean)^2)) part

# the mixed model for one visit column: its coefficient table
visit_model <- function(column) {
  fit <- suppressMessages(suppressWarnings(
    lme4::glmer.nb(reformulate("Treatment * d0 + (1 | Tank)", response = column), data = daily)))
  summary(fit)$coefficients
}
model_label <- function(coef) sprintf("TadBot %s\nDay %s\nDay \u00d7 TadBot %s", sigstr(coef["TreatmentTadBot", 4]),
                                      sigstr(coef["d0", 4]), sigstr(coef["TreatmentTadBot:d0", 4]))
cumulative_panel <- function(curve, coef, ylab, title) {
  ggplot(curve, aes(TrialDay, Cumulative, color = Treatment, fill = Treatment)) +
    geom_ribbon(aes(ymin = Cumulative - SE, ymax = Cumulative + SE), color = NA, alpha = 0.2) +
    geom_line(linewidth = 0.9) +
    annotate("text", x = 0.3, y = Inf, label = model_label(coef), hjust = 0, vjust = 1.2,
             size = 3.3, fontface = "bold", lineheight = 1.2) +
    scale_color_manual(values = TREAT_COLORS, guide = "none") +
    scale_fill_manual(values = TREAT_COLORS, guide = "none") +
    scale_x_continuous(breaks = c(0, 3, 7, 10, 14)) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.12))) +
    labs(x = "Trial day", y = ylab, title = title)
}

D_curve <- cumulative_curve("AllVisits")
D_coef  <- visit_model("AllVisits")
D_cross <- 1 - D_coef["TreatmentTadBot", "Estimate"] / D_coef["TreatmentTadBot:d0", "Estimate"]
D_top   <- max(D_curve$Cumulative)
panelD <- cumulative_panel(D_curve, D_coef, "Cumulative visits",
                           sprintf("D. Cumulative visits, TadBot off vs on   (n = %d tanks)", length(tanks_DE))) +
  geom_vline(xintercept = D_cross, linetype = "dashed", color = "grey35", linewidth = 0.45) +
  annotate("text", x = D_cross + 0.25, y = D_top * 0.02, label = sprintf("on advantage ends\n(day %.0f)", D_cross),
           hjust = 0, vjust = 0, size = 3.0, color = "grey25", lineheight = 0.95) +
  annotate("text", x = 14.2, y = D_top * c(0.37, 0.27), label = c("15 Hz (on)", "0 Hz (off)"),
           color = TREAT_COLORS[c("TadBot", "NegCtrl")], hjust = 1, size = 3.3, fontface = "bold")

E_curve <- cumulative_curve("Visits_Family")
E_coef  <- visit_model("Visits_Family")
panelE <- cumulative_panel(E_curve, E_coef, "Cumulative coordinated visits",
                           sprintf("E. Cumulative coordinated visits   (n = %d tanks)", length(tanks_DE)))

first3_visits <- t.test(tb_pvn$TadBot, tb_pvn$NegCtrl, paired = TRUE)
first3_visits

#	Paired t-test
#data:  tb_pvn$TadBot and tb_pvn$NegCtrl
#t = 3.1538, df = 5, p-value = 0.02527
#alternative hypothesis: true mean difference is not equal to 0
#95 percent confidence interval:
#  0.4212465 4.1343091
#sample estimates:
#  mean difference 
#2.277778 

care_wide <- tadbot |>
  select(Tank, Treatment, CarePC1) |>
  pivot_wider(names_from = Treatment, values_from = CarePC1) |>
  drop_na(NegCtrl, TadBot)
first3_care <- t.test(care_wide$TadBot, care_wide$NegCtrl, paired = TRUE)
first3_care

#Paired t-test

#data:  care_wide$TadBot and care_wide$NegCtrl
#t = -0.18205, df = 5, p-value = 0.8627
#alternative hypothesis: true mean difference is not equal to 0
#95 percent confidence interval:
#  -1.091113  0.946785
#sample estimates:
#  mean difference 
#-0.0721638 



# ---- Printed statistics
day_points <- function(curve, treatment) paste(round(curve$Cumulative[curve$Treatment == treatment & curve$TrialDay %in% c(3, 7, 14)], 1), collapse = " / ")
model_line <- function(coef) sprintf("TadBot on/off on day 1 = %.2fx, p=%.3f | day (off) %+.1f%% per day, p=%.3f | day x TadBot p=%.4f",
  exp(coef["TreatmentTadBot", 1]), coef["TreatmentTadBot", 4], 100 * (exp(coef["d0", 1]) - 1), coef["d0", 4], coef["TreatmentTadBot:d0", 4])
cat(sprintf("Panel C | Begging vs wiggling: n=%d individuals, paired t-test p=%.3f | Begging vs TadBot: Mann-Whitney p=%.3f\n",
    nrow(beg_wig), C_bw_p, C_br_p))
#Panel C | Begging vs wiggling: n=11 individuals, paired t-test p=0.018 | Begging vs TadBot: Mann-Whitney p=0.122
cat(sprintf("Panel D | All visits, 14 days, NB mixed model (n=%d tanks, %d tank-days): %s | rates equal on trial day %.1f\n",
    length(tanks_DE), nrow(daily), model_line(D_coef), D_cross))
#Panel D | All visits, 14 days, NB mixed model (n=6 tanks, 156 tank-days): TadBot on/off on day 1 = 1.50x, p=0.017 | day (off) +3.4% per day, p=0.027 | day x TadBot p=0.0065 | rates equal on trial day 7.7
cat(sprintf("        | Cumulative visits (mean of tanks) on day 3 / 7 / 14: off %s; on %s\n", day_points(D_curve, "NegCtrl"), day_points(D_curve, "TadBot")))
#        | Cumulative visits (mean of tanks) on day 3 / 7 / 14: off 20 / 47.5 / 108.4; on 25.4 / 53 / 108.4
cat(sprintf("Panel E | Coordinated visits, 14 days, NB mixed model: %s\n", model_line(E_coef)))
#Panel E | Coordinated visits, 14 days, NB mixed model: TadBot on/off on day 1 = 1.45x, p=0.385 | day (off) +2.8% per day, p=0.462 | day x TadBot p=0.0752
cat(sprintf("        | Cumulative coordinated visits on day 3 / 7 / 14: off %s; on %s\n", day_points(E_curve, "NegCtrl"), day_points(E_curve, "TadBot")))
#        | Cumulative coordinated visits on day 3 / 7 / 14: off 2.5 / 5 / 11.9; on 2.4 / 4.8 / 9.8
cat(sprintf("First three trial days | all visits/day, off vs on: n=%d tanks, means %.2f vs %.2f, paired t=%+.2f, p=%.3f\n",
    nrow(tb_pvn), mean(tb_pvn$NegCtrl), mean(tb_pvn$TadBot), first3_visits$statistic, first3_visits$p.value))
#First three trial days | all visits/day, off vs on: n=6 tanks, means 6.67 vs 8.94, paired t=+3.15, p=0.025
cat(sprintf("First three trial days | CarePC1, off vs on: n=%d tanks, means %.2f vs %.2f, paired t=%+.2f, p=%.3f\n",
    nrow(care_wide), mean(care_wide$NegCtrl), mean(care_wide$TadBot), first3_care$statistic, first3_care$p.value))
#First three trial days | CarePC1, off vs on: n=6 tanks, means -0.70 vs -0.77, paired t=-0.18, p=0.863


fig3 <- (panelA | panelB | panelC | panelD | panelE) +
  plot_layout(widths = c(2.1, 1.3, 1.1, 1.25, 1.25)) &
  theme(plot.margin = margin(3, 5, 3, 5))

fig3

FIG_W <- 14.791564 * 7.0 / 6.3; FIG_H <- 4.807258
ggsave(file.path(fig_dir, "Figure3_11112026.png"), fig3,
       width = FIG_W, height = FIG_H, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "Figure3_11112026.pdf"), fig3,
       width = FIG_W, height = FIG_H)

cat("Saved Figure3_11112026 to", fig_dir, "\n")
