# FIGURE 1
# Input : Data/event_logs/*.xlsx   one spreadsheet per recording; each row
#                                      is one "In" or "Out" event with a date,
#                                      a clock time and a visitor (Male,
#                                      Female, or Family = both parents)
#         Data/recording_metadata.csv   one row per recording
# Output: Supplementary/derived_tables/VisitDurations_byType_derived_v4.csv
#         Figures/Figure1_ethogram_panels_11102026.pdf / .png
#
# Rules for reading the event logs (Part 1)
#   1. Recordings = the 26 ObservationIDs in recording_metadata.csv. A spreadsheet is used if its name starts with that ID (after an optional "done-"). If two spreadsheets match, the "done-" one is used.
#   2. Columns are found by what they contain (see find_columns).
#   3. Visitor typos are read by their first letters (Mlae -> Male, Femae -> Female, Famly -> Family).
#   4. Typing slips in dates and clock times are corrected, and every correction is printed.
#   6. A visit is kept if it lasts at most 24 h. Visits with In and Out in the same minute are kept with a duration of 0: they count as visits, but analyses duration leave them out.
#   7. The recording day comes from the spreadsheet's Day column. Where that cell is blank, it is filled in from the date, using the samedate-to-day offset as the rest of that recording.
# ============================================================================

library(readxl)
library(ggplot2)
library(patchwork)

# ---- 1. Where things live ---------------------------------------------------
repo_dir <- here::here()

visit_file <- file.path(repo_dir, "Supplementary", "derived_tables", "VisitDurations_byType_derived_v4.csv")
meta_file  <- file.path(repo_dir, "Data", "recording_metadata.csv")
fig_dir    <- file.path(repo_dir, "Figures")
out_name   <- "Figure1_ethogram_panels_11102026"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(visit_file), showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# PART 1 -- Transform raw event logs
# ============================================================================
local({
  raw_dir <- file.path(repo_dir, "Data", "event_logs")
  id_file <- meta_file

  # ---- Which spreadsheet belongs to which recording ----------------------------
  recordings <- read.csv(id_file, stringsAsFactors = FALSE)$ObservationID

  files     <- list.files(raw_dir, pattern = "\\.xlsx$")
  file_tank <- sub(" .*$", "", sub("^done-", "", files))        # "done-A7-1 3_24_22..." -> "A7-1"

  pick_file <- function(tank) {
    hits <- files[file_tank == tank]
    if (length(hits) == 0) stop("No spreadsheet found for ", tank)
    done <- hits[startsWith(hits, "done-")]
    if (length(done) > 0) done[1] else hits[1]
  }
  source_files <- vapply(recordings, pick_file, character(1))
  # ---- Read each spreadsheet and pair In with Out -------------------------------
  # The reading and pairing functions live in shared_functions.R, shared with the
  # RNAscope event logs of Figure 5.
  source(file.path(repo_dir, "Functions", "shared_functions.R"), local = TRUE)


  # ---- Build table ---------------------------------------------------------
  build_one <- function(tank) {
    data.frame(tank = tank, read_visits(source_files[[tank]], raw_dir), source_file = source_files[[tank]])
  }
  visit_table <- do.call(rbind, lapply(recordings, build_one))
  visit_table <- visit_table[, c("tank", "visitor", "day", "duration_min", "source_file")]
  rownames(visit_table) <- NULL

  write.csv(visit_table, visit_file, row.names = FALSE)

  cat(sprintf("Wrote %d visits (%d of them zero-length) from %d recordings to %s\n",
              nrow(visit_table), sum(visit_table$duration_min == 0), length(recordings), visit_file))
  print(table(visit_table$visitor))
})

#you should have 5805 visits total at this point - 101 are zero length and therefore ignored

# ============================================================================
# PART 2 -- Panel B and C
# ============================================================================
# ---- 2. Read cleaned up data -------------------------------------------------------
visits <- read.csv(visit_file)                                   # one row per visit
visits <- visits[visits$duration_min > 0, ]                      # In and Out in the same minute: no duration to compare
meta   <- read.csv(meta_file, stringsAsFactors = FALSE)          # one row per recording
eggs   <- data.frame(ObservationID = meta$ObservationID,
                     Eggs = meta$EggFeedingEvents,                # number of egg-feeding events
                     NFD  = meta$NumberofFullDays)                # number of full recording days

visits$type <- factor(ifelse(visits$visitor == "Family", "Coord", "Single"),
                      levels = c("Single", "Coord"))

# Panel C
eggs <- eggs[!is.na(eggs$NFD) & eggs$NFD > 0 & !is.na(eggs$Eggs), ]

# Panel C x axis: coordinated minutes per day
coord_minutes <- tapply(visits$duration_min[visits$type == "Coord"],
                        visits$tank[visits$type == "Coord"], sum)
eggs$coord_min_per_day <- as.numeric(coord_minutes[eggs$ObservationID]) / eggs$NFD

# ---- 3. Stats ----------------------------------------------------------
stars <- function(p) if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "n.s."

# B mean log10 duration per recording and visit type; their difference is the recording's log10 fold difference (Coord / Single, as a ratio of geometric means). Wilcoxon signed-rank across recordings.

tank_means <- tapply(log10(visits$duration_min), list(visits$tank, visits$type), mean) #transform the duration of each visit 
tank_means <- tank_means[complete.cases(tank_means), ] #only look for complete cases
tank_means
fold       <- data.frame(recording  = rownames(tank_means),
                         log10_fold = tank_means[, "Coord"] - tank_means[, "Single"]) # coordinated minus single-parent: each recording's mean log10 duration of coordinated visits, minus that of its single-parent visits
paired     <- wilcox.test(tank_means[, "Coord"], tank_means[, "Single"], paired = TRUE)
paired #V=305 and p-value is 0.0005246
ranks      <- rank(abs(fold$log10_fold))
rank_pos   <- sum(ranks[fold$log10_fold > 0])
rank_neg   <- sum(ranks[fold$log10_fold < 0])
rank_biserial <- (rank_pos - rank_neg) / sum(ranks) #generates the effect size
rank_biserial #0.7378917


# C: Spearman correlation
rho   <- cor(eggs$coord_min_per_day, eggs$Eggs, method = "spearman")
rho
rho_p <- suppressWarnings(cor.test(eggs$coord_min_per_day, eggs$Eggs,
                                   method = "spearman", exact = FALSE)$p.value)
rho_p

cat(sprintf("B per recording: n = %d, W = %.0f (smaller rank sum; V = %.0f), p = %.3g, rank-biserial r = %.3f, Coord > Single in %d/%d, median fold = %.2fx\n",
            nrow(fold), min(rank_pos, rank_neg), paired$statistic, paired$p.value, rank_biserial,
            sum(fold$log10_fold > 0), nrow(fold), 10^median(fold$log10_fold)))
cat(sprintf("C:          n = %d recordings, Spearman rho = %.3f, p = %.3g\n", nrow(eggs), rho, rho_p))

# text only results

#all visits pooled (Mann-Whitney test, Cliff's delta).
# (The Wilcoxon signed-rank test with matched-pairs rank-biserial correlation
# is the block above: `paired` and `rank_biserial`.)
mw      <- wilcox.test(duration_min ~ type, data = visits)
mw
n_type  <- table(visits$type)
cliff_d <- 1 - 2 * mw$statistic / prod(n_type)          # > 0 means Coord visits are longer
cliff_d
cat(sprintf("B pooled visits: n = %d (Single %d, Coord %d), medians %.0f vs %.0f min, Mann-Whitney p = %.3g, Cliff's delta = %.2f\n",
            nrow(visits), n_type["Single"], n_type["Coord"],
            median(visits$duration_min[visits$type == "Single"]),
            median(visits$duration_min[visits$type == "Coord"]), mw$p.value, cliff_d))

###figures###

theme_set(
  theme_classic(base_size = 11) +
    theme(axis.line    = element_line(linewidth = 0.48),
          axis.ticks   = element_line(linewidth = 0.48),
          axis.ticks.length = unit(1.5, "pt"),
          axis.text    = element_text(size = 11, color = "black"),
          axis.title   = element_text(size = 12),
          axis.title.x = element_text(margin = margin(t = 3)),
          axis.title.y = element_text(margin = margin(r = 3)),
          plot.title   = element_text(face = "bold", size = 15, hjust = 0, margin = margin(b = 2)),  # panel letter
          plot.title.position = "plot",
          plot.margin  = margin(t = 1, r = 4, b = 1, l = 2))
)


fold_ticks <- c(0.5, 1, 2, 4, 8, 16) #fold difference
set.seed(1)                                               # fixes the sideways jitter
panel_B <- ggplot(fold, aes(x = "", y = log10_fold)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey40") +
  geom_boxplot(width = 0.5, fill = "#B3B3B3", color = "#747474", alpha = 0.6,
               outlier.shape = NA, linewidth = 0.48) +
  geom_jitter(width = 0.16, height = 0, shape = 21, fill = "#B3B3B3", color = "black",
              size = 2.7, stroke = 0.48) +
  annotate("text", x = 1, y = Inf, label = stars(paired$p.value), vjust = 1.1,
           size = 14 / .pt, fontface = "bold") +
  scale_x_discrete(expand = expansion(add = 0.6)) +
  scale_y_continuous(breaks = log10(fold_ticks), labels = paste0(fold_ticks, "\u00d7"),
                     expand = expansion(mult = c(0.05, 0.12))) +
  labs(title = "B", y = "Coordinated / single-parent\nvisit duration (fold)",
       x = sprintf("n = %d recordings", nrow(fold))) +
  theme(axis.title.x = element_text(size = 11, vjust = 1),   # vjust: sits under the axis, level with C's first title line
        axis.ticks.x = element_blank())

callout <- data.frame(
  x = 74.4,
  y = 7.05 - 0.6545 * (0:2),
  label = c(sprintf("rho == \"%.2f\"", rho),
            sprintf("italic(p) == \"%.3f\"", rho_p),
            sprintf("n == %d", nrow(eggs))))

panel_C <- ggplot(eggs, aes(coord_min_per_day, Eggs)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              color = "black", linetype = "dashed", linewidth = 0.56) +
  geom_point(shape = 21, fill = "grey70", color = "black", size = 2.7, stroke = 0.48) +
  geom_text(data = callout, aes(x, y, label = label), parse = TRUE,
            hjust = 0, vjust = 1, size = 10 / .pt) +
  scale_x_continuous(breaks = c(0, 50, 100)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6)) +
  coord_cartesian(xlim = c(0, 120)) +
  labs(title = "C", x = "Coordinated visit duration\n(min/day)", y = "Eggs delivered")

width_B <- 124
width_C <- 222    
height  <- 245
figure  <- (panel_B | panel_C) + plot_layout(widths = c(width_B, width_C))

ggsave(file.path(fig_dir, paste0(out_name, ".pdf")), figure,
       width = (width_B + width_C) / 72, height = height / 72)
ggsave(file.path(fig_dir, paste0(out_name, ".png")), figure,
       width = (width_B + width_C) / 72, height = height / 72, dpi = 300, bg = "white")
