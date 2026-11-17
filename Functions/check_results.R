# ============================================================================
# Check that a rebuild reproduces the published numbers
# November 11, 2026
#
# Reads the tables in Supplementary/derived_tables/ (written by the five figure scripts) and
# compares the headline numbers with the values in the reproducible report. Stops with a message at
# the first number that does not match. The mixed models of Figure 2 are not
# checked here because their p-values depend on the version of lme4.
#
# Run after the five figure scripts.
# ============================================================================

repo_dir <- here::here()
visits <- read.csv(file.path(repo_dir, "Supplementary", "derived_tables", "VisitDurations_byType_derived_v4.csv"))
care   <- read.csv(file.path(repo_dir, "Supplementary", "derived_tables", "CarePC1_classification.csv"))
rule   <- read.csv(file.path(repo_dir, "Supplementary", "derived_tables", "CarePC1_rule.csv"))
meta   <- read.csv(file.path(repo_dir, "Data", "recording_metadata.csv"))

check <- function(what, got, expected, tol = 0) {
  ok <- all(abs(got - expected) <= tol)
  cat(sprintf("  %-55s %s  (got %s)\n", what, if (ok) "ok" else "DOES NOT MATCH",
              paste(signif(got, 4), collapse = " / ")))
  if (!ok) stop("Rebuild does not match the published value for: ", what, call. = FALSE)
}

cat("Visit table\n")
check("recordings",                         length(unique(visits$tank)), 26)
check("visits",                             nrow(visits), 5805)
check("father-alone / mother-alone / coordinated",
      as.numeric(table(factor(visits$visitor, c("Male", "Female", "Family")))), c(3097, 1849, 859))
check("zero-length visits",                 sum(visits$duration_min == 0), 101)

cat("Six-metric PCA\n")
pca6 <- read.csv(file.path(repo_dir, "Supplementary", "derived_tables", "CarePCA_sixMetrics_loadings.csv"))
check("PC1 / PC2 variance explained (%)",   100 * c(pca6$VarExplained_PC1[1], pca6$VarExplained_PC2[1]), c(37.7, 24.6), 0.05)
check("loading of coordinated visits on PC1", pca6$PC1[pca6$Variable == "Visits_Family_perDay"], 0.535, 0.0005)
check("loading of father visits on PC1",    pca6$PC1[pca6$Variable == "Visits_Male_perDay"], 0.517, 0.0005)

cat("CarePC1\n")
check("PC1 variance explained (%)",         100 * rule$VarExplained_PC1[1], 61.2, 0.05)
check("cutoff",                             rule$Threshold[1], 0.156, 0.0005)
check("High / Low recordings",              c(sum(care$CareLevel == "High"), sum(care$CareLevel == "Low")), c(11, 15))

cat("Figure 1B\n")
pos   <- visits[visits$duration_min > 0, ]
type  <- ifelse(pos$visitor == "Family", "Coord", "Single")
m     <- tapply(log10(pos$duration_min), list(pos$tank, type), mean)
fold  <- m[, "Coord"] - m[, "Single"]
ranks <- rank(abs(fold))
check("visits with a duration",             nrow(pos), 5704)
check("Wilcoxon signed-rank W",             min(sum(ranks[fold > 0]), sum(ranks[fold < 0])), 46)
check("p",                                  wilcox.test(m[, "Coord"], m[, "Single"], paired = TRUE)$p.value, 0.000525, 5e-6)
check("recordings with longer coordinated visits", sum(fold > 0), 19)
check("median fold",                        10^median(fold), 1.98, 0.005)

cat("Figure 1C\n")
coord_min <- tapply(visits$duration_min[visits$visitor == "Family"],
                    visits$tank[visits$visitor == "Family"], sum)
eggs <- meta[!is.na(meta$EggFeedingEvents), ]
x    <- as.numeric(coord_min[eggs$ObservationID]) / eggs$NumberofFullDays
check("recordings scored for egg feeding",  nrow(eggs), 15)
check("Spearman rho",                       cor(x, eggs$EggFeedingEvents, method = "spearman"), 0.730, 0.0005)

cat("Figure 2C\n")
bl <- care[!is.na(care$bodylength), ]
check("recordings with body length",        nrow(bl), 23)
check("Spearman rho",                       cor(bl$Minutes_Family_perDay, bl$bodylength, method = "spearman"), 0.579, 0.0005)

cat("Figure 3 (TadBot)\n")
tb  <- read.csv(file.path(repo_dir, "Supplementary", "derived_tables", "TadBot_first3days.csv"))
off <- tb[tb$Treatment == "NegCtrl", ]; on <- tb[tb$Treatment == "TadBot", ]
stopifnot(identical(off$Tank, on$Tank))
check("tanks",                              nrow(off), 6)
check("first 3 days, all visits/day, off vs on: p", t.test(on$AllVisits_perDay, off$AllVisits_perDay, paired = TRUE)$p.value, 0.025, 0.0005)
td  <- read.csv(file.path(repo_dir, "Supplementary", "derived_tables", "TadBot_dailyVisits.csv"))
check("tank-days, trial days 1 to 14",      nrow(td), 156)
check("coordinated visits over those days, off / on", as.numeric(tapply(td$Visits_Family, td$Treatment, sum)[c("NegCtrl", "TadBot")]), c(71, 52))
check("visits over those days, off / on",   as.numeric(tapply(td$AllVisits, td$Treatment, sum)[c("NegCtrl", "TadBot")]), c(643, 574))
check("first 3 days, CarePC1, off vs on: p",     t.test(on$CarePC1, off$CarePC1, paired = TRUE)$p.value, 0.863, 0.0005)

cat("Figure 5 (naltrexone)\n")
ph   <- read.csv(file.path(repo_dir, "Supplementary", "derived_tables", "Pharmacology_perTank.csv"))
tank <- unique(ph[, c("Tank", "CareLevel")])
check("High / Low tanks",                   c(sum(tank$CareLevel == "High"), sum(tank$CareLevel == "Low")), c(3, 10))
hf  <- ph[ph$CareLevel == "High" & ph$Visitor == "Male", ]
sal <- hf[hf$Treatment == "Saline", ]; ntx <- hf[hf$Treatment == "Naltrexone", ]
stopifnot(identical(sal$Tank, ntx$Tank))
check("High-care fathers, saline vs naltrexone: p", t.test(ntx$MinutesPerDay, sal$MinutesPerDay, paired = TRUE)$p.value, 0.0364, 0.00005)

cat("Figure 5 (RNAscope pairs)\n")
rp <- read.csv(file.path(repo_dir, "Supplementary", "derived_tables", "RNAscope_CarePC1_perPair.csv"))
check("breeding pairs",                     nrow(rp), 18)
check("High / Low pairs",                   c(sum(rp$CareLevel == "High"), sum(rp$CareLevel == "Low")), c(14, 4))
check("coordinated visits, all pairs",      sum(rp$Visits_Family), 680)
check("coordinated minutes, all pairs",     sum(rp$Minutes_Family), 29748)

cat("All checks passed.\n")
