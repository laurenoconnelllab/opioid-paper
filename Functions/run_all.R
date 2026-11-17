# Runs the five figure scripts in the order of the figures, from the raw data
# to all five figures, and then the checks. shared_functions.R holds the
# functions the figure scripts share.
# Open the .Rproj file (or set the working directory to the project folder), then:  source("Functions/run_all.R")
# Packages: readxl, ggplot2 (3.5 or newer for Figures 4 and 5), patchwork, png,
#           dplyr, tidyr, tibble, lme4, lmerTest, cowplot, scales, tidyverse
#           (Figure 5), here

steps <- c("Figure1_ethogram_panels_11162026.R",      # Figure 1: event logs -> one row per nest visit -> panels B and C
           "Figure2_11162026.R",                      # Figure 2: visits -> CarePC1 and High/Low care -> figure
           "Figure3_11162026.R",                      # Figure 3: TadBot tables -> panels B to E
           "Figure4_11162026.R",                      # Figure 4: mate-pair test (B), scatter (C), bubble plot (D), snRNA-seq images (E to G), assembly
           "Figure5_11162026.R",                      # Figure 5: CarePC1 of the RNAscope pairs and naltrexone tanks, micrograph image, RNAscope and naltrexone panels
           "check_results.R")                         # stops if the rebuilt tables do not give the published numbers

scripts_dir <- here::here("Functions")
for (step in steps) {
  cat("\n==== ", step, " ====\n", sep = "")
  source(file.path(scripts_dir, step), local = new.env())
}

cat("\n==== R and package versions used for this run ====\n")
print(sessionInfo())
