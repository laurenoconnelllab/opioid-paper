## Code, scripts, data for the statistical analyses of the paper:

**"Opioidergic neurons modulate parental attentiveness in biparental poison frogs"**

# Billie Goolsby - Nov 16 2026

Access the interactive file here: https://laurenoconnelllab.github.io/opioid-paper/Reproducible-report.html


## What is in this repository

| Path | What it is |
|---|---|
| `Data/` | Raw data -- nest-visit event logs, recording metadata, begging-frequency measurements, TadBot and naltrexone visit sheets, processed phosphoTRAP tables, RNAscope cell counts, and the finished single-nucleus RNA-seq panel images with their Python notebooks. Each file is described in the report. |
| `Functions/` | The R scripts the report runs, in the order of the figures: one script per figure (each goes from the raw data to the finished figure), the functions they share, a check script, and `run_all.R`. |
| `Figures/` | Figures written by the scripts. `Figures/artwork/` holds the images the scripts place into the figures: the Figure 2 photo strip, the Figure 4A drawing, the Figure 5 micrographs and the Figure 5 naltrexone timeline. |
| `Supplementary/derived_tables/` | Tables written by the scripts (visit table, CarePC1 index and care level, per-cohort tables, statistics tables). Included so each figure script can be run on its own. |
| `opioid-parental-care.Rproj` | RStudio project file. Open this first; the scripts find the folder from it. |

The report covers the quantitative panels of Figures 1 to 5. It does not produce the panels drawn by hand (Figure 1A, Figure 3A). The Figure 4A drawing, the single-nucleus RNA-seq panels (second row of Figure 4, made in Python), the Figure 5 micrographs and the Figure 5 naltrexone timeline are finished images that are seperately available to the reader.

## How to run

1. Open `opioid-parental-care.Rproj` in RStudio.
2. Run:

   ```r
   source("Functions/run_all.R")
   ```

It will rebuild everything in `Supplementary/derived_tables/` and `Figures/` from `Data/`.

Packages: `here`, `knitr`, `readxl`, `ggplot2` (3.5 or newer), `patchwork`, `png`, `dplyr`, `tidyr`, `tibble`, `lme4`, `lmerTest`, `cowplot`, `scales` and `tidyverse`. Rendering the report also needs [Quarto](https://quarto.org).

Tested with R 4.3.3 (Ubuntu 24.04), Quarto 1.6.42, ggplot2 4.0.3, patchwork 1.3.2, readxl 1.4.3, dplyr 1.1.4, tidyr 1.3.1, lme4 1.1-35.1. The versions used for the published page are printed at the end of the report.
