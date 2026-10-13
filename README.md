# tarpon_SIA_move_eco

Code and data for:

Griffin LP, Shipley ON, Adams AJ, Brownscombe JW, Ceriani SA, Cooke SJ, Dello Russo JJ, Newsome SD, Power M, Wilson JK, Danylchuk AJ (2026) Discrete foraging landscapes support large scale migrations of a marine fish. *Movement Ecology* 14:63. https://doi.org/10.1186/s40462-026-00677-3

## Data

`supplementary_data/` holds the article's Data S1 (inshore lizardfish δ13C and δ15N values used to build the isoscapes) and Data S2 (Atlantic tarpon fin and blood plasma values).

The inputs to the analysis are in `data/`:

- `data_cleaning/`: tarpon isotope values (ethanol-stored fin clips are adjusted by `script_1_tar_sia_f.R`)
- `isoscape_data/`: lizardfish values and the δ13C and δ15N isoscapes (means and variances)
- `assignment_data/`: depth, monthly space use from acoustic telemetry, and habitat maps from satellite telemetry, which set each fish's search area
- `move_data/`: predicted latitude by date, which narrows the search area

Tarpon capture locations are rounded to 0.1°, because they are guides' and anglers' fishing locations, and the outputs here were run from the rounded locations. The values in the article come from the exact locations, which the authors hold. The rounded run gives the same highest-probability cells, regions and counts, with distances a few km different (a mean fin distance of 298.6 rather than 299.3 km).

## Running

Open `tarpon_SIA_move_eco.Rproj` and run `source("run_all.R")` (a few minutes on 16 cores; on Windows it runs on one core and takes longer). Packages: dplyr, tidyr, readr, raster, sf, mvnmle, lubridate, here, ggplot2, cowplot, ggmap, pacman.

| Script (run in this order) | Does |
|---|---|
| script_1_tar_sia_f.R | paired ethanol vs frozen tests and ethanol-adjusted fin values |
| script_6_assignment_models_fin.R | fin assignments (5-month window, 95% latitude interval, TDF 3.0/3.5 ‰) |
| script_7_assignment_models_blood.R | plasma assignments (3-month window, 75% latitude interval, TDF 1.1/1.4 ‰) |
| script_6_7_assignment_models_output_summary.R | values reported in the Results, regional averages (Figure 2), fit to the isoscape |
| script_6_7_sensitivity.R | 121 TDF combinations for the 30 fish with fin and plasma samples |
| manuscript_figures.R | Figures 2–4 and S11–S14 |

`source_script_assignment_functions.R` holds the assignment model functions. `scripts/isoscapes/` has the INLA isoscape models as they were run; they need environmental layers and lizardfish files that are not included here, and their outputs are the isoscapes in `data/isoscape_data/`. The telemetry analyses behind the search areas are not included.

## Outputs

| Folder | Contents |
|---|---|
| output/ethanol | paired ethanol vs frozen tests, with and without T-2078 |
| output/fin, output/blood | highest-probability cell and distance for each sample, variance terms, and the Figure 2 regional surfaces |
| output/summary | values reported in the Results, regional averages, fit to the isoscape |
| output/sensitivity | TDF grids and counts of consistent assignments |
| figures | Figures 2–4 and S11–S14 as in the article |

## Notes

- **Seed.** A seed is set before the trophic discrimination draws, so reruns give the same results.
- **Consistent assignments** are those whose distance varies by less than 100 km across the 121 combinations. That gives 21 of 30 fin and 16 of 30 plasma samples with the 75% interval, and 12 of 30 fin samples with the 95% interval used for the main fin run. The rule measures distance, so for a few samples the highest-probability cell moves to another cell at a similar distance (`consistent_peak_moved` in `output/sensitivity/stability_summary.csv`).
- **T-2078.** This fish's ethanol clip is left out of the paired comparison (δ13C 12.3 ‰ below its frozen clip). With it included, the t-test is not significant (p = 0.13) but a Wilcoxon test is (p < 0.001).
- **Figures.** The files in `figures/` are the versions in the article, drawn from the exact locations. Rerunning `manuscript_figures.R` redraws them from the rounded locations. Figures S11 and S14 use a Google basemap. Set `GGMAP_GOOGLE_API_KEY` in `~/.Renviron`; without a key they are drawn on the coastline layer.

Contact: Lucas P. Griffin, University of South Florida (lpgriffin@usf.edu)
