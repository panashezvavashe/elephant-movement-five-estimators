# Home-range estimator comparison for African savannah elephants, Mana Pools National Park

Code and derived results for:

> Zvavashe, P., Mukomberanwa, N. T., Kunedzimwe, F., Taru, P. & Ngorima, P. *Methodological Choices Shape Inferences of African Savannah Elephant (Loxodonta africana) Spatial Ecology Across Temporal Scales in Mana Pools National Park, Zimbabwe.*

Archived version:

The study applies five home-range estimators to GPS tracks from ten collared elephants (85,727 screened fixes). The estimators are:

- autocorrelated kernel density estimation (AKDE);
- the dynamic Brownian bridge movement model (dBBMM);
- the step-anchored elliptical kernel (SAEK);
- time-geographic density (TGD);
- the adaptive local convex hull (aLoCoH).

It then asks whether estimators that disagree on absolute area still rank individuals the same way. The analysis also fits season-year AKDE models and combines them with `ctmm::meta()`.

## Data availability

**The GPS locations are not in this repository.** Elephant locations are sensitive because of poaching risk. They are held by the Zimbabwe Parks and Wildlife Management Authority (ZPWMA) and are available on reasonable request to ZPWMA (Scientific Services Unit, Marongora) and the corresponding author, subject to a data-sharing agreement.

The repository therefore contains:

- **all analysis code**;
- **derived results that hold no coordinates**: areas, model-selection tables, season-year fits, comparison statistics and screening counts. These are the files cited as Supporting Information in the paper.

No fitted model objects (`.rds`), rasters (`.tif`), polygons (`.shp`) or maps are included, because all of them encode locations.

## Requirements

- R 4.6.0 with: ctmm 1.3.0, move 4.2.7, adehabitatHR 0.4.22, adehabitatLT 0.3.29, sf 1.1.3, raster 3.6-32, terra, sp, ggplot2 4.0.3, ggrepel, cowplot, magick, dplyr, tidyr, reshape2, viridis, circular, geosphere, rstatix, vegan. The full session information is in `final/comparison/tables/SessionInfo.txt`.
- Python 3 with pandas and numpy, for `code/01_screening/table1_centroid_v6.py` only.
- On Windows, unset `PROJ_LIB`, `PROJ_DATA` and `GDAL_DATA` before running R if another GIS installation sets them. Most scripts do this themselves with `Sys.unsetenv()`.

## Repository layout

```
.
├── README.md
├── .gitignore                  excludes location data, model objects, rasters, maps and keys
├── data/
│   ├── README.md               how to obtain the GPS data
│   ├── raw/                    AWT collar exports go here (not in the repository)
│   └── processed/              screened fixes are written here (not in the repository)
├── code/
│   ├── 01_screening/           raw exports → screened fixes; Table 1; Appendix 14
│   ├── 02_estimators/          the five estimators and movement metrics (combined2.R + run helpers)
│   ├── 03_akde_fulltrack/      AKDE on every screened fix (Table 3)
│   ├── 04_seasonal/            season-year AKDE + ctmm::meta() (Table 5)
│   ├── 05_comparison/          rank, scale and overlap statistics (Table 6)
│   └── 06_figures/             figures and appendices
└── final/                      derived results (outputs of the scripts)
    ├── AKDE/
    ├── comparison/
    ├── data_quality/
    ├── dBBMM/
    ├── ETD/                    SAEK results (see "Naming" below)
    └── movement/
```

## Input data

`code/01_screening/build_screened_v6.R` reads the raw Africa Wildlife Tracking (AWT) collar exports from `data/raw/`, one file per collar named `Elephant <ID>.csv`. It writes `data/processed/all_elephants_screened_v6.csv`, and every later script reads that file. Collar 4557 is in the raw data but is excluded from this analysis.

**Run every script from the repository root**, for example `Rscript code/01_screening/build_screened_v6.R`. All paths are relative to the root, and outputs are written under `final/`.

## Run order

| Step | Script | Produces | Paper |
|---|---|---|---|
| 1 | `code/01_screening/build_screened_v6.R` | Screened fixes; `final/data_quality/tables/Screening_log_v6.csv` | Section 2.2; Appendix 14 |
| 2 | `code/01_screening/filter_counts_v6.R`, `fix_success_v6.R` | `Filter_counts_v6.csv`, `Fix_success_v6.csv` | Table 1; Appendix 14 |
| 3 | `code/01_screening/table1_centroid_v6.py` | `Table1_centroid_v6.csv` (data quality, wet/dry centroid shift) | Tables 1 and 5 |
| 4 | `code/02_estimators/combined2.R` | All five estimators, movement metrics, diurnal profile, area comparison table | Sections 2.3–2.8; Tables 2 and 4; Figures 2 and 5; Appendices 1–10 |
| 5 | `code/03_akde_fulltrack/full_track_ctmm_v6.R` | AKDE model selection and fits on every screened fix | Section 2.3 |
| 6 | `code/03_akde_fulltrack/adopt_fulltrack_akde_v6.R` | Writes the full-track AKDE as the AKDE results; updates `MethodComparison_Wide.csv` | Table 3; Figure 2; Appendices 2–3 |
| 7 | `code/04_seasonal/rerun_seasonal_v6.R` | Season-year AKDE fits and `ctmm::meta()` wet:dry ratios | Section 3.4; Table 5; Figures 3–4 |
| 8 | `code/04_seasonal/meta_neff10_v6.R` | Wet:dry meta-analysis without fits with N_eff < 10 | Section 3.4 |
| 9 | `code/05_comparison/phase3_stats_v6.R` | Kendall's W, Spearman ρ with Holm correction, log-scale Bland–Altman, agreement after removing the spread of fixes | Section 3.5; Table 6; Figure 7 |
| 10 | `code/05_comparison/kendallw_groups_v6.R` | Kendall's W within the occurrence and UD estimator groups | Section 3.5 |
| 11 | `code/05_comparison/phase3_overlap_v6.R` | Bhattacharyya coefficient of UDs, Jaccard index, centroid distance, UD mass inside the TGD/aLoCoH polygons | Section 3.5; Table 6; Figure 6 |
| 12 | `code/05_comparison/dbbmm_artefact_check_v6.R` | Share of each dBBMM isopleth more than 2 km from any fix | Section 3.3; Appendix 4 |
| 13 | `code/06_figures/diurnal_v6.R`, `saek_panels_v6.R` | Appendix 1 (local time) and Appendix 9 | Appendices 1 and 9 |
| 14 | `code/06_figures/figures_v6.R A`, then `… B` after step 7 | Figure 5 and Appendices 5 and 10 (part A); Figure 4 and Appendices 11 and 12 (part B) | Figures and appendices |

**Running `combined2.R` in parts.** The TGD and SAEK sections are slow. For the published run they were executed per animal and in parallel with `code/02_estimators/run_section_one_id.R` (`Rscript --max-ppsize=500000 code/02_estimators/run_section_one_id.R <section> <ID>`). `resume_v6_from_dBBMM.R` and `finish_v6_HtoL.R` then completed the remaining sections from the saved outputs. They run the same section code as `combined2.R`, so running `combined2.R` from start to finish gives the same results.

**Thinning.** In `combined2.R`, AKDE is fitted to tracks thinned to 1,500 fixes. That thinning hid velocity autocorrelation, so steps 5–6 replace those fits with fits to every screened fix, and the paper reports only the full-track fits. aLoCoH is still fitted to tracks thinned to at most 1,500 fixes, for computational tractability (Section 2.5). Season-year fits are never thinned.

## Naming: ETD = SAEK

The estimator called **SAEK** (step-anchored elliptical kernel) in the paper appears as **ETD** in script section names, folder names (`final/ETD/`) and column names. It is a simplified kernel inspired by, but not equivalent to, the Elliptical Time-Density model of Wall et al. (2014), and was renamed during revision. aLoCoH appears as `aLOCOH` in folder names.

## Derived results included

| File | Contents |
|---|---|
| `final/AKDE/tables/AKDE_Parameters.csv` | Per animal: selected model, n fixes, τ_p, τ_v, N_eff, 95% area with CI, ΔAIC to the best isotropic model |
| `final/AKDE/tables/ModelSelection_<ID>.csv`, `AKDE_Tests_<ID>.txt` | Candidate models with ΔAIC and Akaike weights; asymptote and isotropy tests |
| `final/AKDE/fulltrack_v6/FullTrack_<ID>.csv`, `FullTrack_vs_thinned.csv` | Full-track fits with τ confidence intervals; comparison with the thinned fits |
| `final/AKDE/seasonal_v6/SeasonYear_AKDE.csv` | Every season-year fit: dates, n fixes, model, N_eff, area with CI |
| `final/AKDE/seasonal_v6/SeasonYear_excluded.csv`, `SeasonYear_neff_lt10.csv` | Excluded season-years; fits with N_eff < 10 |
| `final/AKDE/seasonal_v6/Meta_population*.csv`, `PerAnimal_meta*.csv` | `ctmm::meta()` wet:dry ratios (all season-years, one per animal, without 4501, N_eff ≥ 10) |
| `final/AKDE/seasonal_v6/FullVsSeasonal_containment.csv` | Share of wet- and dry-season fixes outside the full-period 95% contour |
| `final/comparison/tables/MethodComparison_Wide.csv` | 95% area of all five estimators per animal |
| `final/comparison/v6/*.csv`, `*.txt` | Kendall's W, Spearman and Holm results, Bland–Altman, overlap statistics per animal and in summary |
| `final/data_quality/tables/*.csv` | Screening log, filter counts, fix success, Table 1 |
| `final/movement/tables/StepLength_Summary.csv`, `Diurnal_population_v6.csv` | Movement metrics (Table 2); population diurnal speed profile |
| `final/ETD/tables/ETD_sigma_sensitivity_all_v6.csv` | SAEK 95% area across σ = 250–2,000 m |
| `final/dBBMM/spatial/dBBMM_Sensitivity_<ID>.csv` | dBBMM 95% area across four margin/window settings (Appendix 5) |
| `final/dBBMM/tables/dBBMM_artefact_check_v6.csv` | dBBMM gap-artefact check |
| `final/comparison/tables/SessionInfo.txt` | R session information |

## Reproducibility

- Bootstrap confidence intervals use 10,000 case resamples with fixed seeds (`set.seed()` in `phase3_stats_v6.R`, `phase3_overlap_v6.R` and `kendallw_groups_v6.R`).
- Holm-adjusted p-values are computed from exact permutation p-values (`cor.test`, n = 10, no ties).
- `ctmm` model selection is deterministic for a given dataset and package version. Results may differ slightly with other ctmm versions.



## Contact

Panashe Zvavashe, Department of Geoinformatics and Environmental Conservation, Chinhoyi University of Technology. panashe.zvavashe@gmail.com
