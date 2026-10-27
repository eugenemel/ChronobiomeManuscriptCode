# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an R-based manuscript analysis codebase for the paper *"Disrupted Daytime Variation in Circulating Factors as a Predictor of Disease and Mortality"* (medRxiv 2026.01.16.26344298, posted as *"Disrupted Diurnal Phase Variation…"*). It models diurnal (time-of-day) variation in ~3,200 circulating biomarkers across ~500,000 UK Biobank participants using Generalized Additive Models (GAMs).

## Running the Analysis

This is an R Markdown (`.Rmd`) notebook project — there is no build/test/lint system. Notebooks are run chunk-by-chunk in RStudio or knitted via:

```r
rmarkdown::render("ChronobiomeUKBB.Rmd")
rmarkdown::render("ChronobiomeSerialOlink.Rmd")
```

To run a single chunk interactively, open the file in RStudio and use Cmd+Shift+Enter on the chunk.

Key R packages required: `dplyr`, `ggplot2`, `mgcv`, `broom`, `parallel`, `doParallel`, `reshape2`, `ggrepel`, `ggpubr`, `survival`. `SPEC.md` lists every package with the version used, and every UK Biobank field the analysis reads; update it when either changes.

## Architecture

### Analysis Pipeline (ChronobiomeUKBB.Rmd)

The main notebook follows a sequential pipeline with checkpointing:

1. **Load & join** — raw preprocessed `.Rds` files from `preprocessed/` are joined into a single dataframe `X` (one row per individual, columns per biomarker)
2. **Z-scale** — all biomarker columns are z-scaled across individuals
3. **Checkpoint** — `X` and `BIOMARKER_NAMES` saved to `outputs/savedX.Rdata`; section 4.2 restores it to skip preprocessing. `params$instance` picks the UKBB visit (0 = baseline, default; 1 = repeat assessment; 2 = imaging visit): instance N > 0 writes `outputs/savedX_instanceN.Rdata` and `outputs/saveOBSERVED_instanceN.Rds` instead; section 7.4 compares instance 1 against baseline.
4. **Fit splines** — `gamwrap()` fits GAMs for each biomarker × 6 population strata (good/poor health × male/female/all-sex-adjusted), results stored in `OBSERVED` → `outputs/saveOBSERVED.Rds`
5. **Downstream analyses** — sex differences (section 8), health differences (section 9), disease associations (section 10+) all load from the saved `.Rds` checkpoints

### Core Functions (ChronoFunctions.R)

- **`COVARIATES` / `SURVIVAL_COVARIATES` / `covariateTerms()`** — single source of truth for model adjustment. Every biomarker model uses `COVARIATES` (age, sex, smoking, town_depr, fastingtime, storagetime, chronotype, sleepduration); survival models use `SURVIVAL_COVARIATES` (adds hourminute, health_rating). Build formulas with `covariateTerms()` rather than hard-coding covariates.

- **`getGamFit(biomarker, rowsubset, ...)`** — fits `mgcv::bam()` with P-spline on `hourminute` (k=24), adjusting for `COVARIATES` (sex dropped unless `adjustforSex`; `dropCovariates` removes others, e.g. chronotype within a chronotype group). Returns either prediction values, plot data, or tidy coefficient table depending on flags.
- **`gamwrap(biomarker)`** — calls all 6 population-specific wrappers (`healthyMales`, `healthyFemales`, `poorMales`, `poorFemales`, `goodHealth`, `poorHealth`) and row-binds results.
- **`testDiffSpline(biomarker, popname, rowsA, rowsB)`** — tests whether two groups have different spline shapes: 100 refits of each group on random 80% subsamples, a GAM on the simulated difference curve, median statistics (~200 GAM fits per call). **`testDiffSplineFast`** has the same interface and output but fits each group once and draws the subsample variability from the fit's posterior (20–40× faster; agrees within the original's run-to-run noise, though it can be less conservative for outlier-prone biomarkers). The wrappers `testDiffSplineSex`, `testDiffSplineHealth`, `testDiffSplineWomen` and `modelPastFutureEvents` take `fast = TRUE` (default) or `fast = FALSE` for the original.
- **`testSplineInteraction(biomarker, ...)`** — direct GAM interaction test using `s(hourminute, by=group)` formula.
- **`plotGamFit / plotBiomarkerSet / plotGamFitDisease`** — visualization functions; `plotBiomarkerSet` returns a named list of ggplot objects for use with `ggpubr::ggarrange`.

### Key Data Objects

| Object | Description |
|---|---|
| `X` | Main modelling dataframe; one row per UKBB participant at baseline |
| `BIOMARKER_NAMES` | `data.frame(biomarker, type)` where type ∈ {CBC, ChemPanel, Protein, Lipids} |
| `OBSERVED` / `ALL` | GAM fit results; one row per biomarker × population; columns include `statistic`, `p.value`, `range`, `beta1`–`beta23` |
| `FDRLINE$fdr01` | Empirical FDR threshold for significance (statistic ≥ 6 ≈ 1% FDR) |

### Population Strata

Six strata used throughout: `good_health_males`, `good_health_females`, `poor_health_males`, `poor_health_females`, `good_health_all` (sex-adjusted), `poor_health_all` (sex-adjusted). Health is self-reported: Good/Excellent vs. Fair/Poor.

### Second Notebook (ChronobiomeSerialOlink.Rmd)

Compares UKBB GAM fits to a serial-sampling Olink proteomics dataset (the individual-level `support_files/2025-09-17_*_NPX-values.txt` files, which are gitignored; `support_files/README.md` describes each file). Uses `getGamFitSerial()` which fits a GAM (REML) with a random intercept for `Subject_id` (`s(Subject_id, bs = "re")`). It sources `ChronoFunctions.R` and loads the baseline UKBB checkpoints itself (`outputs/savedX.Rdata` unless `X` is already in the session, and `outputs/saveOBSERVED.Rds`). The `compareFits()` / `correlationSerialUKBB()` functions overlay UKBB and serial curves on the same plot.

## Data Files

Individual-level data is **not included** in the repository. Required files in `preprocessed/`:
- `UKBB_BaselineCovariates.Rds`, `UKBB_PrevCond.Rds`, `UKBB_HealthRating.Rds`
- `UKBB_BLOODSAMPLES.Rds`, `UKBB_CBC.Rds`, `UKBB_Biomarkers.Rds`, `UKBB_Lipids.Rds`
- `PPP_PLATE_ADJUSTED_IMPUTE1_3K.Rds` (2,919 Olink proteins)
- `UKBB_SLEEP.Rds` (chronotype and sleep duration covariates), `age_at_menopause.Rds` (section 15.5)
- `ALLEVENTS.Rds` (disease event outcomes for section 10+)

Preprocessed files are generated by `preprocess_UKBB.Rmd`, which requires the `ukbb` R package and institutional cloud credentials (`~/.ukbb_v5_auth.json`). The exceptions are `PPP_PLATE_ADJUSTED_IMPUTE1_3K.Rds`, `UKBB_SLEEP.Rds` and `age_at_menopause.Rds`, which were made outside this repository (`SPEC.md` section 2.1).

Intermediate `.Rds`/`.Rdata` checkpoints are saved in `outputs/` (`outputs/saveOBSERVED.Rds`, `outputs/savedX.Rdata`, etc.) and allow resuming the pipeline without rerunning expensive GAM fitting steps. Section 11.1 reads every `outputs/*_splinefit.Rds` as a disease fit, so give other checkpoints names that do not end that way.

## Parallelism

Computationally intensive `foreach` loops (sections 7–9 and 15) use forked workers via `registerDoParallel(cores = 6)`: they share `X` (~13 GB at baseline) with the session, so do not use `makeCluster()` with `.export = ls(globalenv())`, which copies everything to every worker. No `stopCluster()` is needed. Forking from RStudio is officially discouraged; if it hangs, run the step as a background job or with `Rscript`. The section 7 spline fit over all ~3,200 biomarkers takes about 8 minutes this way. Sections 14.1, 14.2 and 16.2 fork the same way with `mclapply(..., mc.cores = 6)`. `modelPastFutureEvents` (section 11) parallelizes over biomarkers internally (`cores = 6`, forked `mclapply`), so call it from a plain loop over diseases, not inside another parallel block; it takes roughly 20 minutes per disease.
