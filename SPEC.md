# Data and Software Specification

This file lists the UK Biobank (UKB) data fields the analysis uses and the R and package versions it was run with. Field IDs and titles follow the [UKB Showcase](https://biobank.ndph.ox.ac.uk/showcase/). The biomarker field lists in section 2.3 were matched to the columns of the preprocessed files and checked against the 3,226 biomarkers fitted in `outputs/saveOBSERVED.Rds`.

## 1. Software

R 4.4.1 (2024-06-14) on macOS, aarch64-apple-darwin20.

| Package | Version | Used for |
|---|---|---|
| dplyr | 1.1.4 | data manipulation throughout |
| ggplot2 | 3.5.1 | figures |
| mgcv | 1.9-1 | GAM fits (`bam`, `gam`); depends on nlme 3.1-166 |
| broom | 1.0.7 | coefficient tables from model fits |
| parallel | (base R) | forked `mclapply` workers |
| doParallel | 1.0.17 | forked `foreach` workers; with foreach 1.5.2 and iterators 1.0.14 |
| reshape2 | 1.4.4 | `dcast` / `melt` reshaping |
| ggrepel | 0.9.6 | point labels |
| ggpubr | 0.6.0 | multi-panel figures |
| survival | 3.8-3 | Cox models, time-to-event outcomes |
| survminer | 0.5.0 | Kaplan-Meier plot of the mortality risk score (section 13.6) |
| gbm3 | 3.0.1 | biomarker selection for the mortality model (section 13.5); <https://github.com/gbm-developers/gbm3> |
| pheatmap | 1.0.12 | biomarker heatmap (section 14.7) |
| gridExtra | 2.3 | FDR plots (section 14.6) |
| scales | 1.3.0 | axis scales |
| readr | 2.1.5 | reading `support_files/ICD10.csv` |
| MASS | 7.3-64 | robust regression line in pairs plots (`rlm`) |
| rmarkdown / knitr | 2.29 / 1.49 | knitting the notebooks |

Preprocessing only (`preprocess_UKBB.Rmd`):

| Package | Version | Used for |
|---|---|---|
| ukbb | 0.0.0.9002 | institutional package for reading UKB fields from cloud storage; not publicly available |
| mice | 3.17.0 | imputing missing baseline covariates |

To install these versions:

```r
install.packages("remotes")
versions <- c(dplyr = "1.1.4", ggplot2 = "3.5.1", mgcv = "1.9-1", nlme = "3.1-166", broom = "1.0.7",
              doParallel = "1.0.17", foreach = "1.5.2", iterators = "1.0.14", reshape2 = "1.4.4",
              ggrepel = "0.9.6", ggpubr = "0.6.0", survival = "3.8-3", survminer = "0.5.0", gbm3 = "3.0.1",
              pheatmap = "1.0.12", gridExtra = "2.3", scales = "1.3.0", readr = "2.1.5", MASS = "7.3-64",
              rmarkdown = "2.29", knitr = "1.49", mice = "3.17.0")
for (p in names(versions)) remotes::install_version(p, versions[[p]], upgrade = "never")
```

## 2. UK Biobank data

The analysis uses one assessment visit at a time, chosen by `params$instance` in `ChronobiomeUKBB.Rmd`: instance 0 (baseline, 2006–2010, the default), 1 (repeat assessment, 2012–2013) or 2 (first imaging visit, from 2014). The data were obtained under UK Biobank application 18448; individual-level data are not included in this repository.

### 2.1. Preprocessed files read by the analysis

| File in `preprocessed/` | UKB source | Made by |
|---|---|---|
| `UKBB_BLOODSAMPLES.Rds` | fields 68, 74, 3166 | `preprocess_UKBB.Rmd` (field 3166 is read from a separate CSV export) |
| `UKBB_CBC.Rds` | category 100081, 31 fields | `preprocess_UKBB.Rmd` |
| `UKBB_Biomarkers.Rds` | category 17518, 30 fields | `preprocess_UKBB.Rmd` |
| `UKBB_Lipids.Rds` | category 220, 251 fields | `preprocess_UKBB.Rmd` |
| `PPP_PLATE_ADJUSTED_IMPUTE1_3K.Rds` | UKB-PPP Olink NPX values, 2,919 proteins | outside this repository |
| `UKBB_BaselineCovariates.Rds` | fields 21003, 31, 189, 1239, 20161, 1558, 21001, 53 | `preprocess_UKBB.Rmd` |
| `UKBB_HealthRating.Rds` | field 2178 | `preprocess_UKBB.Rmd` |
| `UKBB_SLEEP.Rds` | fields 1160, 1170, 1180, 1190, 1200, 1210, 1220 | outside this repository |
| `age_at_menopause.Rds` | field 3581 | outside this repository |
| `UKBB_PrevCond.Rds` | fields 2443, 2453, 2966, 3627, 3761, 3786, 3894, 4012, 4022, 4056 (instance 0) | `preprocess_UKBB.Rmd` |
| `ALLEVENTS.Rds` | first occurrences, cancer and death registries, hospital inpatient records (section 2.4) | `preprocess_UKBB.Rmd` |

`preprocess_UKBB.Rmd` also extracts spirometry, grip strength, imaging availability, heel bone density, cognition, cardiac and MRI fat measures, genetic principal components, family history and blood pressure. The analysis notebooks do not read any of those files.

### 2.2. Time of day, covariates and strata

| Field | Title | Variable | Use |
|---|---|---|---|
| [3166](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=3166) | Time blood sample collected | `hourminute`, `storagetime` | time of day of the blood draw in decimal hours; only draws between 09:00 and 20:00 are kept. `storagetime` is the time from the draw to the latest draw in the file, in units of 10^8 s (about 3.2 years); covariate |
| [74](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=74) | Fasting time | `fastingtime` | covariate |
| [21003](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21003) | Age when attended assessment centre | `age` | covariate |
| [31](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=31) | Sex | `sex` | covariate in sex-adjusted fits; stratifier |
| [189](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=189) | Townsend deprivation index at recruitment | `town_depr` | covariate |
| [1239](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1239) | Current tobacco smoking | `smoking` | covariate, recoded as current smoker yes/no |
| [1180](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1180) | Morning/evening person (chronotype) | `chronotype` | covariate; exposure in section 15 |
| [1160](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1160) | Sleep duration | `sleepduration` | covariate |
| [2178](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=2178) | Overall health rating | `health_rating` | strata: good health (Excellent, Good) vs poor health (Fair, Poor); covariate in survival models |
| [53](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=53) | Date of attending assessment centre | `visit_date`, `first_visit_date` | start of follow-up |
| [3581](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=3581) | Age at menopause (last menstrual period) | `age_at_menopause` | post-menopause status (section 15.5) |

The model covariates are set once, in `COVARIATES` and `SURVIVAL_COVARIATES` at the top of `ChronoFunctions.R`.

These fields are extracted but not used in any model:

- [68](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=68) (number of blood samples taken), stored with the draw time.
- [20161](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20161) (pack years of smoking), [1558](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1558) (alcohol intake frequency) and [21001](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21001) (BMI). `preprocess_UKBB.Rmd` imputes missing values in `UKBB_BaselineCovariates.Rds` with `mice` (random forest, one imputation) using every column of that file, so these three fields do affect the imputed values of age, sex, Townsend index and smoking.
- [1170](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1170), [1190](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1190), [1200](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1200), [1210](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1210), [1220](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1220) (getting up in the morning, daytime napping, sleeplessness, snoring, daytime dozing), stored in `UKBB_SLEEP.Rds`.
- The self-reported prior conditions in `UKBB_PrevCond.Rds`: [2443](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=2443) diabetes, [2453](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=2453) cancer, and age diagnosed with high blood pressure ([2966](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=2966)), angina ([3627](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=3627)), hay fever/rhinitis/eczema ([3761](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=3761)), asthma ([3786](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=3786)), heart attack ([3894](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=3894)), deep-vein thrombosis ([4012](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=4012)), pulmonary embolism ([4022](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=4022)) and stroke ([4056](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=4056)). They are joined into the time-to-event tables but no model uses them.

### 2.3. Biomarkers

| Panel | UKB source | Fields | Biomarkers analysed | Type in `BIOMARKER_NAMES` |
|---|---|---|---|---|
| Blood count | [category 100081](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=100081) | 31 | 31 | `CBC` |
| Blood biochemistry | [category 17518](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=17518) | 30 | 30 | `ChemPanel` |
| NMR metabolomics (Nightingale) | [category 220](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220) | 251 | 246 | `Lipids` |
| Olink proteomics (UKB-PPP, Olink Explore 3072) | [category 1839](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=1839) | — | 2,919 | `Protein` |
| **Total** | | | **3,226** | |

- **NMR metabolomics:** albumin, apolipoprotein B, creatinine, glucose and HDL cholesterol are dropped from the NMR panel in section 1.1 because they have the same column names as the blood biochemistry measurements.
- **Olink:** NPX values are not individual Showcase fields. UKB releases them as separate record-level tables; category 1839 holds the per-sample fields 30900–30902 (proteins measured, plate, well). `PPP_PLATE_ADJUSTED_IMPUTE1_3K.Rds` holds 2,919 proteins at instance 0 (48,581 participants) and instance 2 (1,160). It was made outside this repository:
  - Quality control followed Sethi, Raj, Wright and Melamud, [*Plasma Proteomic Determinants of Common Causes of Mortality*](https://doi.org/10.21203/rs.3.rs-2626017/v1) (Research Square preprint, 2023), which follows Sun et al. 2022 ([bioRxiv 2022.06.17.496443](https://doi.org/10.1101/2022.06.17.496443)): NPX values as reported by Olink, with measurements that failed plate QC or sample QC removed, and proteins missing in most samples removed. That paper used the first release of 1,472 proteins; this file covers the full 2,919-protein release.
  - Plate effects were removed for each protein by regressing NPX on Olink plate ID.
  - Missing values were then imputed using mice lasso regression; the resulting dataset contains no missing values.
- **Scaling:** all biomarkers are z-scaled across participants before fitting (section 3).

<details>
<summary>Blood count fields (category 100081)</summary>

| Field | Title | Biomarker name in the analysis |
|---|---|---|
| [30000](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30000) | White blood cell (leukocyte) count | `white_blood_cell_leukocyte_count` |
| [30010](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30010) | Red blood cell (erythrocyte) count | `red_blood_cell_erythrocyte_count` |
| [30020](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30020) | Haemoglobin concentration | `haemoglobin_concentration` |
| [30030](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30030) | Haematocrit percentage | `haematocrit_percentage` |
| [30040](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30040) | Mean corpuscular volume | `mean_corpuscular_volume` |
| [30050](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30050) | Mean corpuscular haemoglobin | `mean_corpuscular_haemoglobin` |
| [30060](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30060) | Mean corpuscular haemoglobin concentration | `mean_corpuscular_haemoglobin_concentration` |
| [30070](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30070) | Red blood cell (erythrocyte) distribution width | `red_blood_cell_erythrocyte_distribution_width` |
| [30080](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30080) | Platelet count | `platelet_count` |
| [30090](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30090) | Platelet crit | `platelet_crit` |
| [30100](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30100) | Mean platelet (thrombocyte) volume | `mean_platelet_thrombocyte_volume` |
| [30110](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30110) | Platelet distribution width | `platelet_distribution_width` |
| [30120](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30120) | Lymphocyte count | `lymphocyte_count` |
| [30130](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30130) | Monocyte count | `monocyte_count` |
| [30140](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30140) | Neutrophill count | `neutrophill_count` |
| [30150](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30150) | Eosinophill count | `eosinophill_count` |
| [30160](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30160) | Basophill count | `basophill_count` |
| [30170](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30170) | Nucleated red blood cell count | `nucleated_red_blood_cell_count` |
| [30180](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30180) | Lymphocyte percentage | `lymphocyte_percentage` |
| [30190](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30190) | Monocyte percentage | `monocyte_percentage` |
| [30200](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30200) | Neutrophill percentage | `neutrophill_percentage` |
| [30210](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30210) | Eosinophill percentage | `eosinophill_percentage` |
| [30220](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30220) | Basophill percentage | `basophill_percentage` |
| [30230](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30230) | Nucleated red blood cell percentage | `nucleated_red_blood_cell_percentage` |
| [30240](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30240) | Reticulocyte percentage | `reticulocyte_percentage` |
| [30250](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30250) | Reticulocyte count | `reticulocyte_count` |
| [30260](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30260) | Mean reticulocyte volume | `mean_reticulocyte_volume` |
| [30270](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30270) | Mean sphered cell volume | `mean_sphered_cell_volume` |
| [30280](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30280) | Immature reticulocyte fraction | `immature_reticulocyte_fraction` |
| [30290](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30290) | High light scatter reticulocyte percentage | `high_light_scatter_reticulocyte_percentage` |
| [30300](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30300) | High light scatter reticulocyte count | `high_light_scatter_reticulocyte_count` |

</details>

<details>
<summary>Blood biochemistry fields (category 17518)</summary>

| Field | Title | Biomarker name in the analysis |
|---|---|---|
| [30600](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30600) | Albumin | `albumin` |
| [30610](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30610) | Alkaline phosphatase | `alkaline_phosphatase` |
| [30620](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30620) | Alanine aminotransferase | `alanine_aminotransferase` |
| [30630](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30630) | Apolipoprotein A | `apolipoprotein_a` |
| [30640](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30640) | Apolipoprotein B | `apolipoprotein_b` |
| [30650](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30650) | Aspartate aminotransferase | `aspartate_aminotransferase` |
| [30660](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30660) | Direct bilirubin | `direct_bilirubin` |
| [30670](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30670) | Urea | `urea` |
| [30680](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30680) | Calcium | `calcium` |
| [30690](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30690) | Cholesterol | `cholesterol` |
| [30700](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30700) | Creatinine | `creatinine` |
| [30710](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30710) | C-reactive protein | `c_reactive_protein` |
| [30720](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30720) | Cystatin C | `cystatin_c` |
| [30730](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30730) | Gamma glutamyltransferase | `gamma_glutamyltransferase` |
| [30740](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30740) | Glucose | `glucose` |
| [30750](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30750) | Glycated haemoglobin (HbA1c) | `glycated_haemoglobin_hba1c` |
| [30760](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30760) | HDL cholesterol | `hdl_cholesterol` |
| [30770](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30770) | IGF-1 | `igf_1` |
| [30780](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30780) | LDL direct | `ldl_direct` |
| [30790](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30790) | Lipoprotein A | `lipoprotein_a` |
| [30800](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30800) | Oestradiol | `oestradiol` |
| [30810](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30810) | Phosphate | `phosphate` |
| [30820](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30820) | Rheumatoid factor | `rheumatoid_factor` |
| [30830](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30830) | SHBG | `shbg` |
| [30840](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30840) | Total bilirubin | `total_bilirubin` |
| [30850](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30850) | Testosterone | `testosterone` |
| [30860](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30860) | Total protein | `total_protein` |
| [30870](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30870) | Triglycerides | `triglycerides` |
| [30880](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30880) | Urate | `urate` |
| [30890](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30890) | Vitamin D | `vitamin_d` |

</details>

<details>
<summary>NMR metabolomics fields (category 220)</summary>

| Field | Title | Biomarker name in the analysis |
|---|---|---|
| [20280](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20280) | Glucose-lactate | `glucose_lactate` |
| [20281](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20281) | Spectrometer-corrected alanine | `spectrometer_corrected_alanine` |
| [23400](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23400) | Total Cholesterol | `total_cholesterol` |
| [23401](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23401) | Total Cholesterol Minus HDL-C | `total_cholesterol_minus_hdl_c` |
| [23402](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23402) | Remnant Cholesterol (Non-HDL, Non-LDL -Cholesterol) | `remnant_cholesterol_non_hdl_non_ldl__cholesterol` |
| [23403](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23403) | VLDL Cholesterol | `vldl_cholesterol` |
| [23404](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23404) | Clinical LDL Cholesterol | `clinical_ldl_cholesterol` |
| [23405](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23405) | LDL Cholesterol | `ldl_cholesterol` |
| [23406](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23406) | HDL Cholesterol | not used: same name as a blood biochemistry field |
| [23407](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23407) | Total Triglycerides | `total_triglycerides` |
| [23408](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23408) | Triglycerides in VLDL | `triglycerides_in_vldl` |
| [23409](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23409) | Triglycerides in LDL | `triglycerides_in_ldl` |
| [23410](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23410) | Triglycerides in HDL | `triglycerides_in_hdl` |
| [23411](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23411) | Total Phospholipids in Lipoprotein Particles | `total_phospholipids_in_lipoprotein_particles` |
| [23412](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23412) | Phospholipids in VLDL | `phospholipids_in_vldl` |
| [23413](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23413) | Phospholipids in LDL | `phospholipids_in_ldl` |
| [23414](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23414) | Phospholipids in HDL | `phospholipids_in_hdl` |
| [23415](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23415) | Total Esterified Cholesterol | `total_esterified_cholesterol` |
| [23416](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23416) | Cholesteryl Esters in VLDL | `cholesteryl_esters_in_vldl` |
| [23417](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23417) | Cholesteryl Esters in LDL | `cholesteryl_esters_in_ldl` |
| [23418](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23418) | Cholesteryl Esters in HDL | `cholesteryl_esters_in_hdl` |
| [23419](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23419) | Total Free Cholesterol | `total_free_cholesterol` |
| [23420](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23420) | Free Cholesterol in VLDL | `free_cholesterol_in_vldl` |
| [23421](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23421) | Free Cholesterol in LDL | `free_cholesterol_in_ldl` |
| [23422](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23422) | Free Cholesterol in HDL | `free_cholesterol_in_hdl` |
| [23423](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23423) | Total Lipids in Lipoprotein Particles | `total_lipids_in_lipoprotein_particles` |
| [23424](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23424) | Total Lipids in VLDL | `total_lipids_in_vldl` |
| [23425](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23425) | Total Lipids in LDL | `total_lipids_in_ldl` |
| [23426](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23426) | Total Lipids in HDL | `total_lipids_in_hdl` |
| [23427](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23427) | Total Concentration of Lipoprotein Particles | `total_concentration_of_lipoprotein_particles` |
| [23428](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23428) | Concentration of VLDL Particles | `concentration_of_vldl_particles` |
| [23429](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23429) | Concentration of LDL Particles | `concentration_of_ldl_particles` |
| [23430](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23430) | Concentration of HDL Particles | `concentration_of_hdl_particles` |
| [23431](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23431) | Average Diameter for VLDL Particles | `average_diameter_for_vldl_particles` |
| [23432](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23432) | Average Diameter for LDL Particles | `average_diameter_for_ldl_particles` |
| [23433](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23433) | Average Diameter for HDL Particles | `average_diameter_for_hdl_particles` |
| [23434](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23434) | Phosphoglycerides | `phosphoglycerides` |
| [23435](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23435) | Triglycerides to Phosphoglycerides ratio | `triglycerides_to_phosphoglycerides_ratio` |
| [23436](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23436) | Total Cholines | `total_cholines` |
| [23437](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23437) | Phosphatidylcholines | `phosphatidylcholines` |
| [23438](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23438) | Sphingomyelins | `sphingomyelins` |
| [23439](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23439) | Apolipoprotein B | not used: same name as a blood biochemistry field |
| [23440](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23440) | Apolipoprotein A1 | `apolipoprotein_a1` |
| [23441](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23441) | Apolipoprotein B to Apolipoprotein A1 ratio | `apolipoprotein_b_to_apolipoprotein_a1_ratio` |
| [23442](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23442) | Total Fatty Acids | `total_fatty_acids` |
| [23443](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23443) | Degree of Unsaturation | `degree_of_unsaturation` |
| [23444](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23444) | Omega-3 Fatty Acids | `omega_3_fatty_acids` |
| [23445](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23445) | Omega-6 Fatty Acids | `omega_6_fatty_acids` |
| [23446](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23446) | Polyunsaturated Fatty Acids | `polyunsaturated_fatty_acids` |
| [23447](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23447) | Monounsaturated Fatty Acids | `monounsaturated_fatty_acids` |
| [23448](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23448) | Saturated Fatty Acids | `saturated_fatty_acids` |
| [23449](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23449) | Linoleic Acid | `linoleic_acid` |
| [23450](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23450) | Docosahexaenoic Acid | `docosahexaenoic_acid` |
| [23451](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23451) | Omega-3 Fatty Acids to Total Fatty Acids percentage | `omega_3_fatty_acids_to_total_fatty_acids_percentage` |
| [23452](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23452) | Omega-6 Fatty Acids to Total Fatty Acids percentage | `omega_6_fatty_acids_to_total_fatty_acids_percentage` |
| [23453](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23453) | Polyunsaturated Fatty Acids to Total Fatty Acids percentage | `polyunsaturated_fatty_acids_to_total_fatty_acids_percentage` |
| [23454](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23454) | Monounsaturated Fatty Acids to Total Fatty Acids percentage | `monounsaturated_fatty_acids_to_total_fatty_acids_percentage` |
| [23455](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23455) | Saturated Fatty Acids to Total Fatty Acids percentage | `saturated_fatty_acids_to_total_fatty_acids_percentage` |
| [23456](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23456) | Linoleic Acid to Total Fatty Acids percentage | `linoleic_acid_to_total_fatty_acids_percentage` |
| [23457](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23457) | Docosahexaenoic Acid to Total Fatty Acids percentage | `docosahexaenoic_acid_to_total_fatty_acids_percentage` |
| [23458](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23458) | Polyunsaturated Fatty Acids to Monounsaturated Fatty Acids ratio | `polyunsaturated_fatty_acids_to_monounsaturated_fatty_acids_ratio` |
| [23459](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23459) | Omega-6 Fatty Acids to Omega-3 Fatty Acids ratio | `omega_6_fatty_acids_to_omega_3_fatty_acids_ratio` |
| [23460](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23460) | Alanine | `alanine` |
| [23461](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23461) | Glutamine | `glutamine` |
| [23462](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23462) | Glycine | `glycine` |
| [23463](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23463) | Histidine | `histidine` |
| [23464](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23464) | Total Concentration of Branched-Chain Amino Acids (Leucine + Isoleucine + Valine) | `total_concentration_of_branched_chain_amino_acids_leucine___isoleucine___valine` |
| [23465](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23465) | Isoleucine | `isoleucine` |
| [23466](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23466) | Leucine | `leucine` |
| [23467](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23467) | Valine | `valine` |
| [23468](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23468) | Phenylalanine | `phenylalanine` |
| [23469](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23469) | Tyrosine | `tyrosine` |
| [23470](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23470) | Glucose | not used: same name as a blood biochemistry field |
| [23471](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23471) | Lactate | `lactate` |
| [23472](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23472) | Pyruvate | `pyruvate` |
| [23473](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23473) | Citrate | `citrate` |
| [23474](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23474) | 3-Hydroxybutyrate | `nmr_3_hydroxybutyrate` |
| [23475](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23475) | Acetate | `acetate` |
| [23476](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23476) | Acetoacetate | `acetoacetate` |
| [23477](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23477) | Acetone | `acetone` |
| [23478](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23478) | Creatinine | not used: same name as a blood biochemistry field |
| [23479](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23479) | Albumin | not used: same name as a blood biochemistry field |
| [23480](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23480) | Glycoprotein Acetyls | `glycoprotein_acetyls` |
| [23481](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23481) | Concentration of Chylomicrons and Extremely Large VLDL Particles | `concentration_of_chylomicrons_and_extremely_large_vldl_particles` |
| [23482](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23482) | Total Lipids in Chylomicrons and Extremely Large VLDL | `total_lipids_in_chylomicrons_and_extremely_large_vldl` |
| [23483](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23483) | Phospholipids in Chylomicrons and Extremely Large VLDL | `phospholipids_in_chylomicrons_and_extremely_large_vldl` |
| [23484](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23484) | Cholesterol in Chylomicrons and Extremely Large VLDL | `cholesterol_in_chylomicrons_and_extremely_large_vldl` |
| [23485](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23485) | Cholesteryl Esters in Chylomicrons and Extremely Large VLDL | `cholesteryl_esters_in_chylomicrons_and_extremely_large_vldl` |
| [23486](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23486) | Free Cholesterol in Chylomicrons and Extremely Large VLDL | `free_cholesterol_in_chylomicrons_and_extremely_large_vldl` |
| [23487](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23487) | Triglycerides in Chylomicrons and Extremely Large VLDL | `triglycerides_in_chylomicrons_and_extremely_large_vldl` |
| [23488](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23488) | Concentration of Very Large VLDL Particles | `concentration_of_very_large_vldl_particles` |
| [23489](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23489) | Total Lipids in Very Large VLDL | `total_lipids_in_very_large_vldl` |
| [23490](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23490) | Phospholipids in Very Large VLDL | `phospholipids_in_very_large_vldl` |
| [23491](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23491) | Cholesterol in Very Large VLDL | `cholesterol_in_very_large_vldl` |
| [23492](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23492) | Cholesteryl Esters in Very Large VLDL | `cholesteryl_esters_in_very_large_vldl` |
| [23493](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23493) | Free Cholesterol in Very Large VLDL | `free_cholesterol_in_very_large_vldl` |
| [23494](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23494) | Triglycerides in Very Large VLDL | `triglycerides_in_very_large_vldl` |
| [23495](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23495) | Concentration of Large VLDL Particles | `concentration_of_large_vldl_particles` |
| [23496](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23496) | Total Lipids in Large VLDL | `total_lipids_in_large_vldl` |
| [23497](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23497) | Phospholipids in Large VLDL | `phospholipids_in_large_vldl` |
| [23498](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23498) | Cholesterol in Large VLDL | `cholesterol_in_large_vldl` |
| [23499](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23499) | Cholesteryl Esters in Large VLDL | `cholesteryl_esters_in_large_vldl` |
| [23500](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23500) | Free Cholesterol in Large VLDL | `free_cholesterol_in_large_vldl` |
| [23501](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23501) | Triglycerides in Large VLDL | `triglycerides_in_large_vldl` |
| [23502](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23502) | Concentration of Medium VLDL Particles | `concentration_of_medium_vldl_particles` |
| [23503](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23503) | Total Lipids in Medium VLDL | `total_lipids_in_medium_vldl` |
| [23504](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23504) | Phospholipids in Medium VLDL | `phospholipids_in_medium_vldl` |
| [23505](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23505) | Cholesterol in Medium VLDL | `cholesterol_in_medium_vldl` |
| [23506](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23506) | Cholesteryl Esters in Medium VLDL | `cholesteryl_esters_in_medium_vldl` |
| [23507](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23507) | Free Cholesterol in Medium VLDL | `free_cholesterol_in_medium_vldl` |
| [23508](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23508) | Triglycerides in Medium VLDL | `triglycerides_in_medium_vldl` |
| [23509](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23509) | Concentration of Small VLDL Particles | `concentration_of_small_vldl_particles` |
| [23510](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23510) | Total Lipids in Small VLDL | `total_lipids_in_small_vldl` |
| [23511](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23511) | Phospholipids in Small VLDL | `phospholipids_in_small_vldl` |
| [23512](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23512) | Cholesterol in Small VLDL | `cholesterol_in_small_vldl` |
| [23513](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23513) | Cholesteryl Esters in Small VLDL | `cholesteryl_esters_in_small_vldl` |
| [23514](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23514) | Free Cholesterol in Small VLDL | `free_cholesterol_in_small_vldl` |
| [23515](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23515) | Triglycerides in Small VLDL | `triglycerides_in_small_vldl` |
| [23516](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23516) | Concentration of Very Small VLDL Particles | `concentration_of_very_small_vldl_particles` |
| [23517](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23517) | Total Lipids in Very Small VLDL | `total_lipids_in_very_small_vldl` |
| [23518](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23518) | Phospholipids in Very Small VLDL | `phospholipids_in_very_small_vldl` |
| [23519](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23519) | Cholesterol in Very Small VLDL | `cholesterol_in_very_small_vldl` |
| [23520](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23520) | Cholesteryl Esters in Very Small VLDL | `cholesteryl_esters_in_very_small_vldl` |
| [23521](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23521) | Free Cholesterol in Very Small VLDL | `free_cholesterol_in_very_small_vldl` |
| [23522](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23522) | Triglycerides in Very Small VLDL | `triglycerides_in_very_small_vldl` |
| [23523](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23523) | Concentration of IDL Particles | `concentration_of_idl_particles` |
| [23524](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23524) | Total Lipids in IDL | `total_lipids_in_idl` |
| [23525](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23525) | Phospholipids in IDL | `phospholipids_in_idl` |
| [23526](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23526) | Cholesterol in IDL | `cholesterol_in_idl` |
| [23527](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23527) | Cholesteryl Esters in IDL | `cholesteryl_esters_in_idl` |
| [23528](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23528) | Free Cholesterol in IDL | `free_cholesterol_in_idl` |
| [23529](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23529) | Triglycerides in IDL | `triglycerides_in_idl` |
| [23530](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23530) | Concentration of Large LDL Particles | `concentration_of_large_ldl_particles` |
| [23531](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23531) | Total Lipids in Large LDL | `total_lipids_in_large_ldl` |
| [23532](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23532) | Phospholipids in Large LDL | `phospholipids_in_large_ldl` |
| [23533](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23533) | Cholesterol in Large LDL | `cholesterol_in_large_ldl` |
| [23534](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23534) | Cholesteryl Esters in Large LDL | `cholesteryl_esters_in_large_ldl` |
| [23535](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23535) | Free Cholesterol in Large LDL | `free_cholesterol_in_large_ldl` |
| [23536](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23536) | Triglycerides in Large LDL | `triglycerides_in_large_ldl` |
| [23537](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23537) | Concentration of Medium LDL Particles | `concentration_of_medium_ldl_particles` |
| [23538](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23538) | Total Lipids in Medium LDL | `total_lipids_in_medium_ldl` |
| [23539](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23539) | Phospholipids in Medium LDL | `phospholipids_in_medium_ldl` |
| [23540](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23540) | Cholesterol in Medium LDL | `cholesterol_in_medium_ldl` |
| [23541](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23541) | Cholesteryl Esters in Medium LDL | `cholesteryl_esters_in_medium_ldl` |
| [23542](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23542) | Free Cholesterol in Medium LDL | `free_cholesterol_in_medium_ldl` |
| [23543](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23543) | Triglycerides in Medium LDL | `triglycerides_in_medium_ldl` |
| [23544](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23544) | Concentration of Small LDL Particles | `concentration_of_small_ldl_particles` |
| [23545](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23545) | Total Lipids in Small LDL | `total_lipids_in_small_ldl` |
| [23546](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23546) | Phospholipids in Small LDL | `phospholipids_in_small_ldl` |
| [23547](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23547) | Cholesterol in Small LDL | `cholesterol_in_small_ldl` |
| [23548](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23548) | Cholesteryl Esters in Small LDL | `cholesteryl_esters_in_small_ldl` |
| [23549](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23549) | Free Cholesterol in Small LDL | `free_cholesterol_in_small_ldl` |
| [23550](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23550) | Triglycerides in Small LDL | `triglycerides_in_small_ldl` |
| [23551](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23551) | Concentration of Very Large HDL Particles | `concentration_of_very_large_hdl_particles` |
| [23552](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23552) | Total Lipids in Very Large HDL | `total_lipids_in_very_large_hdl` |
| [23553](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23553) | Phospholipids in Very Large HDL | `phospholipids_in_very_large_hdl` |
| [23554](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23554) | Cholesterol in Very Large HDL | `cholesterol_in_very_large_hdl` |
| [23555](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23555) | Cholesteryl Esters in Very Large HDL | `cholesteryl_esters_in_very_large_hdl` |
| [23556](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23556) | Free Cholesterol in Very Large HDL | `free_cholesterol_in_very_large_hdl` |
| [23557](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23557) | Triglycerides in Very Large HDL | `triglycerides_in_very_large_hdl` |
| [23558](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23558) | Concentration of Large HDL Particles | `concentration_of_large_hdl_particles` |
| [23559](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23559) | Total Lipids in Large HDL | `total_lipids_in_large_hdl` |
| [23560](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23560) | Phospholipids in Large HDL | `phospholipids_in_large_hdl` |
| [23561](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23561) | Cholesterol in Large HDL | `cholesterol_in_large_hdl` |
| [23562](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23562) | Cholesteryl Esters in Large HDL | `cholesteryl_esters_in_large_hdl` |
| [23563](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23563) | Free Cholesterol in Large HDL | `free_cholesterol_in_large_hdl` |
| [23564](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23564) | Triglycerides in Large HDL | `triglycerides_in_large_hdl` |
| [23565](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23565) | Concentration of Medium HDL Particles | `concentration_of_medium_hdl_particles` |
| [23566](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23566) | Total Lipids in Medium HDL | `total_lipids_in_medium_hdl` |
| [23567](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23567) | Phospholipids in Medium HDL | `phospholipids_in_medium_hdl` |
| [23568](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23568) | Cholesterol in Medium HDL | `cholesterol_in_medium_hdl` |
| [23569](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23569) | Cholesteryl Esters in Medium HDL | `cholesteryl_esters_in_medium_hdl` |
| [23570](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23570) | Free Cholesterol in Medium HDL | `free_cholesterol_in_medium_hdl` |
| [23571](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23571) | Triglycerides in Medium HDL | `triglycerides_in_medium_hdl` |
| [23572](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23572) | Concentration of Small HDL Particles | `concentration_of_small_hdl_particles` |
| [23573](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23573) | Total Lipids in Small HDL | `total_lipids_in_small_hdl` |
| [23574](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23574) | Phospholipids in Small HDL | `phospholipids_in_small_hdl` |
| [23575](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23575) | Cholesterol in Small HDL | `cholesterol_in_small_hdl` |
| [23576](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23576) | Cholesteryl Esters in Small HDL | `cholesteryl_esters_in_small_hdl` |
| [23577](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23577) | Free Cholesterol in Small HDL | `free_cholesterol_in_small_hdl` |
| [23578](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23578) | Triglycerides in Small HDL | `triglycerides_in_small_hdl` |
| [23579](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23579) | Phospholipids to Total Lipids in Chylomicrons and Extremely Large VLDL percentage | `phospholipids_to_total_lipids_in_chylomicrons_and_extremely_large_vldl_percentage` |
| [23580](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23580) | Cholesterol to Total Lipids in Chylomicrons and Extremely Large VLDL percentage | `cholesterol_to_total_lipids_in_chylomicrons_and_extremely_large_vldl_percentage` |
| [23581](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23581) | Cholesteryl Esters to Total Lipids in Chylomicrons and Extremely Large VLDL percentage | `cholesteryl_esters_to_total_lipids_in_chylomicrons_and_extremely_large_vldl_percentage` |
| [23582](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23582) | Free Cholesterol to Total Lipids in Chylomicrons and Extremely Large VLDL percentage | `free_cholesterol_to_total_lipids_in_chylomicrons_and_extremely_large_vldl_percentage` |
| [23583](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23583) | Triglycerides to Total Lipids in Chylomicrons and Extremely Large VLDL percentage | `triglycerides_to_total_lipids_in_chylomicrons_and_extremely_large_vldl_percentage` |
| [23584](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23584) | Phospholipids to Total Lipids in Very Large VLDL percentage | `phospholipids_to_total_lipids_in_very_large_vldl_percentage` |
| [23585](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23585) | Cholesterol to Total Lipids in Very Large VLDL percentage | `cholesterol_to_total_lipids_in_very_large_vldl_percentage` |
| [23586](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23586) | Cholesteryl Esters to Total Lipids in Very Large VLDL percentage | `cholesteryl_esters_to_total_lipids_in_very_large_vldl_percentage` |
| [23587](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23587) | Free Cholesterol to Total Lipids in Very Large VLDL percentage | `free_cholesterol_to_total_lipids_in_very_large_vldl_percentage` |
| [23588](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23588) | Triglycerides to Total Lipids in Very Large VLDL percentage | `triglycerides_to_total_lipids_in_very_large_vldl_percentage` |
| [23589](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23589) | Phospholipids to Total Lipids in Large VLDL percentage | `phospholipids_to_total_lipids_in_large_vldl_percentage` |
| [23590](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23590) | Cholesterol to Total Lipids in Large VLDL percentage | `cholesterol_to_total_lipids_in_large_vldl_percentage` |
| [23591](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23591) | Cholesteryl Esters to Total Lipids in Large VLDL percentage | `cholesteryl_esters_to_total_lipids_in_large_vldl_percentage` |
| [23592](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23592) | Free Cholesterol to Total Lipids in Large VLDL percentage | `free_cholesterol_to_total_lipids_in_large_vldl_percentage` |
| [23593](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23593) | Triglycerides to Total Lipids in Large VLDL percentage | `triglycerides_to_total_lipids_in_large_vldl_percentage` |
| [23594](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23594) | Phospholipids to Total Lipids in Medium VLDL percentage | `phospholipids_to_total_lipids_in_medium_vldl_percentage` |
| [23595](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23595) | Cholesterol to Total Lipids in Medium VLDL percentage | `cholesterol_to_total_lipids_in_medium_vldl_percentage` |
| [23596](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23596) | Cholesteryl Esters to Total Lipids in Medium VLDL percentage | `cholesteryl_esters_to_total_lipids_in_medium_vldl_percentage` |
| [23597](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23597) | Free Cholesterol to Total Lipids in Medium VLDL percentage | `free_cholesterol_to_total_lipids_in_medium_vldl_percentage` |
| [23598](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23598) | Triglycerides to Total Lipids in Medium VLDL percentage | `triglycerides_to_total_lipids_in_medium_vldl_percentage` |
| [23599](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23599) | Phospholipids to Total Lipids in Small VLDL percentage | `phospholipids_to_total_lipids_in_small_vldl_percentage` |
| [23600](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23600) | Cholesterol to Total Lipids in Small VLDL percentage | `cholesterol_to_total_lipids_in_small_vldl_percentage` |
| [23601](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23601) | Cholesteryl Esters to Total Lipids in Small VLDL percentage | `cholesteryl_esters_to_total_lipids_in_small_vldl_percentage` |
| [23602](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23602) | Free Cholesterol to Total Lipids in Small VLDL percentage | `free_cholesterol_to_total_lipids_in_small_vldl_percentage` |
| [23603](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23603) | Triglycerides to Total Lipids in Small VLDL percentage | `triglycerides_to_total_lipids_in_small_vldl_percentage` |
| [23604](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23604) | Phospholipids to Total Lipids in Very Small VLDL percentage | `phospholipids_to_total_lipids_in_very_small_vldl_percentage` |
| [23605](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23605) | Cholesterol to Total Lipids in Very Small VLDL percentage | `cholesterol_to_total_lipids_in_very_small_vldl_percentage` |
| [23606](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23606) | Cholesteryl Esters to Total Lipids in Very Small VLDL percentage | `cholesteryl_esters_to_total_lipids_in_very_small_vldl_percentage` |
| [23607](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23607) | Free Cholesterol to Total Lipids in Very Small VLDL percentage | `free_cholesterol_to_total_lipids_in_very_small_vldl_percentage` |
| [23608](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23608) | Triglycerides to Total Lipids in Very Small VLDL percentage | `triglycerides_to_total_lipids_in_very_small_vldl_percentage` |
| [23609](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23609) | Phospholipids to Total Lipids in IDL percentage | `phospholipids_to_total_lipids_in_idl_percentage` |
| [23610](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23610) | Cholesterol to Total Lipids in IDL percentage | `cholesterol_to_total_lipids_in_idl_percentage` |
| [23611](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23611) | Cholesteryl Esters to Total Lipids in IDL percentage | `cholesteryl_esters_to_total_lipids_in_idl_percentage` |
| [23612](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23612) | Free Cholesterol to Total Lipids in IDL percentage | `free_cholesterol_to_total_lipids_in_idl_percentage` |
| [23613](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23613) | Triglycerides to Total Lipids in IDL percentage | `triglycerides_to_total_lipids_in_idl_percentage` |
| [23614](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23614) | Phospholipids to Total Lipids in Large LDL percentage | `phospholipids_to_total_lipids_in_large_ldl_percentage` |
| [23615](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23615) | Cholesterol to Total Lipids in Large LDL percentage | `cholesterol_to_total_lipids_in_large_ldl_percentage` |
| [23616](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23616) | Cholesteryl Esters to Total Lipids in Large LDL percentage | `cholesteryl_esters_to_total_lipids_in_large_ldl_percentage` |
| [23617](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23617) | Free Cholesterol to Total Lipids in Large LDL percentage | `free_cholesterol_to_total_lipids_in_large_ldl_percentage` |
| [23618](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23618) | Triglycerides to Total Lipids in Large LDL percentage | `triglycerides_to_total_lipids_in_large_ldl_percentage` |
| [23619](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23619) | Phospholipids to Total Lipids in Medium LDL percentage | `phospholipids_to_total_lipids_in_medium_ldl_percentage` |
| [23620](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23620) | Cholesterol to Total Lipids in Medium LDL percentage | `cholesterol_to_total_lipids_in_medium_ldl_percentage` |
| [23621](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23621) | Cholesteryl Esters to Total Lipids in Medium LDL percentage | `cholesteryl_esters_to_total_lipids_in_medium_ldl_percentage` |
| [23622](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23622) | Free Cholesterol to Total Lipids in Medium LDL percentage | `free_cholesterol_to_total_lipids_in_medium_ldl_percentage` |
| [23623](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23623) | Triglycerides to Total Lipids in Medium LDL percentage | `triglycerides_to_total_lipids_in_medium_ldl_percentage` |
| [23624](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23624) | Phospholipids to Total Lipids in Small LDL percentage | `phospholipids_to_total_lipids_in_small_ldl_percentage` |
| [23625](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23625) | Cholesterol to Total Lipids in Small LDL percentage | `cholesterol_to_total_lipids_in_small_ldl_percentage` |
| [23626](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23626) | Cholesteryl Esters to Total Lipids in Small LDL percentage | `cholesteryl_esters_to_total_lipids_in_small_ldl_percentage` |
| [23627](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23627) | Free Cholesterol to Total Lipids in Small LDL percentage | `free_cholesterol_to_total_lipids_in_small_ldl_percentage` |
| [23628](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23628) | Triglycerides to Total Lipids in Small LDL percentage | `triglycerides_to_total_lipids_in_small_ldl_percentage` |
| [23629](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23629) | Phospholipids to Total Lipids in Very Large HDL percentage | `phospholipids_to_total_lipids_in_very_large_hdl_percentage` |
| [23630](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23630) | Cholesterol to Total Lipids in Very Large HDL percentage | `cholesterol_to_total_lipids_in_very_large_hdl_percentage` |
| [23631](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23631) | Cholesteryl Esters to Total Lipids in Very Large HDL percentage | `cholesteryl_esters_to_total_lipids_in_very_large_hdl_percentage` |
| [23632](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23632) | Free Cholesterol to Total Lipids in Very Large HDL percentage | `free_cholesterol_to_total_lipids_in_very_large_hdl_percentage` |
| [23633](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23633) | Triglycerides to Total Lipids in Very Large HDL percentage | `triglycerides_to_total_lipids_in_very_large_hdl_percentage` |
| [23634](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23634) | Phospholipids to Total Lipids in Large HDL percentage | `phospholipids_to_total_lipids_in_large_hdl_percentage` |
| [23635](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23635) | Cholesterol to Total Lipids in Large HDL percentage | `cholesterol_to_total_lipids_in_large_hdl_percentage` |
| [23636](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23636) | Cholesteryl Esters to Total Lipids in Large HDL percentage | `cholesteryl_esters_to_total_lipids_in_large_hdl_percentage` |
| [23637](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23637) | Free Cholesterol to Total Lipids in Large HDL percentage | `free_cholesterol_to_total_lipids_in_large_hdl_percentage` |
| [23638](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23638) | Triglycerides to Total Lipids in Large HDL percentage | `triglycerides_to_total_lipids_in_large_hdl_percentage` |
| [23639](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23639) | Phospholipids to Total Lipids in Medium HDL percentage | `phospholipids_to_total_lipids_in_medium_hdl_percentage` |
| [23640](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23640) | Cholesterol to Total Lipids in Medium HDL percentage | `cholesterol_to_total_lipids_in_medium_hdl_percentage` |
| [23641](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23641) | Cholesteryl Esters to Total Lipids in Medium HDL percentage | `cholesteryl_esters_to_total_lipids_in_medium_hdl_percentage` |
| [23642](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23642) | Free Cholesterol to Total Lipids in Medium HDL percentage | `free_cholesterol_to_total_lipids_in_medium_hdl_percentage` |
| [23643](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23643) | Triglycerides to Total Lipids in Medium HDL percentage | `triglycerides_to_total_lipids_in_medium_hdl_percentage` |
| [23644](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23644) | Phospholipids to Total Lipids in Small HDL percentage | `phospholipids_to_total_lipids_in_small_hdl_percentage` |
| [23645](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23645) | Cholesterol to Total Lipids in Small HDL percentage | `cholesterol_to_total_lipids_in_small_hdl_percentage` |
| [23646](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23646) | Cholesteryl Esters to Total Lipids in Small HDL percentage | `cholesteryl_esters_to_total_lipids_in_small_hdl_percentage` |
| [23647](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23647) | Free Cholesterol to Total Lipids in Small HDL percentage | `free_cholesterol_to_total_lipids_in_small_hdl_percentage` |
| [23648](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=23648) | Triglycerides to Total Lipids in Small HDL percentage | `triglycerides_to_total_lipids_in_small_hdl_percentage` |

</details>

### 2.4. Disease and death outcomes

`ALLEVENTS.Rds` holds one row per participant per diagnosis: the earliest date of each 3-character ICD-10 code across these sources.

| Source | UKB fields or tables |
|---|---|
| First occurrences | the "Date … first reported" fields in categories 2401–2417 (field IDs 130000–132604) |
| Cancer registry | [40005](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=40005) date of cancer diagnosis, [40006](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=40006) type of cancer (ICD-10) |
| Death registry | [40000](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=40000) date of death, plus `date_of_death` from the record-level `death` table; coded as `Death` |
| Hospital inpatient records (HES) | record-level tables `hesin` (`epistart`) and `hesin_diag` (`diag_icd10`) |

Follow-up starts at the assessment visit (field 53). Participants without an event are censored at death or at `params$censor_date` (2023-09-01, the latest event date in the file). `preprocess_UKBB.Rmd` also extracts causes of death and age at death ([40001](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=40001), [40002](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=40002), [40007](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=40007)), but the analysis does not use them.

## 3. Other data

| File | Contents | In repository |
|---|---|---|
| `support_files/2025-09-17_*_NPX-values.txt` (3 files) | serial-sampling Olink study: sample metadata and NPX values for all proteins, read by `ChronobiomeSerialOlink.Rmd`, and NPX values for a 201-protein subset, not read | no (individual-level) |
| `support_files/Cosinor-lmer_Fit-across-all-samples_Olink_Plasma_Proteomics.txt` | cosinor fits per protein from the serial-sampling study, read by `ChronobiomeSerialOlink.Rmd` and by sections 13.4 and 14.9 of `ChronobiomeUKBB.Rmd` | yes |
| `support_files/ICD10.csv` | ICD-10 code titles, read in sections 10.2 and 11.1 | yes |
