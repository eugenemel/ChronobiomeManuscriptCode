# Supporting files

## Reference and summary files

- `ICD10.csv` - ICD-10 code titles: 12,131 codes with columns `id` (code), `text` (title) and `parent` (parent code). Used to label diseases by their 3-character code in `ChronobiomeUKBB.Rmd`, sections 10.2 and 11.1.
- `Cosinor-lmer_Fit-across-all-samples_Olink_Plasma_Proteomics.txt` - summary results from the serial-sampling Olink study (Chronobiome [1]): a cosinor fit (mixed model across all subjects) for each of 2,925 proteins, one row per Olink assay (2,943 rows; 6 proteins, such as IL6 and TNF, are measured on all four Explore panels and appear four times). Columns: `OlinkID`, `UniProt`, `Assay` (protein symbol), `Panel`, `Cosinor_Pvalue`, `Cosinor_Qvalue`, `Test_type`, `Amplitude`, `Phi` (phase in radians), `Phase` (peak time, hours). Used in `ChronobiomeUKBB.Rmd` section 14.9 to compare UKBB and Chronobiome cyclers. Also read in section 13.4 and in `ChronobiomeSerialOlink.Rmd`, where none of its columns are used.

## Serial-sampling study - individual-level data

Not included in the repository; can be requested from the corresponding authors of the Chronobiome study [1]. Plasma from 20 healthy volunteers (10 old, 10 young; 9 women, 11 men) sampled every 4 hours over about two days (9 samples per subject, 179 samples), measured with Olink Explore (NPX values).

- `2025-09-17_Common-chronobiome-and-ukbb-cyclers_Sample-metadata_AB-proteomics_Plasma_NPX-values.txt` - one row per sample: `Sample_id`, `Subject_id`, `Sex`, `Age_group` (Old, Young), `Day`, `Clock_time` (hhmm) and `Cont_time` (hours since the start of sampling). Read by `ChronobiomeSerialOlink.Rmd`.
- `2025-09-17_All-protein-assays_Feature-matrix_AB-proteomics_Plasma_NPX-values.txt` - NPX values for 2,925 proteins: one row per protein (`Feature_id` = OlinkID, `Symbol`), one column per sample. Read by `ChronobiomeSerialOlink.Rmd`.

## Reference

1. Skarke C, Lahens NF, Mrcela A, Olsson N, Glessner J, Naik A, Theken KN, Chivily C, Douglas D, Hess KN, Gaun A, Moustafa AM, Daniel SG, Bittinger K, Teegarden S, Lordan R, El Jamal N, Meng H, Hollingsworth T, Woolfork A, Sengupta A, Brooks TG, Weljie A, Grant G, O'Brien J, McAllister F, Melamud E, FitzGerald GA. *Age and the Diurnal Oscillatory Features of the Human Chronobiome.* medRxiv 2026.01.21.26344528 (2026). doi: [10.64898/2026.01.21.26344528](https://doi.org/10.64898/2026.01.21.26344528)
