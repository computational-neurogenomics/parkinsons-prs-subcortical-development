# Analysis code for the manuscript: "Genetic liability for Parkinson's disease shapes subcortical brain structure in youth cohorts aged 8-30 years"

This repository contains the custom analytical pipeline used to investigate the association between Parkinson's disease polygenic risk scores (PRS) and subcortical brain volumes in the ABCD and QTIM youth cohorts.

## Software Requirements
* **R**: dplyr, purrr, neuroCombat, lme4, performance, broom.mixed, tidyr, openxlsx
* **Python 3**: pandas, numpy, seaborn, matplotlib, statsmodels

## Pipeline Overview
The analysis is divided into cohort-specific pipelines. Users must set their local working directory to the root of this repository and ensure their own phenotypic/genotypic data is stored in the respective ./data/ folders before running the scripts.

### ABCD Cohort (/ABCD_analysis/)
1. 1_abcd_preprocessing_pipeline_and_harmonisation.R: Merges phenotypic and PRS data, performs neuroCombat harmonisation for scanner effects, and standardises variables.
2. 2_abcd_ols_sensitivity_analyses_and_results_visualisation.py: Runs Ordinary Least-Squares (OLS) regressions, performs ICV/Scanner sensitivity analyses, and generates heatmap visualisations.

### QTIM Cohort (/QTIM_analysis/)
1. 1_qtim_preprocessing_pipeline.R: Merges twin cohort phenotypic and PRS data, formats the FID column for downstream models, and standardises variables.
2. 2_qtim_lme_and_icv_sens_analysis.R: Runs Linear Mixed-Effects (LME) models controlling for familial relatedness, followed by ICV sensitivity analyses.
3. 3_qtim_results_visualisation.py: Extracts LME statistical outputs to generate final heatmap visualisations.

## Data Privacy
In compliance with strict Data Use Agreements and ethical approvals, no raw participant-level phenotypic, imaging, or genetic data is hosted in this repository.
