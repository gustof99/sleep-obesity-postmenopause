# Sleep behaviours and incident obesity in postmenopausal women

This repository contains the R analysis code for the study examining the prospective associations between multidimensional sleep behaviours and ICD-coded incident obesity among postmenopausal women in the UK Biobank.

## Data availability

Individual-level UK Biobank data cannot be publicly shared because they are available only to approved researchers through the UK Biobank application process. This repository therefore provides the analysis code, but not the raw data, derived individual-level datasets, or participant identifiers.

Researchers who obtain access to UK Biobank data can use this code to reproduce the cohort derivation, variable construction, statistical analyses, and generation of tables and figures.

## Code

`Y001_sleep_obesity.R` contains the analysis code for cohort derivation, variable construction, Cox regression models, restricted cubic spline analysis, sensitivity analyses, subgroup and interaction analyses, missing-data analyses, proportional hazards diagnostics, exploratory indirect-association analyses, and generation of tables and figures.

## Software

The analyses were conducted using R. Package versions and session information should be generated using the session information section in the R script.

## Random seed

The random seed was fixed at `20250611` for multiple imputation and bootstrap-based exploratory indirect-association analyses.

## Citation

If using this code, please cite the associated article and the archived version of this repository once the DOI is available.
