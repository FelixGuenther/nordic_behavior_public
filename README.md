# nordic_behavior_public

This is the public code repository for the manuscript "Quantifying the impact of social activities on SARS-CoV-2 transmission using Google mobility reports"
It contains data and a targets-Pipeline for estimating the proposed dynamic transmission model to data from two German regions using cmdstanr. 

Content:
- Folder R contains custom R functions for data preprocessing, postprocessing of the fitted model (MCMC samples), and figures
- Folder data contains raw data files
- Folder stan contains the dynamic transmission model in stan code
- _targets contains files related to the analysis pipeline based on the Rpackage targets (e.g., storage of intermediate results)

The Rscript _targets.R specifies the analysis pipeline using target syntax, it can be evaluated in parallel by evaluting run_parallel.R

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10514858.svg)](https://doi.org/10.5281/zenodo.10514858)

