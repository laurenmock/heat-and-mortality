# Estimating excess mortality during heat waves using a causal inference framework and non-randomized time series
This repository contains all code for this project for reproducibility.

### Data
Data were colleccted as part of the National Morbidity, Mortality, and Air Pollution Study (NMMAPS).

### Code
[analysis/1_processing.R](analysis/1_processing.R): prepare data for analysis, EDA  
[analysis/2_matching.R](analysis/2_matching.R): match treated and control days  
[analysis/3_assess_matching.R](analysis/3_assess_matching.R): assess success of matching  
[analysis/4_inference.R](analysis/4_inference.R): compare mortality on treated and control days, create figures
