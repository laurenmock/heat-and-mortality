# Estimating excess mortality during heat waves using a causal inference framework and non-randomized time series
This repository contains all code for this project for reproducibility.

### Data
Data were colleccted as part of the National Morbidity, Mortality, and Air Pollution Study (NMMAPS).

### Code
[R_code/1_processing.R](R_code/1_processing.R): prepare data for analysis, EDA  
[R_code/2_matching.R](R_code/2_matching.R): match treated and control days  
[R_code/3_assess_matching.R](3_assess_matching.R): assess success of matching  
[R_code/4_inference.R](4_inference.R): compare mortality on treated and control days, create figures
