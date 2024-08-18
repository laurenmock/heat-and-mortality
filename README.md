# Estimating excess mortality during heat waves using a causal inference framework and non-randomized time series
This repository contains all code for this project for reproducibility.

### Data
Data were colleccted as part of the National Morbidity, Mortality, and Air Pollution Study (NMMAPS).

### Code
[analysis/1_preprocess_data.R](analysis/1_preprocess_data.R): prepare data for analysis  
[analysis/2_EDA.R](analysis/2_EDA.R): exploratory data analysis
[analysis/3_matching.R](analysis/3_matching.R): match control and treated days
[analysis/4_assess_matching.R](analysis/4_assess_matching.R): assess success of matching  
[analysis/5_inference.R](analysis/5_inference.R): compare mortality on treated and control days, create figures
[analysis/6_create_tables.R](analysis/6_create_tables.R):
