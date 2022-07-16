
###############################################################
#### Causal Inference for Time Series (Heat and Mortality) ####
###############################################################

### Functions to Define Treatment (treated vs. control) ###

library(Rcpp)


# These functions create the is_treated column.
# The is_treated column is TRUE for treated, FALSE for control, and NA for unusable (doesn't meet the criteria of either).


# FUNCTION INPUTS
#----------------
# data: data frame on which to add the is_treated column
# lag: number of preceding days to include in treatment 
# exposure: "is_high" column being used for treatment


#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#

### comparing low, low, low vs. high, high, high

treatment_LLL_vs_HHH <- function(data, lag, exposure){
  
  # use inputted exposure column name to get vector of exposure values
  selected_column <- eval(substitute(exposure), data)
  
  # initialize column for control, treatment, or NA
  data$is_treated <- rep(NA, nrow(data))
  
  # loop through rows of data frame, starting with lag+1
  for(l in (lag + 1):nrow(data)){
    
    # only look at days and lag days with non-missing exposure values
    if(all(!is.na(selected_column[(l - lag):l]))){
      
      # assign control day if current day and previous lag days were all low
      if(all((selected_column[(l - lag):l]) == FALSE)){
        data$is_treated[l] <- FALSE
      }
      
      # assign treated day if current day and previous lag days were all high
      if(all((selected_column[(l - lag):l]) == TRUE)){
        data$is_treated[l] <- TRUE
      }
    }
    
  }
  return(data)
}


#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#


###########################################################
################## Functions for Matching #################
###########################################################


# pairwise difference between real-univariate covariate of treated VS control group
cppFunction('NumericMatrix pairDistCpp(NumericVector treated, NumericVector control) {
            NumericMatrix D(treated.size(), control.size());
            for (int t = 0; t < treated.size(); t++) {
            for (int c = 0; c < control.size(); c++) {
            D(t,c) = treated[t] - control[c];
            }
            }
            return D;
            }')
# pairwise absolute difference between real-univariate covariates of treated VS control group
cppFunction('NumericMatrix pairAbsDistCpp(NumericVector treated, NumericVector control) {
            NumericMatrix D(treated.size(), control.size());
            for (int t = 0; t < treated.size(); t++) {
            for (int c = 0; c < control.size(); c++) {
            D(t,c) = abs(treated[t] - control[c]);
            }
            }
            return D;
}')
# pairwise difference between factor-valued (i.e. bounded integer-valued) covariates 
# (e.g. day of the week, month, ...) of treated VS control group, assuming the facotr levels are cyclic
# and only the shortest difference modulo nb_levels matters.
pairModuloDist = function(factors_treated, factors_control, nb_levels) {
  return (pmin(pairDistCpp(as.integer(factors_treated),as.integer(factors_control))%%nb_levels,
               t(pairDistCpp(as.integer(factors_control),as.integer(factors_treated))%%nb_levels)))
}
# pairwise difference between covariates of treated VS control group
# Inputs: treated/control are of covariate vectors (one entry per unit, for a given covariate)
# Outputs: pairwise difference matrix
pairdifference = function(treated, control){
  if(is.factor(treated[1])){
    # if factor-valued, use shortest difference modulo number of levels
    return (pairModuloDist(treated, control, length(levels(treated[1]))))
  } else {
    return (pairAbsDistCpp(treated,control))
  }
}
# pairwise discrepancy between treated VS control group
# Inputs: treated/control are lists or dataframes (with one column per covariate, one row per unit), 
# thresholds is a LIST of values for each covariate to be matched (match is admissible
# if and only if the differences between covariates are all less than the associated thresholds), 
# standard_deviations is an optional LIST of values to standardize/reweight the differences
discrepancyMatrix = function(treated, control, thresholds, scaling = NULL){
  nb_covariates = ncol(treated)
  # matrix of pairwise discrepancies (computed as standardized L1-distance)
  D = matrix(0, nrow = nrow(treated), ncol = nrow(control))
  # keep track of pairs that are non-admissible matches
  non_admissible = matrix(FALSE, nrow = nrow(treated), ncol = nrow(control))
  for (i in 1:nb_covariates){
    # only compute the distances for covariates that are matched on (i.e. finite thresholds)
    if (thresholds[[i]]<Inf){
      differences = pairdifference(treated[[i]], control[[i]])
      D = D + differences*scaling[[i]]
      # The user is responsible for inputing complete data (i.e. impute missing data beforhand if needed).
      # In the undesirable case where some covariates that are matched on are NA, we exclude the corresponding
      # unit from the matching (default behavior for convenience, but NOT for statistical validity, especially
      # if the missing-data mechanism is non-ignorable)
      differences[is.na(differences)] = Inf
      non_admissible = non_admissible | (differences > thresholds[[i]])
    }
  }
  D = D/nb_covariates # "standardize" the discrepancies (just for convenience, doesn't change the matching at all)
  D[non_admissible] = Inf # give infinite penalty to non-admissible pairs
  return (D)
}


