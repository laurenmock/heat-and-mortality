
###############################################################
#### Causal Inference for Time Series (Heat and Mortality) ####
###############################################################

### Script 4 ###
# Inputs: Processed data from script 1 and matched data from script 2
# Causal inference!

#------------------------------------------------------------------------#

library(ggplot2)
library(tidyverse)
library(pander)
library(MASS)
library(gridExtra)
library(kableExtra)

#------------------------------------------------------------------------#
#--------- set working directory/paths ---------#
#------------------------------------------------------------------------#

# set path for processed data from script 1 and matched data from script 2
processed_data_path <- "data/processed/"

# set path for temporary Fisher data (will delete later!)
Fisher_path <- "data/processed/Fisher_delete_later/"

# set path to write out PDF with figures (use "" for general path)
figures_path <- "figures/"

# open PDF (will automatically write in figures after running the whole script)
pdf(file = paste0(figures_path, "inference.pdf"))


# load matched data
matched <- read.csv(paste0(processed_data_path, "matched.csv"))
# change row names to I:N
row.names(matched) <- seq(1:nrow(matched))


#------------------------------------------------------------------------#

# select only the columns I will need
columns <- c("date", "city", "dow", "month", "year", "week",
             "tmax", "tmax_lag_1", "tmax_lag_2",
             "tmax_lag_3", "tmax_lag_4", "tmax_lag_5",
             "pm10_lag_3", "pm10_lag_4", 
             "pm25_lag_3", "pm25_lag_4",
             "o3_lag_3", "o3_lag_4", 
             "no2_lag_3", "no2_lag_4", 
             "so2_lag_3", "so2_lag_4", 
             "co_lag_3", "co_lag_4",
             "is_treated", "death_sum3", "id", "pair")

matched <- matched %>% dplyr::select(all_of(columns))



##############################################################################################
#---------------------------------------  FREQUENTIST  --------------------------------------#
##############################################################################################

# city names
cities <- unique(matched$city)

# calculate sample means by city and treatment group
mean_table <- matched %>% 
  group_by(city, is_treated) %>%
  summarise(mean = mean(death_sum3),
            n = n()*2) %>%
  spread(key = "is_treated", value = "mean") %>%
  rename(mean.control = `FALSE`,
         mean.treated = `TRUE`) %>%
  mutate(tau = mean.treated - mean.control)


# set N and K for each city
N <- mean_table$n
K <- mean_table$n/2

# --------------- ignoring matched pairs to get variance/CIs (block) ---------------#

# calculate sample variances by city and treatment group
var_table <- matched %>%
  group_by(city, is_treated) %>%
  summarise(variance = var(death_sum3)) %>%
  spread(key = "is_treated", value = "variance") %>%
  rename(var.control = `FALSE`,
         var.treated = `TRUE`)

# join into one table
stats_table <- left_join(mean_table, var_table, by = "city")

# add column for variance of tau (unpaired)
stats_table$var.tau.unpaired <- (stats_table$var.control + stats_table$var.treated) / (stats_table$n/2)

# confidence intervals
stats_table$lower.CI.tau.unpaired <- stats_table$tau - (1.96*sqrt(stats_table$var.tau.unpaired))
stats_table$upper.CI.tau.unpaired <- stats_table$tau + (1.96*sqrt(stats_table$var.tau.unpaired))


# --------------- using matched pairs to get variance/CIs (paired) ---------------#

# DR chapter 3, pages 37, 38, 39

# order by pair (with each FALSE day first)
matched <- matched %>%
  arrange(pair, is_treated)


# initialize columns for variance of tau and 95% CIs (using pairs)
stats_table$lower.CI.tau.paired <- NA
stats_table$upper.CI.tau.paired <- NA

# loop through 5 cities
for(i in 1:length(cities)){
  
  # filter to one city
  one_city <- matched %>% filter(city == cities[i])
  
  # vector of differences in deaths on control vs. treated days for each pair (keep only odd diffs between rows)
  D_k <- diff(one_city$death_sum3, lag = 1)[seq_len(stats_table$n[i]) %% 2 == 1]
  
  # mean difference between control vs. treated (same as tau in the stats_table)
  D_bar <- mean(D_k)
  
  # (s2_D)/k --> unbiased estimate of var(tau)
  var.tau <- sum((D_k-D_bar)^2) / (K[i]*(K[i]-1))
  
  # paired t-statistic
  t_k_0.025 <- qt(0.025, df = nrow(one_city)/2) %>% abs()
  
  # 95% confidence intervals
  stats_table$lower.CI.tau.paired[i] <- D_bar - t_k_0.025*sqrt(var.tau)
  stats_table$upper.CI.tau.paired[i] <- D_bar + t_k_0.025*sqrt(var.tau)
  
}


# --------------- add Young's results (Bayesian) to table ---------------#

# will need to update these!!!

stats_table$tau.bayes <- c(9.36, 7.72, 25.06, 1.91, 3.88)
stats_table$var.tau.bayes <- c(2.91, 2.62, 5.60, 1.02, 1.08)^2

stats_table$lower.CI.tau.bayes <- stats_table$tau.bayes - (1.96*sqrt(stats_table$var.tau.bayes))
stats_table$upper.CI.tau.bayes <- stats_table$tau.bayes + (1.96*sqrt(stats_table$var.tau.bayes))


#------------------------------ REGRESSION ------------------------------#

# histograms of deaths
# all cities
matched %>%
  ggplot() +
  geom_density(aes(death_sum3, fill = city), col = "white", alpha = 0.7) +
  ggtitle("Distribution of Death Counts") +
  xlab("Deaths over 3 Day Experiment") +
  theme_classic()
# split by city
matched %>%
  ggplot() +
  geom_histogram(aes(death_sum3), bins = 12, col = "white") +
  facet_wrap(~city, scales= "free") +
  xlab("Deaths over 3 Day Experiment") +
  theme_bw()


# which model?
city_means <- tapply(matched$death_sum3, matched$city, mean)
city_vars <- tapply(matched$death_sum3, matched$city, var)
data.frame(city_means, city_vars) %>% pander()



#---------- negative binomial without covariates ----------#


### without pairs ###


# model for each city
negB <- list()
IRR <- vector()
lower.CI <- vector()
upper.CI <- vector()

for(i in 1:length(cities)){
  negB[[i]] <- glm.nb(death_sum3 ~ is_treated, data = matched %>% filter(city == cities[i]))
  
  # calculate incidence rate ratio for treated vs. untreated
  IRR[i] <- negB[[i]]$coefficients[2] %>% exp()
  
  # IRR confidence interval
  lower.CI[[i]] <- confint(negB[[i]])[2,1] %>% exp()
  upper.CI[[i]] <- confint(negB[[i]])[2,2] %>% exp()
  
}

regression_table <- data.frame(city = cities, IRR = IRR, lower.CI = lower.CI, upper.CI = upper.CI)


### with pairs ###

# model for each city
negB.pairs <- list()
IRR.pairs <- vector()
lower.CI.pairs <- vector()
upper.CI.pairs <- vector()

for(i in 1:length(cities)){
  negB.pairs[[i]] <- glm.nb(death_sum3 ~ is_treated + as.factor(pair), data = matched %>% filter(city == cities[i]))
  
  # calculate incidence rate ratio for treated vs. untreated
  IRR.pairs[i] <- negB.pairs[[i]]$coefficients[2] %>% exp()
  
  # IRR confidence interval
  lower.CI.pairs[[i]] <- confint(negB.pairs[[i]])[2,1] %>% exp()
  upper.CI.pairs[[i]] <- confint(negB.pairs[[i]])[2,2] %>% exp()
  
}


regression_table_pairs <- data.frame(city = cities, IRR.pairs = IRR.pairs, 
                                     lower.CI.pairs = lower.CI.pairs, upper.CI.pairs = upper.CI.pairs)


regression_table <- left_join(regression_table, regression_table_pairs, by = "city")




###  with covariates ###

# model for each city
negB.covs <- list() 
IRR.covs <- vector()
lower.CI.covs <- vector()
upper.CI.covs <- vector()

for(i in 1:length(cities)){
  negB.covs[[i]] <- glm.nb(death_sum3 ~ is_treated + as.factor(pair) + 
                             dow + month + o3_lag_3 + no2_lag_3, data = matched %>% filter(city == cities[i]))
  
  # calculate incidence rate ratio for treated vs. untreated
  IRR.covs[i] <- negB.covs[[i]]$coefficients[2] %>% exp()
  
  # IRR confidence interval
  lower.CI.covs[[i]] <- confint(negB.covs[[i]])[2,1] %>% exp()
  upper.CI.covs[[i]] <- confint(negB.covs[[i]])[2,2] %>% exp()
  
}


regression_table_covs <- data.frame(city = cities, IRR.covs = IRR.covs, 
                                     lower.CI.covs = lower.CI.covs, upper.CI.covs = upper.CI.covs)


regression_table <- left_join(regression_table, regression_table_covs, by = "city")



##############################################################################################
#-----------------------  FISHER P-VALUES  ----------------------#
##############################################################################################

#---------- Fisher p-values ----------#


# initialize tau
# tau_null <- list()
# 
# # loop through each city
# for(i in 1:length(cities)){
#   
#   # make each element of tau_null (a list) a vector to store simulated taus under the null
#   tau_null[[i]] <- vector()
# 
#   # for each city, loop through simulations
#   for(j in 1:10000){
#     
#     one_city <- matched %>% filter(city == cities[i])
#     
#     # randomly select one ID from each pair
#     treated_IDs <- one_city %>%
#       group_by(pair) %>%
#       slice_sample(n = 1) %>%
#       pull(id)
#     
#     # new column is_treated_sim to indicate treated vs. control in this simulation
#     one_city$is_treated_sim <- ifelse(one_city$id %in% treated_IDs, 1, 0)
#     
#     # calculate tau (Neyman inference)
#     tau_null[[i]][j] <- tapply(one_city$death_sum3, one_city$is_treated_sim, mean) %>%
#       diff()
#     
#   }
# }
# 
# # bind results into a table
# names(tau_null) <- cities
# tau_null_df <- bind_rows(tau_null) %>%
#   gather(key = "city", value = "tau")

#save(tau_null_df, file = "tau_null.Rda")
load(file = paste0(Fisher_path, "tau_null.Rda"))

# p <- list()
# for(i in 1:length(cities)){
#   p[[i]] <- tau_null_df %>%
#     filter(city == cities[i]) %>%
#     ggplot() +
#     geom_histogram(aes(tau), col = "white") +
#     labs(title = cities[i],
#          x = "Tau",
#          y = "Count") +
#     theme_minimal() +
#     geom_vline(xintercept = stats_table$tau[i], col = "red", size = 2)
# }
# 
# 
# grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], nrow = 5)


# calculate prob. of values under the null that were more extreme than the observed tau
# loop through 5 cities
p_val_tau <- vector()
for(i in 1:length(cities)){
  
  # filter to one city
  one_city <- tau_null_df %>% filter(city == cities[i])
  
  # observed tau
  tau_obs <- stats_table$tau[i]
  
  # calculate % of tau_null that is more extreme than tau_null
  p_val_tau[i] <- mean(abs(one_city$tau) >= abs(tau_obs))
  
}

tau_table <- data.frame(city = cities, tau = stats_table$tau, fisher.p.value.tau = p_val_tau)
tau_table %>% pander()


#-------------------- now try with regression betas --------------------#

# function to apply to each city


# IRR_null <- list()
# 
# # loop through each city
# for(i in 1:length(cities)){
# 
#   # make each element of IRR_null (a list) a vector to store simulated IRRs under the null
#   IRR_null[[i]] <- vector()
# 
#   # for each city, loop through simulations
#   for(j in 1:10000){
# 
#     one_city <- matched %>% filter(city == cities[i])
# 
#     # randomly select one ID from each pair
#     treated_IDs <- one_city %>%
#       group_by(pair) %>%
#       slice_sample(n = 1) %>%
#       pull(id)
# 
#     # new column is_treated_sim to indicate treated vs. control in this simulation
#     one_city$is_treated_sim <- ifelse(one_city$id %in% treated_IDs, 1, 0)
# 
#     # fit regression model
#     negB <- glm.nb(death_sum3 ~ is_treated_sim, data = one_city)
#     
#     # calculate incidence rate ratio for treated vs. untreated
#     IRR_null[[i]][j] <- negB$coefficients[2] %>% exp()
# 
#   }
# }
# 
# # bind results into a table
# names(IRR_null) <- cities
# IRR_null_df <- bind_rows(IRR_null) %>%
#   gather(key = "city", value = "IRR")
# 
# save(IRR_null_df, "IRR_null.Rda")
load(file = paste0(Fisher_path, "IRR_null.Rda"))



# plot IRRs
# p <- list()
# for(i in 1:length(cities)){
#   p[[i]] <- IRR_null_df %>%
#     filter(city == cities[i]) %>%
#     ggplot() +
#     geom_histogram(aes(IRR), col = "white") +
#     labs(title = cities[i],
#          x = "IRR",
#          y = "Count") +
#     theme_minimal() +
#     geom_vline(xintercept = regression_table$IRR[i], col = "red", size = 2)
# }
# 
# 
# grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], nrow = 5)


# calculate prob. of values under the null that were more extreme than the observed IRR
# loop through 5 cities
p_val_IRR <- vector()
for(i in 1:length(cities)){
  
  # filter to one city
  one_city <- IRR_null_df %>% filter(city == cities[i])
  
  # observed tau
  IRR_obs <- regression_table$IRR[i]
  
  # calculate % of tau_null that is more extreme than tau_null
  p_val_IRR[i] <- mean(abs(one_city$IRR) >= abs(IRR_obs))
  
}

# table of p-vals
IRR_table <- data.frame(city = cities, IRR = regression_table$IRR, fisher.p.value.IRR = p_val_IRR)

# manually change fisher p-value for NY
IRR_table$fisher.p.value.IRR[3] <- "< 0.0001"

IRR_table %>% pander()



# ---------------------------------------------------------------------------------------------------- #

# bind all results into one table
big_tau <- left_join(stats_table, tau_table, by = c("city", "tau")) %>%
  dplyr::select(-c(mean.control, mean.treated, var.control, var.treated, var.tau.unpaired)) #%>%
  #mutate_if(is.numeric, round, 4) 

# manually change fisher p-value for NY
big_tau$fisher.p.value.tau[3] <- "< 0.0001"

# hi <- big_tau %>%
#   unite_ci(col = "tau.unpaired", tau, lower.CI.tau.unpaired, upper.CI.tau.unpaired, m100 = FALSE) %>%
#   unite_ci(col = "tau.paired", tau, lower.CI.tau.paired, upper.CI.tau.paired, m100 = FALSE)

big_IRR <- left_join(regression_table, IRR_table, by = c("city", "IRR"))

# try making a separate table for each city


##############################################################################################
#-------------------------  Fisher/Fiducial 95% Confidence Intervals  -----------------------#
##############################################################################################

# fill in columns for known and missing outcomes
matched$Y.0 <- NA
matched$Y.1 <- NA

for(i in 1:nrow(matched)){

  # if control, assign deaths to Y.0
  if(!matched$is_treated[i]){
    matched$Y.0[i] <- matched$death_sum3[i]
  }

  # if treated, assign deaths to Y.1
  if(matched$is_treated[i]){
    matched$Y.1[i] <- matched$death_sum3[i]
  }
}

#----------------------- ATE -----------------------#

# # set a sequence of a values --> will test if each is in the 95% CI
# a_seq <- seq(-10, 40, by = 0.1)
# 
# # initialize
# tau_a <- list()
# tau_obs <- vector()
# p_val_tau_a <- vector()
# 
# # function to test each value of a
# get_tau_a <- function(a){
#   
#   # fill in all the "missing values" given value of a
#   matched$Y.0 <- ifelse(is.na(matched$Y.0), matched$Y.1 - a, matched$Y.0)
#   matched$Y.1 <- ifelse(is.na(matched$Y.1), matched$Y.0 + a, matched$Y.1)
#   
#   # loop through each city
#   for(i in 1:length(cities)){
#     
#     one_city <- matched %>% filter(city == cities[i])
#     
#     # make each element of tau_a a vector to store simulated taus under a
#     tau_a[[i]] <- vector()
#     
#     # loop through simulations
#     for(j in 1:1000){
#       
#       # randomly select one ID from each pair
#       treated_IDs <- one_city %>%
#         group_by(pair) %>%
#         slice_sample(n = 1) %>%
#         pull(id)
#       
#       # new column is_treated_sim to indicate treated vs. control in this simulation
#       one_city$is_treated_sim <- ifelse(one_city$id %in% treated_IDs, TRUE, FALSE)
#       
#       # new column death_sum3_sim with the number of deaths under the simulation
#       one_city$death_sum3_sim <- ifelse(one_city$is_treated_sim, one_city$Y.1, one_city$Y.0)
#       
#       # calculate tau under this value of a
#       tau_a[[i]][j] <- tapply(one_city$death_sum3_sim, one_city$is_treated_sim, mean) %>% diff()
#       
#     }
#     
#     # observed tau for city i
#     tau_obs <- stats_table$tau[i]
#     
#     # calculate % of tau_a that are more extreme than tau_obs
#     p_val_tau_a[i] <- mean(tau_a[[i]] >= tau_obs)
#     
#   }
#   
#   # return all p-values for the given value of a
#   return(p_val_tau_a)
#   
# }
# 
# # # apply function to the sequence of a values
# # a_table <- lapply(a_seq, get_tau_a)
# # 
# # # bind into a single data frame
# # a_table <- do.call(rbind.data.frame, a_table)
# # rownames(a_table) <- a_seq
# # colnames(a_table) <- cities

#save(a_table, file = "a_table.Rda")
load(file = paste0(Fisher_path, "a_table.Rda"))

# the p-values when a=0 matches the ones I got before! (in tau_table)

# chic: [3.5, 15.5]
# la: [1, 13]
# ny: [18, 35]
# pitt: [-2.5, 7]
# seat: [1, 6]

# add to stats table
stats_table$lower.CI.tau.fid <- c(3.5, 1, 18, -2.5, 1)
stats_table$upper.CI.tau.fid <- c(15.5, 13, 35, 7, 6)


#----------------------- RATE RATIO -----------------------#

# # set a sequence of a values --> will test if each is in the 95% CI
# a_seq_RR <- seq(0.94, 1.14, by = 0.002)
# 
# # initialize
# RR_a <- list()
# RR_obs <- vector()
# p_val_RR_a <- vector()
# 
# # function to test each value of a
# get_RR_a <- function(a){
#   
#   # fill in all the "missing values" given value of a
#   # round to integer because these are counts!
#   matched$Y.0 <- ifelse(is.na(matched$Y.0), matched$Y.1 * (1/a), matched$Y.0) %>% round()
#   matched$Y.1 <- ifelse(is.na(matched$Y.1), matched$Y.0 * a, matched$Y.1) %>% round()
#   
#   # loop through each city
#   for(i in 1:length(cities)){
#     
#     one_city <- matched %>% filter(city == cities[i])
#     
#     # make each element of RR_a a vector to store simulated RRs under a
#     RR_a[[i]] <- vector()
#     
#     # loop through simulations
#     for(j in 1:100){
#       
#       # randomly select one ID from each pair
#       treated_IDs <- one_city %>%
#         group_by(pair) %>%
#         slice_sample(n = 1) %>%
#         pull(id)
#       
#       # new column is_treated_sim to indicate treated vs. control in this simulation
#       one_city$is_treated_sim <- ifelse(one_city$id %in% treated_IDs, TRUE, FALSE)
#       
#       # new column death_sum3_sim with the number of deaths under the simulation
#       one_city$death_sum3_sim <- ifelse(one_city$is_treated_sim, one_city$Y.1, one_city$Y.0)
#       
#       # fit regression model
#       negB <- glm.nb(death_sum3_sim ~ is_treated_sim, data = one_city)
# 
#       # calculate incidence rate ratio for treated vs. control
#       RR_a[[i]][j] <- negB$coefficients[2] %>% exp()
#       
#     }
#     
#     # observed RR for city i
#     RR_obs <- regression_table$IRR[i]
#     
#     # calculate % of RR_a that are more extreme than RR_obs
#     p_val_RR_a[i] <- mean(RR_a[[i]] >= RR_obs)
#     
#   }
#   
#   # return all p-values for the given value of a
#   return(p_val_RR_a)
#   
# }
# 
# # # apply function to the sequence of a values
# a_table_RR <- lapply(a_seq_RR, get_RR_a)
# # 
# # # bind into a single data frame
# a_table_RR <- do.call(rbind.data.frame, a_table_RR)
# rownames(a_table_RR) <- a_seq_RR
# colnames(a_table_RR) <- cities

#save(a_table_RR, file = "a_table_RR.Rda")
load(file = paste0(Fisher_path, "a_table_RR.Rda"))


# chic: [1.008, 1.07]
# la: [1.002, 1.044]
# ny: [1.052, 1.102]
# pitt: [0.97, 1.076]
# seat: [1.002, 1.114]

# add to regression table
regression_table$lower.CI.fid <- c(1.008, 1.002, 1.052, 0.97, 1.002)
regression_table$upper.CI.fid <- c(1.07, 1.044, 1.102, 1.076, 1.114)



##############################################################################################
#---------------  UNMATCHED DATA--traditional adjusted analysis  -------------#
##############################################################################################

unmatched <- read.csv(paste0(processed_data_path, "all_processed.csv"))
unmatched <- unmatched %>%
  filter(!is.na(is_treated))

# new column for sum of deaths on 3 exp. days
unmatched$death_sum3 <- unmatched$death + unmatched$death_lag_1 + unmatched$death_lag_2

# histograms of deaths
# all cities
unmatched %>%
  ggplot() +
  geom_density(aes(death_sum3, fill = city), col = "white", alpha = 0.7) +
  ggtitle("Distribution of Death Counts") +
  xlab("Deaths over 3 Day Experiment") +
  theme_classic()
# split by city
unmatched %>%
  ggplot() +
  geom_histogram(aes(death_sum3), bins = 12, col = "white") +
  facet_wrap(~city, scales= "free") +
  xlab("Deaths over 3 Day Experiment") +
  theme_bw()


# which model?
city_means <- tapply(unmatched$death_sum3, unmatched$city, mean)
city_vars <- tapply(unmatched$death_sum3, unmatched$city, var)

data.frame(city_means, city_vars) %>% pander()


#---------- adjusted negative binomial model ----------#

# model for each city
negB.covs <- list()
IRR.covs.unmatch <- vector()
lower.CI.covs.unmatch <- vector()
upper.CI.covs.unmatch <- vector()

for(i in 1:length(cities)){
  negB.covs[[i]] <- glm.nb(death_sum3 ~ is_treated + 
                             dow + month + o3_lag_3 + no2_lag_3, data = unmatched %>% filter(city == cities[i]))
  
  # calculate incidence rate ratio for treated vs. untreated
  IRR.covs.unmatch[i] <- negB.covs[[i]]$coefficients[2] %>% exp()
  
  # IRR confidence interval
  lower.CI.covs.unmatch[[i]] <- confint(negB.covs[[i]])[2,1] %>% exp()
  upper.CI.covs.unmatch[[i]] <- confint(negB.covs[[i]])[2,2] %>% exp()
  
}


regression_table_unmatched <- data.frame(city = cities, IRR.covs.unmatch = IRR.covs.unmatch, 
                                    lower.CI.covs.unmatch = lower.CI.covs.unmatch, 
                                    upper.CI.covs.unmatch = upper.CI.covs.unmatch)

# join to regression table
regression_table <- left_join(regression_table, regression_table_unmatched, by = "city")



##############################################################################################
#-------------------------  TABLES  -----------------------#
##############################################################################################

# table with number of treated/matched units in each city

# table for unmatched data
n_tab <- unmatched %>%
  group_by(city, is_treated) %>%
  summarise(n = n()) %>%
  spread(key = "is_treated", value = "n") %>%
  rename(n.control = `FALSE`,
         n.treated = `TRUE`) %>%
  mutate(n.total = n.control + n.treated)

# add column for matched data
n_tab$n.pair <- stats_table$n

# reformat
n_tab <- n_tab[,c(1,4,2,3,5)]
n_tab <- n_tab %>% t() %>% as.data.frame()
names(n_tab) <- cities
n_tab <- n_tab[-1,]

kbl(n_tab) %>%
  kable_classic()

#----------------------------------------------#

# table for matching criteria

# matched on:
matched_on <- c("City", "Year", "Week of the year", 
                "Day of week category (Monday, Tuesday/Wednesday/Thursday, Friday, Saturday/Sunday)", 
                "NO2 day before experiment (ppb)", "O3 day before experiment (ppb)")
matched_by <- c("exact", "0", "3", "exact", "5", "5")

matched_tab <- data.frame(Covariate = matched_on,
                          `Discrepancy Allowed` = matched_by)

kbl(matched_tab) %>%
  kable_classic()

#----------------------------------------------#

# table for ATE

big_tau$city <- c("Chicago", "Los Angeles", "New York", "Pittsburgh", "Seattle")
big_tau$tau <- big_tau$tau %>% round(1)
big_tau$pairs <- big_tau$n / 2

big_tau %>%
  dplyr::select(city, pairs, tau, fisher.p.value.tau) %>%
  rename(` ` = city,
         `N Pairs` = pairs,
         ATE = tau,
         `Fisher p-value` = fisher.p.value.tau) %>% 
  kbl() %>%
  kable_paper(full_width = F, "striped")
  #kable_material(c("striped", full_width = F))

#----------------------------------------------#

# table for RR

IRR_table$city <- c("Chicago", "Los Angeles", "New York", "Pittsburgh", "Seattle")
IRR_table$IRR <- IRR_table$IRR %>% round(4)
IRR_table$pairs <- big_tau$pairs
IRR_table <- IRR_table[,c(1,4,2,3)]

IRR_table %>%
  dplyr::select(city, pairs, IRR, fisher.p.value.IRR) %>%
  rename(` ` = city,
         `N Pairs` = pairs,
         `Rate Ratio (Unadjusted)` = IRR,
         `Fisher p-value (Unadjusted)` = fisher.p.value.IRR) %>% 
  kbl() %>%
  kable_material(full_width = F)


##############################################################################################
#-------------------------  FIGURES  -----------------------#
##############################################################################################


# visualization for average causal effect
stats_table %>%
  ggplot(aes(x = tau, y = city)) +
  # unpaired error bars
  geom_errorbar(aes(xmin = lower.CI.tau.unpaired, xmax = upper.CI.tau.unpaired, col = "Unpaired Variance"), 
                size = 1, width = 0.15, position = position_nudge(y = 0.2)) +
  # paired error bars
  geom_errorbar(aes(xmin = lower.CI.tau.paired, xmax = upper.CI.tau.paired, col = "Paired Variance"), 
                size = 1, width = 0.15) +
  # Fiducial error bars
  geom_errorbar(aes(xmin = lower.CI.tau.fid, xmax = upper.CI.tau.fid, col = "Fisherian"), 
                size = 1, width = 0.15, position = position_nudge(y = -0.2)) +
  # Bayesian error bars
  # geom_errorbar(aes(xmin = lower.CI.tau.bayes, xmax = upper.CI.tau.bayes, col = "Bayesian"), 
  #               size = 1, width = 0.15, position = position_nudge(y = -0.3), lty = 2) +
  # points for tau
  geom_point(aes(x = tau, y = city, color = "ATE"), size = 3) +
  #geom_point(aes(x = tau, y = city), size = 3, position = position_nudge(y = 0.15)) +
  #geom_point(aes(x = tau, y = city), size = 3, position = position_nudge(y = -0.15)) +
  # points for Bayes
  # geom_point(aes(x = tau.bayes, y = city), size = 3, col = "#fd7f6f", position = position_nudge(y = -0.3)) +
  # labels
  #geom_text(aes(label = round(tau, 1)), nudge_y = .3) + 
  #geom_text(aes(x = tau.bayes, label = round(tau.bayes, 1)), nudge_y = -.45, col = "#fd7f6f") + 
  labs(title = "",
       x = "Average Treatment Effect",
       y = "",
       color = "Uncertainty Interval") +
  # manual legend
  
  # scale_color_manual(values = c("Unpaired" = "#52B2CF", "Paired" = "#084887", 
  #                               "Fiducial" = "#F58A07", "Neymanian Inference" = "black",
  #                               #"Bayesian" = "#fd7f6f",
  #                               guide = guide_legend(override.aes = list(
  #                                 linetype = c("solid", "solid", "blank"),
  #                                 shape = c(NA, NA, 1)))
  #                               )) +
  
  scale_color_manual(values = c("ATE" = "black",
                                "Unpaired Variance" = "#52B2CF", "Paired Variance" = "#084887", 
                                "Fisherian" = "#F58A07"#,
                                #"Bayesian" = "#fd7f6f",
                                )) +
  
  guides(color = guide_legend(override.aes = list(
                         linetype = c("blank", "solid", "solid", "solid"),
                         shape = c(16, NA, NA, NA)))) +
                         
  # reverse order of cities
  scale_y_discrete(labels = c("Seattle", "Pittsburgh", "New York", "Los Angeles", "Chicago"), limits = rev) +
  theme_minimal() +
  geom_vline(aes(xintercept = 0), lty = 2) +
  theme(axis.line.x = element_line(), 
        legend.position = c(.95, .95), 
        legend.justification = c("right", "top"))


# visualization for rate ratio
regression_table %>%
  ggplot(aes(y = city)) +
  # without pairs
  geom_errorbar(aes(xmin = lower.CI, xmax = upper.CI, col = "Unadjusted"), 
                size = 1, width = 0.15, position = position_nudge(y = 0.15)) +
  # Fiducial (for IRR without pairs)
  # geom_errorbar(aes(xmin = lower.CI.fid, xmax = upper.CI.fid, col = "Fiducial"), 
  #               size = 1, width = 0.2, position = position_nudge(y = 0.2), lty = 2) +
  # with pairs
  geom_errorbar(aes(xmin = lower.CI.pairs, xmax = upper.CI.pairs, col = "Adjusted by Pair"), 
                size = 1, width = 0.15) +
  # with pairs and covariates
  geom_errorbar(aes(xmin = lower.CI.covs, xmax = upper.CI.covs, col = "Adjusted"), 
                size = 1, width = 0.15, position = position_nudge(y = -0.15)) +
  # unmatched data (adjusted)
  geom_errorbar(aes(xmin = lower.CI.covs.unmatch, xmax = upper.CI.covs.unmatch, col = "Original Data (Adjusted)"), 
                size = 1, width = 0.15, position = position_nudge(y = -0.3)) +
  # points for tau
  geom_point(aes(x = IRR, y = city), size = 2, col = "#b881b1", position = position_nudge(y = 0.15)) +
  geom_point(aes(x = IRR.pairs, y = city), size = 2, col = "#e54a50") +
  geom_point(aes(x = IRR.covs, y = city), size = 2, col = "#ffa600", position = position_nudge(y = -0.15)) +
  geom_point(aes(x = IRR.covs.unmatch, y = city), size = 2, col = "#0d88e6", position = position_nudge(y = -0.3)) +
  #geom_point(aes(x = IRR, y = city), size = 2, col = "#b30000", position = position_nudge(y = 0.4)) +
  labs(title = "",
       x = "Rate Ratio",
       y = "",
       color = "Regression Model") +
  # manual legend
  scale_color_manual(values = c("Unadjusted" = "#b881b1", 
                                #"Fiducial" = "#b30000", 
                                "Adjusted by Pair" = "#e54a50", 
                                "Adjusted" = "#ffa600", 
                                "Original Data (Adjusted)" = "#0d88e6")) +
  # reverse order of cities
  scale_y_discrete(labels = c("Seattle", "Pittsburgh", "New York", "Los Angeles", "Chicago"), limits = rev) +
  # manually change legend shape
  guides(color = guide_legend(override.aes = list(linetype = c(1)))) +
  theme_minimal() +
  geom_vline(aes(xintercept = 1), lty = 2) +
  theme(axis.line.x = element_line(), 
        legend.position = c(1, 1), 
        legend.justification = c("right", "top"))

dev.off()

