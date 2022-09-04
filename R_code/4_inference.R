
###############################################################
#### Causal Inference for Time Series (Heat and Mortality) ####
###############################################################

### Script 4 ###
# Inputs: Processed data from script 1 and matched data from script 2
# Causal inference!

#------------------------------------------------------------------------#

# install package for split violin plot
# see stack overflow: https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
# devtools::install_github("psyteachr/introdataviz")

library(ggplot2)
library(tidyverse)
library(MASS)
library(gridExtra)
library(kableExtra)
library(DataCombine)
library(readr)


#------------------------------------------------------------------------#
#--------- set working directory/paths ---------#
#------------------------------------------------------------------------#

# set path for processed data from script 1 and matched data from script 2
processed_data_path <- "data/processed/"

# set path for Fisher matrices
Fisher_data_path <- "data/Fisher_matrices/"

# set path for tables
table_path <- "figures/tables/"

# set path for figures with results
fig_path <- "figures/results/"

#------------------------------------------------------------------------#
#--------- load data ---------#
#------------------------------------------------------------------------#

# load original data
unmatched <- read.csv(paste0(processed_data_path, "all_processed.csv"))

# load matched data
matched <- read.csv(paste0(processed_data_path, "matched.csv"))

#------------------------------------------------------------------------#

# select only the columns I will need from the matched data
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
neyman_table <- matched %>% 
  group_by(city, is_treated) %>%
  summarise(mean = mean(death_sum3),
            n = n()*2) %>%
  spread(key = "is_treated", value = "mean") %>%
  rename(mean.control = `FALSE`,
         mean.treated = `TRUE`) %>%
  mutate(tau = mean.treated - mean.control)


# set N and K for each city
N <- neyman_table$n # total number of treated + control units
K <- neyman_table$n/2 # number of treated/control pairs

# --------------- ignoring matched pairs to get variance/CIs (block) ---------------#

# calculate sample variances by city and treatment group
variances <- matched %>%
  group_by(city, is_treated) %>%
  summarise(variance = var(death_sum3)) %>%
  spread(key = "is_treated", value = "variance") %>%
  rename(var.control = `FALSE`,
         var.treated = `TRUE`)

# join into one table
neyman_table <- left_join(neyman_table, variances, by = "city")

# add column for variance of tau (unpaired)
neyman_table$var.tau.unpaired <- (neyman_table$var.control + neyman_table$var.treated) / (neyman_table$n/2)

# confidence intervals
neyman_table$lower.CI.tau.unpaired <- neyman_table$tau - (1.96*sqrt(neyman_table$var.tau.unpaired))
neyman_table$upper.CI.tau.unpaired <- neyman_table$tau + (1.96*sqrt(neyman_table$var.tau.unpaired))


# --------------- using matched pairs to get variance/CIs (paired) ---------------#

# DR chapter 3, pages 37, 38, 39

# order by pair (with each FALSE day first)
matched <- matched %>%
  arrange(pair, is_treated)

# initialize columns for variance of tau and 95% CIs (using pairs)
neyman_table$lower.CI.tau.paired <- NA
neyman_table$upper.CI.tau.paired <- NA

# loop through 5 cities
for(i in 1:length(cities)){
  
  # filter to one city
  one_city <- matched %>% filter(city == cities[i])
  
  # vector of differences in deaths on control vs. treated days for each pair 
  # (keep only odd diffs between rows)
  D_k <- diff(one_city$death_sum3, lag = 1)[seq_len(neyman_table$n[i]) %% 2 == 1]
  
  # mean difference between control vs. treated (same as tau in the stats_table)
  D_bar <- mean(D_k)
  
  # (s2_D)/k --> unbiased estimate of var(tau)
  var.tau <- sum((D_k-D_bar)^2) / (K[i]*(K[i]-1))
  
  # paired t-statistic
  t_k_0.025 <- qt(0.025, df = nrow(one_city)/2) %>% abs()
  
  # 95% confidence intervals
  neyman_table$lower.CI.tau.paired[i] <- D_bar - t_k_0.025*sqrt(var.tau)
  neyman_table$upper.CI.tau.paired[i] <- D_bar + t_k_0.025*sqrt(var.tau)
  
}


#------------------------------ REGRESSION ------------------------------#

# Poisson or negative binomial? check means and variances for each city
data.frame(city_means = tapply(matched$death_sum3, matched$city, mean), 
           city_vars = tapply(matched$death_sum3, matched$city, var))

# overdispersion in every city except Seattle --> use the negative binomial model
# produces warning because Seattle is underdispersed

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
                             month + o3_lag_3 + no2_lag_3, data = matched %>% filter(city == cities[i]))
  
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
#-----------------------  BAYESIAN  ----------------------#
##############################################################################################

# will need to update these!!! read in from csv from Young

bayes_table <- data.frame(city = cities,
                          est.bayes = c(9.36, 7.72, 25.06, 1.91, 3.88),
                          var.bayes = c(2.91, 2.62, 5.60, 1.02, 1.08)^2)

bayes_table$lower.CI.bayes <- bayes_table$est.bayes - (1.96*sqrt(bayes_table$var.bayes))
bayes_table$upper.CI.bayes <- bayes_table$est.bayes + (1.96*sqrt(bayes_table$var.bayes))


##############################################################################################
#-----------------------  FISHERIAN / FIDUCIAL  ----------------------#
##############################################################################################

# fidiciual_matrix<-function(rand_mat,a_sequence,city_index,pairs){
#   
#   a_seq <- a_sequence ## values to test for interval
#   
#   ## initialize
#   tau_a <- list()
#   tau_obs <- vector()
#   p_val_tau_a <- vector()
#   
#   ## function to test each value of a
#   inc<-1
#   for(a in a_seq){
#     
#     one_city <- matched %>% filter(city == cities[city_index])
#     #     
#     #     # make each element of tau_a a vector to store simulated taus under a
#     tau_a <- vector()
#     #     
#     #     # loop through simulations
#     for(j in 1:100000){
#       
#       tau_a[j]<- ( ( (1/pairs) * (sum(one_city$death_sum3[rand_mat[,j]==1 & one_city$is_treated == 1]) +
#                                     sum((one_city$death_sum3[rand_mat[,j]==1 & one_city$is_treated == 0]) + a))) -
#                      
#                      ((1/pairs) * (sum(one_city$death_sum3[rand_mat[,j]==0 & one_city$is_treated == 0])
#                                    + sum( ((one_city$death_sum3[rand_mat[,j]==0 & one_city$is_treated == 1]) - a) ))) )
#     }
#   
#     #     # observed tau for city i
#     tau_obs <- neyman_table$tau[city_index]
#     #     
#     #     # calculate % of tau_a that are more extreme than tau_obs
#     p_val_tau_a[inc] <- signif(mean((tau_a) >= (tau_obs)),4)
#     if(p_val_tau_a[inc] == .025){
#       print(a)
#     }
#     inc<-inc+1
#     
#   }
#   
#   ret<-data.frame(p_val_tau_a)
#   row.names(ret)<-a_seq
#   
#   return(invisible(ret))
# }
# 
# 
# ### Functions taken from updated (not yet updated on CRAN) Counternull package https://github.com/ymabene/Counternull
# 
# find_test_stat_diff_means<-function(sample_data,variable){
#   # mean for experimental group (exposed)
#   on_mean <-mean((variable)[sample_data[,1]=="1"])
#   # mean for control group (non exposure)
#   off_mean <-mean((variable)[sample_data[,1]=="0"])
#   # difference
#   test_stat<-on_mean - off_mean
#   return(invisible(test_stat))
# }
# 
# 
# permutation_null_diff_means<-function(rand_matrix,variable,iterations){
#   # permutation vector with differences of means
#   perm_samples<-matrix(ncol=1,nrow=iterations)
#   # creates distribution
#   for(k in 1:iterations)
#   {
#     on<-mean(variable[rand_matrix[,k]==1]) # exposed
#     off<-mean(variable[rand_matrix[,k]==0]) # not exposed
#     perm_samples[k]<-on-off
#   }
#   return(invisible(perm_samples))
# }
# 
# 
# create_null_distribution<-function(sample_data, extreme, rand_matrix,
#                                    permutation_null_function,test_stat,
#                                    variable,iterations){
#   # Creates permutation vector
#   perm_samples<-permutation_null_function(rand_matrix,variable,iterations)
#   
#   # creates histogram and prints p-value
#   null_hist<-hist(perm_samples,breaks=100,col = "gold",
#                   main=paste("Null Distribution"), xlab="Test Statistics")
#   abline(v=test_stat,col="black",lty=2, lwd=5)
#   if (extreme==0){ # smaller test statistics are more extreme
#     pvalue<-sum(perm_samples<=(test_stat))/iterations
#     
#   } else { # larger test statistics are more extreme
#     pvalue<-sum(perm_samples>=(test_stat))/iterations
#     
#   }
#   print(paste("Test Statistic =",test_stat))
#   print(paste("Pvalue =",pvalue))
#   return(invisible(perm_samples))
# }
# 
# #----------------------- Fisher exact p-values -----------------------#
# 
# # read in matrices from Yasmine
# load(paste0(Fisher_data_path, "rand_mat_c.rda"))
# load(paste0(Fisher_data_path, "rand_mat_la.rda"))
# load(paste0(Fisher_data_path, "rand_mat_ny.rda"))
# load(paste0(Fisher_data_path, "rand_mat_p.rda"))
# load(paste0(Fisher_data_path, "rand_mat_s.rda"))
# 
# 
# # Chicago
# c <- matched %>% filter(city == cities[1])
# c<-c[,c(25,1:24,26:28)] # reorder columns
# c$is_treated<-as.numeric(c$is_treated) # make is_treated numeric
# dm_c<-find_test_stat_diff_means(c,c$death_sum3) # ATE: 12.4737
# perm_c<-create_null_distribution(c,1,rand_mat_c,permutation_null_diff_means,dm_c,c$death_sum3,100000)
# # p-value: .00018
# 
# 
# # LA
# la <- matched %>% filter(city == cities[2])
# la<-la[,c(25,1:24,26:28)] # reorder columns
# la$is_treated<-as.numeric(la$is_treated)
# dm_la<-find_test_stat_diff_means(la,la$death_sum3) #ATE: 8.458
# perm_la<-create_null_distribution(la,1,rand_mat_la,permutation_null_diff_means,dm_la,la$death_sum3,100000)
# # p-value: .00764
# 
# # NY
# ny <- matched %>% filter(city == cities[3])
# ny<-ny[,c(25,1:24,26:28)] # reorder columns
# ny$is_treated<-as.numeric(ny$is_treated)
# dm_ny<-find_test_stat_diff_means(ny,ny$death_sum3) #ATE: 20.313
# perm_ny<-create_null_distribution(ny,1,rand_mat_ny,permutation_null_diff_means,dm_ny,ny$death_sum3,100000)
# # p-value: .00012
# 
# 
# # Pittsburgh
# p <- matched %>% filter(city == cities[4])
# p<-p[,c(25,1:24,26:28)] # reorder columns
# p$is_treated<-as.numeric(p$is_treated)
# dm_p<-find_test_stat_diff_means(p,p$death_sum3) #ATE: .463
# perm_p<-create_null_distribution(p,1,rand_mat_p,permutation_null_diff_means,dm_p,p$death_sum3,100000)
# #p-value= .4192
# 
# 
# # Seattle
# s <- matched %>% filter(city == cities[5])
# s<-s[,c(25,1:24,26:28)] # reorder columns
# s$is_treated<-as.numeric(s$is_treated)
# dm_s<-find_test_stat_diff_means(s,s$death_sum3) #ATE: 3.727
# perm_s<-create_null_distribution(s,1,rand_mat_s,permutation_null_diff_means,dm_s,s$death_sum3,100000)
# # p-value: .0072
# 
# 
# #----------------------- adjusted p-values -----------------------#
# 
# all_p_values<-function(data){ # obtain all p-values or each city and iterations
#   
#   all_p_values<-data.frame(matrix(NA,nrow=100000,ncol=5))
#   
#   for(i in 1:5){
#     
#     for(j in 1:100000){
#       
#       all_p_values[j,i]<-sum(data[,i]>=(data[j,i]))/100000
#       
#     }
#   }
#   return(all_p_values)
# }
# 
# 
# min_p<-function(p_values){ # obtain minimum p-value for each iteration
#   
#   minimum<-vector()
#   
#   for(i in 1:100000){
#     
#     minimum[i]<-min(p_values[i,])
#   }
#   
#   return(minimum)
#   
# }
# 
# all_test_stat <- data.frame("chi"=perm_c,"la"=perm_la,"ny"=perm_ny,"pitt"=perm_p,"seat"=perm_s)
# #p_values <- all_p_values(all_test_stat)
# #write.csv(p_values, file = "data/all_p_values.csv")
# p_values <- read.csv(file = "data/all_p_values.csv")
# 
# minimum_p_values<-min_p(p_values)
# 
# adjusted_p<-function(minimum, obs){
#   
#   return(sum(minimum<=obs)/100000)
#   
# }
# 
# 
# adjusted_p(minimum_p_values,.0018) # chicago
# # .00866
# 
# adjusted_p(minimum_p_values,.00764) # la
# # .03689
# 
# adjusted_p(minimum_p_values,.0012) # ny
# # .00585
# 
# adjusted_p(minimum_p_values,.41945) # pitt
# # .93279
# 
# adjusted_p(minimum_p_values,.00717) # seat
# # .03514S
# 
# 
# #----------------------- Fiducial intervals -----------------------#
# 
# # Chicago 
# 
# # [6.05,18.88]
# # tau: 12.47
# # p-value: .00018
# 
# c_interval<-fidiciual_matrix(rand_mat_c,seq(6.0,6.1,.01),1,38) 
# c_interval<-fidiciual_matrix(rand_mat_c,seq(18.85,18.9,.01),1,38) 
# View(c_interval)
# 
# 
# # LA 
# 
# # [1.65,15.22]
# # tau: 8.46
# # p-value: .00764
# la_interval<-fidiciual_matrix(rand_mat_la,seq(1.6,1.7,.01),2,59) 
# la_interval<-fidiciual_matrix(rand_mat_la,seq(15.2,15.3,.01),2,59)
# View(la_interval)
# 
# 
# # New York 
# 
# # [10.24,30.41]
# # tau: 20.3125
# # p-value: .00012
# ny_interval<-fidiciual_matrix(rand_mat_ny,seq(10.2,10.3,.01),3,32) 
# ny_interval<-fidiciual_matrix(rand_mat_ny,seq(30.4,30.5,.01),3,32)
# View(ny_interval)
# 
# 
# # Pittsburgh 
# 
# # [-3.77,4.68]
# # tau: .463
# # p-value: .419
# # counternull value: .75
# 
# p_interval<-fidiciual_matrix(rand_mat_p,seq(-3.8,-3.7,.01),4,41)
# p_interval<-fidiciual_matrix(rand_mat_p,seq(4.6,4.7,.01),4,41)
# View(p_interval)
# 
# 
# # Seattle
# 
# # [.84,6.63]
# # tau: 3.727
# # p-value: .00717
# s_interval<-fidiciual_matrix(rand_mat_s,seq(.8,.9,.01),5,33)
# s_interval<-fidiciual_matrix(rand_mat_s,seq(6.6,6.7,.01),5,33)
# View(s_interval)

#------------------------------------------------------------#

# read in results from Yasmine

library(readxl)
fisher_table <- read_excel("data/results.xlsx")


##############################################################################################
#---------------  UNMATCHED DATA--traditional adjusted analysis  -------------#
##############################################################################################

unmatched <- read.csv(paste0(processed_data_path, "all_processed.csv"))

# count number of observations per city
total_obs <- unmatched %>%
  group_by(city) %>%
  summarise(n())

unmatched <- unmatched %>%
  filter(!is.na(is_treated))

# new column for sum of deaths on 3 exp. days
unmatched$death_sum3 <- unmatched$death + unmatched$death_lag_1 + unmatched$death_lag_2

# Poisson or negative binomial? check means and variances for each city
data.frame(city_means = tapply(unmatched$death_sum3, unmatched$city, mean), 
           city_vars = tapply(unmatched$death_sum3, unmatched$city, var))

# overdispersion in every city except Seattle --> use the negative binomial model
# produces warning because Seattle is underdispersed

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

# get all results in one table

all_results <- left_join(neyman_table, regression_table, by = "city") %>%
  left_join(., bayes_table, by = "city") %>%
  left_join(., fisher_table, by = "city")


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
  as.data.frame()

# add column for total number of days with data for each city
n_tab$total <- total_obs$`n()`

# add column for matched data
n_tab$n.pair <- neyman_table$n/2

# reformat
n_tab <- n_tab[,c(1,4,2,3,5)]
#row.names(n_tab) <- c(n_tab[,1])
n_tab <- n_tab %>%
  `row.names<-`(c("Chicago", "Los Angeles", "New York", "Pittsburgh", "Seattle")) %>%
  dplyr::select(-c(city))

# add row for % of treated units ("heat waves") that were had a control match
n_tab$pct_matched <- (100*(n_tab$n.pair/n_tab$n.treated)) %>% round(1)

# kable
n_tab %>%
  setNames(c("N total", "N control", 
             "N treated", "N pair", 
             "% Treated Units Matched")) %>%
  kable(format = "html") %>%
  kable_styling(bootstrap_options = c("striped")) # %>%
  # save_kable(file = paste0(table_path, "matching_counts.png"))


# make the column names subscript! (may not be possible in R script)

#----------------------------------------------#

# table for matching criteria

# matched on:
matched_on <- c("City", "Year", "Week of the year",
                "Day of the week",
                "NO2 day before experiment (ppb)", "O3 day before experiment (ppb)")
matched_by <- c("exact", "0", "3", "exact", "5", "5")

matched_tab <- data.frame(cov = matched_on,
                          disc = matched_by)

matched_tab %>%
  setNames(c("Covariate", "Discrepancy Allowed")) %>%
  kable(format = "html") %>%
  kable_styling(bootstrap_options = c("striped"))


#----------------------------------------------#

# NEED TO MAKE NEW TABLES WITH ALL RESULTS (for appendix)



# table for ATE

# ATE_tab <- all_results
# 
# ATE_tab$city <- c("Chicago", "Los Angeles", "New York", "Pittsburgh", "Seattle")
# ATE_tab$tau <- ATE_tab$tau %>% round(1)
# ATE_tab$pairs <- ATE_tab$n / 2
# 
# ATE_tab %>%
#   #dplyr::select(city, pairs, tau, fisher.p.value.tau) %>%
#   rename(` ` = city,
#          `N Pairs` = pairs,
#          A.T.E. = tau,
#          `Fisher adjusted p-value` = `adjusted p_value`) %>% 
#   kbl() %>%
#   kable_paper(full_width = F, "striped")
#   #kable_material(c("striped", full_width = F))

# table for RR

# IRR_table$city <- c("Chicago", "Los Angeles", "New York", "Pittsburgh", "Seattle")
# IRR_table$IRR <- IRR_table$IRR %>% round(4)
# IRR_table$pairs <- big_tau$pairs
# IRR_table <- IRR_table[,c(1,4,2,3)]
# 
# IRR_table %>%
#   dplyr::select(city, pairs, IRR, fisher.p.value.IRR) %>%
#   rename(` ` = city,
#          `N Pairs` = pairs,
#          `Rate Ratio (Unadjusted)` = IRR,
#          `Fisher p-value (Unadjusted)` = fisher.p.value.IRR) %>% 
#   kbl() %>%
#   kable_material(full_width = F)

##############################################################################################
#-------------------------  FIGURES  -----------------------#
##############################################################################################

#----- Split Violin Plot -----#

city_names <- c(
  `chic` = "Chicago",
  `la` = "Los Angeles",
  `ny` = "New York",
  `pitt` = "Pittsburgh",
  `seat` = "Seattle"
)

# library(introdataviz)
# ggplot(matched, aes(x = 1, y = death_sum3, fill = is_treated)) +
#   geom_split_violin(col = "white") +
#   facet_wrap(~city, nrow = 1, scales = "free", labeller = as_labeller(city_names)) +
#   labs(x = "Density", 
#        y = "Deaths per 3 Day Period",
#        fill = "") +
#   scale_fill_manual(labels = c("Control", "Treated"),
#                     values = c("orange", "firebrick3")) +
#   theme_classic() +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         strip.background = element_blank(),
#         axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))

# try boxplots instead
png(file = paste0(fig_path, "mortality_boxplots.png"), width = 700, height = 300)

matched %>%
  ggplot(aes(y = death_sum3, fill = is_treated)) +
  geom_boxplot() +
  facet_wrap(~city, nrow = 1, scales = "free", labeller = as_labeller(city_names)) +
  labs(x = "", 
       y = "Deaths per 3 Day Period",
       fill = "") +
  scale_fill_manual(labels = c("Control", "Treated"),
                    values = c("orange", "firebrick3")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))

dev.off()


#----- ATE -----#

png(file = paste0(fig_path, "ATE_results.png"), width = 400, height = 400)

all_results %>%
  ggplot(aes(x = tau, y = city)) +
  # unpaired error bars
  geom_errorbar(aes(xmin = lower.CI.tau.unpaired, xmax = upper.CI.tau.unpaired, col = "Unpaired Variance"), 
                size = 1, width = 0.15, position = position_nudge(y = 0.15)) +
  # paired error bars
  geom_errorbar(aes(xmin = lower.CI.tau.paired, xmax = upper.CI.tau.paired, col = "Paired Variance"), 
                size = 1, width = 0.15) +
  # Fiducial error bars
  geom_errorbar(aes(xmin = fisher_low, xmax = fisher_high, col = "Fisherian"), 
                size = 1, width = 0.15, position = position_nudge(y = -0.15)) +
  # Bayesian error bars
  geom_errorbar(aes(xmin = lower.CI.bayes, xmax = upper.CI.bayes, col = "Bayesian"),
                size = 1, width = 0.15, position = position_nudge(y = -0.3)) +
  # points for tau
  geom_point(aes(x = tau, y = city, color = "ATE"), size = 3) +
  geom_point(aes(x = tau, y = city), size = 3, position = position_nudge(y = 0.15)) +
  geom_point(aes(x = tau, y = city), size = 3, position = position_nudge(y = -0.15)) +
  # points for Bayes
  geom_point(aes(x = est.bayes, y = city), size = 3, col = "gold", position = position_nudge(y = -0.3)) +
  # labels
  #geom_text(aes(label = round(tau, 1)), nudge_y = .3) + 
  #geom_text(aes(x = tau.bayes, label = round(tau.bayes, 1)), nudge_y = -.45, col = "#fd7f6f") + 
  labs(title = "",
       x = "Estimated Average Treatment Effect",
       y = "",
       color = "Uncertainty Interval") +
  # manual legend
  scale_color_manual(values = c("ATE" = "black",
                                "Unpaired Variance" = "#52B2CF", 
                                "Paired Variance" = "#084887", 
                                "Fisherian" = "#F58A07",
                                "Bayesian" = "gold"
                                )) +
  guides(color = guide_legend(override.aes = list(
                         linetype = c("blank", "solid", "solid", "solid", "solid"),
                         shape = c(16, NA, NA, NA, NA)))) +
  # reverse order of cities
  scale_y_discrete(labels = c("Seattle", "Pittsburgh", "New York", "Los Angeles", "Chicago"), limits = rev) +
  theme_minimal() +
  geom_vline(aes(xintercept = 0), lty = 2) +
  theme(axis.line.x = element_line(), 
        legend.position = c(.95, .95), 
        legend.justification = c("right", "top"))

dev.off()


#----- rate ratio -----# 

png(file = paste0(fig_path, "rate_ratio_results.png"), width = 500, height = 400)


all_results %>%
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

