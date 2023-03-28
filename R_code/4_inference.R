
###############################################################
#### Causal Inference for Time Series (Heat and Mortality) ####
###############################################################

### Script 4 ###
# Inputs: Processed data from script 1 and matched data from script 2
# Causal inference!

#------------------------------------------------------------------------#

library(ggplot2)
library(tidyverse)
library(MASS)
library(gridExtra)
library(kableExtra)
library(DataCombine)
library(readr)
library(AER)


#------------------------------------------------------------------------#
#--------- set working directory/paths ---------#
#------------------------------------------------------------------------#

# set path for processed data from script 1 and matched data from script 2
processed_data_path <- "data/processed/"

# set path for Fisher results
Fisher_path <- "data/Fisher_results/"

# set path for Fisher results

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


##############################################################################################
##############################################################################################

#---------- get population adjusted ATE ----------#

# get 1990 census data
# https://www.census.gov/data/tables/time-series/demo/popest/2000-subcounties-eval-estimates.html

# Chicago
# 2783485

# Los Angeles
# 3485567

# NY
# 7322564

# Pittsburgh
# 369962

# Seattle
# 516262

neyman_table$pop = c(2783485, 3485567, 7322564, 369962, 516262)
# ATE per 100,000 people
neyman_table$ate_adj <- (neyman_table$tau / neyman_table$pop) * 100000

##############################################################################################
##############################################################################################

#------------------------------ POISSON REGRESSION ------------------------------#

### without pairs ###

# model for each city
pois <- list()
IRR <- vector()
lower.CI <- vector()
upper.CI <- vector()

for(i in 1:length(cities)){
  pois[[i]] <- glm(death_sum3 ~ is_treated, data = matched %>% filter(city == cities[i]),
                   famil = poisson())
  
  # calculate incidence rate ratio for treated vs. untreated
  IRR[i] <- pois[[i]]$coefficients[2] %>% exp()
  
  # IRR confidence interval
  lower.CI[[i]] <- confint(pois[[i]])[2,1] %>% exp()
  upper.CI[[i]] <- confint(pois[[i]])[2,2] %>% exp()
  
}

regression_table <- data.frame(city = cities, IRR = IRR, lower.CI = lower.CI, upper.CI = upper.CI)


### with pairs ###

# model for each city
pois.pairs <- list()
IRR.pairs <- vector()
lower.CI.pairs <- vector()
upper.CI.pairs <- vector()

for(i in 1:length(cities)){
  pois.pairs[[i]] <- glm(death_sum3 ~ is_treated + as.factor(pair), 
                         data = matched %>% filter(city == cities[i]),
                         family = poisson())
  
  # calculate incidence rate ratio for treated vs. untreated
  IRR.pairs[i] <- pois.pairs[[i]]$coefficients[2] %>% exp()
  
  # IRR confidence interval
  lower.CI.pairs[[i]] <- confint(pois.pairs[[i]])[2,1] %>% exp()
  upper.CI.pairs[[i]] <- confint(pois.pairs[[i]])[2,2] %>% exp()
  
}

regression_table_pairs <- data.frame(city = cities, IRR.pairs = IRR.pairs, 
                                     lower.CI.pairs = lower.CI.pairs, upper.CI.pairs = upper.CI.pairs)

regression_table <- left_join(regression_table, regression_table_pairs, by = "city")


###  with covariates ###

# model for each city
pois.covs <- list() 
IRR.covs <- vector()
lower.CI.covs <- vector()
upper.CI.covs <- vector()

for(i in 1:length(cities)){
  pois.covs[[i]] <- glm(death_sum3 ~ is_treated + as.factor(pair) + 
                             month + o3_lag_3 + no2_lag_3, 
                        data = matched %>% filter(city == cities[i]),
                        family = poisson())
  
  # calculate incidence rate ratio for treated vs. untreated
  IRR.covs[i] <- pois.covs[[i]]$coefficients[2] %>% exp()
  
  # IRR confidence interval
  lower.CI.covs[[i]] <- confint(pois.covs[[i]])[2,1] %>% exp()
  upper.CI.covs[[i]] <- confint(pois.covs[[i]])[2,2] %>% exp()
  
}


regression_table_covs <- data.frame(city = cities, IRR.covs = IRR.covs, 
                                     lower.CI.covs = lower.CI.covs, upper.CI.covs = upper.CI.covs)


regression_table <- left_join(regression_table, regression_table_covs, by = "city")


#----- confirm Poisson fits well (not over-dispersed) -----#

for(i in 1:5){
  # test with AER package
  dispersiontest(pois[[i]])
}

for(i in 1:5){
  # test with AER package
  dispersiontest(pois.pairs[[i]])
}

for(i in 1:5){
  # test with AER package
  dispersiontest(pois.covs[[i]])
}

##############################################################################################
#-----------------------  BAYESIAN  ----------------------#
##############################################################################################

# updated 3/6

bayes_table <- data.frame(city = cities,
                          est.bayes = c(12.495369410363244, 8.569217644037433, 20.53897429451777, 
                                        0.4228827702335256, 3.6558323694724923),
                          var.bayes = c(3.096168660929421, 2.557833564322541, 5.79524887493058, 
                                        1.0471630267477539, 0.9556067978051604)^2)

bayes_table$lower.CI.bayes <- bayes_table$est.bayes - (1.96*sqrt(bayes_table$var.bayes))
bayes_table$upper.CI.bayes <- bayes_table$est.bayes + (1.96*sqrt(bayes_table$var.bayes))


##############################################################################################
#-----------------------  FISHERIAN / FIDUCIAL  ----------------------#
##############################################################################################

# read in results from Yasmine

library(readxl)
fisher_diff <- read_excel(paste0(Fisher_path, "Fisher_results.xlsx"), sheet = "ATE")
fisher_ratio <- read_excel(paste0(Fisher_path, "Fisher_results.xlsx"), sheet = "RR")


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

#---------- adjusted Poisson model ----------#

# model for each city
pois.covs.unmatch <- list()
IRR.covs.unmatch <- vector()
lower.CI.covs.unmatch <- vector()
upper.CI.covs.unmatch <- vector()

for(i in 1:length(cities)){
  pois.covs.unmatch[[i]] <- glm(death_sum3 ~ is_treated + 
                             dow + month + o3_lag_3 + no2_lag_3,
                           data = unmatched %>% filter(city == cities[i]),
                        family = poisson())
  
  # calculate incidence rate ratio for treated vs. untreated
  IRR.covs.unmatch[i] <- pois.covs.unmatch[[i]]$coefficients[2] %>% exp()
  
  # IRR confidence interval
  lower.CI.covs.unmatch[[i]] <- confint(pois.covs.unmatch[[i]])[2,1] %>% exp()
  upper.CI.covs.unmatch[[i]] <- confint(pois.covs.unmatch[[i]])[2,2] %>% exp()
  
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
  left_join(., fisher_diff, by = "city") %>%
  left_join(., fisher_ratio, by = "city")

write.csv(all_results, file = paste0(processed_data_path, "inference_results.csv"), row.names = FALSE)


##############################################################################################
#-------------------------  FIGURES  -----------------------#
##############################################################################################

# order by population
matched$city <- factor(matched$city, levels = c("pitt", "seat", "chic", "la", "ny"))
cities_ord <- c("pitt", "seat", "chic", "la", "ny")

city_names <- c(
  `chic` = "Chicago",
  `la` = "Los Angeles",
  `ny` = "New York",
  `pitt` = "Pittsburgh",
  `seat` = "Seattle"
)

# boxplots
pdf(file = paste0(fig_path, "mortality_boxplots.pdf"), width = 11, height = 5)

matched %>%
  ggplot(aes(x = is_treated, y = death_sum3, fill = is_treated, col = is_treated)) +
  geom_boxplot(aes(col = is_treated)) +
  #geom_jitter(position = position_jitter(0.2)) + 
  facet_wrap(~city, nrow = 1, 
             #scales = "free", 
             labeller = as_labeller(city_names)) +
  labs(x = "", 
       y = "Deaths per Three-Day Period",
       fill = "") +
  scale_fill_manual(name = "",
                    labels = c("Control (Warm)", "Treated (Hot)"),
                    values = c("orange", "tomato3")) +
  scale_color_manual(name = "",
                     labels = c("Control (Warm)", "Treated (Hot)"),
                    values = c("orange4", "firebrick4")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        #axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 18),
        legend.key.size = unit(1.3, 'cm')
        )

dev.off()


# theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#----- ATE -----#

# make a data frame for ATE plot
ate_df <- data.frame(city = rep(cities, 4),
                     ate.method = c(rep("Crude Mean", 5*3), rep("Posterior Mean", 5)),
                     ate.est = c(rep(all_results$tau, 3), 
                                 rep(all_results$est.bayes)),
                     ci.method = rep(c("Neyman (Unpaired Variance)", "Neyman (Paired Variance)", 
                                       "Fisherian", "Bayesian"), each = 5),
                     lower.ci = c(all_results$lower.CI.tau.unpaired, all_results$lower.CI.tau.paired, 
                                  all_results$fisher_ATE_low, all_results$lower.CI.bayes),
                     upper.ci = c(all_results$upper.CI.tau.unpaired, all_results$upper.CI.tau.paired, 
                                  all_results$fisher_ATE_high, all_results$upper.CI.bayes))

# order by population
ate_df$city <- factor(ate_df$city, levels = c("pitt", "seat", "chic", "la", "ny"))
ate_df$ci.method <- factor(ate_df$ci.method, levels = c("Neyman (Unpaired Variance)", 
                                                        "Neyman (Paired Variance)", 
                                                        "Fisherian", "Bayesian"))

# set nudge values
nudge_ate <- rep(c(0.24, 0.08, -0.08, -0.24), each = 5)

# plot ATE results
pdf(file = paste0(fig_path, "ATE_results.pdf"), width = 12, height = 6)

ate_df |>
  ggplot(aes(x = ate.est, y = city)) +
  geom_errorbar(aes(xmin = lower.ci, xmax = upper.ci, color = ci.method),
                width = 0.1, size = 1.3,
                position = position_nudge(y = nudge_ate)) +
  scale_color_manual(name = "Estimated 95% Uncertainty Interval",
                     values = c("Neyman (Unpaired Variance)" = "#52B2CF",
                                "Neyman (Paired Variance)" = "#084887",
                                "Fisherian" = "darkgoldenrod2",
                                "Bayesian" = "indianred3")) +
  geom_vline(aes(xintercept = 0), lty = 2) +
  geom_point(aes(x = ate.est, y = city, fill = ate.method),
             position = position_nudge(y = nudge_ate),
             size = 2.5, pch = 21, stroke = 1.3) +
  scale_fill_manual(name = "Estimated ATE",
                    values = c("Crude Mean" = "gray80",
                               "Posterior Mean" = "indianred3")) +
  scale_y_discrete(labels = c("New York", "Los Angeles", "Chicago", 
                              "Seattle", "Pittsburgh"),
                   limits = rev) +
  labs(title = "",
       x = "Estimated Average Treatment Effect",
       y = "") +
  theme_minimal() +
  theme(axis.line.x = element_line(),
        text = element_text(size = 14),
        axis.text.y = element_text(size = 16)) +
  guides(fill = guide_legend(order = 1), 
         color = guide_legend(order = 2))

dev.off()

#----- rate ratio -----# 

# make a data frame for RR plot
rr_df <- data.frame(city = rep(cities, 3),
                    rr.data = c(rep("Before Matching", 5), rep("After Matching", 5*2)),
                    rr.est = c(all_results$IRR.covs.unmatch, 
                               rep(all_results$IRR, 2)),
                    ci.method = rep(c("Poisson Regression", 
                                      "Poisson Regression", "Fisherian"), each = 5),
                    #estimand = c(rep("Pop. Level", 5*2), rep("Indiv. Level", 5)),
                    lower.ci = c(all_results$lower.CI.covs.unmatch, all_results$lower.CI, 
                                 all_results$fisher_RR_low),
                    upper.ci = c(all_results$upper.CI.covs.unmatch, all_results$upper.CI, 
                                 all_results$fisher_RR_high))

# order by population
rr_df$city <- factor(rr_df$city, levels = c("pitt", "seat", "chic", "la", "ny"))
rr_df$ci.method <- factor(rr_df$ci.method, levels = c("Poisson Regression", "Fisherian"))

# set nudge values
nudge_rr <- rep(c(0.2, 0, -0.2), each = 5)

pdf(file = paste0(fig_path, "rate_ratio_results.pdf"), width = 12, height = 6)

# plot RR results
rr_df |>
  ggplot(aes(x = rr.est, y = city)) +
  geom_errorbar(aes(xmin = lower.ci, xmax = upper.ci, color = ci.method),
                width = 0.1, size = 1.3,
                position = position_nudge(y = nudge_rr)) +
  scale_color_manual(name = "Estimated 95% Uncertainty Interval",
                     values = c("Poisson Regression" = "mediumslateblue",
                                "Fisherian" = "darkgoldenrod2")) +
  geom_vline(aes(xintercept = 1), lty = 2) +
  geom_point(aes(x = rr.est, y = city, fill = rr.data),
             position = position_nudge(y = nudge_rr),
             size = 2.5, pch = 21, stroke = 1.3) +
  # geom_point(aes(x = rr.est, y = city, fill = rr.data, shape = estimand),
  #            position = position_nudge(y = nudge_rr),
  #            size = 2.5, stroke = 1.3) +
  scale_fill_manual(name = "",
                    values = c("Before Matching" = "white",
                               "After Matching" = "black")) +
  scale_y_discrete(labels = c("New York", "Los Angeles", "Chicago", 
                              "Seattle", "Pittsburgh"), 
                   limits = rev) +
  # scale_shape_manual(name = "Estimand (?)",
  #                    values = c("Pop. Level" = 21,
  #                               "Indiv. Level" = 24)) +
  labs(title = "",
       x = "Estimated Rate Ratio",
       y = "") +
  theme_minimal() +
  theme(axis.line.x = element_line(),
        text = element_text(size = 14),
        axis.text.y = element_text(size = 16)) +
  guides(fill = guide_legend(order = 1), 
         color = guide_legend(order = 2))

dev.off()

