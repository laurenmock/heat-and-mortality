
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
                          est.bayes = c(12.21, 9.93, 21.6, 0.19, 3.13),
                          var.bayes = c(3.52, 2.79, 6.98, 1.17, 1.11)^2)

bayes_table$lower.CI.bayes <- bayes_table$est.bayes - (1.96*sqrt(bayes_table$var.bayes))
bayes_table$upper.CI.bayes <- bayes_table$est.bayes + (1.96*sqrt(bayes_table$var.bayes))


##############################################################################################
#-----------------------  FISHERIAN / FIDUCIAL  ----------------------#
##############################################################################################

# read in results from Yasmine

library(readxl)
fisher_diff <- read_excel("data/Fisher_results/Fisher_results.xlsx", sheet = "ATE")
fisher_ratio <- read_excel("data/Fisher_results/Fisher_results.xlsx", sheet = "RR")


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
                             dow + month + o3_lag_3 + no2_lag_3,
                           data = unmatched %>% filter(city == cities[i]))
  
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
  left_join(., fisher_diff, by = "city") %>%
  left_join(., fisher_ratio, by = "city")

write.csv(all_results, file = paste0(processed_data_path, "inference_results.csv"), row.names = FALSE)


##############################################################################################
#-------------------------  FIGURES  -----------------------#
##############################################################################################


city_names <- c(
  `chic` = "Chicago",
  `la` = "Los Angeles",
  `ny` = "New York",
  `pitt` = "Pittsburgh",
  `seat` = "Seattle"
)

# boxplots
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
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        )

dev.off()


theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

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
  geom_errorbar(aes(xmin = fisher_ATE_low, xmax = fisher_ATE_high, col = "Fisherian"), 
                size = 1, width = 0.15, position = position_nudge(y = -0.15)) +
  # Bayesian error bars
  geom_errorbar(aes(xmin = lower.CI.bayes, xmax = upper.CI.bayes, col = "Bayesian"),
                size = 1, width = 0.15, position = position_nudge(y = -0.3)) +
  # points for tau
  geom_point(aes(x = tau, y = city, color = "ATE"), size = 2) +
  geom_point(aes(x = tau, y = city), size = 2, position = position_nudge(y = 0.15)) +
  geom_point(aes(x = tau, y = city), size = 2, position = position_nudge(y = -0.15)) +
  # points for Bayes
  geom_point(aes(x = est.bayes, y = city), size = 2, col = "gold", position = position_nudge(y = -0.3)) +
  # labels
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
                size = 1, width = 0.15, position = position_nudge(y = 0.3)) +
  # Fiducial (for IRR without pairs)
  geom_errorbar(aes(xmin = fisher_RR_low, xmax = fisher_RR_high, col = "Fiducial"),
                size = 1, width = 0.15, position = position_nudge(y = 0.15)) +
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
  geom_point(aes(x = IRR, y = city), size = 2, col = "#b881b1", position = position_nudge(y = 0.3)) +
  geom_point(aes(x = IRR, y = city), size = 2, col = "#b881b1", position = position_nudge(y = 0.15)) +
  geom_point(aes(x = IRR.pairs, y = city), size = 2, col = "#e54a50") +
  geom_point(aes(x = IRR.covs, y = city), size = 2, col = "#ffa600", position = position_nudge(y = -0.15)) +
  geom_point(aes(x = IRR.covs.unmatch, y = city), size = 2, col = "#0d88e6", position = position_nudge(y = -0.3)) +
  labs(title = "",
       x = "Rate Ratio",
       y = "",
       color = "Regression Model") +
  # manual legend
  scale_color_manual(values = c("Unadjusted" = "#b881b1", 
                                "Fiducial" = "#b30000", 
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

