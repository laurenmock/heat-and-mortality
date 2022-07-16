
###############################################################
#### Causal Inference for Time Series (Heat and Mortality) ####
###############################################################

### Script 3 ###
# Inputs: Processed data from script 1 and matched data from script 2
# Outputs: Matching analysis, PDF with all figures


#------------------------------------------------------------------------#

library(ggplot2)
library(tidyverse)
library(pander)
library(gridExtra)
library(ggpubr)


#------------------------------------------------------------------------#
#--------- set working directory/paths ---------#
#------------------------------------------------------------------------#

# set path for processed data from script 1 and matched data from script 2
processed_data_path <- "data/processed/"

# set path to write out PDF with figures (use "" for general path)
figures_path <- "figures/"

#------------------------------------------------------------------------#
#--------- load data ---------#
#------------------------------------------------------------------------#

# open PDF (will automatically write in figures after running the whole script)
pdf(file = paste0(figures_path, "balance_diagnostics.pdf"))


# load treated data (before matching)
dat_before <- read.csv(paste0(processed_data_path, "all_processed.csv"))

#load matched data (after matching)
dat_after <- read.csv(paste0(processed_data_path, "matched.csv"))

# create column in matched data to indicate that these rows were in the matched subset (after join below)
dat_after$matched <- 1


#------------------------------------------------------------------------------------------------#

# remove missing values for treatment (but save NA treatment for later!)
dat_before_wNA <- dat_before

dat_before <- dat_before %>%
  filter(!is.na(is_treated))

#------------------------------------------------------------------------------------------------#

# create one data frame with before/after matching data together
dat_all <- left_join(dat_before, dat_after)

# matching column indicates 0 for rows that were not matched, 1 for rows that were matched
dat_all$matched <- ifelse(is.na(dat_all$matched), 0, 1)


#------------------------------------------------------------------------------------------------#
#------------------- Love Plots (code from Alice!) -------------------#
#------------------------------------------------------------------------------------------------#

# create a standardized difference in means function (stdif)
stdif <- function(X, W){
  mean_diff = mean(X[W == TRUE], na.rm = TRUE) - mean(X[W == FALSE], na.rm = TRUE)
  SE = sqrt(0.5 * (var(X[W == TRUE], na.rm = TRUE) + var(X[W == FALSE], na.rm = TRUE)))
  return(mean_diff/SE)
}

# covariates of interest for love plot
# covs <- c("no2_lag_3", "no2_lag_4", "no2_lag_5",
#           "o3_lag_3", "o3_lag_4", "o3_lag_5", 
#           "tmax_lag_3", "tmax_lag_4", "tmax_lag_5",
#           "year", "month_Sep", "month_Aug", "month_Jul", "month_Jun",
#           "dow_Sunday", "dow_Saturday", "dow_Friday",
#           "dow_Thursday", "dow_Wednesday", "dow_Tuesday", "dow_Monday")

covs <- c("no2_lag_3", "no2_lag_4", "no2_lag_5",
          "o3_lag_3", "o3_lag_4", "o3_lag_5", 
          "tmax_lag_3", "tmax_lag_4", "tmax_lag_5",
          "month_Sep", "month_Aug", "month_Jul", "month_Jun", "year")

# covariate labels
# cov_labs <- c(expression('NO'[2]*' Lag 1'), expression('NO'[2]*' Lag 2'), expression('NO'[2]*' Lag 3'),
#               expression('O'[2]*' Lag 1'), expression('O'[2]*' Lag 2'), expression('O'[2]*' Lag 3'),
#               "Temp. Lag 1", "Temp. Lag 2", "Temp. Lag 3",
#               "Year", "September", "August", "July", "June",
#               "Sunday", "Saturday", "Friday",
#               "Thursday", "Wednesday", "Tuesday", "Monday")

cov_labs <- c(expression('NO'[2]*' Lag 1'), expression('NO'[2]*' Lag 2'), expression('NO'[2]*' Lag 3'),
              expression('O'[3]*' Lag 1'), expression('O'[3]*' Lag 2'), expression('O'[3]*' Lag 3'),
              "Temp. Lag 1", "Temp. Lag 2", "Temp. Lag 3",
              "September", "August", "July", "June", "Year")

#-----------------------------------------#

# standardized diff. before matching
stdif_before <- sapply(dat_before %>% dplyr::select(covs),
                           function(x) stdif(x, dat_before %>% pull(is_treated)))

# standardized diff. after matching
stdif_after <- sapply(dat_after %>% dplyr::select(covs),
                       function(x) stdif(x, dat_after %>% pull(is_treated)))


# love plot
# set margins so we can see the whole plot
par(mar=c(6,6,4,4))
plot(stdif_before, 1:length(covs), 
     axes = FALSE, 
     main = "Love Plot",
     ylab = "",
     xlab = "Standardized Difference in Covariate Means", 
     xlim = c(-1,2))
axis(1)
axis(2, at = 1:length(covs), las = 1, lab = cov_labs)
abline(v = 0)
abline(v = 0.1, lty = 3)
abline(v = -0.1, lty = 3)
abline(h = 1:length(covs), col = "cornsilk3", lty=6)
points(stdif_after, 1:length(covs), pch = 17, col = "darkorange2")
legend("topright", legend = c("Initial Data", "Matched Data"), 
       col = c("black", "darkorange2"), pch = c(1, 17), box.col = "cornsilk3", cex = 0.8)
# reset margins to default
par(mar=c(4,4,4,4))


#------------------------------------------------------------------------------------------------#
#------------------- Covariate Balance in Treated vs. Control Days  -------------------#
#------------------------------------------------------------------------------------------------#

# re-level so TRUE is the first level of treatment
dat_all$is_treated <- as.character(dat_all$is_treated)
dat_all$is_treated <- factor(dat_all$is_treated, levels = c("TRUE", "FALSE"))

# week of the year density

# initial data
week1 <- dat_all %>%
  ggplot() +
  geom_density(aes(week, fill = is_treated), alpha = 0.5, color = "white") +
  labs(title = "Initial Data", x = "Week of the Year", y = "Density", fill = "") +
  scale_fill_manual(labels = c("Treated (Hot)", "Control (Warm)"),
                     values = c("tomato", "deepskyblue3")) +
  ylim(0,.1) +
  theme_classic() + 
  theme(plot.title = element_text(size = 12),
        axis.title = element_text(size = 10))

# matched data
week2 <- dat_all %>%
  filter(matched == 1) %>%
  ggplot() +
  geom_density(aes(week, fill = is_treated), alpha = 0.5, color = "white") +
  labs(title = "Matched Data", x = "Week of the Year", y = "Density", fill = "") +
  scale_fill_manual(labels = c("Treated (Hot)", "Control (Warm)"),
                     values = c("tomato", "deepskyblue3")) +
  ylim(0,.1) +
  theme_classic() +
  theme(plot.title = element_text(size = 12),
        axis.title = element_text(size = 10))

ggarrange(week1, week2, nrow = 1, common.legend = TRUE)


#--------------- VIOLIN PLOTS for pollutants ---------------#

# o3
dat_after %>%
  ggplot(aes(x = is_treated, y = o3, fill = is_treated)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  ggtitle("Matched Data") +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(~city)

# no2
dat_after %>%
  ggplot(aes(x = is_treated, y = no2, fill = is_treated)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  ggtitle("Matched Data") +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(~city)

# pm10
dat_after %>%
  ggplot(aes(x = is_treated, y = pm10, fill = is_treated)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  ggtitle("Matched Data") +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(~city)

# pm25
dat_after %>%
  ggplot(aes(x = is_treated, y = pm25, fill = is_treated)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  ggtitle("Matched Data") +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(~city)


#--------------- look at lags ---------------#

### tmax ###

tmax_means <- dat_after %>%
  gather(key = "lag", value = "lag_tmax", 
         c("tmax", "tmax_lag_1", "tmax_lag_2", "tmax_lag_3", "tmax_lag_4", "tmax_lag_5")) %>%
  group_by(is_treated, lag) %>%
  summarise(mean_tmax = mean(lag_tmax, na.rm = TRUE))

# re-level
tmax_means$lag <- factor(tmax_means$lag, levels = c("tmax_lag_5", "tmax_lag_4", "tmax_lag_3", 
                                                    "tmax_lag_2", "tmax_lag_1", "tmax"))
tmax_means$is_treated <- as.character(tmax_means$is_treated)
tmax_means$is_treated <- factor(tmax_means$is_treated, levels = c("TRUE", "FALSE"))

tmax_lag <- tmax_means %>%
  ggplot(aes(x = as.factor(lag), y = mean_tmax, group = is_treated, color = is_treated)) +
  geom_line() +
  geom_point(size = 3, shape = c(rep(16,3), rep(1,3), rep(16,3), rep(1,3))) +
  labs(title = "High Temperature", x = "", y = expression(paste("Mean High Temp. (", ~degree, "F)")), color = "") + 
  scale_x_discrete(labels = c("Lag 3", "Lag 2", "Lag 1", 
                              "Day 1", "Day 2", "Day 3")) +
  scale_color_manual(labels = c("Treated (Hot)", "Control (Warm)"),
                     values = c("red", "blue")) +
  theme_bw() + 
  theme(legend.text = element_text(size = 17),
        legend.key.size = unit(2.2, 'cm'),
        plot.title = element_text(size = 15))


### OZONE ###

o3_means <- dat_after %>%
  gather(key = "lag", value = "lag_o3", 
         c("o3_lag_5", "o3_lag_4", "o3_lag_3", "o3_lag_2", "o3_lag_1", "o3")) %>%
  group_by(is_treated, lag) %>%
  summarise(mean_o3 = mean(lag_o3, na.rm = TRUE))

# re-level
o3_means$lag <- factor(o3_means$lag, levels = c("o3_lag_5", "o3_lag_4", "o3_lag_3", 
                                                "o3_lag_2", "o3_lag_1", "o3"))
o3_means$is_treated <- as.character(o3_means$is_treated)
o3_means$is_treated <- factor(o3_means$is_treated, levels = c("TRUE", "FALSE"))

o3_lag <- o3_means %>%
  ggplot(aes(x = as.factor(lag), y = mean_o3, group = is_treated, color = is_treated)) +
  geom_line() +
  geom_point(size = 3, shape = c(rep(16,3), rep(1,3), rep(16,3), rep(1,3))) +
  labs(title = expression('O'[3]), x = "", y = expression('Mean O'[3]*' (ppb)'), color = "") + 
  scale_x_discrete(labels = c("Lag 3", "Lag 2", "Lag 1", 
                              "Day 1", "Day 2", "Day 3")) +
  scale_color_manual(labels = c("Treated", "Control"),
                     values = c("red", "blue")) +
  theme_bw() #+
  #theme(legend.position = "none")


### NO2 ###

no2_means <- dat_after %>%
  gather(key = "lag", value = "lag_no2", 
         c("no2_lag_5", "no2_lag_4", "no2_lag_3", "no2_lag_2", "no2_lag_1", "no2")) %>%
  group_by(is_treated, lag) %>%
  summarise(mean_no2 = mean(lag_no2, na.rm = TRUE))

# re-level
no2_means$lag <- factor(no2_means$lag, levels = c("no2_lag_5", "no2_lag_4", "no2_lag_3", 
                                                  "no2_lag_2", "no2_lag_1", "no2"))
no2_means$is_treated <- as.character(no2_means$is_treated)
no2_means$is_treated <- factor(no2_means$is_treated, levels = c("TRUE", "FALSE"))

no2_lag <- no2_means %>%
  ggplot(aes(x = as.factor(lag), y = mean_no2, group = is_treated, color = is_treated)) +
  geom_line() +
  geom_point(size = 3, shape = c(rep(16,3), rep(1,3), rep(16,3), rep(1,3))) +
  labs(title = expression('NO'[2]), x = "", y = expression('Mean NO'[2]*' (ppb)'), color = "") +
  scale_x_discrete(labels = c("Lag 3", "Lag 2", "Lag 1", 
                              "Day 1", "Day 2", "Day 3")) +
  scale_color_manual(labels = c("Treated", "Control"),
                     values = c("red", "blue")) +
  theme_bw() #+
  #theme(legend.position = "none")


### PM10 ###

pm10_means <- dat_after %>%
  gather(key = "lag", value = "lag_pm10", 
         c("pm10_lag_5", "pm10_lag_4", "pm10_lag_3", "pm10_lag_2", "pm10_lag_1", "pm10")) %>%
  group_by(is_treated, lag) %>%
  summarise(mean_pm10 = mean(lag_pm10, na.rm = TRUE))

# re-level
pm10_means$lag <- factor(pm10_means$lag, levels = c("pm10_lag_5", "pm10_lag_4", "pm10_lag_3", 
                                                    "pm10_lag_2", "pm10_lag_1", "pm10"))
pm10_means$is_treated <- as.character(pm10_means$is_treated)
pm10_means$is_treated <- factor(pm10_means$is_treated, levels = c("TRUE", "FALSE"))

pm10_lag <- pm10_means %>%
  ggplot(aes(x = as.factor(lag), y = mean_pm10, group = is_treated, color = is_treated)) +
  geom_line() +
  geom_point(size = 3, shape = c(rep(16,3), rep(1,3), rep(16,3), rep(1,3))) +
  labs(title = expression('PM'[10]), x = "", y = expression('Mean PM'[10]*' (\u00b5g/m\u00b3)'), color = "") +
  scale_x_discrete(labels = c("Lag 3", "Lag 2", "Lag 1", 
                              "Day 1", "Day 2", "Day 3")) +
  scale_color_manual(labels = c("Treated", "Control"),
                     values = c("red", "blue")) +
  theme_bw() #+
  #theme(legend.position = "none")


### PM25 ###

pm25_means <- dat_after %>%
  gather(key = "lag", value = "lag_pm25", 
         c("pm25_lag_5", "pm25_lag_4", "pm25_lag_3", "pm25_lag_2", "pm25_lag_1", "pm25")) %>%
  group_by(is_treated, lag) %>%
  summarise(mean_pm25 = mean(lag_pm25, na.rm = TRUE))

# re-level
pm25_means$lag <- factor(pm25_means$lag, levels = c("pm25_lag_5", "pm25_lag_4", "pm25_lag_3", 
                                                    "pm25_lag_2", "pm25_lag_1", "pm25"))
pm25_means$is_treated <- as.character(pm25_means$is_treated)
pm25_means$is_treated <- factor(pm25_means$is_treated, levels = c("TRUE", "FALSE"))

pm25_lag <- pm25_means %>%
  ggplot(aes(x = as.factor(lag), y = mean_pm25, group = is_treated, color = is_treated)) +
  geom_line() +
  geom_point(size = 3, shape = c(rep(16,3), rep(1,3), rep(16,3), rep(1,3))) +
  labs(title = expression('PM'[2.5]), x = "", y = expression('Mean PM'[2.5]*' (\u00b5g/m\u00b3)'), color = "") +
  scale_x_discrete(labels = c("Lag 3", "Lag 2", "Lag 1", 
                              "Day 1", "Day 2", "Day 3")) +
  scale_color_manual(labels = c("Treated", "Control"),
                     values = c("red", "blue")) +
  theme_bw() #+
  #theme(legend.position = "none")


### RHUM ###

rhum_means <- dat_after %>%
  gather(key = "lag", value = "lag_rhum", 
         c("rhum_lag_5", "rhum_lag_4", "rhum_lag_3", "rhum_lag_2", "rhum_lag_1", "rhum")) %>%
  group_by(is_treated, lag) %>%
  summarise(mean_rhum = mean(lag_rhum, na.rm = TRUE))

# re-level
rhum_means$lag <- factor(rhum_means$lag, levels = c("rhum_lag_5", "rhum_lag_4", "rhum_lag_3", 
                                                    "rhum_lag_2", "rhum_lag_1", "rhum"))
rhum_means$is_treated <- as.character(rhum_means$is_treated)
rhum_means$is_treated <- factor(rhum_means$is_treated, levels = c("TRUE", "FALSE"))

rhum_lag <- rhum_means %>%
  ggplot(aes(x = as.factor(lag), y = mean_rhum, group = is_treated, color = is_treated)) +
  geom_line() +
  geom_point(size = 3, shape = c(rep(16,3), rep(1,3), rep(16,3), rep(1,3))) +
  labs(title = "Relative Humidity", x = "", y = "Relative Humidity (%)", color = "") +
  scale_x_discrete(labels = c("Lag 3", "Lag 2", "Lag 1", 
                              "Day 1", "Day 2", "Day 3")) +
  scale_color_manual(labels = c("Treated", "Control"),
                     values = c("red", "blue")) +
  theme_bw() #+
#theme(legend.position = "none")



# all lags in one plot
ggarrange(tmax_lag, o3_lag, no2_lag, pm10_lag, pm25_lag, rhum_lag, common.legend = TRUE)


#-------------- example for tmax with just one pair ----------------------#

tmax_ex <- dat_after %>%
  #filter(pair == 70) %>%
  filter(pair == 53) %>%
  gather(key = "lag", value = "tmax", 
         c("tmax", "tmax_lag_1", "tmax_lag_2", "tmax_lag_3", "tmax_lag_4", "tmax_lag_5")) %>%
  group_by(is_treated, lag) %>%
  dplyr::select(is_treated, lag, tmax)

# re-level
tmax_ex$lag <- factor(tmax_ex$lag, levels = c("tmax_lag_5", "tmax_lag_4", "tmax_lag_3", 
                                                    "tmax_lag_2", "tmax_lag_1", "tmax"))
tmax_ex$is_treated <- as.character(tmax_ex$is_treated)
tmax_ex$is_treated <- factor(tmax_ex$is_treated, levels = c("TRUE", "FALSE"))

tmax_ex %>%
  ggplot(aes(x = as.factor(lag), y = tmax, group = is_treated, color = is_treated)) +
  geom_line() +
  # manually set different shapes for lag vs. experiment
  geom_point(size = 3, shape = c(rep(16,6), rep(1,6))) +
  #geom_point(size = 3, shape = c(rep(16,3), rep(1,3), rep(16,3), rep(1,3))) +
  labs(title = "Daily High Temperature in Los Angeles in 1990: Matched Pair Example", 
       x = "", y = expression(paste(~degree, "F")), color = "") + 
  scale_x_discrete(labels = c("Lag 3", "Lag 2", "Lag 1",
                              "Day 1", "Day 2", "Day 3")) +
  scale_color_manual(labels = c("Aug. 30th \u2013 Sep. 4th", "Sep. 13th \u2013 18th"),
                     values = c("red", "blue")) +
  # line for summer median
  geom_hline(yintercept = 80, lty = 2) +
  # label for summer median
  # geom_text(aes(4, 80, label = "Summer Median", vjust = -1, hjust = -1), 
  #           show.legend = FALSE) + 
  #geom_text(aes(0, 80, label = 80, vjust = -1)) + 
  theme_minimal() +
  theme(legend.position = c(.2,.9),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.text = element_text(size = 9))


#-------------- hot, warm, neither examples ----------------------#

# make mini dataframe
hot_ex <- c(82, 85, 91)
warm_ex <- c(75, 77, 76)
# using Sep 25 to Sep 27
neither_ex <- c(82, 80, 74)
days <- rep(c("Day 1", "Day 2", "Day 3"), 3)
trts <- c(rep("Hot", 3), rep ("Warm", 3), rep("Undefined", 3))

# check <- dat_before_wNA %>%
#   filter(city == "la" & year == 1990 & month == "Sep") %>%
#   dplyr::select(date, tmax, is_treated)

mini_ex <- data.frame(date = days,
                      is_treated = trts,
                      tmax = c(hot_ex, warm_ex, neither_ex))
#levels(mini_ex$is_treated) <- c("Hot", "Warm", "Neither")

mini_ex$is_treated <- factor(mini_ex$is_treated, levels = c("Hot", "Warm", "Undefined"))

# plot treatment examples
mini_ex %>%
  ggplot(aes(x = as.factor(date), y = tmax, group = is_treated, color = is_treated)) +
  geom_line() +
  geom_point(size = 3) +
  labs(title = "Daily High Temperature in Los Angeles in September 1990", 
       x = "", y = expression(paste(~degree, "F")), color = "Time Period") +
  ylim(73, 92) +
  geom_hline(aes(yintercept = 80, col = "median"), lty = 2) +
  scale_color_manual(labels = c("Sep. 2 \u2013 Sep. 4", "Sep. 16 \u2013 Sep. 18", 
                                "Sep. 25 \u2013 Sep. 27", "LA Summer Median"),
                     values = c("Hot" = "red", "Warm" = "blue", "Undefined" = "#C77CFF", "median" = "black")) +
  theme_bw() +
  theme(#legend.position = "NONE",
        axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  guides(color = guide_legend(override.aes = list(linetype = c(1,1,1,2), shape = c(16,16,16,NA)))) +
  facet_grid(~ is_treated)

#------------------------------------------------------------------------------------------------#
#------------------- Covariate Balance in Initial vs. Matched Data  -------------------#
#------------------------------------------------------------------------------------------------#

# month
month_comp <- dat_all %>%
  group_by(month) %>%
  summarise(count = n()) %>%
  mutate(percent = 100*count/sum(count)) %>%
  ggplot() +
  # all data
  geom_line(aes(x = month, y = percent, group = 1, color = "Initial Data")) +
  # matched data
  geom_line(aes(x = month, y = percent, group = 1, color = "Matched Data"), 
            data = dat_after %>%
              group_by(month) %>%
              summarise(count = n()) %>%
              mutate(percent = 100*count/sum(count))) +
  #ylim(0, 30) +
  # manual legend
  scale_color_manual(values = c("Initial Data" = "gold3", 
                                "Matched Data" = "slateblue2")) +
  labs(color = "") +
  theme_minimal() #+
  #theme(axis.title.y = element_text(angle=0))


# dow
dow_comp <- dat_all %>%
  group_by(dow) %>%
  summarise(count = n()) %>%
  mutate(percent = 100*count/sum(count)) %>%
  ggplot() +
  # all data
  geom_line(aes(x = dow, y = percent, group = 1, color = "Initial Data")) +
  # matched data
  geom_line(aes(x = dow, y = percent, group = 1, color = "Matched Data"), 
            data = dat_after %>%
              group_by(dow) %>%
              summarise(count = n()) %>%
              mutate(percent = 100*count/sum(count))) +
  ylim(0, 20) +
  # manual legend
  scale_color_manual(values = c("Initial Data" = "gold3", 
                                "Matched Data" = "slateblue2")) +
  labs(color = "") +
  theme_minimal()


# year
year_comp <- dat_all %>%
  group_by(year) %>%
  summarise(count = n()) %>%
  mutate(percent = 100*count/sum(count)) %>%
  ggplot() +
  # all data
  geom_line(aes(x = year, y = percent, group = 1, color = "Initial Data")) +
  # matched data
  geom_line(aes(x = year, y = percent, group = 1, color = "Matched Data"), 
            data = dat_after %>%
              group_by(year) %>%
              summarise(count = n()) %>%
              mutate(percent = 100*count/sum(count))) +
  ylim(0, 13) +
  # manual legend
  scale_color_manual(values = c("Initial Data" = "gold3", 
                                "Matched Data" = "slateblue2")) +
  labs(color = "") +
  theme_minimal()


# week
week_comp <- dat_all %>%
  group_by(week) %>%
  summarise(count = n()) %>%
  mutate(percent = 100*count/sum(count)) %>%
  ggplot() +
  # all data
  geom_line(aes(x = week, y = percent, group = 1, color = "Initial Data")) +
  # matched data
  geom_line(aes(x = week, y = percent, group = 1, color = "Matched Data"), 
            data = dat_after %>%
              group_by(week) %>%
              summarise(count = n()) %>%
              mutate(percent = 100*count/sum(count))) +
  #ylim(0, 8) +
  # manual legend
  scale_color_manual(values = c("Initial Data" = "gold3", 
                                "Matched Data" = "slateblue2")) +
  labs(color = "") +
  theme_minimal()

# all comparisons in one plot
grid.arrange(month_comp, dow_comp, year_comp, week_comp)


dev.off()

