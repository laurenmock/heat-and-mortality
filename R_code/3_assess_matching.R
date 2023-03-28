
###############################################################
#### Causal Inference for Time Series (Heat and Mortality) ####
###############################################################

### Script 3 ###
# Inputs: Processed data from script 1 and matched data from script 2
# Outputs: All figures saved in figures folder

#------------------------------------------------------------------------#

library(ggplot2)
library(tidyverse)
library(pander)
library(gridExtra)
library(ggpubr)
library(scales)


#------------------------------------------------------------------------#
#--------- set working directory/paths ---------#
#------------------------------------------------------------------------#

# set path for processed data from script 1 and matched data from script 2
processed_data_path <- "data/processed/"

# set path for figures
fig_path <- "figures/covariate_balance/"

# set path for examples of treatment
trt_ex_path <- "figures/treatment_examples/"

#------------------------------------------------------------------------#
#--------- load data ---------#
#------------------------------------------------------------------------#

# load processed data (before matching)
dat_before <- read.csv(paste0(processed_data_path, "all_processed.csv"))

# load matched data (after matching)
dat_after <- read.csv(paste0(processed_data_path, "matched.csv"))

# create column in matched data to indicate that these rows were in the matched sub (after join below)
dat_after$matched <- 1

#------------------------------------------------------------------------------------------------#

# save full df (includes days with undefined treatment)
dat_before_undef <- dat_before

# now remove days with undefined treatment
dat_before <- dat_before %>%
  filter(!is.na(is_treated))

#------------------------------------------------------------------------------------------------#

# create one data frame with before/after matching data together
dat_all <- left_join(dat_before, dat_after)

# matching column indicates 0 for rows that were not matched, 1 for rows that were matched
dat_all$matched <- ifelse(is.na(dat_all$matched), 0, 1)


#------------------------------------------------------------------------------------------------#

# get a new dataframe that can be faceted to original and matched
# where matched data is included in original

dat_all$sub <- "Original Data"
dat_after$sub <- "Matched Data"
init_v_match <- rbind(dat_all, dat_after)
init_v_match$sub <- factor(init_v_match$sub, levels = c("Original Data", "Matched Data"))

# has duplicate rows for matched data (with "sub" column to differentiate)

#------------------------------------------------------------------------------------------------#
#------------------- missing covariates in matched data -------------------#
#------------------------------------------------------------------------------------------------#

relevant_variables <- c("pm10_lag_3","pm25_lag_3","o3_lag_3","no2_lag_3")

for(i in 1:length(relevant_variables)){
  num_miss = dat_after |>
    pull(relevant_variables[i]) |>
    is.na() |>
    sum()
  cat(relevant_variables[i], num_miss, "\n")
}

#------------------------------------------------------------------------------------------------#
#------------------- Love Plots (code from Alice) -------------------#
#------------------------------------------------------------------------------------------------#

# create a standardized difference in means function (stdif)
stdif <- function(X, W){
  mean_diff = mean(X[W == TRUE], na.rm = TRUE) - mean(X[W == FALSE], na.rm = TRUE)
  SE = sqrt(0.5 * (var(X[W == TRUE], na.rm = TRUE) + var(X[W == FALSE], na.rm = TRUE)))
  return(mean_diff/SE)
}

covs <- c("no2_lag_5", "no2_lag_4", "no2_lag_3",
          "o3_lag_5", "o3_lag_4", "o3_lag_3", 
          "tmax_lag_5", "tmax_lag_4", "tmax_lag_3",
          "month_Sep", "month_Aug", "month_Jul", "month_Jun", "year")

cov_labs <- c(expression('NO'[2]*' Lag 3'), expression('NO'[2]*' Lag 2'), expression('NO'[2]*' Lag 1'),
              expression('O'[3]*' Lag 3'), expression('O'[3]*' Lag 2'), expression('O'[3]*' Lag 1'),
              "High Temp. Lag 3", "High Temp. Lag 2", "High Temp. Lag 1",
              "September", "August", "July", "June", "Year")

#-----------------------------------------#

# standardized diff. before matching
stdif_before <- sapply(dat_before %>% dplyr::select(covs),
                           function(x) stdif(x, dat_before %>% pull(is_treated)))

# standardized diff. after matching
stdif_after <- sapply(dat_after %>% dplyr::select(covs),
                       function(x) stdif(x, dat_after %>% pull(is_treated)))

# love plot

pdf(file = paste0(fig_path, "love_plot.pdf"), width = 6, height = 6)

# set margins so we can see the whole plot
par(mar=c(4,8,4,4))
plot(stdif_before, 1:length(covs), 
     axes = FALSE, 
     ylab = "",
     xlab = "Standardized Difference in Covariate Means", 
     xlim = c(-1,2),
     cex = 1.5)
axis(1)
axis(2, at = 1:length(covs), las = 1, lab = cov_labs)
abline(v = 0)
abline(v = 0.1, lty = 3)
abline(v = -0.1, lty = 3)
abline(h = 1:length(covs), col = "cornsilk3", lty=6)
points(stdif_after, 1:length(covs), pch = 17, col = "steelblue1", cex = 1.5)
legend("topright", legend = c("Initial Data", "Matched Data"), 
       col = c("black", "steelblue1"), pch = c(1, 17), box.col = "cornsilk3",
       bg = "white", pt.cex = 1, cex = 1
       )
# reset margins to default
par(mar=c(4,4,4,4))

dev.off()

#------------------------------------------------------------------------------------------------#
#------------------- Covariate Balance in Treated vs. Control Days  -------------------#
#------------------------------------------------------------------------------------------------#

# color palette
pal <- c("orange", "firebrick3")

#----------------- week of the year density -----------------#

# re-level so TRUE is the first level of treatment
dat_all$is_treated <- as.character(dat_all$is_treated)
dat_all$is_treated <- factor(dat_all$is_treated, levels = c("TRUE", "FALSE"))

# ggplot layers for the treated and control week density plots
# density_gg <- list(
#   geom_density(aes(week, fill = is_treated), alpha = 0.7, color = "white"),
#   xlab("Week of the Year"),
#   ylab("Density"),
#   labs(fill = ""),
#   scale_fill_manual(labels = c("Treated (Hot)", "Control (Warm)"),
#                     values = pal),
#   scale_x_continuous(expand = c(0, 0), breaks = pretty_breaks()),
#   scale_y_continuous(limits = c(0, 0.1), expand = c(0, 0)),
#   theme_classic(),
#   theme(plot.title = element_text(size = 12),
#         axis.title = element_text(size = 10),
#         legend.text = element_text(size = 10),
#         legend.key.size = unit(0.8, 'cm'))
# )

#---------- week of the year -----------#

pdf(file = paste0(fig_path, "week_density.pdf"), width = 8, height = 4)

init_v_match %>%
  ggplot() +
  geom_density(aes(week, col = is_treated, fill = is_treated), alpha = 0.5, size = 1.5) +
  xlab("Week of the Year") +
  ylab("Density") +
  labs(fill = "") +
  scale_color_manual(name = "",
                     labels = c("Control (Warm)", "Treated (Hot)"),
                  values = pal) +
  scale_fill_manual(name = "",
                    labels = c("Control (Warm)", "Treated (Hot)"),
                    values = pal) +
  scale_x_continuous(expand = c(0, 0), breaks = pretty_breaks()) + 
  scale_y_continuous(limits = c(0, 0.1), expand = c(0, 0)) +
  theme_bw() +
  theme(plot.title = element_text(size = 12),
      axis.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.key.size = unit(0.8, 'cm'),
      legend.position = "bottom",
      panel.border = element_rect(colour = "gray", fill = NA, size = 0.5),
      strip.background = element_rect(fill = "white", color = "gray"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      strip.text.x = element_text(size = 12)) + 
  facet_grid(~sub)

dev.off()


#--------------- time series ---------------#

# ggplot layers for all time series
time_series_gg = list(
  geom_line(size = 1.2),
  # geom_point(size = 2.5, 
  #            shape = c(rep(16,3), rep(1,3), rep(16,3), rep(1,3)),
  #            stroke = 2,
  #            # shape = 16,
  #            #alpha = c(rep(1,3), rep(0.5,3), rep(1,3), rep(0.5,3))
  #            ),
  geom_point(size = 2.5,
             shape = 16),
  scale_x_discrete(labels = c("Lag 3", "Lag 2", "Lag 1", 
                              "Day 1", "Day 2", "Day 3")),
  scale_color_manual(labels = c("Control (Warm)", "Treated (Hot)"),
                     values = rev(pal)),
  theme_classic(),
  theme(legend.text = element_text(size = 12),
        legend.key.size = unit(1.5, 'cm'),
        #line = element_blank(), # remove axes
        #panel.border = element_rect(colour = "black", fill = NA, size = 0.5), # add border
        plot.title = element_text(size = 12))
)


### tmax ###

tmax_means <- dat_after %>%
  gather(key = "lag", value = "lag_tmax", 
         c("tmax", "tmax_lag_1", "tmax_lag_2", "tmax_lag_3", "tmax_lag_4", "tmax_lag_5")) %>%
  group_by(is_treated, lag) %>%
  summarise(mean_tmax = mean(lag_tmax, na.rm = TRUE))

# re-level
tmax_means$lag <- factor(tmax_means$lag, levels = c("tmax_lag_5", "tmax_lag_4", "tmax_lag_3", 
                                                    "tmax_lag_2", "tmax_lag_1", "tmax"))

# tmax plot
tmax_lag <- tmax_means %>%
  ggplot(aes(x = as.factor(lag), y = mean_tmax, group = is_treated, color = is_treated)) +
  labs(title = "Max. Temperature", 
       x = "", 
       y = expression(paste("Mean Max. Temp. (", ~degree, "F)")), 
       color = "") + 
  ylim(62, 88) + 
  time_series_gg


### OZONE ###

o3_means <- dat_after %>%
  gather(key = "lag", value = "lag_o3", 
         c("o3_lag_5", "o3_lag_4", "o3_lag_3", "o3_lag_2", "o3_lag_1", "o3")) %>%
  group_by(is_treated, lag) %>%
  summarise(mean_o3 = mean(lag_o3, na.rm = TRUE))

# re-level
o3_means$lag <- factor(o3_means$lag, levels = c("o3_lag_5", "o3_lag_4", "o3_lag_3", 
                                                "o3_lag_2", "o3_lag_1", "o3"))

# o3 plot
o3_lag <- o3_means %>%
  ggplot(aes(x = as.factor(lag), y = mean_o3, group = is_treated, color = is_treated)) +
  labs(title = expression('O'[3]), 
       x = "", 
       y = expression('Mean O'[3]*' (\u00b5g/m\u00b3)'), 
       color = "") + 
  ylim(22, 34) + 
  time_series_gg


### NO2 ###

no2_means <- dat_after %>%
  gather(key = "lag", value = "lag_no2", 
         c("no2_lag_5", "no2_lag_4", "no2_lag_3", "no2_lag_2", "no2_lag_1", "no2")) %>%
  group_by(is_treated, lag) %>%
  summarise(mean_no2 = mean(lag_no2, na.rm = TRUE))

# re-level
no2_means$lag <- factor(no2_means$lag, levels = c("no2_lag_5", "no2_lag_4", "no2_lag_3", 
                                                  "no2_lag_2", "no2_lag_1", "no2"))

# no2 plot
no2_lag <- no2_means %>%
  ggplot(aes(x = as.factor(lag), y = mean_no2, group = is_treated, color = is_treated)) +
  labs(title = expression('NO'[2]), 
       x = "", 
       y = expression('Mean NO'[2]*' (\u00b5g/m\u00b3)'), 
       color = "") +
  ylim(22, 34) + 
  time_series_gg

### PM10 ###

pm10_means <- dat_after %>%
  gather(key = "lag", value = "lag_pm10", 
         c("pm10_lag_5", "pm10_lag_4", "pm10_lag_3", "pm10_lag_2", "pm10_lag_1", "pm10")) %>%
  group_by(is_treated, lag) %>%
  summarise(mean_pm10 = mean(lag_pm10, na.rm = TRUE))

# re-level
pm10_means$lag <- factor(pm10_means$lag, levels = c("pm10_lag_5", "pm10_lag_4", "pm10_lag_3", 
                                                    "pm10_lag_2", "pm10_lag_1", "pm10"))

# pm10 plot
pm10_lag <- pm10_means %>%
  ggplot(aes(x = as.factor(lag), y = mean_pm10, group = is_treated, color = is_treated)) +
  labs(title = expression('PM'[10]), 
       x = "", 
       y = expression('Mean PM'[10]*' (\u00b5g/m\u00b3)'), 
       color = "") +
  ylim(10, 50) + 
  time_series_gg


### PM25 ###

pm25_means <- dat_after %>%
  gather(key = "lag", value = "lag_pm25", 
         c("pm25_lag_5", "pm25_lag_4", "pm25_lag_3", "pm25_lag_2", "pm25_lag_1", "pm25")) %>%
  group_by(is_treated, lag) %>%
  summarise(mean_pm25 = mean(lag_pm25, na.rm = TRUE))

# re-level
pm25_means$lag <- factor(pm25_means$lag, levels = c("pm25_lag_5", "pm25_lag_4", "pm25_lag_3", 
                                                    "pm25_lag_2", "pm25_lag_1", "pm25"))

# pm25 plot
pm25_lag <- pm25_means %>%
  ggplot(aes(x = as.factor(lag), y = mean_pm25, group = is_treated, color = is_treated)) +
  labs(title = expression('PM'[2.5]), 
       x = "", 
       y = expression('Mean PM'[2.5]*' (\u00b5g/m\u00b3)'), 
       color = "") +
  ylim(10, 50) + 
  time_series_gg


### relative humidity ###

rhum_means <- dat_after %>%
  gather(key = "lag", value = "lag_rhum", 
         c("rhum_lag_5", "rhum_lag_4", "rhum_lag_3", "rhum_lag_2", "rhum_lag_1", "rhum")) %>%
  group_by(is_treated, lag) %>%
  summarise(mean_rhum = mean(lag_rhum, na.rm = TRUE))

# re-level
rhum_means$lag <- factor(rhum_means$lag, levels = c("rhum_lag_5", "rhum_lag_4", "rhum_lag_3", 
                                                    "rhum_lag_2", "rhum_lag_1", "rhum"))

# relative humidity plot
rhum_lag <- rhum_means %>%
  ggplot(aes(x = as.factor(lag), y = mean_rhum, group = is_treated, color = is_treated)) +
  labs(title = "Relative Humidity", 
       x = "", 
       y = "Relative Humidity (%)", 
       color = "") +
  ylim(62, 88) + 
  time_series_gg

###

# all lags in one plot
pdf(file = paste0(fig_path, "time_series.pdf"), width = 10, height = 5)

ggarrange(tmax_lag, o3_lag, no2_lag, rhum_lag, pm10_lag, pm25_lag, 
          common.legend = TRUE, legend = "bottom")

dev.off()


#------------------------------------------------------------------------------------------------#
#------------------- Covariate Balance in Initial vs. Matched Data  -------------------#
#------------------------------------------------------------------------------------------------#

#---------- year -----------#

# plot
year_plot <- init_v_match |>
  ggplot() +
  geom_density(aes(year, fill = as.factor(sub)), alpha = 0.7, color = "white") +
  labs(x = "Year", y = "Density", fill = "") +
  scale_fill_manual(labels = c("Original Data", "Matched Data"),
                    values = c("#7570B3", "#66A61E")) +
  scale_x_continuous(expand = c(0, 0), breaks = pretty_breaks()) +
  scale_y_continuous(limits = c(0, 0.1), expand = c(0, 0)) +
  theme_classic() +
  theme(legend.position = "bottom")


#---------- week of the year ----------#

# plot
week_plot <- init_v_match |>
  ggplot() +
  geom_density(aes(week, fill = as.factor(sub)), alpha = 0.5, color = "white") +
  labs(x = "Week of the Year", y = "Density", fill = "") +
  scale_fill_manual(labels = c("Original Data", "Matched Data"),
                    values = c("#377EB8", "#E7298A")) +
  scale_x_continuous(expand = c(0, 0), breaks = pretty_breaks()) +
  scale_y_continuous(limits = c(0, 0.1), expand = c(0, 0)) +
  theme_classic() +
  theme(legend.position = "bottom")

# continuous comparisons
pdf(file = paste0(fig_path, "initial_vs_matched_cont.pdf"), width = 10, height = 4)
ggarrange(year_plot, week_plot, ncol = 2)
dev.off()

#---------- month ----------#

# get data frames with percentages
month_before <- dat_all %>%
  group_by(month) %>%
  summarise(count = n()) %>%
  mutate(percent = 100*count/sum(count),
         type = "Original Data")
month_after <- dat_after %>%
  group_by(month) %>%
  summarise(count = n()) %>%
  mutate(percent = 100*count/sum(count),
         type = "Matched Data")
month_comp <- rbind(month_before, month_after)

# fix month labels
month_comp$month <- factor(month_comp$month, levels = c("Jun", "Jul", "Aug", "Sep"))
levels(month_comp$month) <- c("June", "July", "Aug.", "Sept.")

# fix data labels
month_comp$type <- factor(month_comp$type, levels = c("Original Data", "Matched Data"))

# plot
month_plot <- month_comp |>
  ggplot() +
  geom_bar(aes(fill = month, x = type, y = percent), 
           position="fill", stat="identity", color = "white", alpha = 0.5) +
  labs(x = "", y = "Proportion", fill = "Month") +
  scale_fill_brewer(palette = "Set1") + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()



#----------- day of week ----------#

# get dataframe with percentages
dow_before <- dat_all %>%
  group_by(dow) %>%
  summarise(count = n()) %>%
  mutate(percent = 100*count/sum(count),
         type = "Original Data")
dow_after <- dat_after %>%
  group_by(dow) %>%
  summarise(count = n()) %>%
  mutate(percent = 100*count/sum(count),
         type = "Matched Data")
dow_comp <- rbind(dow_before, dow_after)

# reorder levels
dow_comp$dow <- factor(dow_comp$dow, levels = c("Monday", "Tuesday", "Wednesday", "Thursday",
                                                    "Friday", "Saturday", "Sunday"))

# fix data labels
dow_comp$type <- factor(dow_comp$type, levels = c("Original Data", "Matched Data"))

# plot
dow_plot <- dow_comp |>
  ggplot() +
  geom_bar(aes(fill = as.factor(dow), x = type, y = percent), 
           position="fill", stat="identity", color = "white", alpha = 0.5) +
  labs(x = "", y = "Proportion", fill = "Day of Week") +
  scale_fill_brewer(palette = "Dark2") + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()


# categorical comparisons
pdf(file = paste0(fig_path, "initial_vs_matched_cat.pdf"), width = 11, height = 4)
ggarrange(month_plot, dow_plot, ncol = 2)
dev.off()


#------------------------------------------------------------------------------------------------#
#------------------- Examples to Visualize Treatment/Time Series  -------------------#
#------------------------------------------------------------------------------------------------#

#-------------- example for tmax with just one pair ----------------------#

tmax_ex <- dat_after %>%
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

png(file = paste0(trt_ex_path, "time_series_example.png"), width = 550, height = 300)

tmax_ex %>%
  ggplot(aes(x = as.factor(lag), y = tmax, group = is_treated)) +
  geom_line(aes(color = is_treated)) +
  # manually set different shapes for lag vs. experiment
  geom_point(aes(color = is_treated), size = 3, shape = c(rep(16,6), rep(1,6))) +
  #geom_point(size = 3, shape = c(rep(16,3), rep(1,3), rep(16,3), rep(1,3))) +
  labs(title = "Daily High Temperature in Los Angeles in 1990: Matched Pair Example", 
       x = "", y = expression(paste(~degree, "F")), color = "") + 
  scale_x_discrete(labels = c("Lag 3", "Lag 2", "Lag 1",
                              "Day 1", "Day 2", "Day 3")) +
  scale_color_manual(labels = c("Aug. 30th \u2013 Sep. 4th", "Sep. 13th \u2013 18th"),
                     values = pal) +
  # line for summer median
  geom_hline(yintercept = 80, lty = 2) +
  # label for summer median
  geom_text(aes(4, 80, label = "LA Summer Median", vjust = -1, hjust = -1),
            show.legend = FALSE) +
  theme_minimal() +
  theme(legend.position = c(.2,.9),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.text = element_text(size = 9))

dev.off()

#-------------- hot, warm, neither examples ----------------------#

# make mini dataframe

find_ex <- dat_before_undef %>%
  filter(city == "la" & year == 1990 & month == "Sep") %>%
  dplyr::select(date, tmax, is_treated)

days <- rep(c("Day 1", "Day 2", "Day 3"), 3)
trts <- c(rep ("Warm", 3), rep("Hot", 3), rep("Undefined", 3))

mini_ex <- data.frame(date = days,
                      is_treated = trts,
                      tmax = c(75, 77, 76, 82, 85, 91, 82, 80, 74))

mini_ex$is_treated <- factor(mini_ex$is_treated, levels = c("Warm", "Hot", "Undefined"))


pdf(file = paste0(trt_ex_path, "trt_cntrl_undef.pdf"), width = 7, height = 2)

# plot treatment examples
mini_ex %>%
  ggplot(aes(x = as.factor(date), y = tmax, group = is_treated, color = is_treated)) +
  geom_hline(aes(yintercept = 80, col = "median"), lty = 2) +
  geom_line() +
  geom_point(size = 3) +
  labs(x = "", y = expression(paste(~degree, "F")), color = "") +
  ylim(73, 92) +
  scale_color_manual(labels = c("Sep. 16 \u2013 Sep. 18", "Sep. 2 \u2013 Sep. 4", 
                                "Sep. 25 \u2013 Sep. 27", "LA Summer Median"),
                     values = c("Warm" = "orange",
                                "Hot" = "firebrick3",
                                "Undefined" = "gray50", 
                                "median" = "black")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5),
        panel.border = element_rect(colour = "gray", fill = NA, size = 0.5),
        strip.background = element_rect(fill = "white", color = "gray"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(linetype = c(1,1,1,2), 
                                                  shape = c(16,16,16,NA), 
                                                  col = c("orange", "firebrick3",
                                                          "gray50", "black")))) +
  facet_grid(~ is_treated)

dev.off()

