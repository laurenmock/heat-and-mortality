
###############################################################
#### Causal Inference for Time Series (Heat and Mortality) ####
###############################################################

### Script 2 ### 
# Input: Raw Data 
# (la, la.mortality, chic, chic.mortality, ny, ny.mortality, 
# pitt, pitt.mortality, seat, seat.mortality)
# Output: Processed data (all_processed.csv)

#------------------------------------------------------------------------#
#--------- load libraries ---------#
#------------------------------------------------------------------------#

library(ggplot2)
library(tidyverse)
library(ggpubr)

#------------------------------------------------------------------------#
#--------- set file paths ---------#
#------------------------------------------------------------------------#

# path for loading processed data
processed_data_path <- "data/intermediate/"

# path for EDA figures
EDA_path <- "results/figures/EDA/"

#------------------------------------------------------------------------#
#--------- import data ---------#
#------------------------------------------------------------------------#

# load processed data from script 1
processed_data <- read.csv(file = paste0(processed_data_path, "all_processed.csv"))


#------------------------------------------------------------------------#
# ----------- EDA ----------- # 
#------------------------------------------------------------------------#


#----- summer trends in pollutants -----#

gathered <- processed_data %>%
  gather(key = "poll_name", value = "poll_value", c(12,16:21)) 

# fix city factor levels
gathered$city <- factor(gathered$city,
                        levels = c("pitt", "seat", "chic", "la", "ny"))
levels(gathered$city) <- c("Pittsburgh", "Seattle", "Chicago", "LA", "NY")

# fix pollutant factor levels
gathered$poll_name <- factor(gathered$poll_name,
                             levels = c("tmax", "o3", "no2", "so2", "co", "pm25", "pm10"))
levels(gathered$poll_name) <- c("Temp.[max.]", "O[3]", "NO[2]", "SO[2]", "CO", "PM[2.5]", "PM[10]")

# fix month factor levels
gathered$month <- factor(gathered$month,
                         levels = c("Jun", "Jul", "Aug", "Sep"))

pdf(file = paste0(EDA_path, "summer_trends.pdf"), width = 5, height = 6)

# # plotted by week
# gathered %>%
#   group_by(week, city, poll_name) %>%
#   summarise(mean_poll = mean(poll_value, na.rm = TRUE)) %>%
#   ggplot() +
#   geom_line(aes(x = week, y = mean_poll, color = poll_name), lwd = 1.3) +
#   scale_color_brewer(palette = "Dark2") +
#   facet_grid(poll_name ~ city, 
#              scales = "free",
#              labeller = label_parsed) +
#   labs(x = "Week of the Year", 
#        y = "Mean Pollutant Level") +
#   theme_bw() +
#   theme(legend.position = "none",
#         panel.border = element_rect(colour = "gray", fill = NA, size = 0.5),
#         strip.background = element_rect(fill = "white", color = "gray"),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# plotted by month
gathered %>%
  group_by(month, city, poll_name) %>%
  summarise(mean_poll = mean(poll_value, na.rm = TRUE)) %>%
  ggplot() +
  geom_line(aes(x = month, y = mean_poll, color = poll_name, group = poll_name), lwd = 1.3) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(poll_name ~ city,
             scales = "free",
             labeller = label_parsed) +
  labs(x = "Week of the Year",
       y = "Mean Pollutant Level") +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "gray", fill = NA, size = 0.5),
        strip.background = element_rect(fill = "white", color = "gray"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# # plotted by date
# # class character (not a great solution)
# gathered$month_day <- format(as.Date(gathered$date), "%m-%d")
# gathered %>%
#   group_by(month_day, city, poll_name) %>%
#   summarise(mean_poll = mean(poll_value, na.rm = TRUE)) %>%
#   ggplot() +
#   geom_line(aes(x = month_day, y = mean_poll, color = poll_name, group = poll_name), lwd = 1.3) +
#   scale_color_brewer(palette = "Dark2") +
#   facet_grid(poll_name ~ city, 
#              scales = "free",
#              labeller = label_parsed) +
#   labs(x = "Week of the Year", 
#        y = "Mean Pollutant Level") +
#   theme_bw() +
#   theme(legend.position = "none",
#         panel.border = element_rect(colour = "gray", fill = NA, size = 0.5),
#         strip.background = element_rect(fill = "white", color = "gray"),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dev.off()


#----- correlation between pollutants and tmax during summer months -----#

cities <- c("pitt", "seat", "chic", "la", "ny")

corr_plots <- list()
for(i in 1:length(cities)){
  corr_plots[[i]] <- processed_data %>%
    gather(key = "poll_name", value = "poll_value", 16:21) %>%
    filter(poll_name == names(processed_data)[i+15]) %>%
    ggplot(aes(x = poll_value, y = tmax)) +
    geom_point(alpha = 0.2, color = "dodgerblue") +
    geom_smooth() +
    labs(x = names(processed_data)[i+15]) +
    facet_grid(~ city, scales = "free")
}

pdf(file = paste0(EDA_path, "tmax_correlations.pdf"), width = 5, height = 7)

ggarrange(corr_plots[[1]], corr_plots[[2]], corr_plots[[3]], 
                 corr_plots[[4]], corr_plots[[5]], nrow = 5)

dev.off()

#----- daily mortality density -----#

processed_data %>%
  ggplot() +
  geom_density(aes(death, fill = city), col = "white", alpha = 0.7) +
  ggtitle("Distribution of Death Counts") +
  xlab("Deaths over 3 Day Experiment") +
  theme_classic()

#----- daily mortality histograms -----#

processed_data %>%
  ggplot() +
  geom_histogram(aes(death), bins = 20, col = "white") +
  facet_wrap(~city, scales= "free") +
  xlab("Deaths over 3 Day Experiment") +
  theme_bw()

#----- missing pollutant measurements -----#

pdf(file = paste0(EDA_path, "missing_pollutants.pdf"), width = 6, height = 2)

processed_data |>
  group_by(city, year) |>
  summarize(pct_measured = (sum(!is.na(o3)) / n())*100) |>
  ggplot() +
  geom_bar(aes(x = year, y = pct_measured), stat = "identity") +
  labs(title = "o3", x = "Year", y = "% of days w/ data") +
  facet_wrap(~city, nrow = 1)

processed_data |>
  group_by(city, year) |>
  summarize(pct_measured = (sum(!is.na(no2)) / n())*100) |>
  ggplot() +
  geom_bar(aes(x = year, y = pct_measured), stat = "identity") +
  labs(title = "no2", x = "Year", y = "% of days w/ data") +
  facet_wrap(~city, nrow = 1)

processed_data |>
  group_by(city, year) |>
  summarize(pct_measured = (sum(!is.na(co)) / n())*100) |>
  ggplot() +
  geom_bar(aes(x = year, y = pct_measured), stat = "identity") +
  labs(title = "co", x = "Year", y = "% of w/ data") +
  facet_wrap(~city, nrow = 1)

processed_data |>
  group_by(city, year) |>
  summarize(pct_measured = (sum(!is.na(so2)) / n())*100) |>
  ggplot() +
  geom_bar(aes(x = year, y = pct_measured), stat = "identity") +
  labs(title = "so2", x = "", y = "% of days w/ data") +
  facet_wrap(~city, nrow = 1)

processed_data |>
  group_by(city, year) |>
  summarize(pct_measured = (sum(!is.na(pm25)) / n())*100) |>
  ggplot() +
  geom_bar(aes(x = year, y = pct_measured), stat = "identity") +
  labs(title = "pm25", x = "", y = "% of days w/ data") +
  facet_wrap(~city, nrow = 1)

processed_data |>
  group_by(city, year) |>
  summarize(pct_measured = (sum(!is.na(pm10)) / n())*100) |>
  ggplot() +
  geom_bar(aes(x = year, y = pct_measured), stat = "identity") +
  labs(title = "pm10", x = "", y = "% of days w/ data") +
  facet_wrap(~city, nrow = 1)

processed_data |>
  group_by(city, year) |>
  summarize(pct_measured = (sum(!is.na(rhum)) / n())*100) |>
  ggplot() +
  geom_bar(aes(x = year, y = pct_measured), stat = "identity") +
  labs(title = "rhum", x = "", y = "% of days w/ data") +
  facet_wrap(~city, nrow = 1)

dev.off()

