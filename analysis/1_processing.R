
###############################################################
#### Causal Inference for Time Series (Heat and Mortality) ####
###############################################################

### Script 1 ### 
# Input: Raw Data 
# (la, la.mortality, chic, chic.mortality, ny, ny.mortality, 
# pitt, pitt.mortality, seat, seat.mortality)
# Output: Processed data (all_processed.csv)

#------------------------------------------------------------------------#
#--------- load libraries ---------#
#------------------------------------------------------------------------#

library(ggplot2)
library(tidyverse)
library(gridExtra)
library(lubridate)
library(fastDummies)


#------------------------------------------------------------------------#
#--------- set file paths ---------#
#------------------------------------------------------------------------#

# path for folder with raw data files (should only contain those files)
raw_data_path <- "data/raw/"

# path for writing processed data
processed_data_path <- "data/processed/"

# path for EDA figures
EDA_path <- "figures/EDA/"

#------------------------------------------------------------------------#
#--------- import data ---------#
#------------------------------------------------------------------------#

# get all file names in the raw data folder

# names of files in the working directory that end in ".rda" and don't contain "mortality"
files <- grep(list.files(path = raw_data_path,
                         pattern = "\\.rda$"), 
              pattern = 'mortality', invert = TRUE, value = TRUE)
# names of files in the working directory that end in ".mortality.rda"
files_mort <- list.files(path = raw_data_path,
                         pattern = "\\.mortality.rda$")


# function to load files into a list
load_files <- function(file) {
  temp <- new.env()
  load(file, envir = temp)
  as.list(temp)
}


### load RDA files into two lists (city data, then mortality by age data)
### then combine each list into a data frame

# city.rda

# load files into list
files_list <- Map(load_files, file.path(folder = raw_data_path, files))
# name list elements
names(files_list) <- tools::file_path_sans_ext(files)
# convert to list of data frames (instead of list of lists)
files_list <- lapply(files_list, function(x) do.call(rbind, x))
# join all cities into one data frame 
raw_data <- bind_rows(files_list)
rownames(raw_data) <- NULL # remove row names


# city.mortality.rda

# load files into list
files_mort_list <- Map(load_files, file.path(folder = raw_data_path, files_mort)) 
# name list elements
names(files_mort_list) <- names(files_list)
# convert to list of data frames (instead of list of lists)
files_mort_list <- lapply(files_mort_list, function(x) do.call(rbind.data.frame, x))
# join all cities into one data frame
raw_data_mortality <- bind_rows(files_mort_list, .id = "city")
rownames(raw_data_mortality) <- NULL # remove row names

#------------------------------------------------------------------------#

# remove the deaths of people under 65
raw_data$death <- subset(raw_data_mortality, agecat == "65to74")$death + 
  subset(raw_data_mortality, agecat == "75p")$death

#------------------------------------------------------------------------#
#--------- check for missing data ---------#
#------------------------------------------------------------------------#

# count total number of unique days in original data
raw_data |>
  filter(month(raw_data$date) %in% c(6,7,8,9)) |>
  nrow()

# count total number of days missing outcome (death)
raw_data |>
  filter(month(raw_data$date) %in% c(6,7,8,9)) |>
  pull(death) |>
  is.na() |>
  sum()

# count total number of days missing exposure (temp)
raw_data |>
  filter(month(raw_data$date) %in% c(6,7,8,9)) |>
  pull(tmax) |>
  is.na() |>
  sum()


# Specify the relevant covariates to include
# relevant_variables <- names(raw_data) # all covariates
relevant_variables = c("dptp","rhum",
                       "pm10","pm25","o3","no2", "so2","co")

for(i in 1:length(relevant_variables)){
  num_miss = raw_data |>
    filter(month(raw_data$date) %in% c(6,7,8,9)) |>
    pull(relevant_variables[i]) |>
    is.na() |>
    sum()
  cat(relevant_variables[i], num_miss, "\n")
}

# check in raw data
if (exists("relevant_variables") && !is.null(relevant_variables)){
  rows_with_NAs = c()
  for (i in 1:length(relevant_variables)){
    is_missing = is.na(raw_data[relevant_variables[i]])
    rows_with_NAs = c(rows_with_NAs,(1:nrow(raw_data))[is_missing])
    nb_missing = sum(is_missing)
    nb_missing_percent = round(100*nb_missing/nrow(raw_data), 2)
    if (nb_missing>0){
      # print(raw_data[is.na(raw_data[relevant_variables[i]]),])
      cat("\n >>>>>>>>> WARNING: MISSING", nb_missing_percent, "%", relevant_variables[i],
          "values <<<<<<<<<\n")
    }
  }
  rows_with_NAs = sort(unique(rows_with_NAs))
  # print(raw_data$date[rows_with_NAs])
}

# check in mortality data
mortality_with_NAs = is.na(raw_data_mortality$death)
print(subset(raw_data_mortality[mortality_with_NAs,], agecat %in% c("65to74","75p")))

# remove any rows with missing outcome (very few rows)
raw_data <- raw_data %>%
  filter(!is.na(death))

#------------------------------------------------------------------------#
#--------- add columns for lags (values on previous days) ---------#
#------------------------------------------------------------------------#

# remove existing lags (so that the new lags will all be consistent)
raw_data <- raw_data %>%
  # remove columns that start with l1, l2, or l3
  dplyr::select(-str_subset(names(raw_data), "^l[1-3]"))


# For each covariate, create new covariate defined as the value of that covariate on 
# (day - 1), (day - 2), ... up to (day - nb_lag)

# city names
cities <- names(files_list)

# make each city a df in a list
cities_list <- list()
for(i in 1:length(cities)){
  cities_list[[i]] <- raw_data %>% filter(city == cities[i])
}

# number of lags to calculate for each covariate
nb_lag <- 5


# function that I will apply to each city to calculate lags
calc_lags <- function(df){
  
  # loop through each column
  for(j in 1:length(names(df))){
    
    # loop through number of lags to calculate
    for(l in 1:nb_lag){
      
      # create new column name
      new_col_name <- paste(names(df)[j], "_lag_", toString(l), sep = "")
      
      # add a lag column to the data frame
      df <- df %>% cbind(lag(df[[j]], n = l))

      # name column
      names(df)[ncol(df)] <- new_col_name
      
    }
  }
  
  # create a column for deaths on the following day
  df$death_nextDay <- lead(df$death)
  
  return(df)
}

# apply function to get lags for each city
cities_list_lags <- lapply(cities_list, calc_lags)

# bind all the cities back together
processed_data <- bind_rows(cities_list_lags)

#------------------------------------------------------------------------#

# add explicit covariates for day of week, month, year (coded as integers)
month = factor(c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
               levels = c("Jan","Feb","Mar","Apr","May","Jun", "Jul","Aug","Sep","Oct","Nov","Dec"))
processed_data$month = month[as.numeric(format(processed_data$date,"%m"))]
processed_data$year = as.numeric(format(processed_data$date,"%Y"))

#------------------------------------------------------------------------#
#-------- new columns for matching --------# 
#------------------------------------------------------------------------#

# weekend vs. weekday
processed_data$weekend <-
  ifelse(processed_data$dow == "Saturday" | processed_data$dow == "Sunday", TRUE, FALSE)

# week of the year
processed_data$week <- isoweek(processed_data$date)

#------------------------------------------------------------------------#
# ------- new columns for dummy variables (need for love plots) -------- # 
#------------------------------------------------------------------------#

# day of week
processed_data <- dummy_cols(processed_data, select_columns = c("dow", "month"))


#------------------------------------------------------------------------#
# ----------- EDA ----------- # 
#------------------------------------------------------------------------#

#----- summer trends in pollutants -----#

gathered <- processed_data %>%
  filter(month %in% c("Jun", "Jul", "Aug", "Sep")) %>%
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

# trends_labels <- c(
#   `no2` = "Chicago"
# )

pdf(file = paste0(EDA_path, "summer_trends.pdf"), width = 5, height = 6)

# plotted by week
gathered %>%
  group_by(week, city, poll_name) %>%
  summarise(mean_poll = mean(poll_value, na.rm = TRUE)) %>%
  ggplot() +
  geom_line(aes(x = week, y = mean_poll, color = poll_name), lwd = 1.3) +
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


# plotted by date
# class character (not a great solution)
gathered$month_day <- format(gathered$date, "%m-%d")
gathered %>%
  group_by(month_day, city, poll_name) %>%
  summarise(mean_poll = mean(poll_value, na.rm = TRUE)) %>%
  ggplot() +
  geom_line(aes(x = month_day, y = mean_poll, color = poll_name, group = poll_name), lwd = 1.3) +
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

dev.off()


#----- correlation between pollutants and tmax during summer months -----#

corr_plots <- list()
for(i in 1:length(cities)){
  corr_plots[[i]] <- processed_data %>%
    filter(month %in% c("Jun", "Jul", "Aug", "Sep")) %>%
    gather(key = "poll_name", value = "poll_value", 16:21) %>%
    filter(poll_name == names(processed_data)[i+15]) %>%
    ggplot(aes(x = poll_value, y = tmax)) +
    geom_point(alpha = 0.2, color = "dodgerblue") +
    geom_smooth() +
    labs(x = names(processed_data)[i+15]) +
    facet_grid(~ city, scales = "free")
}

g <- arrangeGrob(corr_plots[[1]], corr_plots[[2]], corr_plots[[3]], 
                 corr_plots[[4]], corr_plots[[5]], nrow = 5)

ggsave(file = paste0(EDA_path, "tmax_correlations.png"), width = 5, height = 8, g)


#----- daily mortality density -----#

png(file = paste0(EDA_path, "mortality_density.png"), width = 500, height = 500)

processed_data %>%
  ggplot() +
  geom_density(aes(death, fill = city), col = "white", alpha = 0.7) +
  ggtitle("Distribution of Death Counts") +
  xlab("Deaths over 3 Day Experiment") +
  theme_classic()

dev.off()

#----- daily mortality histograms -----#

png(file = paste0(EDA_path, "mortality_histograms.png"), width = 500, height = 500)

processed_data %>%
  ggplot() +
  geom_histogram(aes(death), bins = 20, col = "white") +
  facet_wrap(~city, scales= "free") +
  xlab("Deaths over 3 Day Experiment") +
  theme_bw()

dev.off()


##########################################################################
#--------------------------- DEFINE TREATMENT ---------------------------#
##########################################################################

#------------------------------------------------------------------------#
#-------- Define single days as high/low temperature --------#
#------------------------------------------------------------------------#

# initialize list
city_df <- list()

# loop through cities
for(i in 1:length(cities)){
  
  # filter to include only one city at a time
  city_df[[i]] <- filter(processed_data, city == cities[i]) # look at one city at a time
  
  # find low/high tmax days
  
  # find median in summer months
  tmax_median <- city_df[[i]] %>%
    filter(month %in% c("Jun", "Jul", "Aug", "Sep")) %>%
    pull(tmax) %>%
    median(na.rm = TRUE)
  
  # find low/high tmax days
  tmax_low <- city_df[[i]]$tmax <= tmax_median
  tmax_high <- city_df[[i]]$tmax > tmax_median
  
  # is_tmax_high?
  city_df[[i]]$is_tmax_high[tmax_low] <- FALSE
  city_df[[i]]$is_tmax_high[tmax_high] <- TRUE
  
  # all non-summer days should be NA
  city_df[[i]]$is_tmax_high <- ifelse(city_df[[i]]$month %in% c("Jun", "Jul", "Aug", "Sep"), 
                                      city_df[[i]]$is_tmax_high, # is summer, nothing changes
                                      NA) # if not summer, make NA
}


# bind all the cities back together
treated_data <- bind_rows(city_df)


#------------------------------------------------------------------------#
#--- Define Treatment (Depends on Previous Days) ---#
#------------------------------------------------------------------------#

# length of time period
time_prd <- 3

# initialize column for treatment
treated_data$is_treated <- NA

# loop through rows of data frame, starting with time_prd
for(l in (time_prd):nrow(treated_data)){

  # only look at days and lag days with non-missing exposure values
  if(all(!is.na(treated_data$is_tmax_high[(l - (time_prd-1)):l]))){

    # assign control day if current day and previous lag days were all low
    if(all((treated_data$is_tmax_high[(l - (time_prd-1)):l]) == FALSE)){
      treated_data$is_treated[l] <- FALSE
    }

    # assign treated day if current day and previous lag days were all high
    if(all((treated_data$is_tmax_high[(l - (time_prd-1)):l]) == TRUE)){
      treated_data$is_treated[l] <- TRUE
    }
  }
}

# write CSV for processed data
write.csv(treated_data, 
          file = paste0(processed_data_path, "all_processed.csv"),
          row.names = FALSE)

#--------------------------------------------------------------------------------------------#
#------------ Check relevant columns to make sure treated vs. control makes sense -----------#
#--------------------------------------------------------------------------------------------#

check_treatment <- treated_data %>%
  dplyr::select(c("date", "month", "city", "death", "tmax", "is_treated"))

