
###############################################################
#### Causal Inference for Time Series (Heat and Mortality) ####
###############################################################

### Script 3 ###
# Input: Processed data from script 1
# Output: Data frame that only contains matched pairs (matched.csv)

#------------------------------------------------------------------------#
#------------------------------------------------------------------------#

library(ggplot2)
library(tidyverse)

#------------------------------------------------------------------------#
#--------- set working directory/paths ---------#
#------------------------------------------------------------------------#

# source matching functions
source("analysis/matching_functions.R")

# path to read in processed data and write out matched data
processed_data_path <- "data/intermediate/"

#------------------------------------------------------------------------#
#------------------------------------------------------------------------#

# load processed data from script 1
dat_processed <- read.csv(file = paste0(processed_data_path, "all_processed.csv"))

# remove rows where is_treated is undefined
matching_data <- dat_processed[!is.na(dat_processed$is_treated),] 
N <- nrow(matching_data)

#------------------------------------------------------------------------#
#--- columns with characters must be factors (if using for matching) ----#
#------------------------------------------------------------------------#

matching_data$city <- as.factor(matching_data$city)
matching_data$dow <- as.factor(matching_data$dow)


#------------------------------------------------------------------------#
#--------- explore treated vs. control units ---------#
#------------------------------------------------------------------------#

treated_units = subset(matching_data, is_treated)
control_units = subset(matching_data, !is_treated)
N_treated = nrow(treated_units)
N_control = nrow(control_units)
cat("Number of treated units:", N_treated,"\nNumber of control units:", N_control)


### Some quick exploration of the (marginal) distributions of the covariates among treated
# table(matching_data$dow[matching_data$is_treated])
# table(matching_data$month[matching_data$is_treated])
# table(matching_data$year[matching_data$is_treated])
# hist(matching_data$death_DayMinus1[matching_data$is_treated],breaks=N_treated)
# hist(matching_data$tmpd[matching_data$is_treated],breaks=N_treated)

#------------------------------------------------------------------------#
#--------- set scaling for matching ---------#
#------------------------------------------------------------------------#

# Optional weights for each covariate when computing the distances

# WARNING: the order of the items in scaling needs to be the same as the order of the covariates (i.e. columns)
scaling =  rep(list(1), ncol(matching_data)) # L: 1 for each column
names(scaling) = colnames(matching_data)

# # The following function makes it easier to set the scaling factors (using covariate name instead of index)
# set_scaling = function(scaling, covariate_name, scaling_factor = NULL){
#   updated_scaling = scaling
#   if (covariate_name %in% colnames(matching_data)){
#     index = which(colnames(matching_data)==covariate_name)
#     if (!is.null(scaling_factor)){
#       updated_scaling[[index]] = scaling_factor
#     } else {
#       s = tryCatch(1/sd(matching_data[[index]][!is.na(matching_data[[index]])]),
#                    error = function(e){1})
#       updated_scaling[[index]] = s
#     }
#   } else {
#     "Error: no covariate corresponds to that name"
#   }
#   return (updated_scaling)
# }


# set the scaling factors for some covariates
# scaling$tmpd = 1/sd(matching_data$tmpd[!is.na(matching_data$tmpd)])
# scaling$death = 1/sd(matching_data$death[!is.na(matching_data$death)])

#------------------------------------------------------------------------#
#--------- set thresholds for matching ---------#
#------------------------------------------------------------------------#

# set the thresholds for each covariate, default is Inf (i.e. doesn't need to match)
thresholds = rep(list(Inf), ncol(matching_data))
names(thresholds) = colnames(matching_data)

# set threshold values (determines which covariates must match and how closely)
thresholds$city = 0
thresholds$year = 0
thresholds$dow = 0
thresholds$week = 3
thresholds$no2_lag_3 = 5
thresholds$o3_lag_3 = 5

#------------------------------------------------------------------------#
#--------- discrepancy matrix ---------#
#------------------------------------------------------------------------#

# Compute the discrepancy matrix
discrepancies = discrepancyMatrix(treated_units, control_units, thresholds, scaling)
cat("Number of prospective matched treated units =", sum(rowSums(discrepancies < Inf) > 0))

# make row and column names date.city

rownames(discrepancies) = paste(matching_data$date[which(matching_data$is_treated)], 
                                as.character(matching_data$city[which(matching_data$is_treated)]), sep = ".")
colnames(discrepancies) = paste(matching_data$date[which(!matching_data$is_treated)], 
                                as.character(matching_data$city[which(!matching_data$is_treated)]), sep = ".")

rownames(matching_data) <- paste(matching_data$date, 
                                 as.character(matching_data$city), sep = ".")

#------------------------------------------------------------------------#
#--------- Perform Matching ---------#
#------------------------------------------------------------------------#

library(optmatch)

# # WARNING: pairmatch sometimes returns empty match even when fullmatch returns possible matches !!! WHY ??!!
# # I think it's because the built-in pairmatch doesn't know how to "break ties"
# matched_groups = pairmatch(discrepancies,data=matching_data,remove.unmatchables = TRUE)

#matched_groups = fullmatch(discrepancies, data = matching_data, remove.unmatchables = TRUE, mean.controls = 1)
matched_groups = fullmatch(discrepancies, data = matching_data, remove.unmatchables = TRUE)


# Get list of groups
groups_labels = unique(matched_groups[!is.na(matched_groups)])
groups_list = list()
for (i in 1:length(groups_labels)){
  IDs = names(matched_groups)[(matched_groups == groups_labels[i])]
  groups_list[[i]] = as.Date(IDs[!is.na(IDs)])
}

# Print diagnostics
print(matched_groups, grouped = TRUE)
summary(matched_groups)


#------------------------------------------------------------------------------------------------#
#------------------ Force pair matching via bipartite maximal weighted matching -----------------#
#------------------------------------------------------------------------------------------------#

### this is how we go from several days that match to 1:1 pairs ###

library(igraph)

# We build a bipartite graph with one layer of treated nodes, and another layer of control nodes.
# The nodes are labeled by integers from 1 to (N_treated + N_control)
# By convention, the first N_treated nodes correspond to the treated units, and the remaining N_control
# nodes correspond to the control units.

#-----------------------------------------------------------------------------

# build pseudo-adjacency matrix: edge if and only if match is admissible
# NB: this matrix is rectangular so it is not per say the adjacendy matrix of the graph
# (for this bipartite graph, the adjacency matrix had four blocks: the upper-left block of size
# N_treated by N_treated filled with 0's, bottom-right block of size N_control by N_control filled with 0's,
# top-right block of size N_treated by N_control corresponding to adj defined below, and bottom-left block
# of size N_control by N_treated corresponding to the transpose of adj)
adj = (discrepancies < Inf)

# extract endpoints of edges
edges_mat = which(adj,arr.ind = TRUE)

# build weights, listed in the same order as the edges (we use a decreasing function x --> 1/(1+x) to
# have weights inversely proportional to the discrepancies, since maximum.bipartite.matching
# maximizes the total weight and we want to minimize the discrepancy)
weights = 1/(1+sapply(1:nrow(edges_mat), function(i)discrepancies[edges_mat[i,1], edges_mat[i,2]]))

# format list of edges (encoded as a vector resulting from concatenating the end points of each edge)
# i.e c(edge1_endpoint1, edge1_endpoint2, edge2_endpoint1, edge2_endpoint1, edge3_endpoint1, etc...)
edges_mat[,"col"] = edges_mat[,"col"] + N_treated
edges_vector = c(t(edges_mat))

# NB: by convention, the first N_treated nodes correspond to the treated units, and the remaining N_control
# nodes correspond to the control units (hence the "+ N_treated" to shift the labels of the control nodes)

#-----------------------------------------------------------------------------

# Build the graph from the list of edges
BG = make_bipartite_graph(c(rep(TRUE, N_treated), rep(FALSE, N_control)), edges = edges_vector)

# Find the maximal weighted matching
MBM = maximum.bipartite.matching(BG, weights = weights)

# List the dates (and city) of the matched pairs
pairs_list = list()
N_matched = 0
for (i in 1:N_treated){ # loop through treated units and look for matching controls
  if (!is.na(MBM$matching[i])){
    N_matched = N_matched + 1
    pairs_list[[N_matched]] = c(as.character(treated_units$date[i]), 
                                as.character(control_units$date[MBM$matching[i] - N_treated]), 
                                as.character(treated_units$city[i]))
  }
}
cat("Number of matched treated units =", N_matched)

# Quick sanity check for matched pairs
#relevant_fields = colnames(matching_data)[which(unlist(thresholds) < Inf)]
# relevant_fields = c("city", "month", "year", "o3", "weekend", "dow_cats", "death")
#
# lag <- 2
#
# for (i in 1:5){ # change to 1:N_matched to look at ALL matches
#   t_i = toString(subset(matching_data, date %in% (pairs_list[[i]]) & is_treated)$date)
#   c_i = toString(subset(matching_data, date %in% (pairs_list[[i]]) &! is_treated)$date)
#   cat("\n-------------------- Matched pair",i,":",t_i,",",c_i,"--------------------\n")
#   print(subset(matching_data,date %in% (pairs_list[[i]]))[relevant_fields])
#   cat("...........\n")
#   for (j in 1:2){ # L: changed this to be 1:2 so it includes dates 1 and 2 and doesn't look for 3 (the city name)
#     print(subset(dat_processed, 
#                  date %in% (as.Date(pairs_list[[i]][j]) - 0:lag) & city == pairs_list[[i]][3])[relevant_fields])
#   }
# }


#------------------------------------------------------------------------#
#----------------------- Matching Data Frames  --------------------------#
#------------------------------------------------------------------------#

### new data frame for treated (matched) ###

# initialize list
treated_matched <- list()

# get rows in matching_data where the date/city is in pairs_list
for(i in 1:length(pairs_list)){
  treated_matched[[i]] <- filter(matching_data, date == pairs_list[[i]][1] & city == pairs_list[[i]][3])
}

# bind all rows together
treated_matched <- bind_rows(treated_matched)

# add pair column to match to controls
treated_matched$pair <- 1:nrow(treated_matched)

#------------------------------------------------------------------------#

### new data frame for control (matched) ###

# initialize list
control_matched <- list()

# get rows in matching_data where the date/city is in pairs_list
for(i in 1:length(pairs_list)){
  control_matched[[i]] <- filter(matching_data, date == pairs_list[[i]][2] & city == pairs_list[[i]][3])
}

# bind all rows together
control_matched <- bind_rows(control_matched)

# add pair column to match to treated
control_matched$pair <- 1:nrow(control_matched)

#------------------------------------------------------------------------#

# full data frame with all paired matches

matched <- rbind(treated_matched, control_matched)

#------------------------------------------------------------------------#

# new column for sum of deaths on 3 exp. days
matched$death_sum3 <- matched$death + matched$death_lag_1 + matched$death_lag_2

# new column for i (index)
matched$id <- 1:nrow(matched)

#------------------------------------------------------------------------#
#------------------ SAVE MATCHED DATA -----------------------#
#------------------------------------------------------------------------#

# write CSV for matched data
write.csv(matched, 
          file = paste0(processed_data_path, "matched.csv"), 
          row.names = FALSE)

