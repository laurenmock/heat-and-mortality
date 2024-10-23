library(ggplot2)
library(tidyverse)
library(pander)
library(MASS)
library(gridExtra)
library(kableExtra)
library(dplyr)
library(tidyr)
library(DataCombine)
library(readr)


## This file computes Fisherian Intervals, Counternull Values, and Randomized-Adjusted P-values


matched<-read.csv("timeseries\\temperature\\matched_data_temp.csv")



################################################################ FUNCTIONS #####################################################################

## Helper function to compute p-value with tau value = a

a_func<-function(x, one_city, a, pairs){
  
  
  
  
  ( ( (1/pairs) * (sum(one_city$death_sum3[x==1 & one_city$is_treated == 1], na.rm = TRUE) +
                     sum((one_city$death_sum3[x==1 & one_city$is_treated == 0]) + a, na.rm = TRUE))) -
      
      ((1/pairs) * (sum(one_city$death_sum3[x==0 & one_city$is_treated == 0], na.rm = TRUE)
                    + sum( ((one_city$death_sum3[x==0 & one_city$is_treated == 1]) - a), na.rm = TRUE ))) )
  
  
  
}

## Helper function to compute p-values given a sequence of tau values

a_vec<-function(x){
  
  rand_mat <-rand_mat_p
  pairs<-41
  tau_obs<-stats_table$tau[4]
  signif(mean((apply(rand_mat, MARGIN  = 2, FUN = a_func, one_city = p, a = x, pairs = pairs)) >= (tau_obs)),4)
}


## This function computes p-values for various tau values for form Fisher interval
fidiciual_matrix<-function(rand_mat,a_sequence,city, city_index, pairs){
  
  a_seq <- a_sequence ## values to test for interval
  
  ## initialize
  p_val_tau_a <- vector()
  
  p_val_tau_a<- lapply(a_seq, FUN = a_vec)
  
  ret<-data.frame(p_val_tau_a)
  colnames(ret)<-a_seq
  
  return(invisible(ret))
  
  
}



#### Functions taken from Counternull package https://github.com/ymabene/Counternull


## Mean Difference test statistic
find_test_stat_diff_means<-function(sample_data,variable){
  # mean for experimental group (exposed)
  on_mean <-mean((variable)[sample_data$is_treated=="1"])
  # mean for control group (non exposure)
  off_mean <-mean((variable)[sample_data$is_treated=="0"])
  # difference
  test_stat<-on_mean - off_mean
  return(invisible(test_stat))
}


## Compute test statistics for numerous assignment permutations
permutation_null_diff_means<-function(rand_matrix,variable,iterations){
  # permutation vector with differences of means
  perm_samples<-matrix(ncol=1,nrow=iterations)
  # creates distribution
  for(k in 1:iterations)
  {
    on<-mean(variable[rand_matrix[,k]==1]) # exposed
    off<-mean(variable[rand_matrix[,k]==0]) # not exposed
    perm_samples[k]<-on-off
  }
  return(invisible(perm_samples))
}


## Creates plot of null randomization distribution
create_null_distribution<-function(sample_data, extreme, rand_matrix,
                                   permutation_null_function,test_stat,
                                   variable,iterations){
  # Creates permutation vector
  perm_samples<-permutation_null_function(rand_matrix,variable,iterations)
  
  # creates histogram and prints p-value
  null_hist<-hist(perm_samples,breaks=100,col = "gold",
                  main=paste("Null Distribution"), xlab="Test Statistics")
  abline(v=test_stat,col="black",lty=2, lwd=5)
  if (extreme==0){ # smaller test statistics are more extreme
    pvalue<-sum(perm_samples<=(test_stat))/iterations
    
  } else { # larger test statistics are more extreme
    pvalue<-sum(perm_samples>=(test_stat))/iterations
    
  }
  print(paste("Test Statistic =",test_stat))
  print(paste("Pvalue =",pvalue))
  return(invisible(perm_samples))
}








########################################################### Process Data ##########################################################




row.names(matched) <- seq(1:nrow(matched))
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



# order by pair (with each FALSE day first)
matched <- matched %>%
  arrange(pair, is_treated)


cities <- unique(matched$city)


############################################################### Fisher Exact P-values ###############################################################

c<-matched %>% filter(city == cities[1])
la<-matched %>% filter(city == cities[2])
ny<-matched %>% filter(city == cities[3])
p<-matched %>% filter(city == cities[4])
s<-matched %>% filter(city == cities[5])


mean((c$death_sum3)[rand_mat_c[,1]=="1"]) - mean((c$death_sum3)[rand_mat_c[,1]=="0"])


# Chicago

c$is_treated<-as.numeric(c$is_treated) # make is_treated numeric
dm_c<-find_test_stat_diff_means(c,c$death_sum3) # ATE: 12.4737
perm_c<-create_null_distribution(c,1,rand_mat_c,permutation_null_diff_means,dm_c,c$death_sum3,100000)
# p-value: .0002


# LA
la<-la[,c(25,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,26,27,28)] # reorder columns
la$is_treated<-as.numeric(la$is_treated)
dm_la<-find_test_stat_diff_means(la,la$death_sum3) #ATE: 8.458
perm_la<-create_null_distribution(la,1,rand_mat_la,permutation_null_diff_means,dm_la,la$death_sum3,100000)
# p-value: .00759

# NY
ny<-ny[,c(25,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,26,27,28)] # reorder columns
ny$is_treated<-as.numeric(ny$is_treated)
dm_ny<-find_test_stat_diff_means(ny,ny$death_sum3) #ATE: 20.313
perm_ny<-create_null_distribution(ny,1,rand_mat_ny,permutation_null_diff_means,dm_ny,ny$death_sum3,100000)
# p-value: .00013


# Pittsburgh
p<-p[,c(25,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,26,27,28)] # reorder columns
p$is_treated<-as.numeric(p$is_treated)
dm_p<-find_test_stat_diff_means(p,p$death_sum3) #ATE: .463
perm_p<-create_null_distribution(p,1,rand_mat_p,permutation_null_diff_means,dm_p,p$death_sum3,100000)
#p-value= .41762


# Seattle
s<-s[,c(25,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,26,27,28)] # reorder columns
s$is_treated<-as.numeric(s$is_treated)
dm_s<-find_test_stat_diff_means(s,s$death_sum3) #ATE: 3.727
perm_s<-create_null_distribution(s,1,rand_mat_s,permutation_null_diff_means,dm_s,s$death_sum3,100000)
# p-value: .00709






##################################################### Fidiciual/Fisher Intervals ##############################################################


# Chicago 

# [6.0,18.94]
# tau: 12.47
# p-value: .0002

library(Counternull)
rand_mat_c<-create_randomization_matrix(c, 100000,76, 1)

write.table(rand_mat_c, file = "timeseries\\temperature\\rand_mat_c.rda")
rand_mat_c<-read.table("timeseries\\temperature\\rand_mat_c.rda")

c_interval<-fidiciual_matrix(rand_mat_c,seq(6.0,6.1,.01),1,38) 

c_interval<-fidiciual_matrix(rand_mat_c,seq(5.96,6,.004),1,38) 

c_interval<-fidiciual_matrix(rand_mat_c,seq(18.9,19.1,.02),c,1,38) 
View(c_interval)



# LA 

# [1.72,15.23]
# tau: 8.46
# p-value: .00759

test<-fidiciual_matrix(rand_mat_la,seq(0,1,1),1,59) 

rand_mat_la<-create_randomization_matrix(la, 100000,118, 1)
write.table(rand_mat_la, file = "timeseries\\temperature\\rand_mat_la.rda")


la_interval<-fidiciual_matrix(rand_mat_la,seq(1.7,1.8,.01),2,59) 
la_interval<-fidiciual_matrix(rand_mat_la,seq(15.2,15.3,.01),2,59)
View(la_interval)

rand_mat_la<-read.table("timeseries\\temperature\\rand_mat_la.rda", skipNul = T)

# New York 

# [10.27,30.36]
# tau: 20.3125
# p-value: .00012

test<-fidiciual_matrix(rand_mat_ny,seq(0,1,1),1,32) 

rand_mat_ny<-create_randomization_matrix(ny, 100000,64, 1)
write.table(rand_mat_ny, file = "timeseries\\temperature\\rand_mat_ny.rda")


ny_interval<-fidiciual_matrix(rand_mat_ny,seq(10.2,10.3,.01),3,32) 
ny_interval<-fidiciual_matrix(rand_mat_ny,seq(30.3,30.4,.01),3,32)
View(ny_interval)



# Pittsburgh 

# [-3.73,4.72]
# tau: .463
# p-value: .419
# counternull value: .75

rand_mat_p<-create_randomization_matrix(c, 100000,82, 1)
write.table(rand_mat_p, file = "timeseries\\temperature\\rand_mat_p.rda")

p_interval<-fidiciual_matrix(rand_mat_p,seq(-3.8,-3.7,.01),4,41)

p_interval<-fidiciual_matrix(rand_mat_p,seq(4.7,4.8,.01),4,41)

View(p_interval)



# Seattle

# [.85,6.61]
# tau: 3.727
# p-value: .00717

rand_mat_s<-create_randomization_matrix(c, 101000,66, 1)
write.table(rand_mat_s, file = "timeseries\\temperature\\rand_mat_s.rda")


s_interval<-fidiciual_matrix(rand_mat_s,seq(.8,.9,.01),5,33)
s_interval<-fidiciual_matrix(rand_mat_s,seq(6.6,6.7,.01),5,33)
View(s_interval)




################################################################### Adjusted p-values ######################################################


# Computes p values treating each observation as "observed"
all_p_values<-function(data){ # obtain all p-values or each city and iterations
  
  
  all_p_values<-data.frame(matrix(NA,nrow=100000,ncol=5))
  
  
  for(i in 1:5){
    
    for(j in 1:100000){
      
      all_p_values[j,i]<-sum(data[,i]>=(data[j,i]))/100000
      
      print(j)
      
    }
    
    
    
  }
  
  return(all_p_values)
  
}


min_p<-function(p_values){ # obtain minimum p-value for each iteration
  
  minimum<-vector()
  
  for(i in 1:100000){
    
    minimum[i]<-min(p_values[i,])
  }
  
  return(minimum)
  
}



all_test_stat<-data.frame("chi"=perm_c,"la"=perm_la,"ny"=perm_ny,"pitt"=perm_p,"seat"=perm_s)

p_values<-all_p_values(all_test_stat)

write.csv(p_values,file="timeseries\\temperature\\all_p_values.csv")


minimum_p_values<-min_p(p_values)


adjusted_p<-function(minimum, obs){
  
  
  return(sum(minimum<=obs)/100000)
  
  
}


adjusted_p(minimum_p_values,.0002) # chicago
# .00093

adjusted_p(minimum_p_values,.00759) # la
#  0.03655

adjusted_p(minimum_p_values,.00013) # ny
# .00585

adjusted_p(minimum_p_values,.41762) # pitt
# .93038


adjusted_p(minimum_p_values,.00709) # seat
# .03461

