rm(list = ls())
.rs.restartR()

rm(list = ls())
gc()
setwd("~/R/plant")
# library(deSolve)
library(tidyverse)

# Use latest version 
RcppR6::install(".")
devtools::load_all()

source('scripts/ES20/test_objects.R')

ts_values <- seq(from=0, to=0.75, by=0.005)
h_0_values <- seq(from=0.5, to=40, by=0.5)
# h_0_values <- seq(from=0.2, to=1, by=0.2)
b_s_values <- seq(from=0.02, to=0.5, by=0.02)
# b_s_values <- seq(from=0.02, to=0.5, by=0.06)

# opt_ts_M <- matrix(0, nrow=length(h_0_values), ncol = length(b_s_values))
# opt_ts_S <- matrix(0, nrow=length(h_0_values), ncol = length(b_s_values))
# 
# for(i in 1:length(h_0_values)){
#   for(j in 1:length(b_s_values)){
#     out_Mass <- vector()
#     out_Store <- vector()
#     valid <- vector()
#     
#     for(k in 1:length(ts_values)){
#       
#       p0_es <- scm_base_parameters("ES20")
#       p0_es$disturbance_mean_interval <- 30.0
#       
#       p1_es <- expand_parameters(trait_matrix(c(ts_values[k], 21.9, h_0_values[i], b_s_values[j]), c("t_s", "a_s", "height_0", "b_s")), p0_es, FALSE)
#       pl_1_es <- ES20_Plant(s = p1_es$strategies[[1]])
#       pl_es <- pl_1_es # ES20_Plant()
#       env_es <- get_constant_environment_ES20(stress=0.75)
#       
#       tt <- seq(0, 1, length.out=100)
#       
#       # Run single year
#       res_es <- grow_plant_to_time(pl_es, tt, env_es)
#       res_es_df <- as.data.frame(res_es$state)
#       
#       final <- res_es_df[nrow(res_es_df),]
#       mass <- final[,'mass_leaf'] + final[, 'mass_sapwood'] + final[, 'mass_bark'] + final[,'mass_root']
#       
#       out_Mass[k] <- mass
#       out_Store[k] <- final[,'mass_storage']
#       
#       if(any(res_es_df$mass_storage < 0, na.rm=TRUE)){
#         valid[k] <- FALSE
#       }
#       else{
#         valid[k] <- TRUE
#       }
#     }
#     
#     print("finished")
#     
#     ## Evaluate the outputs
#     
#     # get rid of results that are invalid
#     
#     out_Mass[!valid] <- -99999
#     out_Store[!valid] <- -99999
#     
#     if(any(valid, na.rm=TRUE)){
#       opt_ts_M[i,j] = ts_values[which.max(out_Mass)]
#       opt_ts_S[i,j] = ts_values[which.max(out_Store)]
#     }
#     else{
#       opt_ts_M[i,j] = NA
#       opt_ts_S[i,j] = NA
#     }
#     
#   }
# }



p0_es <- scm_base_parameters("ES20")
p0_es$disturbance_mean_interval <- 30.0

p1_es <- expand_parameters(trait_matrix(c(ts_values[18], 21.9, h_0_values[27], b_s_values[23]), c("t_s", "a_s", "height_0", "b_s")), p0_es, FALSE)
pl_1_es <- ES20_Plant(s = p1_es$strategies[[1]])
pl_es <- pl_1_es # ES20_Plant()
env_es <- get_constant_environment_ES20(stress=0.75)

tt <- seq(0, 0.991, length.out=100)

# Run single year
res_es <- grow_plant_to_time(pl_es, tt, env_es)
res_es_df <- as.data.frame(res_es$state)