rm(list = ls())
.rs.restartR()

# rm(list = ls())
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

print("Problematic Values")
print("height0")
print(h_0_values[27])
print("b_s")
print(b_s_values[23])
print("ts")
print(paste(c(ts_values[18] * 365), "d"), collapse=" ") #NoteL for ts_values <i=18 and 19 and 20 this has actually worked. it's so weird

p0_es <- scm_base_parameters("ES20")
p0_es$disturbance_mean_interval <- 30.0

p1_es <- expand_parameters(trait_matrix(c(ts_values[18], 21.9, h_0_values[27], b_s_values[23]), c("t_s", "a_s", "height_0", "b_s")), p0_es, FALSE)
pl_1_es <- ES20_Plant(s = p1_es$strategies[[1]])
pl_es <- pl_1_es # ES20_Plant()
env_es <- get_constant_environment_ES20(stress=0.75)

tt <- seq(0, 1.1, length.out=100)

# Run single year
res_es <- grow_plant_to_time(pl_es, tt, env_es)
res_es_df <- as.data.frame(res_es$state)