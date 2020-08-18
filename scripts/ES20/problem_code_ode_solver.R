##LOAD LATEST PACKAGE

# restart R - clear the C++ output files cache (can't find a different way of doing this. 
# Makes sure the latest compiled version is used)
rm(list = ls())
# .rs.restartR()

gc()
setwd("~/R/plant")
# library(deSolve)
library(tidyverse)

# Use latest version 
RcppR6::install(".")
devtools::load_all()

##PREPARE PLANT

# Normal regular stressed plant

p0_es <- scm_base_parameters("ES20")
p0_es$disturbance_mean_interval <- 30.0
p1_es <- expand_parameters(trait_matrix(c(0.41, 21.9), c("t_s", "a_s")), p0_es, FALSE)
pl_1_es <- ES20_Plant(s = p1_es$strategies[[1]])

pl_es <- pl_1_es # ES20_Plant()
env_es <- ES20_fixed_environment(1.0)


# SINGLE YEAR TESTS
tt <- seq(0, 1, length.out=30)


# Run single year
res_es <- grow_plant_to_time(pl_es, tt, env_es)
