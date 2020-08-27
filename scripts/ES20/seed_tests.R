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


##PREPARE PLANTS 

# Normal regular stressed plant -> height from DF seeds

p0_es <- scm_base_parameters("ES20")
p0_es$disturbance_mean_interval <- 30.0

p1_es <- expand_parameters(trait_matrix(c(0.6, 21.9), c("t_s", "a_s")), p0_es, FALSE)
pl_1_es <- ES20_Plant(s = p1_es$strategies[[1]])
pl_es_1 <- pl_1_es # ES20_Plant()
# env_es <- ES20_fixed_environment(1.0)

# Normal regular stressed plant -> height 2x DF seeds

p0_es <- scm_base_parameters("ES20")
p0_es$disturbance_mean_interval <- 30.0

p1_es <- expand_parameters(trait_matrix(c(0.6, 21.9, 2*3.8e-5), c("t_s", "a_s", "omega")), p0_es, FALSE)
pl_1_es <- ES20_Plant(s = p1_es$strategies[[1]])
pl_es_2 <- pl_1_es # ES20_Plant()
# env_es <- ES20_fixed_environment(1.0)

# Normal regular stressed plant -> initial height of 0.4 m

p0_es <- scm_base_parameters("ES20")
p0_es$disturbance_mean_interval <- 30.0

p1_es <- expand_parameters(trait_matrix(c(0.6, 21.9, 0.4), c("t_s", "a_s", "height_0")), p0_es, FALSE)
pl_1_es <- ES20_Plant(s = p1_es$strategies[[1]])
pl_es_3 <- pl_1_es # ES20_Plant()
env_es <- ES20_fixed_environment(1.0)


env_es <- get_constant_environment_ES20(stress = 0.75)

# SINGLE YEAR TESTS
tt <- seq(0, 4, length.out=100)

# Run single year
res_es_1 <- grow_plant_to_time(pl_es_1, tt, env_es)
res_es_2 <- grow_plant_to_time(pl_es_2, tt, env_es)
res_es_3 <- grow_plant_to_time(pl_es_3, tt, env_es)

res_es_df_1 <- as.data.frame(res_es_1$state)
res_es_df_2 <- as.data.frame(res_es_2$state)
res_es_df_3 <- as.data.frame(res_es_3$state)

res_es_df_1 <- cbind(res_es_df_1, res_es_1$aux_size)
res_es_df_2 <- cbind(res_es_df_2, res_es_2$aux_size)
res_es_df_3 <- cbind(res_es_df_3, res_es_3$aux_size)

res_es_df_1$model <- rep("orig", times=nrow(res_es_df_1))
res_es_df_2$model <- rep("2 * orig", times=nrow(res_es_df_2))
res_es_df_3$model <- rep("h0 0.4m", times=nrow(res_es_df_3))

res_es_df_1$time <- tt
res_es_df_2$time <- tt
res_es_df_3$time <- tt

res_all <- plyr::rbind.fill(res_es_df_1, res_es_df_2, res_es_df_3)

# res_all <- res_es_df
# Plot some results


p <- ggplot(res_all, aes(x = time, y = height)) + 
  geom_line(aes(color = model)) +
  scale_x_continuous(name = "Time (yr)") +
  scale_y_continuous(name = "Plant Height (m)") + 
  theme_linedraw()
p





