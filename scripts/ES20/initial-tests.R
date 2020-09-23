rm(list = ls())
gc()
setwd("~/R/plant")
library(deSolve)
library(ggplot2)
library(tidyverse)

# Use latest version 
RcppR6::install(".")
# devtools::document(".")
devtools::load_all()




## Simple running tests


#Test strategy
s2 <- ES20_Strategy()

p <- scm_base_parameters("ES20")
s1 <- ES20_Strategy(trait_matrix(0.0825, "lma"), p)

# This bit does not work yet. So it must be a problem on my part 
p3 <- ES20_Parameters()

# p0 <- scm_base_parameters("ES20")
# 
# class(p3)
# 
# out  <- run_stochastic_collect(p3)
# 
# ## Test individual plant FF16
# s <- FF16_Strategy()
# pl <- FF16_Plant(s)
# 
# env <- FF16_fixed_environment(1.0)
# 
# pl$ode_rates
# pl$compute_rates(env)
# 

#### ES20 (individual plant)

pl <- ES20_Plant(s2)

env <- ES20_fixed_environment(1.0)

pl$ode_rates
# pl$aux_names
pl$ode_names
# pl$ode_size
pl$ode_state

pl$compute_rates(env)

pl$ode_rates
pl$ode_state <- pl$ode_rates + pl$ode_state
pl$ode_state


## Common settings

tt <- seq(0, 1, length.out=10001)

## grow a plant - FF16

pl_ff <- FF16_Plant()
env_ff <- FF16_fixed_environment(1.0)

## grow a plant - not quite through the system but getting there


pl_es <- ES20_Plant()
env_es <- ES20_fixed_environment(1.0)


res_ff <- grow_plant_to_time(pl_ff, tt, env_ff)
res_es <- grow_plant_to_time(pl_es, tt, env_es)







plot(height ~ tt, res_es$state, type="l", las=1,
     xlab="Time (years)", ylab="Height (m)", col="#669966", ylim=c(0, max(res_ff$state[,'height'])))
lines(height ~ tt, res_ff$state)


res_ff$state['height', 10001]

plot(net_mass_production_dt ~ tt, res_es$aux_size, type="l", las=1,
     xlab="Time (years)", ylab="NMP", col="#669966")
lines(mortality ~ tt, res_ff$state)


plot(fecundity ~ tt, res_ff$state, type="l", las=1,
     xlab="Time (years)", ylab="Fecundity (unknown - not sure)")
lines(fecundity ~ tt, res_es$state, col="#669966")

plot(mass_storage ~ tt, res_es$state, type="l", las=1,
     xlab="Time (years)", ylab="Mass Storage (kg)")

plot(mass_heartwood ~ tt, res_ff$state, type="l", las=1,
     xlab="Time (years)", ylab="Mass HW (kg)")

pl#mass of components

plot(mass_sapwood ~ tt, res_es$state, type="l", las=1,
     xlab="Time (years)", ylab="Mass (kg)", col="#00cccc",
     ylim = c(0, max(mass_leaf, mass_sapwood, mass_heartwood, mass_bark, mass_root)))
lines(mass_leaf ~ tt, res_es$state, type="l", las=1, col="#669966")
lines(mass_heartwood ~ tt, res_es$state, type="l", las=1, col="#993300")
lines(mass_bark ~ tt, res_es$state, type="l", las=1, col="#FFCC00")
lines(mass_root ~ tt, res_es$state, type="l", las=1, col="#CC3399")


#proportion of storage

plot(mass_storage/(mass_sapwood+mass_leaf+mass_bark+mass_root) ~ tt, res_es$state, type="l", las=1,
     xlab="Time (years)", ylab="Proportion of storage to live mass")

plot(dbiomass_dt/res_es$state[,'mass_storage'] ~ tt, res_es$aux_size, type="l", las=1,
     xlab="Time (years)", ylab="Proportion of storage to live mass")



### From Daniel's code:


p0_ff <- scm_base_parameters("FF16")
p0_ff$disturbance_mean_interval <- 30.0

p1_ff <- expand_parameters(trait_matrix(0.0825, "lma"), p0_ff, FALSE)

p0_es <- scm_base_parameters("ES20")
p0_es$disturbance_mean_interval <- 30.0

p1_es <- expand_parameters(trait_matrix(0.0825, "lma"), p0_es, FALSE)
p1_es <- expand_parameters(trait_matrix(0.0325, "lma"), p1_es, FALSE)
p1_es <- expand_parameters(trait_matrix(0.9, "a_s"), p1_es, FALSE)


pl_1_es <- ES20_Plant(s = p1_es$strategies[[1]])
pl_2_es <- ES20_Plant(s = p1_es$strategies[[3]])
pl_3_es <- ES20_Plant(s = p1_es$strategies[[4]])


pl_1_es$strategy$lma
pl_2_es$strategy$lma
pl_2_es$strategy$a_s
pl_3_es$strategy$a_s



p1$seed_rain <- 20






#p1_eq <- equilibrium_seed_rain(p1)

p2 <- build_schedule(p1)
data1 <- run_scm_collect(p2)

matplot(data1$time, data1$species[[1]]["height", , ], lty=1, col=make_transparent("black", 0.25), type="l",
        las=1, xlab="Time (years)", ylab="Height (m)")


seed_rain_out <- attr(p_new, "seed_rain_out", exact=TRUE)