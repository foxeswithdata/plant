##LOAD LATEST PACKAGE

# restart R - clear the C++ output files cache (can't find a different way of doing this. 
# Makes sure the latest compiled version is used)
# .rs.restartR()

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

# Non-stressed plant

env_ns <- get_constant_environment_ES20()
pl_ns <- get_non_stressed_plant_ES20()

# Normal regular stressed plant

pl_es <- ES20_Plant()
env_es <- ES20_fixed_environment(1.0)

# FF16 plant 

pl_ff <- FF16_Plant()
env_ff <- FF16_fixed_environment(1.0)

# SINGLE YEAR TESTS

tt <- seq(0, 1, length.out=11)


# Run single year
res_ns <- grow_plant_to_time(pl_ns, tt, env_ns)
res_es <- grow_plant_to_time(pl_es, tt, env_es)
res_ff <- grow_plant_to_time(pl_ff, tt, env_ff)


# Plot some results

plot(height ~ tt, res_es$state, type="b", las=1,
     xlab="Time (years)", ylab="Height (m)", col="#669966", ylim=c(0, max(c(res_ns$state[,'height'],res_ff$state[,'height']))))
lines(height ~ tt, res_ns$state, pch = 18, type="b")
lines(height ~ tt, res_ff$state, col="red", pch = 18, type="b")

plot(mass_storage ~ tt, res_es$state, type="l", las=1,
     xlab="Time (years)", ylab="Mass Storage (kg)", col="#669966", ylim=c(0, max(c(res_ns$state[,'mass_storage']))))
lines(mass_storage ~ tt, res_ns$state)

plot(fecundity ~ tt, res_es$state, type="l", las=1,
     xlab="Time (years)", ylab="Mass Storage (kg)", col="#669966", ylim=c(0, max(c(res_ns$state[,'fecundity']))))
lines(fecundity ~ tt, res_ns$state)

plot(mass_storage/(mass_sapwood+mass_leaf+mass_bark+mass_root) ~ tt, res_es$state, type="l", las=1, col="#669966",
     xlab="Time (years)", ylab="Proportion of storage to live mass")
lines(mass_storage/(mass_sapwood+mass_leaf+mass_bark+mass_root) ~ tt, res_ns$state)

plot(dbiomass_dt/res_es$state[,'mass_storage'] ~ tt, res_es$aux_size, type="l", las=1, col="#669966",
     xlab="Time (years)", ylab="Biomass dt / mass storage", ylim = c(0, 1.1))
lines(dbiomass_dt/res_ns$state[,'mass_storage'] ~ tt, res_ns$aux_size)


plot(dbiomass_dt ~ tt, res_es$aux_size, type="l", las=1, col="#669966",
     xlab="Time (years)", ylab="Biomass dt", ylim=c(0, max(c(res_ns$state[,'mass_storage']))))
lines(dbiomass_dt ~ tt, res_ns$aux_size)
lines(mass_storage ~ tt, res_es$state, col="#669966", lty=2)
lines(mass_storage ~ tt, res_ns$state,  lty=2)

plot((net_mass_production_dt + pl_es$strategy$a_y * pl_es$strategy$a_bio *  respiration_dt) ~ tt, res_es$aux_size, type="l", las=1, col="#669966",
     xlab="Time (years)", ylab="Net Mass Production")
lines(net_mass_production_dt ~ tt, res_es$aux_size, lty=2, col="#669966")
lines(pl_es$strategy$a_y * pl_es$strategy$a_bio * respiration_dt ~ tt, res_es$aux_size, lty=3, col="#669966")


lines((net_mass_production_dt + respiration_dt) ~ tt, res_ns$aux_size)
lines(net_mass_production_dt ~ tt, res_ff$aux_size, col="red")



plot((net_mass_production_dt + respiration_dt) ~ tt, res_ns$aux_size, type="l", las=1, col="#669966",
     xlab="Time (years)", ylab="Net Mass Production", ylim=c(0, 0.02))
lines(net_mass_production_dt ~ tt, res_ns$aux_size, lty=2, col="#669966")
lines(respiration_dt ~ tt, res_ns$aux_size, lty=3, col="#669966")



plot(mass_sapwood ~ tt, res_es$state, type="l", las=1,
     xlab="Time (years)", ylab="Mass (kg)", col="#00cccc",
     ylim = c(0, max(mass_leaf, mass_sapwood, mass_heartwood, mass_bark, mass_root)))
lines(mass_leaf ~ tt, res_es$state, type="l", las=1, col="#669966")
lines(mass_heartwood ~ tt, res_es$state, type="l", las=1, col="#993300")
lines(mass_bark ~ tt, res_es$state, type="l", las=1, col="#FFCC00")
lines(mass_root ~ tt, res_es$state, type="l", las=1, col="#CC3399")






plot(competition_effect ~ tt, res_ff$aux_size, type="l", las=1,
     xlab="Time (years)", ylab="area (m2)", col="red")
lines(area_leaf ~ tt, res_es$state, col="#669966")
lines(area_leaf ~ tt, res_ns$state)


plot(mass_leaf ~ tt, res_ff$aux_size, type="l", las=1,
     xlab="Time (years)", ylab="leaf mass (kgC)", col="red", ylim=c(0, 0.17))
lines(mass_leaf ~ tt, res_es$state, col="#669966")
lines(mass_leaf ~ tt, res_ns$state)

plot(mass_sapwood ~ tt, res_ff$aux_size, type="l", las=1,
     xlab="Time (years)", ylab="leaf mass (kgC)", col="red")
lines(mass_sapwood ~ tt, res_es$state, col="#669966")
lines(mass_sapwood ~ tt, res_ns$state)

plot(mass_root ~ tt, res_ff$aux_size, type="l", las=1,
     xlab="Time (years)", ylab="leaf mass (kgC)", col="red")
lines(mass_root ~ tt, res_es$state, col="#669966")
lines(mass_root ~ tt, res_ns$state)

plot(mass_bark ~ tt, res_ff$aux_size, type="l", las=1,
     xlab="Time (years)", ylab="leaf mass (kgC)", col="red")
lines(mass_bark ~ tt, res_es$state, col="#669966")
lines(mass_bark ~ tt, res_ns$state)
