##LOAD LATEST PACKAGE

# restart R - clear the C++ output files cache (can't find a different way of doing this. 
# Makes sure the latest compiled version is used)
rm(list = ls())
.rs.restartR()

rm(list = ls())
gc()
setwd("~/R/plant") 
# library(deSolve)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)

# Use latest version 
RcppR6::install(".")
devtools::load_all()

source('scripts/ES20/test_objects.R')


##PREPARE PLANTS 

# Non-stressed plant
# 
# env_ns <- get_constant_environment_ES20()
# pl_ns <- get_non_stressed_plant_ES20()


# Normal regular stressed plant

# p0_es <- scm_base_parameters("ES20")
# p0_es$disturbance_mean_interval <- 30.0
# 
# 
# 
# p1_es <- expand_parameters(trait_matrix(c(0.045, 21.6, 0.02, 2), c("t_s", "a_s", "b_s", "height_0")), p0_es, mutant=FALSE)

species_list <- data.frame(a_s = c(0.10, 0.10, 0.30, 0.30)*365, t_s = c(0.66, 0.33, 0.66, 0.33))

p0 <- scm_base_parameters("ES20")
p0$control <- get_random_environment(p0$control, ceiling(p0$cohort_schedule_max_time), 0.75, 30/365)
p0$disturbance_mean_interval <- 10.0
p0$patch_area <- 1

p1 <- expand_parameters(trait_matrix(c(species_list$t_s[1], species_list$a_s[1], 0.1), c("t_s", "a_s","omega")), p0, mutant=FALSE)


pl_1_es <- ES20_Individual(s = p1$strategies[[1]])

pl_es <- pl_1_es # ES20_Plant()
# env_es <- ES20_fixed_environment(1.0)
env_es <- p1$environment

# FF16 plant 

# p0_ff <- scm_base_parameters("FF16")
# p0_ff$disturbance_mean_interval <- 30.0


# p1_ff <- expand_parameters(trait_matrix(c(0.1242302, 505), c("lma", "rho")), p0_ff)
# 
# pl_1_ff <- FF16_Individual(s = p1_ff$strategies[[1]])

# pl_ff <- pl_1_ff # ES20_Plant()
# env_ff <- FF16_fixed_environment(1.0)


# SINGLE YEAR TESTS
tt <- seq(0, 5, length.out=50)
 

# Run single year
# res_ns <- grow_plant_to_time_expanded(pl_ns, tt, env_ns)
res_es <- grow_plant_to_time_expanded(pl_es, tt, env_es)
# res_ff <- grow_plant_to_time_expanded(pl_ff, tt, env_ff)


# p1_es <- expand_parameters(trait_matrix(c(0.195, 21.6, 0.02, 1.5), c("t_s", "a_s", "b_s1", "height_0")), p0_es, mutant=FALSE)
# p1_es_2 <- expand_parameters(trait_matrix(c(ts_values[k], 21.9, h_0_values[i], b_s_values[j]), c("t_s", "a_s", "height_0", "b_s1")), p0_es, mutant=FALSE)
# pl_1_es <- ES20_Individual(s = p1_es_2$strategies[[1]])
# pl_es <- pl_1_es # ES20_Plant()
# env_es <- get_constant_environment_ES20(stress=0.75)
# tt <- seq(0, 1.1, length.out=50)
# res_es <- grow_plant_to_time_expanded(pl_es, tt, env_es)


# res_ns_df <- as.data.frame(res_ns$state)
res_es_df <- as.data.frame(res_es$state)
# res_ff_df <- as.data.frame(res_ff$state)


# res_ns_df <- cbind(res_ns_df, res_ns$aux_size)
res_es_df <- cbind(res_es_df, res_es$aux_size)


# res_ff_df <- cbind(res_ff_df, res_ff$aux_size)

# Fix the leaf area in original model 
# colnames(res_ff_df)[colnames(res_ff_df) %in% c("competition_effect")] <- "area_leaf"



# res_ns_df$model <- rep("NS", times=nrow(res_ns_df))
res_es_df$model <- rep("ES", times=nrow(res_es_df))
# res_ff_df$model <- rep("FF", times=nrow(res_ff_df))

# res_ns_df$time <- tt
res_es_df$time <- tt
# res_ff_df$time <- tt


res_mass_long_es <- pivot_longer(res_es_df, 
                                 c('mass_leaf', 'mass_sapwood', 'mass_heartwood', 'mass_bark', 'mass_storage', 'mass_root'), names_to = "mass_type", values_to = "mass_value" )
# res_mass_long_ns <- pivot_longer(res_ns_df,
                                 # c('mass_leaf', 'mass_sapwood', 'mass_heartwood', 'mass_bark', 'mass_storage', 'mass_root'), names_to = "mass_type", values_to = "mass_value" )


rates_mass_long_es <- pivot_longer(res_es_df, 
                                 c('dmass_leaf_darea_leaf', 'dmass_sapwood_darea_leaf', 'dmass_bark_darea_leaf', 'dmass_root_darea_leaf'), names_to = "mass_type", values_to = "mass_value" )

growth_rates_mass_long_es <- pivot_longer(res_es_df, 
                                   c('dmass_leaf_dt', 'dmass_sapwood_dt', 'dmass_bark_dt', 'dmass_root_dt'), names_to = "mass_type", values_to = "mass_value" )


growth_rates_mass_long_es_2 <- growth_rates_mass_long_es

growth_rates_mass_long_es_2$mass_value <- growth_rates_mass_long_es_2$mass_value / growth_rates_mass_long_es_2$dbiomass_dt
growth_rates_mass_long_es_2$mass_value[growth_rates_mass_long_es_2$mass_value == Inf] = 0

percent <- res_mass_long_es  %>%
  subset(mass_type %in% c('mass_leaf', 'mass_sapwood', 'mass_bark', 'mass_root')) %>%
  group_by(time, mass_type) %>%
  summarise(n = sum(mass_value)) %>%
  mutate(percentage = n / sum(n))


percent_2 <- res_mass_long_es  %>%
  subset(mass_type %in% c('dmass_leaf_dt', 'dmass_sapwood_dt', 'dmass_bark_dt', 'dmass_root_dt')) %>%
  group_by(time, mass_type) %>%
  summarise(n = sum(mass_value)) %>%
  mutate(percentage = n / sum(n))


percent_hd <- res_mass_long_es  %>%
  subset(mass_type %in% c('mass_leaf', 'mass_sapwood', 'mass_bark', 'mass_root', 'mass_heartwood')) %>%
  group_by(time, mass_type) %>%
  summarise(n = sum(mass_value)) %>%
  mutate(percentage = n / sum(n))

res_all <- rbind(res_ns_df, res_es_df)
res_all <- plyr::rbind.fill(res_all, res_ff_df)

res_all <- res_es_df
# Plot some results

## Firstly height

p <- ggplot(res_all, aes(x = time, y = height)) + 
  geom_line(aes(color = model)) +
  # geom_hline(yintercept = pl_1_es$strategy$hmat, linetype = "dashed") +
  scale_x_continuous(name = "Time (yr)") +
  scale_y_continuous(name = "Plant Height (m)") + 
  theme_linedraw()
p

## adjustments

p <- ggplot(res_all, aes(x = time, y = height_adjustment)) +
  geom_line(aes(color = model)) +
  scale_x_continuous(name = "Time (yr)") +
  scale_y_continuous(name = "Plant Height Adjustment (m)") +
  theme_linedraw()
p

# p <- ggplot(res_all, aes(x = time, y = mass_sapwood_adjustment)) + 
#   geom_line(aes(color = model)) +
#   scale_x_continuous(name = "Time (yr)") +
#   scale_y_continuous(name = "Plant Height (m)") + 
#   theme_linedraw()
# p
# 
# p <- ggplot(res_all, aes(x = time, y = mass_sapwood_difference)) + 
#   geom_line(aes(color = model)) +
#   scale_x_continuous(name = "Time (yr)") +
#   scale_y_continuous(name = "Plant Height (m)") + 
#   theme_linedraw()
# p
# 
# p <- ggplot(res_all, aes(x = time, y = height_difference)) + 
#   geom_line(aes(color = model)) +
#   scale_x_continuous(name = "Time (yr)") +
#   scale_y_continuous(name = "Plant Height (m)") + 
#   theme_linedraw()
# p
# 
# p <- ggplot(res_all, aes(x = mass_sapwood_difference, y = mass_sapwood_adjustment)) + 
#   geom_point(aes(color = model)) +
#   scale_x_continuous(name = "Time (yr)", limit = c(-2.5, 0)) +
#   scale_y_continuous(name = "Plant Height (m)") + 
#   theme_linedraw()
# p
# 
# p <- ggplot(res_all, aes(x = time, y = fecundity_dt_abs * (pl_1_es$strategy$omega + pl_1_es$strategy$a_f3))) + 
#   geom_line(aes(color = model)) +
#   geom_line(aes(y = mass_storage)) +
#   scale_x_continuous(name = "Time (yr)") +
#   scale_y_continuous(name = "Plant Height Adjustment (m)") + 
#   theme_linedraw()
# p

## Now the storage mass and proportion (should be on same graph at a later stage)


p <- ggplot(res_all, aes(x = time)) + 
  geom_line(aes(y = mass_storage, color = model)) +
  # geom_line(aes(y = (mass_storage/(mass_sapwood+mass_leaf+mass_bark+mass_root))/coeff, color=model), linetype = "dashed") +
  scale_x_continuous(name = "Time (yr)") +
  scale_y_continuous(name = "Mass Storage(kg)")+ 
  theme_linedraw()
p

p <- ggplot(res_all, aes(x = time)) + 
  geom_line(aes(y = (storage_portion), color=model), linetype = "dashed") +
  scale_x_continuous(name = "Time (yr)") +
  scale_y_continuous(name = "Storage Concentration")+
  theme_linedraw()
p


p <- ggplot(res_all, aes(x = time)) + 
  geom_line(aes(y = (mass_sapwood), color=model)) +
  scale_x_continuous(name = "Time (yr)") +
  scale_y_continuous(name = "Sapwood Mass kgC")+
  theme_linedraw()
p



## Now for the net biomass production

p <- ggplot(res_all, aes(x = time, y = net_mass_production_dt)) + 
  geom_line(aes(color = model)) +
  geom_line(aes(y = respiration_dt, color=model), linetype = "dashed") +
  scale_x_continuous(name = "Time (yr)") +
  scale_y_continuous(name = "Net Biomass Production + Respiration Rate (kgCyr-1)") + 
  theme_linedraw()

p



p <- ggplot(res_all, aes(x = time, y = net_mass_production_dt + respiration_dt, linetype=model)) + 
  geom_line() +
  geom_line(aes(y = area_leaf, linetype = model), color = "blue") +
  geom_line(aes(y = mass_sapwood, linetype =model), color = "green") +
  geom_line(aes(y = respiration_dt,linetype = model), color = "red") +
  # geom_line(aes(y = net_mass_production_dt), linetype = "dashed", color="blue") +
  # scale_x_continuous(name = "Time (yr)", limit = c(0,10)) +
  scale_y_continuous(name = "GPP (kgCyr-1)") + 
  theme_linedraw()

p


p <- ggplot(res_all, aes(x = time, y = dbiomass_dt, color = model)) + 
  geom_line() +
  scale_x_continuous(name = "Time (yr)", limits = c(15, 20)) +
  scale_y_continuous(name = "dBiomass_dt") + 
  theme_linedraw()
p


p <- ggplot(res_all, aes(x = time, y = dmass_leaf_darea_leaf, linetype=model)) + 
  geom_line() +
  scale_x_continuous(name = "Time (yr)") +
  scale_y_continuous(name = "dmass_leaf_darea_leaf") + 
  theme_linedraw()
p


max(growth_rates_mass_long_es_2$mass_value)



p <- ggplot(growth_rates_mass_long_es_2, aes(x = time, y = mass_value, fill = mass_type)) + 
  geom_area(alpha=0.5) + 
  scale_x_continuous(name = "Time (yr)") +
  scale_y_continuous(name = "Proportion Allocation") + 
  theme_linedraw()
p






## Let's do some stacked mass graphs

res_mass_long_es_sub <- subset(res_mass_long_es, mass_type != "mass_storage")

p <- ggplot(res_mass_long_es_sub, aes(x = time, y = mass_value, fill = mass_type)) + 
  geom_area(alpha=0.5) + 
  scale_x_continuous(name = "Time (yr)") +
  scale_y_continuous(name = "Plant Mass Pools kgC") + 
  theme_linedraw()
p

p <- ggplot(percent, aes(x = time, y = percentage * 100, fill = mass_type)) + 
  geom_area(alpha=0.5) + 
  scale_x_continuous(name = "Time (yr)") +
  scale_y_continuous(name = "Percentage total live biomass") + 
  theme_linedraw()
p

p <- ggplot(percent_hd, aes(x = time, y = percentage * 100, fill = mass_type)) + 
  geom_area(alpha=0.5) + 
  scale_x_continuous(name = "Time (yr)") +
  scale_y_continuous(name = "Percentage total biomass") + 
  theme_linedraw()
p

p <- ggplot(res_mass_long_es, aes(x = time, y = mass_value, fill = mass_type)) + 
  geom_area(alpha=0.5) + 
  scale_x_continuous(name = "Time (yr)") +
  scale_y_continuous(name = "Plant Mass Pools kgC") + 
  theme_linedraw()
p


res_mass_long_ns_sub <- subset(res_mass_long_ns, mass_type != "mass_storage")

p <- ggplot(res_mass_long_ns_sub, aes(x = time, y = mass_value, fill = mass_type)) + 
  geom_area(alpha=0.5) + 
  scale_x_continuous(name = "Time (yr)") +
  scale_y_continuous(name = "Plant Mass Pools kgC") + 
  theme_linedraw()
p


## Leaf area! 

p <- ggplot(res_all, aes(x=time, y = area_leaf, color = model)) + 
  geom_line() + 
  scale_x_continuous(name = "Time (yr)", breaks = 0:20) +
  scale_y_continuous(name = "Leaf Area (m2)") + 
  theme_linedraw()
p


## The components I'm struggling with: Fecundity and Mortality!

offset = 60
coeff = 0.2


storage_portion = diff(res_all$fecundity * (pl_1_es$strategy$omega+pl_1_es$strategy$a_f3))
time = res_all$time[1: (length(res_all$time) - 1)]

plot(time, storage_portion)


p <- ggplot(res_all, aes(x=time, y = fecundity * (pl_1_es$strategy$omega+pl_1_es$strategy$a_f3), color = model)) + 
  geom_line() + 
  # geom_line(aes(y = mortality-offset), linetype = "dashed") + 
  scale_x_continuous(name = "Time (yr)") +
  scale_y_continuous(name = "Fecundity") +
  theme_linedraw()
p









p <- ggplot(res_all, aes(x=time, y = mortality, linetype = model)) + 
  geom_line() + 
  scale_x_continuous(name = "Time (yr)") +
  scale_y_continuous(name = "Mortality (yr^-1)")+ 
  theme_linedraw()
p


p <- ggplot(res_es_df, aes(x=storage_portion, y = mortality)) + 
  geom_point(aes(color = time)) + 
  scale_x_continuous(name = "Storage Proportion") +
  scale_y_continuous(name = "Mortality (yr^-1)")+ 
  theme_linedraw()
p

p <- ggplot(res_all, aes(x=time, y = fecundity, color = model)) + 
  geom_line() + 
  scale_x_continuous(name = "Time (yr)") +
  scale_y_continuous(name = "Fecundity (?)")+ 
  theme_linedraw()
p

## stem diameter


p <- ggplot(res_all, aes(x=time, y = diameter_stem, color = model)) + 
  geom_line() + 
  scale_x_continuous(name = "Time (yr)") +
  scale_y_continuous(name = "Diameter (m)") + 
  theme_linedraw()
p


p <- ggplot(res_all, aes(x=height, y = diameter_stem, shape = model, color = model)) + 
  geom_point() + 
  scale_x_continuous(name = "Height (m)") +
  scale_y_continuous(name = "Diameter (m)") + 
  theme_linedraw()
p



p <- ggplot(res_all, aes(y=area_leaf, x = height, color = model)) + 
  geom_point() + 
  scale_y_continuous(name = "Area Leaf (m2)") +
  scale_x_continuous(name = "Height (m)") + 
  theme_linedraw()
p



## from BAAD DB static data (cause I'm lazy)

dbh <- c(0.109, 0.249, 0.235, 0.129, 0.128, 0.145, 0.128, 0.594, 0.400, 0.290,
         0.217, 0.384, 0.105, 0.123, 0.207, 0.127, 0.135, 0.192, 0.158, 1.690, 0.240,
         0.134, 0.276, 0.178, 0.110, 0.520, 0.238, 0.608, 0.400, 0.284, 0.185, 0.429,
         0.197, 0.181, 0.197, 0.135, 0.127, 0.670, 0.228, 0.173, 0.109, 0.111, 0.116,
         0.249, 0.123, 0.148, 0.184, 0.269, 0.102, 0.633, 1.183, 0.137, 0.258, 0.234)

height <- c(15.960, 28.558, 26.855, 18.980,
            13.590, 18.080, 14.620, 36.380, 34.533, 24.950, 21.730, 34.818, 15.980, 14.770,
            22.200, 10.940, 14.360, 16.910, 17.950, 49.117, 27.060, 15.940, 30.220, 17.320,
            18.160, 35.160, 26.270, 42.300, 36.087, 26.323, 22.710, 32.538, 24.120, 22.550,
            22.920, 22.650, 16.750, 36.068,     NA, 17.550, 13.560, 15.050, 16.680, 31.340,
            13.170, 15.850, 23.820, 23.850, 14.120, 36.303, 43.372, 18.460, 24.850, 26.100)

dba <- c(0.00835, 0.00840, 0.00961, 0.00784, 0.00935, 0.01013, 0.01162, 0.01081,
         0.01125, 0.01011, 0.01095, 0.01156, 0.01010, 0.01198, 0.01380, 0.00937,
         0.01167, 0.00886, 0.01408, 0.01317, 0.01209, 0.01390, 0.01372, 0.01512,
         0.01163, 0.01429, 0.01378, 0.01070, 0.01219, 0.01009, 0.00743, 0.01227,
         0.00955, 0.01187, 0.01160, 0.01057, 0.01159, 0.01395, 0.01324, 0.01031,
         0.01128, 0.01252, 0.01230, 0.01155, 0.00937, 0.01272, 0.01462, 0.01355,
         0.01387, 0.01551, 0.01535, 0.01573, 0.01452, 0.01255, 0.00967, 0.01049,
         0.01002, 0.00999, 0.00997, 0.00954, 0.00952, 0.00988, 0.00976, 0.01025,
         0.01093, 0.01046, 0.01287, 0.01477, 0.01060, 0.01304, 0.01154, 0.01046,
         0.01240, 0.01305, 0.01297, 0.01446, 0.01302, 0.01302, 0.01258, 0.01345,
         0.01569, 0.01111, 0.01116, 0.01088, 0.01057, 0.01111, 0.00995, 0.01016,
         0.00803, 0.01113, 0.01225, 0.01099, 0.01071, 0.01137, 0.01325, 0.01077,
         0.00984, 0.01183, 0.01303, 0.01300, 0.01378, 0.01314, 0.01355, 0.01442,
         0.01378, 0.01259, 0.01497, 0.01366, 0.00297, 0.00307, 0.00311, 0.00251,
         0.00338, 0.00257, 0.00339, 0.00239, 0.00266, 0.00158, 0.00232, 0.00268,
         0.00111, 0.00203, 0.00138, 0.00136, 0.00159, 0.00077, 0.00565, 0.00324,
         0.00495, 0.00417, 0.00406, 0.00339, 0.00319, 0.00370, 0.00286, 0.00283,
         0.00136, 0.00199, 0.00175, 0.00267, 0.00316, 0.00231, 0.00178, 0.00204,
         0.00394, 0.00516, 0.00378, 0.00398, 0.00362, 0.00412, 0.00359, 0.00239,
         0.00397, 0.00299, 0.00228, 0.00331, 0.00278, 0.00130, 0.00291, 0.00236,
         0.00189, 0.00213, 0.00594, 0.00455, 0.00709, 0.00538, 0.00516, 0.00472,
         0.00444, 0.00600, 0.00550, 0.00270, 0.00288, 0.00310, 0.00300, 0.00311,
         0.00171, 0.00328, 0.00273, 0.00176, 0.00191, 0.00392, 0.00490, 0.00315,
         0.00497, 0.00390, 0.00500, 0.00406, 0.00436, 0.00244, 0.00268, 0.00287,
         0.00198, 0.00221, 0.00260, 0.00183, 0.00166, 0.00218, 0.00715, 0.00739,
         0.00617, 0.00634, 0.00545, 0.00917, 0.00577, 0.00764, 0.00390, 0.00447,
         0.00388, 0.00533, 0.00678, 0.00422, 0.00366, 0.00416, 0.00407, 0.00508)

height_2 <- c(0.970,  1.050,  1.150,  1.120,  1.040,  1.030,  1.240,  1.090,  1.200,  0.990,
              1.120,  1.310,  1.060,  1.000,  1.210,  1.150,  1.130,  1.020, 1.400 , 1.260,
              1.200,  1.410,  1.180,  1.160,  1.270,  1.160 , 1.330,  1.140,  1.450 , 1.380,
              0.740,  1.290,  1.280,  1.430,  1.310,  1.240,  1.430,  1.610,  1.630 , 1.290,
              1.420,  1.610,  1.590,  1.490,  1.230,  1.300,  1.340,  1.520,  1.230 , 1.660,
              1.380,  1.400,  1.500,  1.500,  1.130,  0.950,  1.120,  1.100,  0.900 , 0.950,
              1.030,  1.060,  1.190,  1.030,  0.930,  1.090,  1.220,  1.240,  0.900 , 1.270,
              1.040,  0.950,  1.320,  1.060,  1.270,  1.370,  1.380,  1.290,  1.300 , 1.100,
              1.070,  1.240,  1.430,  1.420,  1.060,  1.350,  1.530,  1.550,  1.320 , 1.300,
              1.340,  1.270,  1.420,  1.380,  1.690,  1.270 , 1.470,  1.510,  1.480 , 1.400,
              1.150,  1.600,  1.630,  1.480,  1.700,  1.620,  1.590,  1.560,  0.233,  0.201,
              0.282,  0.176,  0.243,  0.164,  0.216,  0.237,  0.240,  0.240,  0.315,  0.184,
              0.084,  0.270,  0.137,  0.205,  0.266,  0.070,  0.349,  0.171,  0.295,  0.267,
              0.232,  0.225,  0.293,  0.235,  0.235,  0.354,  0.189,  0.330,  0.302,  0.335,
              0.371,  0.365,  0.271,  0.403,  0.276,  0.403,  0.280,  0.222,  0.270,  0.383,
              0.276,  0.277,  0.283,  0.480,  0.252,  0.415,  0.450,  0.071,  0.471,  0.304,
              0.300,  0.330,  0.539,  0.509,  0.546,  0.488,  0.480,  0.405,  0.330,  0.520,
              0.420,  0.375,  0.515,  0.400,  0.595,  0.640,  0.370,  0.610,  0.037,  0.230,
              0.235,  0.310,  0.390,  0.340,  0.395,  0.450,  0.460,  0.380,  0.420,  0.490,
              0.320,  0.450,  0.430,  0.390,  0.340,  0.320,  0.370 , 0.460 , 0.388,  0.548,
              0.525,  0.529,  0.423,  0.476,  0.480,  0.478,  0.620,  0.543 , 0.401,  0.624,
              0.575,  0.375,  0.485,  0.449,  0.653,  0.490)

euc_data <- data.frame(height = c(height, height_2), diameter_stem = c(dbh, dba))

p <- ggplot(res_all, aes(x=height, y = diameter_stem)) + 
  geom_point(aes(  color = model)) + 
  geom_point(aes(x = height, y=diameter_stem), data = euc_data) +
  scale_x_continuous(name = "Height (m)", limits = c(0, 2)) +
  scale_y_continuous(name = "Diameter (m)", limits = c(0, 0.02)) + 
  theme_linedraw()
p
