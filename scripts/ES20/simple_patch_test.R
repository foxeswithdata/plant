rm(list=ls())


devtools::load_all()


p0_ff <- scm_base_parameters("FF16")
p1_ff <- expand_parameters(trait_matrix(0.0825, "lma"),  p0_ff, mutant = FALSE)
p1_ff$seed_rain <- 1
out_ff <- run_scm(p1_ff)
out_ff$seed_rains 

p0_es <- scm_base_parameters("ES20")
p0_es$disturbance_mean_interval <- 3
p1_es <- expand_parameters(trait_matrix(0.0825, "lma"),  p0_es, mutant = FALSE)
p1_es$seed_rain <- 1
out <- run_scm(p1_es)
out$seed_rains 


cohort <- ES20_Cohort()

species <- ES20_Species()

patch <- ES20_Patch(p1_es)
