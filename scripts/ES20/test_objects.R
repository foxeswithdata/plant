get_constant_environment_ES20 <- function(stress=1){
  ctrl <- equilibrium_verbose(fast_control())
  stress_regime <- rep(stress, times=3000) + 1:3000 -1
  
  ctrl$stress_regime <- stress_regime
  ctrl$generate_stress <- FALSE
  
  env_es <- ES20_fixed_environment(1.0)
  
  ES20_Environment__reset_stress(obj_ = env_es, new_stress_regime = ctrl$stress_regime)
  
  return(env_es)
}



get_non_stressed_plant_ES20 <- function(){
  p0_es <- scm_base_parameters("ES20")
  p0_es$disturbance_mean_interval <- 30.0
  
  p1_es <- expand_parameters(trait_matrix(c(365,1, 0.2), c("a_s", "t_s", "height_0")), p0_es)
  
  return(pl_1_es <- ES20_Individual(s = p1_es$strategies[[1]]))
}


get_random_environment <- function(control, max_time, stress_mean, stress_sd){
  stress_regime <- get_stress_regime(max_time, stress_mean, stress_sd)
  
  return(assign_random_environment(control, stress_regime))
}

assign_random_environment <- function(control, stress_regime){
  control$stress_regime <- stress_regime
  control$generate_stress <- FALSE
  return(control)
}


get_stress_regime <- function(max_time, stress_mean, stress_sd){
  stress_regime <- rnorm(max_time, stress_mean, stress_sd) + 1:max_time -1
  return(stress_regime)
}

