get_constant_environment_ES20 <- function(){
  ctrl <- equilibrium_verbose(fast_control())
  ctrl$stress_regime <- rep(1, times=3000)
  ctrl$generate_stress <- FALSE
  
  env_es <- ES20_fixed_environment(1.0)
  
  ES20_Environment__reset_stress(obj_ = env_es, new_stress_regime = ctrl$stress_regime)
  
  return(env_es)
}



get_non_stressed_plant_ES20 <- function(){
  p0_es <- scm_base_parameters("ES20")
  p0_es$disturbance_mean_interval <- 30.0
  
  p1_es <- expand_parameters(trait_matrix(c(1,1), c("a_s", "t_s")), p0_es, FALSE)
  
  return(pl_1_es <- ES20_Plant(s = p1_es$strategies[[1]]))
}

