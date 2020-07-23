devtools::load_all()

stochastic_parameter_base <- function(){
  ctrl <- equilibrium_verbose(fast_control())
  ctrl$schedule_eps <- 0.005
  ctrl$equilibrium_eps <- 1e-3
  
  p3 <- Parameters("ES20", "ES20_Env")(patch_area=1.0, control=ctrl, hyperpar=hyperpar("ES20"))
}

