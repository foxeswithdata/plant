rm(list = ls())
gc()
setwd("~/R/plant")

# Use latest version 
RcppR6::install(".")
# devtools::document(".")
devtools::load_all()

#Test strategy
s2 <- ES20_Strategy()

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



