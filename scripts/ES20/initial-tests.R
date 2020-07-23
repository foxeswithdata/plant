.rs.restartR()
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

p0 <- scm_base_parameters("ES20")

class(p3)

out  <- run_stochastic_collect(p3)
