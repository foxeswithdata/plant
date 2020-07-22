rm(list = ls())
setwd("~/R/plant")

# Use latest version 
RcppR6::install(".")
devtools::document(".")
devtools::load_all()

#Test strategy
s <- FF16_Strategy()

s2 <- ES20_Strategy()

s3 <- ES20r_Strategy()

# Initialise parameters. 
p <- FF16_Parameters()
p2 <- FF16r_Parameters()


# This bit does not work yet. So it must be a problem on my part 
p3 <- ES20_Parameters()




out  <- run_stochastic_collect(p)
