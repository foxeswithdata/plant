rm(list = ls())
setwd("~/R/plant")

# Use latest version 
devtools::load_all()

#Test strategy
s <- FF16_Strategy()


# Initialise parameters. This part doesnt work yet
p <- FF16_Parameters()

