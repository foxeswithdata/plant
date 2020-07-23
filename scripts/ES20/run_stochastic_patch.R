devtools::load_all()

ctrl <- equilibrium_verbose(fast_control())
p <- Parameters("FF16", "FF16_Env")(patch_area=1.0, control=ctrl, hyperpar=hyperpar("FF16"))

p <- FF16_Parameters()


out  <- run_stochastic_collect(p)


