devtools::load_all()

# testing passing control stress details. 

# pass the information about ctrl stress regime (mean etc) compare if new details are generated. 

#original 

ctrl <- equilibrium_verbose(fast_control())
p <- Parameters("ES20", "ES20_Env")(patch_area=1.0, control=ctrl, hyperpar=hyperpar("ES20"))
p$environment


print("ORIGINAL")
head(p$environment$stress_regime)

# change the mean of the stress to 0.1 and sd to 0.000001

ctrl <- equilibrium_verbose(fast_control())
ctrl$stress_mean <- 0.1

p <- Parameters("ES20", "ES20_Env")(patch_area=1.0, control=ctrl, hyperpar=hyperpar("ES20"))

ES20_Environment__reset_stress_random(obj_=p$environment, new_mean=0.1, new_sd = 0.000001)

print("CHANGING MEAN")
head(p$environment$stress_regime)



# change the actual stress regime

ctrl <- equilibrium_verbose(fast_control())
ctrl$stress_regime <- rep(5, times=3000)
ctrl$generate_stress <- FALSE

p <- Parameters("ES20", "ES20_Env")(patch_area=1.0, control=ctrl, hyperpar=hyperpar("ES20"))

ES20_Environment__reset_stress(obj_ = p$environment, new_stress_regime = ctrl$stress_regime)

print("ADDING STRESS REGIME MANUALLY")
head(p$environment$stress_regime)
