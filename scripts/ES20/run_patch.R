rm(list=ls())

devtools::load_all()
# source("scripts/ES20/test_objects.R")


# DETERMINISTIC 

# Look at FF16 (Code example)

p0 <- scm_base_parameters("FF16")
p0$patch_area <- 1
# p0$control$equilibrium_nsteps <- 30
# p0$control$equilibrium_solver_name <- "hybrid"
# p0$control$equilibrium_verbose <- TRUE
# p0$disturbance_mean_interval <- 30.0


p1 <- expand_parameters(trait_matrix(c(0.1242302, 0.1), c("lma", "omega")), p0, mutant=FALSE)
p1 <- expand_parameters(trait_matrix(c(0.1242304, 0.1), c("lma", "omega")), p1, mutant=FALSE)
p1 <- expand_parameters(trait_matrix(c(0.1242306, 0.1), c("lma", "omega")), p1, mutant=FALSE)
p1 <- expand_parameters(trait_matrix(c(0.1242308, 0.1), c("lma", "omega")), p1, mutant=FALSE)

# why is this not working?

p1_eq <- equilibrium_seed_rain(p1)


p1_eq$seed_rain

scm <-  FF16_SCM(p1)


for(i in 1:length(p1$cohort_schedule_times_default)){
  cat(i)
  scm$run_next()
}

p1_eq$seed_rain

# p1$seed_rain <- 3.0

data1 <- run_scm_collect(p1_eq)

# Plot (from github instructions)

names(data1)

data1$species[[1]]



data1$seed_rain

data1$time
t <- data1$time
h <- data1$species[[1]]["height", , ]

dim(h)

# h[123, ]

# data1 <- out

matplot(t, h, lty=1, col=make_transparent("black", 0.25), type="l",
        las=1, xlab="Time (years)", ylab="Height (m)")

xlim <- c(0, 1.1)
ylim <- range(data1$env[[length(data1$env)]][, "height"])
plot(NA, xlim=xlim, ylim=ylim, las=1,
     xlab="Canopy openness", ylab="Height (m)")
for (i in data1$env) {
  lines(i[, "canopy_openness"], i[, "height"], col="grey")
}

blues <- c("#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6",
           "#4292C6", "#2171B5", "#08519C", "#08306B")
times <- c(5, 10, 20, 40, data1$time[[length(data1$time)]])
cols <- colorRampPalette(blues[-(1:2)])(length(times))
for (i in seq_along(times)) {
  x <- data1$env[[which.min(abs(times[[i]] - data1$time))]]
  lines(x[, "canopy_openness"], x[, "height"], col=cols[[i]])
  y <- x[nrow(x), "height"]
  points(1, y, pch=19, col=cols[[i]])
  text(1 + strwidth("x"), y, paste(round(times[[i]]), "years"),
       adj=c(0, 0))
}


# look at patches

patches1 <- lapply(seq_along(data1$time), scm_patch, data1)


lai1 <- sapply(patches1, function(x) x$area_leaf_above(0.0))
plot(data1$time, lai1, type="l", las=1, xlab="Time (years)",
     ylab="Leaf area index")

# multiple species

p2 <- expand_parameters(trait_matrix(0.2625, "lma"), p1, FALSE)
p2_eq <- equilibrium_seed_rain(p2)


data2 <- run_scm_collect(p2_eq)

t2 <- data2$time
h1 <- data2$species[[1]]["height", , ]
h2 <- data2$species[[2]]["height", , ]

cols <- c("#e34a33", "#045a8d")
matplot(t2, h1, lty=1, col=make_transparent(cols[[1]], .25), type="l",
        las=1, xlab="Time (years)", ylab="Height (m)")
matlines(t2, h2, lty=1, col=make_transparent(cols[[2]], .25))





## Look at ES20

p0 <- scm_base_parameters("ES20")
p0$disturbance_mean_interval <- 30.0

p1 <- expand_parameters(trait_matrix(c(0.66, 0.30 * 365, 0.1), c("t_s", "a_s", "omega")), p0, mutant=FALSE)

# p1$cohort_schedule_max_time
# p1$cohort_schedule_times_default
# p1$cohort_schedule_times_default[82]



p1_eq <- equilibrium_seed_rain(p1)

p1_eq$seed_rain

p1$seed_rain <- 20
scm <-  ES20_SCM(p1)


for(i in 1:length(p1$cohort_schedule_times_default)){
  cat(i)
  scm$run_next()
}

scm$time


scm$patch$environment$environment_interpolator$xy

scm$patch$species[[1]]$cohorts[[1]]$ode_rates

scm$ode_times

scm$run_next()





ode_times_before <- scm$ode_times
scm$ode_times



scm$run()

scm$
  scm$seed_rains


ode_times_after <- scm$ode_times
length(ode_times_after)

p1$seed_rain

data1 <- run_scm_collect(p1)




## Next we want to get the stochastic patch working
p0 <- scm_base_parameters("FF16")
p0$control$equilibrium_nsteps <- 30
p0$control$equilibrium_solver_name <- "iteration"
p0$control$equilibrium_verbose <- TRUE
p0$disturbance_mean_interval <- 30.0

p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0, FALSE)

p1$seed_rain <- 3
p1$patch_area <- 3

out  <- run_stochastic_collect(p1)

t <- out$time
h <- out$species[[1]]["height", , ]


matplot(t, h, lty=1, col=make_transparent("black", 0.25), type="l",
        las=1, xlab="Time (years)", ylab="Height (m)")



# what about ES20?

p0 <- scm_base_parameters("ES20")
p0$control$equilibrium_nsteps <- 30
# p0$control$stress_sd <- 0
p0$control$generate_stress <- TRUE
p0$control$equilibrium_solver_name <- "iteration"
p0$disturbance_mean_interval <- 1.0


p1 <- expand_parameters(trait_matrix(c(0.4, 21.9,0.4, 0.1), c("t_s", "a_s", "height_0", "b_s")), p0, FALSE)

p1$seed_rain <- 1

p2 <- expand_parameters(trait_matrix(c(0.5, 21.9,0.4, 0.1), c("t_s", "a_s", "height_0", "b_s")), p1, FALSE)

p2$patch_area <- 10

out  <- run_stochastic_collect(p2)

t <- out$time
h <- out$species[[1]]["height", , ]

matplot(t, h, lty=1, col=make_transparent("black", 0.5), type="l",
        las=1, xlab="Time (years)", ylab="Height (m)")
