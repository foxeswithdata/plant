devtools::load_all()

x <- seq(0, 20, length.out =1000)
plot(x, 0.5 + 0.5* sin(x/2/3.14))

p0 <- scm_base_parameters("ES20")
p0$disturbance_mean_interval <- 30.0

p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0, FALSE)
p1$seed_rain <- 20
#p1_eq <- equilibrium_seed_rain(p1)

p2 <- build_schedule(p1)
data1 <- run_scm_collect(p2)

matplot(data1$time, data1$species[[1]]["height", , ], lty=1, col=make_transparent("black", 0.25), type="l",
        las=1, xlab="Time (years)", ylab="Height (m)")


    seed_rain_out <- attr(p_new, "seed_rain_out", exact=TRUE)


Adjust relative accuracy alogorithm:

p0$control$ode_tol_rel <- 1e-5
p0$control$ode_tol_abs <- 1e-5

Adjust cohort spacing timestamp

p0$cohort_schedule_times_default <- 
