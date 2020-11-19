rm(list=ls())
devtools::load_all()

species_list <- data.frame(a_s = c(0.05, 0.05, 0.15, 0.15)*365, t_s = c(0.66, 0.33, 0.66, 0.33))

p0 <- scm_base_parameters("ES20")
p0$disturbance_mean_interval <- 10.0

p0$cohort_schedule_max_time

p1 <- expand_parameters(trait_matrix(c(species_list$t_s[1], species_list$a_s[1]), c("t_s", "a_s")), p0, mutant=FALSE)

# p1$cohort_schedule_max_time
# p1$cohort_schedule_times_default
# p1$cohort_schedule_times_default[82]

p1$seed_rain <- 20
scm <-  ES20_SCM(p1)

for(i in 1:length(p1$cohort_schedule_times_default)){
  cat(i)
  scm$run_next()
}

scm$time
