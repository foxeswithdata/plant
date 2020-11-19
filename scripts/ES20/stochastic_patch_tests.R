rm(list=ls())
devtools::load_all()

library(ggplot2)

source('scripts/ES20/test_objects.R')
source('scripts/ES20/helper.R')
species_list <- data.frame(a_s = c(0.05, 0.05, 0.15, 0.15)*365, t_s = c(0.66, 0.33, 0.66, 0.33))


runs <- data.frame(total_time =c(), sim_length = c())

for(i in 1:1000){
  begin <- Sys.time()
  
  p0 <- scm_base_parameters("ES20")
  p0$control <- get_random_environment(p0$control, ceiling(p0$cohort_schedule_max_time), 0.75, 30/365)
  p0$disturbance_mean_interval <- 30.0
  
  p1 <- expand_parameters(trait_matrix(c(species_list$t_s[1], species_list$a_s[1], 0.000114 * 10), c("t_s", "a_s","omega")), p0, mutant=FALSE)
  
  p1$seed_rain <- 20
  
  ## RUN STOCHASTIC ONE BY ONE
  
  types <- extract_RcppR6_template_types(p1, "Parameters")
  obj <- do.call('StochasticPatchRunner', types)(p1)
  obj$schedule <- stochastic_schedule(p1)
  res <- list(collect_state(obj))
  
  tryCatch({
    while (!obj$complete) {
      obj$run_next()
      res <- c(res, list(collect_all(obj)))
    }
  },
  error = function(err){print("error")}, 
  finally = {
    end <- Sys.time()
    runs = rbind(runs, data.frame(total_time=(end-begin), sim_length = c(obj$time)))
  })
}

send_to_elisa("finished time estimation")





rm(list=ls())
devtools::load_all()

library(ggplot2)

source('scripts/ES20/test_objects.R')
source('scripts/ES20/helper.R')
species_list <- data.frame(a_s = c(0.05, 0.05, 0.15, 0.15)*365, t_s = c(0.66, 0.33, 0.66, 0.33))

begin <- Sys.time()
species_list <- data.frame(a_s = c(0.10, 0.10, 0.30, 0.30)*365, t_s = c(0.50, 0.33, 0.66, 0.50))

p0 <- scm_base_parameters("ES20")
p0$control$environment_light_max_depth <-16
p0$control <- get_random_environment(p0$control, ceiling(p0$cohort_schedule_max_time), 0.75, 30/365)
p0$control$environment_light_tol <- 0.01
p0$disturbance_mean_interval <- 2
p0$patch_area <- 100

p1 <- expand_parameters(trait_matrix(c(species_list$t_s[1], species_list$a_s[1], 0.1), c("t_s", "a_s","omega")), p0, mutant=FALSE)

p1$seed_rain <- 2

# out  <- run_stochastic_collect(p1)


## RUN STOCHASTIC ONE BY ONE

library(log4r)

optts_logfile = "scripts/ES20/competition_patch_test_logfile.txt"

optts_console_appender = console_appender(layout = default_log_layout())
optts_file_appender = file_appender(optts_logfile, append = TRUE, 
                                    layout = default_log_layout())

optts__logger <- log4r::logger(threshold = "INFO", 
                               appenders= list(optts_console_appender,optts_file_appender))

types <- extract_RcppR6_template_types(p1, "Parameters")
obj <- do.call('StochasticPatchRunner', types)(p1)
obj$schedule <- stochastic_schedule(p1)
res <- list(collect_state(obj))


i = 0

while (!obj$complete) {
  tryCatch({
  log4r::info(optts__logger, "\t\t\t******New Time")
  log4r::info(optts__logger, 
              paste0("time: ", obj$time))
  print(obj$time)
  obj$run_next()
  res <- c(res, list(collect_all(obj)))
  }, error=function(e){
    log4r::error(optts__logger, "Error_message")
    log4r::error(optts__logger, e)
    break;
  })
  # invisible(readline(prompt="Press [enter] to continue"))
}

obj$time
end <- Sys.time()
print("Total time: ")
print(end-begin)



out_final <- res[[i+1]]
out_final_2 <- res[[i]]




p <- ggplot(as.data.frame(out_final$env), aes(height, canopy_openness)) +
  geom_line() + 
  geom_vline(xintercept = out_final$species[[1]][1,], alpha=0.5, color = "blue")
p




save(res, file="scripts/ES20/results_stochastic_test_4_years_100m2.RData")


save(runs, file="scripts/ES20/time_est_single_species.RData")




