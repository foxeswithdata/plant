
rm(list=ls())
devtools::load_all()

library(ggplot2)

source('scripts/ES20/test_objects.R')
source('scripts/ES20/helper.R')

begin <- Sys.time()
species_list <- data.frame(a_s = c(0.10, 0.10, 0.30, 0.30)*365, t_s = c(0.50, 0.33, 0.66, 0.50))

p0 <- scm_base_parameters("ES20")
p0$control$environment_light_max_depth <-16
p0$control <- get_random_environment(p0$control, ceiling(p0$cohort_schedule_max_time), 0.75, 30/365)
p0$control$environment_light_tol <- 0.01
p0$disturbance_mean_interval <- 2 #short interval only
p0$patch_area <- 100

p1 <- expand_parameters(trait_matrix(c(species_list$t_s[1], species_list$a_s[1], 0.1, 0.2), c("t_s", "a_s","omega", "a_dG1")), p0, mutant=FALSE)

p1$seed_rain <- 2

## RUN STOCHASTIC ONE BY ONE


types <- extract_RcppR6_template_types(p1, "Parameters")
obj <- do.call('StochasticPatchRunner', types)(p1)
obj$schedule <- stochastic_schedule(p1)
res <- list(collect_state(obj))


while (!obj$complete) {
  print(obj$time)
  obj$run_next()
  res <- c(res, list(collect_all(obj)))
  # invisible(readline(prompt="Press [enter] to continue"))
}

obj$time
end <- Sys.time()
print("Total time: ")
print(end-begin)


##### For the figures

# 1) Environment

final_env <- as.data.frame(res[[length(res)]]$env)

plot(final_env$height, final_env$canopy_openness, type = "l")
abline(v=0.8375843, lty = 3) #add height 0 


# height of live trees
final_trees <- as.data.frame(t(res[[length(res)]]$species[[1]]))
final_trees_live <- subset(final_trees, is_alive == TRUE)
final_trees_live$height


# 2) Height growth throughout simulation

tree_data <- extract_trees(res)

live <- subset(tree_data, is_alive==TRUE)
live$tree_id <- as.factor(live$tree_id)

p <- ggplot(live, aes(x = time, y = height, color = tree_id)) +
  geom_line() + 
  scale_x_continuous("Time [yr]") + 
  theme(legend.position = "none")
p