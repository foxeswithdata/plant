library(ggplot2)
library(RColorBrewer)

# load("scripts/ES20/results_stochastic_test_100_years.RData")

tree_data <- extract_trees(res)

trees_final <- subset(tree_data, time == max(tree_data$time))

live_v2 <- subset(tree_data, is_alive==TRUE)
live <- subset(live_v2, dead_flag != 1)

live$tree_id <- as.factor(live$tree_id)
live_v2$tree_id <- as.factor(live_v2$tree_id)

for(i in 1:ncol(tree_data)){
  if(any(is.na(tree_data[,i]))){
    print(i)
  }
}


# any

p <- ggplot(live, aes(x = time, y = height, color = tree_id)) +
  geom_line() + 
  scale_x_continuous("Time [yr]") + 
  theme(legend.position = "none")
p

p <- ggplot(live_v2, aes(x = time, y = height, color = tree_id)) +
  geom_line() + 
  # scale_x_continuous("Time [yr]", breaks = 0:35, limit = c(1, 1.2)) + 
  theme(legend.position = "none")
p


p <- ggplot(tree_data, aes(x = time, y = mortality, color = as.factor(tree_id))) +
  geom_point() + 
  # scale_x_continuous("Time [yr]", breaks = 0:35, limit = c(1, 1.2)) + 
  theme(legend.position = "none")
p




p <- ggplot(tree_data, aes(x = productivity_area, y = (mortality_growth_dependent_dt), color=stress)) +
  geom_point() +
  scale_color_gradientn(colours = brewer.pal(n = 8, name = 'PiYG'))
  # scale_x_continuous("Time [yr]") + 
  # theme(legend.position = "none")
p

any(tree_data$net_mass_production_dt < 0)
min(tree_data$net_mass_production_dt)


p <- ggplot(tree_data, aes(x = storage_portion, y = mortality_storage_dependent_dt, color = area_leaf)) +
  geom_point() +
  scale_color_gradientn(colours = brewer.pal(n = 8, name = 'PiYG'))
# scale_x_continuous("Time [yr]") + 
# theme(legend.position = "none")
p


p <- ggplot(tree_data, aes(x = productivity_area, y = mortality, color = storage_portion)) +
  geom_point() +
  scale_color_gradientn(colours = brewer.pal(n = 8, name = 'PiYG'))
# scale_x_continuous("Time [yr]") + 
# theme(legend.position = "none")
p



p <- ggplot(live_v2, aes(x = time, y = stress)) +
  geom_point() 
# scale_x_continuous("Time [yr]") + 
# theme(legend.position = "none")
p




p <- ggplot(live_v2, aes(x = time, y = mass_storage, color = tree_id)) +
  geom_line() + 
  scale_x_continuous("Time [yr]", breaks = 0:35) + 
  theme(legend.position = "none")
p


any(live$mass_storage < 0)
any(tree_data$mass_storage < 0, na.rm = TRUE)


min(tree_data$height, na.rm = TRUE)

max(as.numeric(live_v2$tree_id))



#### Examine last elements and whats going on 

trees_final <- subset(tree_data, time == max(tree_data$time))

trees_final_is_alive <- subset(trees_final, is_alive ==TRUE)

trees_at_end <- subset(tree_data, tree_id %in% trees_final_is_alive$tree_id)
trees_at_end$tree_id <- as.factor(trees_at_end$tree_id)

nrow(trees_at_end)


p <- ggplot(trees_at_end, aes(x = time, y = height, color = tree_id)) +
  geom_line() + 
  geom_line(data = data.frame(times = times, stress = (0.84 + stress/1000)), mapping = aes(times, stress), color = "black", alpha = 0.5) +
  scale_x_continuous("Time [yr]") + 
  theme(legend.position = "none")
p


p <- ggplot(trees_at_end, aes(x = time, y = dbiomass_dt, color = tree_id)) +
  geom_line() +
  geom_line(data = data.frame(times = times, stress = (stress/10)), mapping = aes(times, stress), color = "black", alpha = 0.5) +
  scale_x_continuous("Time [yr]") + 
  theme(legend.position = "none")
p

p <- ggplot(trees_at_end, aes(x = time, y = mass_storage, color = tree_id)) +
  geom_line() +
  scale_x_continuous("Time [yr]") + 
  theme(legend.position = "none")
p

p <- ggplot(trees_at_end, aes(x = time, y = storage_portion, color = tree_id)) +
  geom_line() +
  scale_x_continuous("Time [yr]") + 
  theme(legend.position = "none")
p

p <- ggplot(trees_at_end, aes(x = time, y = area_leaf, color = tree_id)) +
  geom_line() +
  scale_x_continuous("Time [yr]") + 
  theme(legend.position = "none")
p

p <- ggplot(trees_at_end, aes(x = time, y = dbiomass_dt/net_mass_production_dt, color = tree_id)) +
  geom_line() +
  scale_x_continuous("Time [yr]") + 
  scale_y_continuous(limits = c(-3,3)) +
  theme(legend.position = "none")
p


p1$environment$stress_regime


