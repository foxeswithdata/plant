# Looking at mortality relationship 

a_dG1 <- pl_es$strategy$a_dG1
a_dG2 <- pl_es$strategy$a_dG2

a_dG1 <- 1
a_dG2 <- 80

prop_storage <- seq(from=0, to=0.6, by=0.01)

mort_rate <- a_dG1 * exp(-a_dG2 * prop_storage)

mortality1 <- data.frame(storage = prop_storage, mortality_rate = mort_rate, a_dg2 = rep(a_dG2, times=length(mort_rate)))

a_dG1 <- 1
a_dG2 <- 40

prop_storage <- seq(from=0, to=0.6, by=0.01)

mort_rate <- a_dG1 * exp(-a_dG2 * prop_storage)

mortality2 <- data.frame(storage = prop_storage, mortality_rate = mort_rate, a_dg2 = rep(a_dG2, times=length(mort_rate)))

a_dG1 <- 1
a_dG2 <- 20

prop_storage <- seq(from=0, to=0.6, by=0.01)

mort_rate <- a_dG1 * exp(-a_dG2 * prop_storage)

mortality3 <- data.frame(storage = prop_storage, mortality_rate = mort_rate, a_dg2 = rep(a_dG2, times=length(mort_rate)))


a_dG1 <- 1
a_dG2 <- 160

prop_storage <- seq(from=0, to=0.6, by=0.01)

mort_rate <- a_dG1 * exp(-a_dG2 * prop_storage)

mortality4 <- data.frame(storage = prop_storage, mortality_rate = mort_rate, a_dg2 = rep(a_dG2, times=length(mort_rate)))


mortality <- rbind(mortality1, mortality2, mortality3, mortality4)

mortality$adg2 <- as.factor(mortality$a_dg2)


p <- ggplot(mortality, aes(x = storage, y=mortality_rate, color=adg2)) + 
  geom_line() +
  geom_vline(xintercept=0, linetype="dashed") +
  scale_x_continuous(name = "Storage Proportion") +
  scale_y_continuous(name = "Mortality Rate (yr-1)") + 
  theme_linedraw()
p




