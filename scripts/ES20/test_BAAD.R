install.packages("devtools")
devtools::install_github("richfitz/datastorr")
devtools::install_github("traitecoevo/baad.data")

rm(list=ls())


## load baad 

baad <- baad.data::baad_data()

baad_dic <- baad$dictionary


all_eucs <- baad$data[grep("eucalypt.*", baad$data$species, ignore.case = TRUE),]

tree_with_sapwood_mass <- subset(baad$data, !is.na(m.ss))
tree_with_heartwood_mass <- subset(baad$data, !is.na(m.sh))

unique(tree_with_sapwood_mass$species)
unique(tree_with_heartwood_mass$species)


euc_saligna <- baad$data[grep("euc.* saligna", baad$data$species, ignore.case=TRUE),]
euc_grandis <- baad$data[grep("euc.* grandis", baad$data$species, ignore.case=TRUE),]
euc_pilularis <- baad$data[grep("euc.* pilularis", baad$data$species, ignore.case=TRUE),]
euc_urophylla <- baad$data[grep("euc.* urophylla", baad$data$species, ignore.case=TRUE),]
euc_nitens <- baad$data[grep("euc.* nitens", baad$data$species, ignore.case=TRUE),]


euc_data <- rbind(euc_saligna, euc_grandis, euc_pilularis, euc_urophylla, euc_nitens)

rm(euc_saligna, euc_grandis, euc_pilularis, euc_urophylla, euc_nitens, baad)




#load packages (tidyverse seems to mask something for baad so doing it later)

library(tidyverse)
# Use latest version 
devtools::load_all()



p0_es <- scm_base_parameters("ES20")
p0_es$disturbance_mean_interval <- 30.0

p1_es <- expand_parameters(trait_matrix(c(0.60, 21.9, 0.2), c("t_s", "a_s", "height_0")), p0_es, FALSE)
pl_1_es <- ES20_Plant(s = p1_es$strategies[[1]])

pl_es <- pl_1_es # ES20_Plant()
env_es <- ES20_fixed_environment(1.0)


# SINGLE YEAR TESTS
tt <- seq(0, 15, length.out=100)


# Run single year
res_es <- grow_plant_to_time(pl_es, tt, env_es)

res_es_df <- as.data.frame(res_es$state)
res_es_df <- cbind(res_es_df, res_es$aux_size)

res_es_df$model <- rep("ES", times=nrow(res_es_df))

res_es_df$time <- tt


res_mass_long_es <- pivot_longer(res_es_df, 
                                 c('mass_leaf', 'mass_sapwood', 'mass_heartwood', 'mass_bark', 'mass_storage', 'mass_root'), names_to = "mass_type", values_to = "mass_value" )

percent <- res_mass_long_es  %>%
  subset(mass_type %in% c('mass_leaf', 'mass_sapwood', 'mass_bark', 'mass_root')) %>%
  group_by(time, mass_type) %>%
  summarise(n = sum(mass_value)) %>%
  mutate(percentage = n / sum(n))

percent_hd <- res_mass_long_es  %>%
  subset(mass_type %in% c('mass_leaf', 'mass_sapwood', 'mass_bark', 'mass_root', 'mass_heartwood')) %>%
  group_by(time, mass_type) %>%
  summarise(n = sum(mass_value)) %>%
  mutate(percentage = n / sum(n))

res_all <- res_es_df



## look at leaf area allometric realationships (i.e. leaf area versus mass)

euc_data$species

## leaf mass 

p <- ggplot(euc_data, aes(x=a.lf, y=m.lf)) + 
  geom_point(alpha=0.25, aes(shape=species)) + 
  geom_point(aes(x=area_leaf, y = mass_leaf), color="green" ,data=res_all) + 
  scale_x_continuous(name = "Leaf Area (m2)") +
  scale_y_continuous(name = "Leaf mass (kg)") + 
  theme_linedraw()
p


p <- ggplot(euc_data, aes(x=a.lf, y=m.lf)) + 
  geom_point(alpha=0.25, aes(shape=species)) + 
  geom_point(aes(x=area_leaf, y = mass_leaf, color=time) , alpha=0.5, data=res_all) + 
  scale_x_continuous(name = "Leaf Area (m2)", limits = c(0, 7.5)) +
  scale_y_continuous(name = "Leaf mass (kg)", limits = c(0, 1.5)) + 
  theme_linedraw()
p



linearMod <- lm(m.lf ~ a.lf, data=euc_data)

summary(linearMod)

## sapwood mass -- no data

p <- ggplot(euc_data, aes(x=a.lf, y=m.ss)) + 
  geom_point(alpha=0.25) + 
  geom_point(aes(x=area_leaf, y = mass_sapwood), color="green" ,data=res_all) + 
  scale_x_continuous(name = "Leaf Area (m2)") +
  scale_y_continuous(name = "Leaf mass (kg)") + 
  theme_linedraw()
p

## heartwood mass -- no data

p <- ggplot(euc_data, aes(x=a.lf, y=m.sh)) + 
  geom_point(alpha=0.25) + 
  geom_point(aes(x=area_leaf, y = mass_heartwood), color="green" ,data=res_all) + 
  scale_x_continuous(name = "Leaf Area (m2)") +
  scale_y_continuous(name = "Leaf mass (kg)") + 
  theme_linedraw()
p

## bark mass 

p <- ggplot(euc_data, aes(x=a.lf, y=m.sb)) + 
  geom_point(alpha=0.25) + 
  geom_point(aes(x=area_leaf, y = mass_bark), color="blue" ,data=res_all) + 
  scale_x_continuous(name = "Leaf Area (m2)") +
  scale_y_continuous(name = "Bark mass (kg)", limits = c(0, 20)) + 
  theme_linedraw()
p


## fine root mass

p <- ggplot(euc_data, aes(x=a.lf, y=m.rf)) + 
  geom_point(alpha=0.25) + 
  geom_point(aes(x=area_leaf, y = mass_root), color="blue", alpha=0.55,data=res_all) + 
  scale_x_continuous(name = "Leaf Area (m2)", limits=c(0, 11)) +
  scale_y_continuous(name = "Fine Root mass (kg)", limits =c(0, 0.4)) + 
  theme_linedraw()
p





##### STEM AREAS

## sapwood area

p <- ggplot(euc_data, aes(x=a.lf, shape=species)) + 
  geom_point(alpha=0.25, aes(y=a.ssba), color = "green") + 
  geom_point(alpha=0.25, aes(y=a.ssbh), color = "blue") + 
  geom_point(alpha=0.25, aes(y=a.ssbc), color = "red") +              
  geom_point(data=res_all, aes(x=area_leaf, y = area_sapwood), shape = 1, color="black", alpha=0.55) + 
  scale_x_continuous(name = "Leaf Area (m2)") +
  scale_y_continuous(name = "Sapwood area (m2)") + 
  theme_linedraw()
p


p <- ggplot(euc_data, aes(x=a.lf, shape=species)) + 
  geom_point(alpha=0.25, aes(y=a.ssba), color = "green") + 
  geom_point(alpha=0.25, aes(y=a.ssbh), color = "blue") + 
  geom_point(alpha=0.25, aes(y=a.ssbc), color = "red") +              
  geom_point(data=res_all, aes(x=area_leaf, y = area_sapwood), shape = 1, color="black", alpha=0.55) + 
  scale_x_continuous(name = "Leaf Area (m2)", limits = c(0, 10)) +
  scale_y_continuous(name = "Sapwood area (m2)", limits = c(0, 0.003)) + 
  theme_linedraw()
p

## Heartwood

p <- ggplot(euc_data, aes(x=a.lf, shape=species)) + 
  geom_point(alpha=0.25, aes(y=a.shba), color = "green") + 
  geom_point(alpha=0.25, aes(y=a.shbh), color = "blue") + 
  geom_point(alpha=0.25, aes(y=a.shbc), color = "red") +              
  geom_point(data=res_all, aes(x=area_leaf, y = area_heartwood), shape = 1, color="black", alpha=0.55) + 
  scale_x_continuous(name = "Leaf Area (m2)") +
  scale_y_continuous(name = "Heartwood area (m2)") + 
  theme_linedraw()
p

## Bark

p <- ggplot(euc_data, aes(x=a.lf, shape=species)) +  
  geom_point(alpha=0.25, aes(y=a.sbbh), color = "blue") +
  geom_point(data=res_all, aes(x=area_leaf, y = area_bark), shape = 1, color="black", alpha=0.55) + 
  scale_x_continuous(name = "Leaf Area (m2)") +
  scale_y_continuous(name = "Bark area (m2)") + 
  theme_linedraw()
p

# Stem

p <- ggplot(euc_data, aes(x=a.lf, shape=species)) + 
  geom_point(alpha=0.25, aes(y=a.stba), color = "green") + 
  geom_point(alpha=0.25, aes(y=a.stbh), color = "blue") + 
  geom_point(alpha=0.25, aes(y=a.stbc), color = "red") +              
  geom_point(data=res_all, aes(x=area_leaf, y = area_heartwood), shape = 1, color="black", alpha=0.55) + 
  scale_x_continuous(name = "Leaf Area (m2)") +
  scale_y_continuous(name = "Stem area (m2)", limits= c(0, 0.05)) + 
  theme_linedraw()
p

# Determining different parameters in the model

mean(euc_data$r.st, na.rm=TRUE)

p <- euc_data %>%
  ggplot( aes(x=r.st)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  # scale_fill_manual(values=c("#69b3a2", "#404080")) +
  # theme_ipsum() +
  labs(fill="")
p












