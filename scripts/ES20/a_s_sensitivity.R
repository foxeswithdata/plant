rm(list = ls())

library(ggplot2)

# Use latest version 
RcppR6::install(".")
devtools::load_all()

source('scripts/ES20/test_objects.R')
source('scripts/ES20/notify_me.R')

ts_values <- c(0.33, 0.66)
tcrit_values <- c(0.75, 0.85)
h_0_values <- c(0.5, 1, 2, 5, 10, 15)
b_s_values <- c(0.1)

a_s_values <- seq(from=0.02, to=0.4, by=0.02)

out <- data.frame()

for(i in ts_values){
  for(j in b_s_values){
    for(k in tcrit_values){
      for(l in h_0_values){
        for(m in a_s_values){
          tryCatch({
            p0_es <- scm_base_parameters("ES20")
            p1_es <- expand_parameters(trait_matrix(c(i, 
                                                      (m * 365), 
                                                      l, 
                                                      j), 
                                                    c("t_s", 
                                                      "a_s", 
                                                      "height_0", 
                                                      "b_s1")), 
                                       p0_es, mutant=FALSE)
            pl_1_es <- ES20_Individual(s = p1_es$strategies[[1]])
            pl_es <- pl_1_es # ES20_Plant()
            env_es <- get_constant_environment_ES20(stress=k)
            
            tt <- seq(0, 1, length.out=20)
            
            # Run single year
            res_es <- grow_plant_to_time_expanded(pl_es, tt, env_es)
            res_es_df <- as.data.frame(res_es$state)
            first <- res_es_df[1,]
            final <- res_es_df[nrow(res_es_df),]
            first_aux <- res_es$aux_size[1,]
            final_aux <- res_es$aux_size[nrow(res_es_df),]
            
            delta_h <- final[,'height'] - first[,'height']
            delta_la <- final[,'area_leaf'] - first[,'area_leaf']
            delta_st <- final[,'mass_storage'] - first[,'mass_storage']
            delta_ma <- ((final[,'mass_leaf'] + final[,'mass_sapwood'] + final[,'mass_bark'] + final[,'mass_root']) 
                         - (first[,'mass_leaf'] + first[,'mass_sapwood'] + first[,'mass_bark'] + first[,'mass_root'])) 
            delta_st_c <- final_aux['storage_portion'] - first_aux['storage_portion']
            
            
            rgr_h <- delta_h/first[,'height']
            rgr_la <- delta_la/first[,'area_leaf']
            rgr_st <- delta_st/first[,'mass_storage']
            rgr_ma <- delta_ma / (first[,'mass_leaf'] + first[,'mass_sapwood'] + first[,'mass_bark'] + first[,'mass_root'])
            rgr_st_c <- delta_st_c/first_aux['storage_portion']
            
            
            dead <- ifelse(any(res_es$aux_size[,'dead_flag']>0), TRUE, FALSE)
            error <- FALSE
            
          }, error = function(e) {
            delta_h <- NA
            delta_la <- NA
            delta_st <- NA
            delta_ma <- NA
            delta_st_c <- NA
            
            rgr_h <- NA
            rgr_la <- NA
            rgr_st <- NA
            rgr_ma <- NA
            rgr_st_c <- NA
            error <- TRUE
            dead <- FALSE
          }, finally = {
            out_temp <- data.frame(
              height_0 = l,
              ts = i,
              bs = j,
              tcrit = k,
              as = m,
              DELTA_H = delta_h, 
              DELTA_LA = delta_la,
              DELTA_MST = delta_st,
              DELTA_MA = delta_ma,
              DELTA_CST = delta_st_c,
              RGR_H = rgr_h, 
              RGR_LA = rgr_la,
              RGR_MST = rgr_st,
              RGR_MA = rgr_ma,
              RGR_CST = rgr_st_c,
              died = dead,
              error = error
            )
            
            out <- rbind(out, out_temp)
          })
          
          
        }
      }
    }
  }
}


save(out, h_0_values, b_s_values, a_s_values, tcrit_values, ts_values, file="scripts/ES20/as_sensitivity.RData")
send_to_elisa("Finished getting a_s sensitivity values :-)")




any(out$died)
out_as_01 <- subset(out, bs == 0.1)


tcrit.labs <- c("tcrit: 0.75y", "tcrit: 0.85y")
names(tcrit.labs) <- c("0.75", "0.85")

# New facet label names for supp variable
ts.labs <- c("ts: 0.33y", "ts: 0.66y")
names(ts.labs) <- c("0.33", "0.66")


p <- ggplot(out_as_01, aes(x = as, y = RGR_H, color = as.factor(height_0))) +
  geom_point() + 
  facet_grid( ts ~ tcrit, 
              labeller = labeller(ts = ts.labs, tcrit = tcrit.labs)) + 
  scale_x_continuous(expression(paste(alpha[s], ' ', kgCkgC^-1, d^-1))) +
  scale_y_continuous(expression(paste(RI[h], ' ', kgCkgC^-1))) +
  scale_color_discrete(name = "Initial Height [m]") +
  theme_bw()
p


p <- ggplot(out_as_01, aes(x = as, y = RGR_LA, color = as.factor(height_0))) +
  geom_point() + 
  facet_grid( ts ~ tcrit, 
              labeller = labeller(ts = ts.labs, tcrit = tcrit.labs)) + 
  scale_x_continuous(expression(paste(alpha[s], ' ', kgCkgC^-1, d^-1))) +
  scale_y_continuous(expression(paste(RI[la], ' ', kgCkgC^-1))) +
  scale_color_discrete(name = "Initial Height [m]") +
  theme_bw()
p

p <- ggplot(out_as_01, aes(x = as, y = RGR_MA, color = as.factor(height_0))) +
  geom_point() + 
  facet_grid( ts ~ tcrit, 
              labeller = labeller(ts = ts.labs, tcrit = tcrit.labs)) + 
  scale_x_continuous(expression(paste(alpha[s], ' ', kgCkgC^-1, d^-1))) +
  scale_y_continuous(expression(paste(RI[ma], ' ', kgCkgC^-1))) +
  scale_color_discrete(name = "Initial Height [m]") +
  theme_bw()
p

p <- ggplot(out_as_01, aes(x = as, y = RGR_MST, color = as.factor(height_0))) +
  geom_point() + 
  facet_grid( ts ~ tcrit, 
              labeller = labeller(ts = ts.labs, tcrit = tcrit.labs)) + 
  scale_x_continuous(expression(paste(alpha[s], ' ', kgCkgC^-1, d^-1))) +
  scale_y_continuous(expression(paste(RI[mst], ' ', kgCkgC^-1))) +
  scale_color_discrete(name = "Initial Height [m]") +
  theme_bw()
p





any(out$RGR_MST < 0)





