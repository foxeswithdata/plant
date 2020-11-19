# rm(list = ls())
# .rs.restartR()

rm(list = ls())
# library(deSolve)
# library(tidyverse)

# Use latest version 
RcppR6::install(".")
devtools::load_all

library(log4r)

source('scripts/ES20/test_objects.R')
source('scripts/ES20/notify_me.R')

begin <- Sys.time()

ts_values <- seq(from=0, to=0.75, by=0.01)
h_0_values_1 <- seq(from=0.2, to=2, by = 0.2)
h_0_values_2 <- seq(from=2.5, to=10, by=0.5)
h_0_values_3 <- seq(from=15, to=25, by=5)
h_0_values <- c(h_0_values_1, h_0_values_2, h_0_values_3)

b_s_values_1 <- seq(from=0.02, to=0.2, by=0.02)
b_s_values_2 <- seq(from=0.25, to=0.5, by=0.05)
b_s_values <- c(b_s_values_1, b_s_values_2)

opt_ts_M <- matrix(0, nrow=length(h_0_values), ncol = length(b_s_values))
opt_ts_S <- matrix(0, nrow=length(h_0_values), ncol = length(b_s_values))



optts_logfile = "scripts/ES20/opt_ts_as_010_tcrit_075_logfile.txt"

optts_console_appender = console_appender(layout = default_log_layout())
optts_file_appender = file_appender(optts_logfile, append = TRUE, 
                                 layout = default_log_layout())

optts__logger <- log4r::logger(threshold = "INFO", 
                           appenders= list(optts_console_appender,optts_file_appender))



for(i in 1:length(h_0_values)){
  for(j in 1:length(b_s_values)){
    out_Mass <- vector()
    out_Store <- vector()
    valid <- vector()


    for(k in 1:length(ts_values)){
      # Run single year
      
      tryCatch({
        p0_es <- scm_base_parameters("ES20")
        p0_es$disturbance_mean_interval <- 30.0
        
        p1_es <- expand_parameters(trait_matrix(c(ts_values[k], (0.10 * 365), h_0_values[i], b_s_values[j]), c("t_s", "a_s", "height_0", "b_s1")), p0_es, mutant=FALSE)
        pl_1_es <- ES20_Individual(s = p1_es$strategies[[1]])
        pl_es <- pl_1_es # ES20_Plant()
        env_es <- get_constant_environment_ES20(stress=0.75)
        
        tt <- seq(0, 1, length.out=50)
        
        # Run single year
        res_es <- grow_plant_to_time(pl_es, tt, env_es)
        res_es_df <- as.data.frame(res_es$state)
        
        final <- res_es_df[nrow(res_es_df),]
        mass <- final[,'mass_leaf'] + final[, 'mass_sapwood'] + final[, 'mass_bark'] + final[,'mass_root']
        
        out_Mass[k] <- mass
        out_Store[k] <- final[,'mass_storage']
        
        if(any(res_es_df$mass_storage < 0, na.rm=TRUE)){
          valid[k] <- FALSE
          log4r::info(optts__logger, "Plant is dead")
          log4r::info(optts__logger, 
                      paste0("ts: ", ts_values[k], " bs: ", b_s_values[j], " h0:", h_0_values[i]))
        }
        else if (any(res_es$aux_size$dead_flag, na.rm = TRUE)){
          valid[k] <- FALSE
          log4r::info(optts__logger, "Plant is dead")
          log4r::info(optts__logger, 
                      paste0("ts: ", ts_values[k], " bs: ", b_s_values[j], " h0:", h_0_values[i]))
        }
        else{
          valid[k] <- TRUE
        }
        
      }, warning = function(w) {
        # print("warning")
      }, error = function(e) {
        
        log4r::error(optts__logger, "Error_message")
        log4r::error(optts__logger, e)
        log4r::info(optts__logger, 
                    paste0("ts: ", ts_values[k], " bs: ", b_s_values[j], " h0:", h_0_values[i]))
        
        out_Mass[k] <- -99999
        out_Store[k] <- -99999
        valid[k] <- FALSE
        
      }, finally = {
        # print('finished')
      })

    }

    # print("finished")

    ## Evaluate the outputs

    # get rid of results that are invalid

    out_Mass[!valid] <- -99999
    out_Store[!valid] <- -99999

    if(any(valid, na.rm=TRUE)){
      opt_ts_M[i,j] = ts_values[which.max(out_Mass)]
      opt_ts_S[i,j] = ts_values[which.max(out_Store)]
    }
    else{
      opt_ts_M[i,j] = NA
      opt_ts_S[i,j] = NA
    }

  }
  save(opt_ts_M, opt_ts_S, h_0_values, b_s_values, file="scripts/ES20/optimal_ts_data_as_010_tcrit_075.RData")
}


save(opt_ts_M, opt_ts_S, h_0_values, b_s_values, file="scripts/ES20/optimal_ts_data_as_010_tcrit_075.RData")

end <- Sys.time()

print(end-begin)

send_to_elisa("finished optimal ts")
