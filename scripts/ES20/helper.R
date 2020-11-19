
extract_trees <- function(all_results){
  species_list <- lapply(1:length(all_results[[length(all_results)]]$species), function(species_ind){
    species_df <- data.frame()
    for(i in 1:length(all_results)){
      new_df <- as.data.frame(t(all_results[[i]]$species[[species_ind]]))
      if(nrow(new_df) > 0){
        new_df$time <- rep(all_results[[i]]$time, times = nrow(new_df))
        new_df$tree_id <- 1:nrow(new_df)
        new_df$species_id <- rep(species_ind, times = nrow(new_df))
      }
      species_df <- rbind(species_df, new_df)
    }
    return(species_df)
  })
  
  species_df <- data.frame()
  
  for(i in length(species_list)){
    species_df <- rbind(species_df, species_list[[i]])
  }
  
  return(species_df)  
  
}

#COLLECT FUNCTIONS

collect_state <- function(obj) {
  obj$state
  
}

collect_aux <- function(obj){
  aux <- lapply(obj$patch$species, function(species){
    sapply(species$plants, function(plant){
      sapply(plant$aux_names, function(name){
        plant$aux(name)
      })
    })
  })
  return(aux)
}

collect_all <- function(obj){
  states <- collect_state(obj)
  aux <- collect_aux(obj)
  states$species <- lapply(1:length(states$species), function(i){
    df <- rbind(rbind(states$species[[i]], attr(states$species[[i]], "is_alive")), aux[[i]])
    rnms <- rownames(df)
    rnms[(nrow(states$species[[i]]) + 1)] <- "is_alive"
    rownames(df) <- rnms
    return(df)
  }) 
  return(states)
}

extract_environment <- function(all_results){
  final_df <- data.frame()
  
  for(i in 1:length(all_results)){
    environment <- as.data.frame(all_results[[i]]$env)
    environment$time <- rep(all_results[[i]]$time, times = nrow(environment))
    final_df <- rbind(final_df, environment)
    print(environment)
    invisible(readline(prompt="Press [enter] to continue"))
  }
  return(final_df)
}

# extract_trees <- function(all_results){
#   species_list <- lapply(1:length(all_results[[length(all_results)]]$species), function(species_ind){
#     species_df <- data.frame()
#     for(i in 1:length(all_results)){
#       new_df <- as.data.frame(t(all_results[[i]]$species[[species_ind]]))
#       new_df$time <- rep(all_results[[i]]$time, times = nrow(new_df))
#       new_df$tree_id <- 1:nrow(new_df)
#       new_df$species_id <- rep(species_Ind, times = nrow(new_df))
#       
#       species_df <- rbind(species_df, new_df)
#     }
#   })
#   
#   return(species_list)  
#   
# }
