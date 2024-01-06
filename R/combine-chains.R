combine_chains <- function(posterior_list) {
  
  combined_list <- list()
  
  # Iterate over each sublist in posterior_list
  for(sublist in posterior_list) {
    # Iterate over each component in the sublist
    for(component_name in names(sublist)) {
      component <- sublist[[component_name]]
      
      # Check the type of the component
      if (is.matrix(component)) {
        if (exists(component_name, combined_list)) {
          combined_list[[component_name]] <- rbind(combined_list[[component_name]], component)
        } else {
          combined_list[[component_name]] <- component
        }
      } else if (is.vector(component) || is.list(component)) {
        if (exists(component_name, combined_list)) {
          combined_list[[component_name]] <- c(combined_list[[component_name]], component)
        } else {
          combined_list[[component_name]] <- component
        }
      }
    }
  }
  
  return(combined_list)
}