mcmc_pdp <- function(data, model, J, K = 40){
  
  # Create grid of variables to explore
  new_grid_vars <- list()
  for(i in 1:length(J)){
    v <- data[, J[i], drop = TRUE]
    n_unique <- min(length(unique(v)), K)
    new_grid_vars[[i]] <- seq(min(v), max(v), length.out = n_unique)
  }
  new_grid <- expand.grid(new_grid_vars) |> setNames(colnames(data)[J])
  
  # Create new dataset to make predictions on (use weights to reduce computation)
  newdata <- new_grid |>
    dplyr::cross_join(data[, !(1:ncol(data) %in% J)]) |>
    summarise(n = n(), .by = everything())
  
  # Make prediction on new dataset
  newdata_w_preds <- lapply(
    model$forests, 
    \(f){
      newdata$tau <- clbart::predict_forest(f, newdata)
      newdata |>
        summarise(tau = weighted.mean(tau, n), .by = colnames(new_grid))
    }) |>
    bind_rows(.id = 'Iteration') |>
    mutate(Iteration = as.numeric(Iteration))
  
  return(newdata_w_preds)
  
}
