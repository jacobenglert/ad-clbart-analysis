# Program Name: ad-hw-results.R
# Author: Jacob Englert
# Date: 04JAN2024
# Description: Process analysis results into interim format used to make
#   tables and figures.

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(clbart)

# Helper Functions --------------------------------------------------------
source(here::here('R','combine-chains.R'))
source(here::here('R','mcmc-pdp.R'))


# Specify Results and Output Directories ----------------------------------
run_date <- '05JAN2024'
res_dir <- here::here('Results','Heatwave', run_date)
out_dir <- here::here('Tables','Heatwave', run_date, 'Data')
dir.create(out_dir, recursive = TRUE)


# Load Parameter Settings -------------------------------------------------

# Original parameter table
params <- read_csv(here::here('Params', 'params.csv'), show_col_types = FALSE)

# Parameter settings to process
param_groups <- params |> 
  filter(num_trees == 100 & beta_rho == 2) |>
  group_split(alpha_rho, num_trees, race_eth_cat, sex)


# Load Data ---------------------------------------------------------------
AD <- read_csv(here::here('Data','Clean','ad-ca-cco.csv'), show_col_types = FALSE) |> 
  na.omit() |>                  # Remove records with missing data
  filter(var(HW) > 0, .by = ID) # Remove individuals who did not experience a heatwave


# Process Results ---------------------------------------------------------

marg_post <- list()
cart_post <- list()
pred_post <- list()

for(group in param_groups){
  
  start <- Sys.time()
  i <- length(marg_post) + 1
  
  # Group Parameters
  alpha_rho <- unique(group$alpha_rho)
  beta_rho <- unique(group$beta_rho)
  num_trees <- unique(group$num_trees)
  
  # Subgroup
  race_eth_cat <- eval(parse(text = unique(group$race_eth_cat)))
  sex <- eval(parse(text = unique(group$sex)))
 
  # Load the model fits and combine them
  group_ids <- sprintf('%04d', group$ID)
  fit <- lapply(group_ids, \(x) read_rds(here::here(res_dir, paste0(x, '.rds')))) |>
    combine_chains()
  
  # Filter based on race/ethnicity
  if(race_eth_cat == 'HISP'){
    AD2 <- filter(AD, HISPANIC == 1)
  } else if(race_eth_cat == 'NHW'){
    AD2 <- filter(AD, WHITE == 1 & HISPANIC == 0)
  } else if(race_eth_cat == 'NHB'){
    AD2 <- filter(AD, BLACK == 1 & HISPANIC == 0)
  } else if(race_eth_cat == 'NHAPI'){
    AD2 <- filter(AD, API == 1 & HISPANIC == 0)
  } else if(race_eth_cat == 'Other'){
    AD2 <- filter(AD, RACEOTH == 1 & HISPANIC == 0)
  }
  
  # Filter based on sex
  if(sex == 'F'){
    AD2 <- filter(AD2, FEMALE == 1)
  } else if(sex == 'M'){
    AD2 <- filter(AD2, FEMALE == 0)
  } 
  
  AD2 <- AD2 |>
    filter(row_number() == 1, .by = ID) |>
    select(colnames(fit$split_props))
  
  # Obtain posterior distribution of CLORs
  tau <- parallel::mclapply(fit$forests, \(f) predict_forest(f, AD2),
                            mc.cores = 5) |>
    do.call(what = rbind)
  
  # Compute ACLOR
  ACLOR <- tau |>
    as.data.frame() |>
    mutate(Iteration = row_number()) |>
    pivot_longer(cols = -Iteration, names_to = 'ID', values_to = 'CLOR') |>
    summarise(est = mean(CLOR), .by = Iteration) |>
    mutate(Variable = 'ACLOR')
  
  # Compute marginal pACLORs
  J <- match(setdiff(colnames(AD2), 'AGE'), colnames(AD2))
  marg_pACLOR <- parallel::mclapply(J, \(j) mcmc_pdp(AD2, fit, J = j, K = 2), mc.cores = 7)
  
  # Compute difference in marginal pACLORs and append ACLOR
  marg_post[[i]] <- marg_pACLOR |>
    bind_rows() |>
    pivot_longer(cols = all_of(setdiff(colnames(AD2), 'AGE')),
                 names_to = 'Variable', values_to = 'Value') |>
    filter(!is.na(Value)) |>
    pivot_wider(id_cols = c(Iteration, Variable), names_from = Value, values_from = tau) |>
    mutate(est = `1` - `0`) |>
    bind_rows(ACLOR) |>
    mutate(race_eth_cat = race_eth_cat, sex = sex,
           alpha_rho = alpha_rho, beta_rho = beta_rho, num_trees = num_trees)
  
  # Obtain lower-dimension posterior summary with CART
  tau_data <- cbind(AD2, tau = colMeans(tau))
  tau_cart <- rpart::rpart(tau ~ ., data = tau_data)
  cart_R2 <- 1 - (sum((tau_data$tau - predict(tau_cart))^2) / sum((tau_data$tau - mean(tau_data$tau))^2))

  # Compute pACLORs from CART summary
  var_list_cart <- na.omit(setdiff(names(tau_cart$variable.importance), 'AGE')[1:3])
  J <- sapply(var_list_cart, \(x) which(colnames(AD2) == x))
  cart_pACLOR <- mcmc_pdp(data = AD2, model = fit, J = J, K = 2)

  # Store CART summary results
  cart_post[[i]] <- list(data = cart_pACLOR |>
                           mutate(race_eth_cat = race_eth_cat, sex = sex,
                                  alpha_rho = alpha_rho, beta_rho = beta_rho, 
                                  num_trees = num_trees, 
                                  R2 = cart_R2), 
                           var_list_cart = var_list_cart)
  
  # Store predictions
  pred_post[[i]] <- tau_data |>
    mutate(race_eth_cat = race_eth_cat, sex = sex,
           alpha_rho = alpha_rho, beta_rho = beta_rho, num_trees = num_trees)

  print(paste0(race_eth_cat, '-', sex, '-', beta_rho, '-', num_trees, ' completed in ', round(difftime(Sys.time(), start, units = 'mins'), 2)))

}

write_rds(marg_post, here::here(out_dir, 'marg-post.rds'))
write_rds(cart_post, here::here(out_dir, 'cart-post.rds'))
write_rds(pred_post, here::here(out_dir, 'pred-post.rds'))



