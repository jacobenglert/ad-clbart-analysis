# Program Name: set-params.R
# Author: Jacob Englert
# Date: 11AUG2023
# Description: Generate hyperparameter tables for real data simulations.


# Load Packages -----------------------------------------------------------
library(tidyverse)


# Specify Paramaters ------------------------------------------------------
alpha_rho <- c(0.95, 0.3)
beta_rho <- c(2, 3)
seed <- 1:5
race_eth_cat <- c("'HISP'","'NHW'","'NHB'","'NHAPI'","'Other'")
sex <- c("'M'","'F'")
num_trees <- c(50, 100)

params <- crossing(race_eth_cat, sex,
                   nesting(alpha_rho, beta_rho), 
                   num_trees, seed) |>
  mutate(iter = 10000, thin = 5, warmup = 5000,
         sigma2_beta = 1, sigma2_beta_update_freq = 50, beta_acc_prob = 0.35,
         n_min = 5,
         alpha_sigma = 1, beta_sigma = 1,
         moves = "c('grow','prune','change')", move_probs = 'c(0.3,0.3,0.4)',
         ID = row_number()) |>
  select(ID, everything())


# Export Parameters -------------------------------------------------------
write_csv(params, here::here('Params', 'params.csv'))
