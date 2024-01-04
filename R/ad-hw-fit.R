# Program Name: ad-hw-fit.R
# Author: Jacob Englert
# Date: 06SEP2023
# Description: Fit CL-BART models to CA AD data using heatwave as exposure

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(clbart)


# Get Parameters ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
names(args) <- c('index')
index <- as.numeric(args['index'])
params <- read_csv(here::here('Params','params.csv'), show_col_types = F)[,-1]
for(i in 1:ncol(params)){
  assign(names(params[index,i]), eval(parse(text = params[index, i, drop = TRUE])))
}


# Load Data ---------------------------------------------------------------
AD <- read_csv(here::here('Data','Clean','ad-ca-cco.csv'), show_col_types = FALSE) |> 
  na.omit() |>                  # Remove records with missing data
  filter(var(HW) > 0, .by = ID) # Remove individuals who did not experience a heatwave

# Filter based on race/ethnicity
if(race_eth_cat == 'HISP'){
  AD <- filter(AD, HISPANIC == 1)
} else if(race_eth_cat == 'NHW'){
  AD <- filter(AD, WHITE == 1 & HISPANIC == 0)
} else if(race_eth_cat == 'NHB'){
  AD <- filter(AD, BLACK == 1 & HISPANIC == 0)
} else if(race_eth_cat == 'NHAPI'){
  AD <- filter(AD, API == 1 & HISPANIC == 0)
} else if(race_eth_cat == 'Other'){
  AD <- filter(AD, RACEOTH == 1 & HISPANIC == 0)
}

# Filter based on sex
if(sex == 'F'){
  AD <- filter(AD, FEMALE == 1)
} else if(sex == 'M'){
  AD <- filter(AD, FEMALE == 0)
} 

x <- cbind(AD$HOLIDAY,
           splines::ns(AD$DP3DMA, 4),
           splines::ns(AD$AVG3DMA, 4))
colnames(x) <- c('HOLIDAY', paste0('DP3DMAns', 1:4), paste0('AVG3DMAns', 1:4))

w <- AD |> select(AGE, CKD, COPD, DEPRESSION, HT, CHF, DIAB, HYPERLIP)
y <- AD$Case
z <- AD$HW
stratum <- match(AD$ID, unique(AD$ID))


# Fit CL-BART Model -------------------------------------------------------
fit <- clbart(w, x, y, z, stratum, 
              num_trees = num_trees, seed = seed,
              iter = iter, thin = thin, warmup = warmup,
              sigma2_beta = sigma2_beta,
              sigma2_beta_update_freq = sigma2_beta_update_freq,
              beta_acc_prob = beta_acc_prob,
              n_min = n_min,
              moves = moves, move_probs = move_probs,
              alpha_rho = alpha_rho, beta_rho = beta_rho,
              alpha_sigma = alpha_sigma, beta_sigma = beta_sigma) 

fit$race_eth_cat <- race_eth_cat
fit$sex <- sex


# Output Results ----------------------------------------------------------
message(paste0("SLURM_ARRAY_TASK_ID is : ", index))
if(!dir.exists(here::here('Results','Heatwave','temp'))) dir.create(here::here('Results','Heatwave','temp'))
saveRDS(fit, here::here('Results','Heatwave','temp', paste0(sprintf("%04d", index), ".rds")))

