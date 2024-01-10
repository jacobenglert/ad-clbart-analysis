#!/usr/bin/env Rscript

library(tidyverse)

in_dir <- here::here('Tables','Data', 'temp')
out_dir <- here::here('Tables','Data')

marg_files <- list.files(in_dir, full.names = TRUE)[str_starts(list.files(in_dir), 'marg')]
cart_files <- list.files(in_dir, full.names = TRUE)[str_starts(list.files(in_dir), 'cart')]
pred_files <- list.files(in_dir, full.names = TRUE)[str_starts(list.files(in_dir), 'pred')]

lapply(marg_files, read_rds) |> write_rds(here::here(out_dir,'marg-post.rds'))
file.remove(marg_files)

lapply(cart_files, read_rds) |> write_rds(here::here(out_dir,'cart-post.rds'))
file.remove(cart_files)

lapply(pred_files, read_rds) |> write_rds(here::here(out_dir,'pred-post.rds'))
file.remove(pred_files)

unlink(in_dir, recursive = TRUE)