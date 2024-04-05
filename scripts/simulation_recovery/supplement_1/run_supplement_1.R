######################################################
## RUN EXTENSIVE SIMULATION-RECOVERY EXPERIMENTS
## Author: James Hay
## Date: 05 Jan 2023
## Summary: runs a suite of simulation-recovery experiments for the first set of sensitivity analyses. First simulates a single set of serosurvey data (script 1.simulate_data.R). We then fit a number of serosolver models to this data, changing some assumptions and inputs to test how model misspecification affects the ability to successfully recover ground-truth parameters. 

library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
library(data.table)
library(foreach)
library(doParallel)
library(parallel)
library(dplyr)
library(tidyverse)

serosolver_wd <- "~/Documents/GitHub/serosolver/"
devtools::load_all(serosolver_wd)
#library(serosolver)

setwd("~/Documents/GitHub/fluscape_infection_histories//scripts/simulation_recovery/supplement_1/")


run_name <- "sim_quarterly_fluscape_supplement_1"
main_wd <- "~/Documents/GitHub/fluscape_infection_histories//"
save_wd <- paste0(main_wd,"/data/simulated/")

antigenic_map_file <- "data/antigenic_maps/antigenic_map_fonville_quarterly_continuous.csv"
antigenic_map_annual_file <- "data/antigenic_maps/antigenic_map_fonville_annual_continuous.csv"

save_wd_chain_plots <- paste0(main_wd,"/figures/supplement_1/chain_plots/")
chain_wd <- paste0(main_wd,"/chains/")

## Save titre data
titre_dat_filename <- paste0(save_wd,"/",run_name,"_titre_data.csv")
titre_dat_filename_annual <- paste0(save_wd,"/",run_name,"_titre_data_annual.csv")
par_tab_filename <- paste0(save_wd,"/",run_name,"_par_tab.csv")
attack_rates_filename <- paste0(save_wd,"/",run_name,"_attack_rates.csv")
infection_histories_filename <- paste0(save_wd,"/",run_name,"_infection_histories.csv")
rho_estimates_file <- paste0(save_wd,"/",run_name,"_rho_estimates.csv")
par_tab_file_blank <- "inputs/par_tab_cr.csv"
rho_file <- "inputs/par_tab_offset_estimates.csv"
antigenic_map_file_wrong <- "data/antigenic_maps/antigenic_map_bedford_quarterly_continuous.csv"



## Run multiple chains in parallel
n_chains <- 5

## If not NULL, set to a vector of integers to choose only a subet of individuals to re-fit to
subset_indivs <- NULL

## MCMC pars
mcmc_pars <- c("save_block"=100,
               "thin"=1000,"thin_hist"=5000,
               "iterations"=2000000,
               "adaptive_period"=500000,
               "burnin"=0,
               "switch_sample"=2,
               "hist_switch_prob"=0,
               "year_swap_propn"=0.8,
               "swap_propn"=0.5,
               "inf_propn"=0.5,
               "hist_sample_prob"=1,
               "move_size"=3, 
               "hist_opt"=1,
               "popt"=0.44,
               "popt_hist"=0.44,
               "opt_freq"=2000,
               propose_from_prior=TRUE)

## Run multiple chains in parallel
cl <- makeCluster(5)
registerDoParallel(cl)

set.seed(1)

## Simulate the serosurvey data
setwd("~/Documents/GitHub/fluscape_infection_histories//scripts/simulation_recovery/supplement_1/")
source("1.simulate_data.R")

## Re-estimate the measurement offset parameters
setwd("~/Documents/GitHub/fluscape_infection_histories//scripts/simulation_recovery/supplement_1/")
source("2.estimate_offsets.R")

## Update some MCMC parameters to work better with the quarterly resolution models
setwd("~/Documents/GitHub/fluscape_infection_histories//scripts/simulation_recovery/supplement_1/")
mcmc_pars["hist_sample_prob"] <- 0.2
mcmc_pars["hist_switch_prob"] <- 0.1
mcmc_pars["year_swap_propn"] <- 1
mcmc_pars["move_size"] <- 4

## Fit the full model at a per-quarter time resolution
source("7.fit_main_model.R")
setwd("~/Documents/GitHub/fluscape_infection_histories//scripts/simulation_recovery/supplement_1/")

## Fit the model using the Bedford map
source("9.fit_wrong_map.R")
setwd("~/Documents/GitHub/fluscape_infection_histories//scripts/simulation_recovery/supplement_1/")

## Fit the correct model but ignoring offsets
source("11.fit_main_model_no_offsets.R")
