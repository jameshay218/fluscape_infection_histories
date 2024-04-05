######################################################
## RUN SUPPLEMENT 1 SCRIPTS
## Author: James Hay
## Date: 09 January 2024
## Summary: runs a suite of simulation-recovery experiments. First simulates a single set of serosurvey data (script 1.simulate_data.R). We then fit a number of serosolver models to this data, changing some assumptions and inputs to test how model misspecification affects the ability to successfully recover ground-truth parameters. This corresponds to results in the first part of the supplementary results

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

top_wd <- "~/Documents/GitHub/fluscape_infection_histories//scripts/simulation_recovery/"
setwd(top_wd)

## ========================================================
## IF YOU'D LIKE TO RUN EACH OF THE AUX SCRIPTS BELOW ALONE,
## RUN THE CHUNK OF CODE BELOW 
## ========================================================
## Run name, filenames and working directory
run_name <- "sim_quarterly_fluscape_supplement_2"
main_wd <- "~/Documents/GitHub/fluscape_infection_histories//"
save_wd <- paste0(main_wd,"/data/simulated/")

source(paste0(main_wd,"/code/main.R"))

antigenic_map_file <- "data/antigenic_maps/antigenic_map_fonville_quarterly_clustered_clustered.csv"
antigenic_map_annual_file <- "data/antigenic_maps/antigenic_map_fonville_annual_clustered_clustered.csv"

save_wd_chain_plots <- paste0(main_wd,"/figures/supplement_2/chain_plots/")
chain_wd <- paste0(main_wd,"/chains/")

## Save titre data
titre_dat_filename <- paste0(save_wd,"/",run_name,"_titre_data.csv")
titre_dat_filename_annual <- paste0(save_wd,"/",run_name,"_titre_data_annual.csv")
par_tab_filename <- paste0(save_wd,"/",run_name,"_par_tab.csv")
attack_rates_filename <- paste0(save_wd,"/",run_name,"_attack_rates.csv")
infection_histories_filename <- paste0(save_wd,"/",run_name,"_infection_histories.csv")
rho_estimates_filename <- paste0(save_wd,"/",run_name,"_rho_estimates.csv")
par_tab_file_blank <- "inputs/par_tab_cr.csv"
rho_file <- "inputs/par_tab_offset_estimates.csv"
antigenic_map_file_wrong <- "data/antigenic_maps/antigenic_map_fonville_quarterly_continuous.csv"

## Run multiple chains in parallel
n_chains <- 5

## If not NULL, set to a vector of integers to choose only a subet of individuals to re-fit to
subset_indivs <- NULL

## MCMC pars
mcmc_pars <- c("save_block"=100,
               "thin"=1000,"thin_hist"=5000,
               "iterations"=5000000,
               "adaptive_period"=1000000,
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
## ========================================================

set.seed(1)

# 1. Simulate the serosurvey data -----------------------------------------------------
setwd(top_wd)
source("1.simulate_data.R")

# 2. Estimate offsets -----------------------------------------------------
## Re-estimate the measurement offset parameters
run_name_tmp <- paste0(run_name, "_offsets")
fluscape_serosolver_main(serosolver_wd,
                         mcmc_pars,
                         subset_indivs=subset_indivs,
                         run_name=run_name_tmp, 
                         main_wd,
                         save_wd=paste0(save_wd_chain_plots,"/",run_name_tmp,"/"), 
                         chain_wd=paste0(chain_wd,"/",run_name_tmp,"/"), 
                         buckets=1, 
                         prior_version=3, 
                         n_chains,
                         titre_dat_filename_annual, 
                         antigenic_map_annual_file,
                         par_tab_file_blank,
                         rho_file,
                         inf_hist_prior=c(2,5),
                         solve_likelihood=TRUE,
                         estimate_measurement_offsets=TRUE,
                         use_measurement_random_effects=FALSE,
                         create_prior_func=NULL,
                         proposal_ratios=NULL)

# Read in and save the estimated offset terms
chains <- load_theta_chains(paste0(chain_wd),burnin=mcmc_pars["adaptive_period"],convert_mcmc=TRUE)

theta_chain <- chains$chain
best_pars <- get_best_pars(theta_chain)

best_rhos <- best_pars[names(best_pars) %like% "rho"]
best_rhos <- best_rhos[!(names(best_rhos) %in% c("rho_mean","rho_sd"))]

write.table(data.frame(names="rho",values=best_rhos,fixed=1,steps=0.1,
                       lower_bound=-3,upper_bound=3,lower_start=-1,upper_start=1,
                       type=3), 
            rho_estimates_file, 
            sep=",",row.names=FALSE)

best_rhos <- best_rhos[best_rhos != 0]
use_names <- par_tab[par_tab$fixed == 0,"names"]

theta_chain %>% as_tibble() %>% select(c(names(best_rhos),"sampno","chain_no")) %>% pivot_longer(-c(sampno, chain_no)) %>%
  ggplot() + geom_line(aes(x=sampno,y=value,col=as.factor(chain_no))) + facet_wrap(~name,scales="free_y")

theta_chain %>% as_tibble() %>% select(c(names(best_rhos),"sampno","chain_no")) %>% pivot_longer(-c(sampno, chain_no)) %>%
  ggplot() + geom_density(aes(x=value,fill=as.factor(chain_no)),alpha=0.25) + facet_wrap(~name,scales="free")


colnames(chains$list[[1]])
chains_list <- chains$list
chains_list <- lapply(chains_list, function(x) x[,c(names(best_rhos),use_names)])

gelman.diag(as.mcmc.list(chains_list))

## Update some MCMC parameters to work better with the quarterly resolution models
setwd("~/Documents/GitHub/fluscape_infection_histories//scripts/simulation_recovery/")
mcmc_pars["hist_sample_prob"] <- 0.2
mcmc_pars["hist_switch_prob"] <- 0.1
mcmc_pars["year_swap_propn"] <- 1
mcmc_pars["move_size"] <- 4

# 6. Fit the full model at a per-quarter time resolution ----------------------------------------------------- 
run_name_tmp <- paste0(run_name, "_main")
fluscape_serosolver_main(serosolver_wd,
                         mcmc_pars,
                         subset_indivs=subset_indivs,
                         run_name=run_name_tmp, 
                         main_wd,
                         save_wd=paste0(save_wd_chain_plots,"/",run_name_tmp,"/"), 
                         chain_wd=paste0(chain_wd,"/",run_name_tmp,"/"), 
                         buckets=4, 
                         prior_version=2, 
                         n_chains,
                         titre_dat_filename, 
                         antigenic_map_file,
                         par_tab_file_blank,
                         rho_file,
                         inf_hist_prior=c(1,1),
                         solve_likelihood=TRUE,
                         estimate_measurement_offsets=FALSE,
                         use_measurement_random_effects=FALSE,
                         append_fixed_rho=TRUE,
                         create_prior_func=NULL,
                         proposal_ratios=NULL)


# 8. Fit the model using the smooth map ----------------------------------------------------- 
run_name_tmp <- paste0(run_name, "_smooth_map")
fluscape_serosolver_main(serosolver_wd,
                         mcmc_pars,
                         subset_indivs=subset_indivs,
                         run_name=run_name_tmp, 
                         main_wd,
                         save_wd=paste0(save_wd_chain_plots,"/",run_name_tmp,"/"), 
                         chain_wd=paste0(chain_wd,"/",run_name_tmp,"/"), 
                         buckets=4, 
                         prior_version=2, 
                         n_chains,
                         titre_dat_filename, 
                         antigenic_map_file_wrong,
                         par_tab_file_blank,
                         rho_file,
                         inf_hist_prior=c(1,1),
                         solve_likelihood=TRUE,
                         estimate_measurement_offsets=FALSE,
                         use_measurement_random_effects=FALSE,
                         append_fixed_rho=TRUE,
                         create_prior_func=NULL,
                         proposal_ratios=NULL)
