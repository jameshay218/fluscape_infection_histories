######################################################
## FIT SEROSOLVER MODEL TO SIMULATED FLUSCAPE DATA TO ESTIMATE OFFSETS
## Author: James Hay
## Date: 21 July 2023
## Summary: fits the annual infection history and antibody kinetics model to the simulated fluscape data
## to estimate the virus-specific measurement offset terms, using prior version 3
run_name_tmp <- paste0(run_name,"_main")
main_wd <- "~/Documents/GitHub/fluscape_infection_histories//"
source(paste0(main_wd,"/code/main.R"))


save_wd_tmp <- paste0(save_wd_chain_plots,"/",run_name_tmp,"/")
chain_wd_tmp <- paste0(chain_wd,"/",run_name_tmp,"/")

buckets <- 4
prior_version <- 2

fluscape_serosolver_main(serosolver_wd,mcmc_pars,
                         subset_indivs=subset_indivs,
                         run_name_tmp, 
                         main_wd,save_wd_tmp, chain_wd_tmp, buckets, prior_version, n_chains,
                         titre_dat_file, antigenic_map_file,par_tab_file_blank,
                         rho_estimates_file,
                         inf_hist_prior=c(1,1),
                         solve_likelihood=TRUE,
                         estimate_measurement_offsets=FALSE,
                         use_measurement_random_effects=FALSE,
                         append_fixed_rho=TRUE,
                         create_prior_func=NULL,
                         proposal_ratios=NULL,
                         mixed_time_cutoff=NULL)
