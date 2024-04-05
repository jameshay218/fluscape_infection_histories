######################################################
## FIT SEROSOLVER MODEL TO SIMULATED FLUSCAPE DATA TO ESTIMATE OFFSETS
## Author: James Hay
## Date: 21 July 2023
## Summary: fits the annual infection history and antibody kinetics model to the simulated fluscape data
## to estimate the virus-specific measurement offset terms, using prior version 3
run_name_tmp <- paste0(run_name, "_offsets")
source(paste0(main_wd,"/code/main.R"))

save_wd_tmp <- paste0(save_wd_chain_plots,"/",run_name_tmp,"/")
chain_wd_tmp <- paste0(chain_wd,"/",run_name_tmp,"/")

buckets <- 1
prior_version <- 3

## Run serosolver
fluscape_serosolver_main(serosolver_wd,
                         mcmc_pars,
                         subset_indivs=subset_indivs,
                         run_name_tmp, 
                         main_wd,
                         save_wd_tmp, chain_wd_tmp, buckets, 
                         prior_version, n_chains,
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

## Read in and save the estimated offset terms
chains <- load_theta_chains(chain_wd_tmp,burnin=mcmc_pars["adaptive_period"],convert_mcmc=TRUE)

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

