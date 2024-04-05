######################################################
## SIMULATION RECOVERY TEST -- VIETNAM-LIKE DATA
## Author: James Hay
## Date: 14 Feb 2023
## Summary: simulates some serosurvey data and fits serosolver

library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
library(data.table)
library(foreach)
library(doParallel)
library(parallel)
library(dplyr)
library(readr)

rerun_chains <- TRUE

serosolver_wd <- "~/Documents/GitHub/serosolver/"
devtools::load_all(serosolver_wd)
#library(serosolver)

run_name <- "sim_vietnam"
main_wd <- "~/Documents/GitHub/fluscape_infection_histories//"
chain_wd <- paste0(main_wd,"/chains/",run_name)
save_wd <- paste0(main_wd,"/figures/chain_plots/")
buckets <- 1
prior_version <- 2
n_chains <- 5

setwd(main_wd)
print(paste0("In directory: ", main_wd))
print(paste0("Saving to: ", save_wd))

set.seed(1)

## Run multiple chains in parallel
cl <- makeCluster(n_chains)
registerDoParallel(cl)

if(!dir.exists(save_wd)) dir.create(save_wd,recursive = TRUE)
if(!dir.exists(chain_wd)) dir.create(chain_wd,recursive = TRUE)

mcmc_pars <- c("save_block"=100,
               "thin"=10,"thin_hist"=100,
               "iterations"=1000000,
               "adaptive_period"=200000,
               "burnin"=0,"switch_sample"=2,"hist_switch_prob"=0.05,
               "year_swap_propn"=0.8,"swap_propn"=0.5,
               "inf_propn"=0.5,"hist_sample_prob"=1,"move_size"=3, "hist_opt"=0,
               "popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,propose_from_prior=FALSE)

## Simulation parameters
n_indivs <- 69
n_groups <- 1
n_samps <- 5
repeats <- 2
samp_min <- 2009
samp_max <- 2015
year_min <- 1968
year_max <- 2015
age_min <- 5
age_max <- 75

sampled_viruses <- seq(year_min, year_max, by=2)
sampling_times <- seq(samp_min, samp_max, by=1)

## Antigenic map
antigenic_map <- read.csv("data/antigenic_maps/antigenic_map_fonville_annual_continuous.csv")
strain_isolation_times <- antigenic_map$inf_times

## Set up parameter table
#par_tab <- read.csv("inputs/par_tab_base.csv",stringsAsFactors=FALSE)
par_tab <- read.csv("inputs/par_tab_cr.csv",stringsAsFactors=FALSE)
par_tab <- par_tab[par_tab$names != "phi",]
par_tab[par_tab$names %in% c("alpha","beta"),c("values")] <- c(1/3,1/3)

## Simulate realistic-ish attack rates
sim_inf_pars=c("mean"=0.15/buckets,"sd"=1.5,large_first_year=TRUE,"bigMean"=0.6/buckets)
attack_rates <- simulate_attack_rates(strain_isolation_times, sim_inf_pars["mean"],sim_inf_pars["sd"],TRUE,sim_inf_pars["bigMean"])
attack_rates[1] <- 0.6
attack_rates[attack_rates > 0.6] <- 0.6
plot(attack_rates)

sim_data <- simulate_data(par_tab=par_tab, group=1, n_indiv=n_indivs, 
                     buckets=buckets,
                     strain_isolation_times=strain_isolation_times,
                     measured_strains=sampled_viruses,
                     sampling_times=sampling_times, nsamps=n_samps, antigenic_map=antigenic_map, 
                     titre_sensoring=0.2, ## Randomly censor 20% of measurements
                     age_min=age_min,age_max=age_max,
                     attack_rates=attack_rates, repeats=repeats,
                     mu_indices = NULL, measurement_indices = NULL,
                     add_noise=TRUE)


plot_data(sim_data$data,sim_data$infection_histories,strain_isolation_times,n_indivs=10)
titre_dat <- sim_data$data %>% tidyr::drop_na()
titre_dat <- titre_dat %>% dplyr::group_by(individual, samples, virus) %>% dplyr::mutate(run=1:n()) %>% ungroup()
titre_dat <- titre_dat %>% dplyr::left_join(sim_data$ages)
titre_dat$DOB <- 1940 ## Pretend that we don't know individual ages
titre_dat <- titre_dat %>% arrange(individual, samples, virus, run) %>% as.data.frame()

## Save titre data
write_csv(titre_dat, file=paste0(save_wd,"/",run_name,"_titre_data.csv"))
## Save parameter table
write_csv(par_tab, file=paste0(save_wd,"/",run_name,"_par_tab.csv"))
## Save attack rates
write_csv(sim_data$attack_rates, file=paste0(save_wd,"/",run_name,"_attack_rates.csv"))
## Save infection histories
write_csv(as.data.frame(sim_data$infection_histories), file=paste0(save_wd,"/",run_name,"_infection_histories.csv"))


f <- create_posterior_func(par_tab,titre_dat,antigenic_map=antigenic_map,strain_isolation_times = strain_isolation_times,
                           version=prior_version,solve_likelihood=TRUE,n_alive=NULL)
if(rerun_chains){
    
## Time runs and use dopar to run multiple chains in parallel
t1 <- Sys.time()
filenames <- paste0(chain_wd, "/",run_name, "_",1:n_chains)
res <- foreach(x = filenames, .packages = c('data.table','plyr',"dplyr")) %dopar% {
    devtools::load_all(serosolver_wd)
    
    index <- 1
    lik <- -Inf
    inf_hist_correct <- 1
    while((!is.finite(lik) || inf_hist_correct > 0) & index < 100){
        start_tab <- generate_start_tab(par_tab)
        start_inf <- setup_infection_histories_total(titre_dat,strain_isolation_times,2,3)
        
        inf_hist_correct <- sum(check_inf_hist(titre_dat, strain_isolation_times, start_inf))
        y <- f(start_tab$values, start_inf)
        lik <- sum(y[[1]])
        index <- index + 1
    }
    
    write.csv(start_tab, paste0(x, "_start_tab.csv"))
    write.csv(start_inf, paste0(x, "_start_inf_hist.csv"))
    write.csv(antigenic_map, paste0(x, "_antigenic_map.csv"))
    write.csv(titre_dat, paste0(x, "_titre_dat.csv"))
    
    res <- serosolver::run_MCMC(start_tab, titre_dat, antigenic_map, start_inf_hist=start_inf,filename=x,
                                CREATE_POSTERIOR_FUNC=create_posterior_func, 
                                CREATE_PRIOR_FUNC = NULL,
                                version=prior_version,
                                mcmc_pars=mcmc_pars,
                                solve_likelihood=TRUE)
}
run_time_fast <- Sys.time() - t1
run_time_fast
}
## Read in chains for trace plot
chains <- load_mcmc_chains(chain_wd,convert_mcmc=TRUE,burnin = mcmc_pars["adaptive_period"],unfixed = TRUE)
pdf(paste0(save_wd,"/",run_name,"_chain.pdf"))
plot(as.mcmc.list(chains$theta_list_chains))
dev.off()

## Read in chains for all other plots
chains <- load_mcmc_chains(chain_wd,convert_mcmc=FALSE,burnin = mcmc_pars["adaptive_period"],unfixed = FALSE)
chain <- as.data.frame(chains$theta_chain)
inf_chain <- chains$inf_chain 

## Plot kinetics parameter estimates and number of infections
p1 <- plot_posteriors_theta(chain, par_tab, burnin=0,samples=100,TRUE,FALSE,FALSE, TRUE,"")
p_total_inf <- plot_total_number_infections(inf_chain,FALSE)
ggsave(paste0(save_wd,"/",run_name,"_total_inf.pdf"),p_total_inf[[1]],height=5,width=7,units="in",dpi=300)

p_cumu_infs <- generate_cumulative_inf_plots(inf_chain,indivs=1:5,real_inf_hist=sim_data$infection_histories,strain_isolation_times = strain_isolation_times)
ggsave(paste0(save_wd,"/",run_name,"_inf_hists.pdf"),p_cumu_infs[[1]],height=8,width=7,units="in",dpi=300)

n_alive <- get_n_alive_group(titre_dat,strain_isolation_times)

## Plot attack rates
n_inf <- sim_data$infection_histories %>% colSums()
true_ar <- data.frame(j=strain_isolation_times,group=1,AR=n_inf/n_alive[1,])
p_ar <- plot_attack_rates(inf_chain, titre_dat, strain_isolation_times,pad_chain=FALSE,plot_den=FALSE,n_alive = n_alive,
                          true_ar=true_ar,
                          prior_pars = c("prior_version"=2,"alpha"=1,"beta"=1))
ggsave(paste0(save_wd,"/",run_name,"_ar.pdf"),p_ar,height=7,width=8,units="in",dpi=300)


## Plot model fits
titre_preds <- get_titre_predictions(chain = chain[chain$chain_no == 1,], 
                                     infection_histories = inf_chain[inf_chain$chain_no == 1,], 
                                     titre_dat = titre_dat, 
                                     individuals = unique(titre_dat$individual),
                                     antigenic_map = antigenic_map, 
                                     par_tab = par_tab,expand_titredat=FALSE)

use_indivs <- sample(unique(titre_dat$individual), 9)
titre_pred_p <- ggplot(titre_preds$predictions[titre_preds$predictions$individual %in% use_indivs,])+
    geom_ribbon(data=titre_preds$predicted_observations[titre_preds$predicted_observations$individual %in% use_indivs,], 
                aes(x=virus,ymin=lower,ymax=upper),fill="grey90") +
    geom_ribbon(aes(x=virus,ymin=lower, ymax=upper),fill="gray70")+
    geom_line(aes(x=virus, y=median))+
    geom_point(aes(x=virus, y=titre))+
    coord_cartesian(ylim=c(0,8))+
    ylab("log titre") +
    xlab("Time of virus circulation") +
    theme_classic() +
    facet_grid(individual~samples)

ggsave(paste0(save_wd,"/",run_name,"_titre_fits.pdf"),titre_pred_p,height=7,width=8,units="in",dpi=300)

