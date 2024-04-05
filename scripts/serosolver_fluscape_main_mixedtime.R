######################################################
## FIT SEROSOLVER MODEL TO FLUSCAPE DATA WITH MEASUREMENT OFFSETS AND MIXED TIME
## Author: James Hay
## Date: 07 JULY 2023
## Summary: fits the quarterly infection history and antibody kinetics model to the Fluscape titer dataset and generates some summary plots. Uses measurement offsets. Up until 2005, fit at annual resolution, then switch to quarterly.

library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
library(data.table)
library(foreach)
library(doParallel)
library(parallel)
library(dplyr)

serosolver_wd <- "~/Documents/GitHub/serosolver/"
devtools::load_all(serosolver_wd)
##library(serosolver)

run_name <- "fluscape_mixed_time"
main_wd <- "~/Documents/GitHub/fluscape_infection_histories//"
chain_wd <- paste0(main_wd,"/chains/",run_name)
save_wd <- paste0(main_wd,"/figures/chain_plots/")
buckets <- 4
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

mcmc_pars <- c("save_block"=100,"thin"=10000,"thin_hist"=35000,"iterations"=30000000,
               "adaptive_period"=5000000,
               "burnin"=0,"switch_sample"=2,"hist_switch_prob"=0.1,
               "year_swap_propn"=1,"swap_propn"=0.5,
               "inf_propn"=1,"hist_sample_prob"=0.2,"move_size"=3*buckets, "hist_opt"=0,
               "popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,propose_from_prior=TRUE)

titre_dat <- read.csv("data/fluscape_data_4_resolution.csv",stringsAsFactors = FALSE)
titre_dat <- titre_dat %>% arrange(individual, samples, virus, run)
titre_dat$group <- 1

## Antigenic map
antigenic_map <- read.csv("data/antigenic_maps/antigenic_map_fonville_quarterly_continuous.csv")
antigenic_map <- antigenic_map[antigenic_map$inf_times >= 1968*buckets & antigenic_map$inf_times <= max(titre_dat$samples),]
strain_isolation_times <- c(seq(1968*4, 2005*4, by=4),8021:max(antigenic_map$inf_times))
antigenic_map <- antigenic_map[antigenic_map$inf_times %in% strain_isolation_times,]

## Set up parameter table
par_tab <- read.csv("inputs/par_tab_base.csv",stringsAsFactors=FALSE)
par_tab <- par_tab[par_tab$names != "phi",]
par_tab[par_tab$names %in% c("alpha","beta"),c("values")] <- c(1,1)
par_tab[par_tab$names == "wane",c("values")] <- 0.8

par_tab_rhos <- read.csv("inputs/par_tab_offset_estimates.csv",stringsAsFactors = FALSE)
par_tab <- bind_rows(par_tab, par_tab_rhos)

## If we were assigning measurement indices for each quarter
strain_isolation_times_full <- c(seq(1968*4, max(antigenic_map$inf_times)))
measurement_indices <- rep(1:48,each=buckets)
measurement_indices <- measurement_indices[1:length(strain_isolation_times_full)]
## Create measurement bias indices, this bit is tricky
measurement_indices_fit <- measurement_indices[match(strain_isolation_times, strain_isolation_times_full)]


f <- create_posterior_func(par_tab,titre_dat,antigenic_map=antigenic_map,strain_isolation_times = strain_isolation_times,
                           version=prior_version,solve_likelihood=TRUE,n_alive=NULL,
                           measurement_indices_by_time=measurement_indices_fit)
print("Starting fits...")

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
        start_inf <- setup_infection_histories_titre(titre_dat,strain_isolation_times,space=5*buckets,titre_cutoff=2,sample_prob=0.9)
        
        inf_hist_correct <- sum(check_inf_hist(titre_dat, strain_isolation_times, start_inf))
        y <- f(start_tab$values, start_inf)
        lik <- sum(y[[1]])
        index <- index + 1
    }
    
    write.csv(start_tab, paste0(x, "_start_tab.csv"))
    write.csv(start_inf, paste0(x, "_start_inf_hist.csv"))
    write.csv(antigenic_map, paste0(x, "_antigenic_map.csv"))
    write.csv(titre_dat, paste0(x, "_titre_dat.csv"))
    
    res <- serosolver::run_MCMC(start_tab, titre_dat, 
                                antigenic_map, start_inf_hist=start_inf,filename=x,
                                CREATE_POSTERIOR_FUNC=create_posterior_func, 
                                CREATE_PRIOR_FUNC = NULL,
                                version=prior_version,
                                measurement_indices=measurement_indices_fit,
                                measurement_random_effects = FALSE,
                                mcmc_pars=mcmc_pars,
                                solve_likelihood=TRUE)
}
run_time_fast <- Sys.time() - t1
run_time_fast

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

## Plot attack rates
p_ar <- plot_attack_rates(inf_chain, titre_dat, strain_isolation_times,pad_chain=FALSE,plot_den=FALSE,n_alive = NULL,
                          prior_pars = c("prior_version"=2,"alpha"=1,"beta"=1))
ggsave(paste0(save_wd,"/",run_name,"_ar.pdf"),p_ar,height=7,width=8,units="in",dpi=300)

## Plot model fits
titre_preds <- get_titre_predictions(chain = chain[chain$chain_no == 1,], 
                                     infection_histories = inf_chain[inf_chain$chain_no == 1,], 
                                     titre_dat = titre_dat, 
                                     measurement_indices_by_time = NULL,
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
