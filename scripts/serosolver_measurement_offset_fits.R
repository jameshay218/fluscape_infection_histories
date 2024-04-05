######################################################
## FIT SEROSOLVER MODEL TO FLUSCAPE DATA WITH MEASUREMENT OFFSETS
## Author: James Hay
## Date: 20 OCTOBER 2023
## Summary: fits the annual infection history and antibody kinetics model to the Fluscape data to
## estimate the virus-specific measurement offset terms, using prior version 3.
## Note the choice of antigenic map -- this version should use the clustered map

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
#library(serosolver)

run_name <- "fluscape_offset_estimation2"
main_wd <- "~/Documents/GitHub/fluscape_infection_histories//"
chain_wd <- paste0(main_wd,"/chains/",run_name)
save_wd <- paste0(main_wd,"/figures/chain_plots/",run_name)
buckets <- 1
prior_version <- 3
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

mcmc_pars <- c("save_block"=100,"thin"=1000,"thin_hist"=10000,
               "iterations"=1000000,
  "adaptive_period"=3000000,
  "burnin"=0,"switch_sample"=2,"hist_switch_prob"=0,
  "year_swap_propn"=0.5,"swap_propn"=0.5,
  "inf_propn"=0.5,"hist_sample_prob"=1,"move_size"=3, "hist_opt"=1,
  "popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,propose_from_prior=TRUE)


titre_dat <- read.csv("data/fluscape_data_1_resolution.csv",stringsAsFactors = FALSE)
titre_dat <- titre_dat %>% arrange(individual, samples, virus, run)
titre_dat$group <- 1

## Antigenic map
antigenic_map <- read.csv("data/antigenic_maps/antigenic_map_fonville_annual_continuous.csv")
antigenic_map <- antigenic_map[antigenic_map$inf_times <= max(titre_dat$samples),]
strain_isolation_times <- antigenic_map$inf_times

## Set up parameter table
par_tab <- read.csv("inputs/par_tab_base.csv",stringsAsFactors=FALSE)
par_tab <- par_tab[par_tab$names != "phi",]
par_tab[par_tab$names %in% c("alpha","beta"),c("values")] <- c(1,1)
par_tab[par_tab$names == "obs_sd",c("names")] <- "error"

## Add entries for each measurement offset
measurement_indices <- match(antigenic_map$inf_times, antigenic_map$inf_times)
for(i in 1:length(measurement_indices)){
    par_tab <- rbind(par_tab, data.frame("names"="rho","values"=rnorm(1,0,1),"fixed"=1,"steps"=0.1,
                                         "lower_bound"=-3,"upper_bound"=3,"lower_start"=-2,"upper_start"=2,"type"=3))
}
fixed <- as.numeric(!(strain_isolation_times %in% unique(titre_dat$virus)))
par_tab[par_tab$names %in% c("rho_mean"),"fixed"] <- 1
par_tab[par_tab$names == "rho_sd","fixed"] <- 0
par_tab[par_tab$names== "rho","fixed"] <- 1
par_tab[par_tab$names == "rho",][which(fixed == 0),"fixed"] <- 0
par_tab[par_tab$names == "rho" & par_tab$fixed == 1, "values"] <- 0


create_prior_func <- function(par_tab){
    prior_function <- function(pars){
        p1 <- dnorm(pars["mu"],0, 100,log=TRUE)
        p2 <- dnorm(pars["mu_short"],0, 100,log=TRUE)
        p3 <- dbeta(pars["wane"],1, 1,log=TRUE)
        p4 <- dbeta(pars["tau"],1, 1,log=TRUE)
        p5 <- dbeta(pars["sigma1"],1, 1,log=TRUE)
        p6 <- dbeta(pars["sigma2"],1, 1,log=TRUE)
        p7 <- dnorm(pars["error"],0, 100,log=TRUE)
        p1 + p2 + p3 + p4 + p5 + p6 + p7
    }
}

f <- create_posterior_func(par_tab,titre_dat,antigenic_map=antigenic_map,strain_isolation_times = strain_isolation_times,
                           version=prior_version,
                           solve_likelihood=TRUE,n_alive=NULL,
                           measurement_indices_by_time=measurement_indices)

## Time runs and use dopar to run multiple chains in parallel
print("Starting fits...")
t1 <- Sys.time()
filenames <- paste0(chain_wd, "/",run_name, "_",1:n_chains)
res <- foreach(x = filenames, .packages = c('data.table','plyr',"dplyr")) %dopar% {
    devtools::load_all(serosolver_wd)
    
    index <- 1
    lik <- -Inf
    inf_hist_correct <- 1
    while((!is.finite(lik) || inf_hist_correct > 0) & index < 100){
        start_tab <- generate_start_tab(par_tab)
        start_inf <- setup_infection_histories_titre(titre_dat,strain_isolation_times,space=5,titre_cutoff=2,sample_prob=0.9)
        
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
                                measurement_indices=measurement_indices,
                                measurement_random_effects = TRUE,
                                mcmc_pars=mcmc_pars,
                                solve_likelihood=TRUE)
}
run_time_fast <- Sys.time() - t1
run_time_fast

## Read in chains for trace plot
chains <- load_mcmc_chains(chain_wd,convert_mcmc=TRUE,
                           par_tab=par_tab,
                           burnin = mcmc_pars["adaptive_period"],unfixed = TRUE,thin=1)
pdf(paste0(save_wd,"/",run_name,"_chain.pdf"))
plot(as.mcmc.list(chains$theta_list_chains))
dev.off()

## Read in chains for all other plots
chains <- load_mcmc_chains(chain_wd,convert_mcmc=FALSE,burnin = mcmc_pars["adaptive_period"],unfixed = FALSE,thin=10)
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
