#library(serosolver)
devtools::load_all("~/Documents/GitHub/serosolver")
library(coda)
library(tidyverse)
library(data.table)
setwd("~/Documents/GitHub/fluscape_infection_histories//chains/sim_quarterly_fluscape_supplement_2_smooth//")

burnin <- 500000

chains <- load_mcmc_chains(convert_mcmc=TRUE,unfixed=TRUE,thin = 1,burnin=0)
start_tab <- load_start_tab()
sim_recovery <- TRUE
est_rhos <- FALSE

if(sim_recovery){
  true_par_tab <- read.csv("~/Documents/GitHub/fluscape_infection_histories//data/simulated/sim_quarterly_fluscape_supplement_2_par_tab.csv")
  true_total_inf <- read_csv("~/Documents/GitHub/fluscape_infection_histories//data/simulated/sim_quarterly_fluscape_supplement_2_infection_histories.csv")
}

use_names <- c("sampno", 
               "mu","wane","mu_short","sigma1","sigma2","error","tau",
               "total_infections","chain_no","lnlike", "likelihood","prior_prob")

if(est_rhos){
  par_tab_rhos <- start_tab[start_tab$names == 'rho',]
  rho_indices <- 1:nrow(par_tab_rhos)
  rho_indices <- rho_indices[which(par_tab_rhos$fixed == 0)] - 1
  rho_indices <- rho_indices[rho_indices > 0]
  use_names <- c(use_names,"rho",paste0("rho.",rho_indices))
  
  if(sim_recovery){
    true_par_tab_rho <- true_par_tab %>% 
      dplyr::filter(names == "rho") %>% 
      dplyr::mutate(i = 1:n()) %>% 
      dplyr::mutate(names = paste0("rho.",i-1)) %>%
      dplyr::mutate(names = if_else(names=="rho.0","rho",names))
    true_par_tab <- bind_rows(true_par_tab %>% filter(names != "rho"), true_par_tab_rho)
  }
}
if(sim_recovery){
  true_par_tab <- true_par_tab[true_par_tab$fixed == 0,]
  true_par_tab <- true_par_tab %>%bind_rows(tibble(names="total_infections",values=sum(true_total_inf),fixed=0,steps=0.1,lower_bound=0,upper_bound=prod(dim(true_total_inf)), lower_start=0,upper_start=1,type=1))
}

max(chains$theta_chain[,"sampno"])
p_trace <- chains$theta_chain %>% as_tibble() %>% 
    select(`use_names`) %>% 
    dplyr::filter(sampno > burnin) %>%
  filter(chain_no != 5) %>%
    pivot_longer(-c(sampno, chain_no)) %>%
     mutate(name = factor(name, levels=use_names)) %>%
    ggplot() + geom_line(aes(x=sampno,y=value,col=as.factor(chain_no))) + 
    facet_wrap(~name,scales="free_y")

p_density <- chains$theta_chain %>% as_tibble() %>% 
    select(`use_names`) %>% 
    filter(sampno > burnin) %>%
  filter(chain_no != 5) %>%
    pivot_longer(-c(sampno, chain_no)) %>%
    mutate(name = factor(name, levels=use_names)) %>%
    ggplot() + geom_density(aes(x=value,fill=as.factor(chain_no)),alpha=0.25) + 
    facet_wrap(~name,scales="free")

if(sim_recovery){
  p_trace <- p_trace + geom_hline(data=true_par_tab[true_par_tab$names %in% use_names,] %>% dplyr::rename(name=names),aes(yintercept=values))
  p_density <- p_density + geom_vline(data=true_par_tab[true_par_tab$names %in% use_names,] %>% dplyr::rename(name=names),aes(xintercept=values),linewidth=1,color="red")
}
p_trace
p_density


chains <- load_theta_chains(convert_mcmc=TRUE,unfixed=TRUE,thin = 1,burnin =burnin)
chains1 <- chains$list[1:4]
chains1 <- lapply(chains1, function(x) as.mcmc(x[,c("mu","mu_short","tau","wane",
                                            "sigma1","sigma2","error")]))
gelman.diag(chains1)

