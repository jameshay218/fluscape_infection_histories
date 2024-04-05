######################################################
## EXTRACT MEASUREMENT OFFSET ESTIMATES
## Author: James Hay
## Date: 21 Feb 2023
## Summary: pulls the virus-specific measurement offset estimates from the previously run chains

library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
library(data.table)
library(foreach)
library(doParallel)
library(parallel)
library(dplyr)
library(tidyr)

#devtools::load_all("~/Documents/GitHub/serosolver")
library(serosolver)

main_wd <- "~/Documents/GitHub/fluscape_infection_histories//"
save_wd <- paste0(main_wd,"/inputs/")

setwd(main_wd)

par_tab <- read.csv("inputs/par_tab_base.csv")
par_tab <- par_tab[par_tab$names != "phi",]
chains <- load_theta_chains("chains/fluscape_offset_estimation2/",burnin=5000000,convert_mcmc=TRUE)

theta_chain <- chains$chain
best_pars <- get_best_pars(theta_chain)


best_rhos <- best_pars[names(best_pars) %like% "rho"]
best_rhos <- best_rhos[!(names(best_rhos) %in% c("rho_mean","rho_sd"))]


write.table(data.frame(names="rho",values=best_rhos,fixed=1,steps=0.1,
                       lower_bound=-3,upper_bound=3,lower_start=-1,upper_start=1,
                       type=3), paste0(save_wd,"par_tab_offset_estimates.csv"), sep=",",row.names=FALSE)

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



