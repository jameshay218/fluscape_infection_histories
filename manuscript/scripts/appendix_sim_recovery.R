## Read and plot against ground truths
## Read in simulated data
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
library(patchwork)
save_plots <- TRUE

serosolver_wd <- "~/Documents/GitHub/serosolver"
devtools::load_all(serosolver_wd)
setwd("~/Documents/GitHub/fluscape_serosolver/")
save_wd <- "~/Documents/GitHub/fluscape_serosolver/manuscript/figures/Appendix/"
theme_main <- theme_bw() + theme(axis.text.x=element_text(size=6),axis.text.y=element_text(size=6),
                                 axis.title.x=element_text(size=7),axis.title.y=element_text(size=7),
                                 strip.text=element_text(size=7))

par_key <- c("error"="ε","mu"="μ_l","mu_short"="μ_s", 
             "sigma1"="σ_l","sigma2"="σ_s",
            "sigma_birth_mod_l"="ν_l","sigma_birth_mod_s"="ν_s",
            "sigma_future_mod_l"="υ_l","sigma_future_mod_s"="υ_s",
            "wane"="ω","tau"="τ","total_infections"="Total infections")

rho_key <- c("ρ", paste0("ρ.",1:200))
names(rho_key) <- c("rho",paste0("rho.",1:200))
par_key <- c(par_key, rho_key)
## Where are ground truths saved?
titre_dat_annual_file <- "data/simulated/sim_quarterly_fluscape_supplement_1_titre_data_annual.csv"
titre_dat_quarterly_file <- "data/simulated/sim_quarterly_fluscape_supplement_1_titre_data.csv"
antigenic_map_annual_file <- "data/antigenic_maps/antigenic_map_fonville_annual_continuous.csv"
antigenic_map_file <- "data/antigenic_maps/antigenic_map_fonville_quarterly_continuous.csv"
par_tab_file <- "data/simulated/sim_quarterly_fluscape_supplement_1_par_tab.csv"

rho_file_true <- "inputs/par_tab_offset_estimates.csv"
rho_file <- "data/simulated/sim_fluscape_offset_estimation_supplement_1_rho_estimates.csv"

## Read in ground-truths
titre_dat <- read.csv(titre_dat_quarterly_file,stringsAsFactors = FALSE)
titre_dat_annual <- read.csv(titre_dat_annual_file,stringsAsFactors = FALSE)
titre_dat <- titre_dat %>% arrange(individual, samples, virus, run)
if(!("group" %in%colnames(titre_dat))){
    titre_dat$group <- 1
}
antigenic_map_annual <- read.csv(antigenic_map_annual_file)
strain_isolation_times_annual <- antigenic_map_annual$inf_times

antigenic_map <- read.csv(antigenic_map_file)
antigenic_map <- antigenic_map[antigenic_map$inf_times >= 1968*4 & antigenic_map$inf_times <= max(titre_dat$samples),]

strain_isolation_times <- antigenic_map$inf_times

true_inf_hist <- read.csv("data/simulated/sim_quarterly_fluscape_supplement_1_infection_histories.csv")

par_tab <- read.csv(par_tab_file)

convert_matrix_to_block_sums <- function(matrix, bucket=4) {
    n <- nrow(matrix)
    m <- ncol(matrix)
    new_matrix <- matrix(0, nrow = n, ncol = ceiling(m/bucket))  # Initialize the new matrix with zeros
    
    for (i in 1:n) {
        for (j in 1:ceiling(m/bucket)) {
            # Calculate the sum of consecutive blocks of 4 entries per row
            start_index <- (bucket * (j - 1) + 1)
            end_index <- min(bucket * j, ncol(matrix))
            block_sum <- sum(matrix[i, start_index:end_index])
            new_matrix[i, j] <- block_sum
        }
    }
    
    return(new_matrix)
}
annual_inf_hist <- convert_matrix_to_block_sums(true_inf_hist)

## Mixed-time true inf_hist is a mixture of the annual and quarterly objects
true_inf_hist_mixed <- cbind(annual_inf_hist[,which(strain_isolation_times_annual < 2005)], true_inf_hist[,which(strain_isolation_times >= 2005*4)])


## Get true AR at quarterly scale
n_infections <- colSums(true_inf_hist)
n_infections_annual <- colSums(annual_inf_hist)
n_alive <- get_n_alive(titre_dat,times=strain_isolation_times)
n_alive_annual <- get_n_alive(titre_dat_annual,times=strain_isolation_times_annual)

true_ar_annual <- data.frame(group=1,j=strain_isolation_times_annual,AR=n_infections_annual/n_alive_annual)
true_ar_quarter <- data.frame(group=1,j=strain_isolation_times,AR=n_infections/n_alive)
true_ar_combined <- bind_rows(true_ar_annual[true_ar_annual$j < 2005,] %>% mutate(j=j*4),
                              true_ar_quarter[true_ar_quarter$j >= 2005*4,])

strain_isolation_times_simulation <- c(seq(1968*4, max(antigenic_map$inf_times)))
measurement_indices <- rep(1:48,each=4)
measurement_indices <- measurement_indices[1:length(strain_isolation_times_simulation)]

fixed <- as.numeric(!(strain_isolation_times_annual %in% unique(titre_dat_annual$virus)))
par_tab[par_tab$names %in% c("rho_mean"),"fixed"] <- 1
par_tab[par_tab$names == "rho","fixed"][fixed == 0] <- 0
#par_tab[par_tab$names == "rho","fixed"] <- 0
par_tab[par_tab$names == "rho","names"] <- paste0("rho.",0:(length(fixed)-1))
par_tab[par_tab$names == "rho.0","names"] <- "rho"
par_tab[par_tab$names == "wane","values"] <- 0.8
par_tab_plot <- par_tab[,c("names","values","fixed","upper_bound","lower_bound")]
par_tab_plot <- bind_rows(par_tab_plot, data.frame(names="total_infections",values=sum(true_inf_hist), fixed=0))


############################################################
## 1. Offset estimation
############################################################
chains <- load_mcmc_chains("~/Documents/GitHub/fluscape_serosolver/chains/sim_fluscape_offset_estimation_supplement_1/",convert_mcmc=FALSE,burnin = 200000,unfixed = FALSE,thin=1)
chain <- as.data.frame(chains$theta_chain)
inf_chain <- chains$inf_chain

true_inf_hist_melted <- reshape2::melt(annual_inf_hist)

p_cumu1 <- generate_cumulative_inf_plots(inf_chain,indivs=1:10,strain_isolation_times=strain_isolation_times_annual,real_inf_hist=annual_inf_hist,number_col=2)
p_cumu1 <- p_cumu1$cumu_infections +
    scale_fill_manual(name="Chain",
                        values=c("1"="green","2"="green",
                                 "3"="green","4"="green","5"="green")) + 
    scale_color_manual(name="Chain",values=c("1"="green","2"="green",
                                             "3"="green","4"="green","5"="green"))

chain <- left_join(chain %>% pivot_longer(-c(sampno,chain_no)) %>% rename(names=name,est_value=value), par_tab_plot)
chain$names <- factor(chain$names, levels=par_tab_plot$names)
chain <- chain %>% filter(fixed == 0)
chain <- chain %>% mutate(set=if_else(names %like% "rho","offsets","parameters"))
chain$names <- as.character(chain$names)
chain <- chain %>% mutate(lower_bound=ifelse(names == "total_infections",0,lower_bound),
                          upper_bound=ifelse(names=="total_infections",prod(dim(true_inf_hist)), upper_bound))
chain$names <- par_key[chain$names]
chain$names <- factor(chain$names, levels=par_key)

p_rho_estimates1 <- ggplot(chain %>% filter(set == "offsets")) + 
    geom_density(aes(x=est_value,fill="Posterior"),alpha=0.25) + 
    geom_vline(aes(xintercept=values,col="True value"))+ 
    geom_vline(xintercept=0,linetype="dashed") +
  scale_x_continuous(limits=c(-1,1)) +
    facet_wrap(~names) +
    ylab("Density") +
    xlab("Value") +
    scale_color_manual(name="Quantity",values=c("True value"="blue")) +
    scale_fill_manual(name="Quantity",values=c("Posterior"="green")) +
    theme_main

p_par_estimates1 <- ggplot(chain %>% filter(set != "offsets")) + 
    #geom_rect(aes(xmin=lower_bound,xmax=upper_bound,ymin=0,ymax=0)) +
    geom_density(aes(x=est_value,fill="Posterior"),alpha=0.25) + 
    geom_vline(aes(xintercept=values,col="True value"))+ 
    facet_wrap(~names,scales="free") +
    ylab("Density") +
    xlab("Value") +
    scale_color_manual(name="Quantity",values=c("True value"="blue")) +
    scale_fill_manual(name="Quantity",values=c("Posterior"="green")) +
    theme_main

if(save_plots){
jahR::save_plots(p_cumu1, save_wd, "rho_estimation_inf_hists",width=6,height=6)
jahR::save_plots(p_rho_estimates1, save_wd, "rho_estimation",width=8,height=5)
jahR::save_plots(p_par_estimates1, save_wd, "rho_estimation_pars",width=8,height=5)
}

############################################################
## 2. Annual estimates without offset terms
############################################################
chains <- load_mcmc_chains("~/Documents/GitHub/fluscape_serosolver/chains/sim_fluscape_annual_no_offsets/",convert_mcmc=FALSE,burnin = 200000,unfixed = FALSE,thin=1)
chain <- as.data.frame(chains$theta_chain)
inf_chain <- chains$inf_chain

chain <- left_join(chain %>% pivot_longer(-c(sampno,chain_no)) %>% rename(names=name,est_value=value), par_tab_plot)
chain$names <- factor(chain$names, levels=par_tab_plot$names)
chain <- chain %>% filter(fixed == 0)
chain$names <- as.character(chain$names)
chain$names <- par_key[chain$names]
chain$names <- factor(chain$names, levels=par_key)
chain_no_offsets <- chain
p_par_estimates2 <- ggplot(chain) + 
    #geom_rect(aes(xmin=lower_bound,xmax=upper_bound,ymin=0,ymax=0)) +
    geom_density(aes(x=est_value,fill="Posterior"),alpha=0.25) + 
    geom_vline(aes(xintercept=values,col="True value"))+ 
    facet_wrap(~names,scales="free") +
    ylab("Density") +
    xlab("Value") +
    scale_color_manual(name="Quantity",values=c("True value"="blue")) +
    scale_fill_manual(name="Quantity",values=c("Posterior"="green")) +
    theme_main+
    labs(tag="A") + ggtitle("Parameter estimates ignoring strain-specific measurement offsets")


p_ar2 <- plot_attack_rates(inf_chain, titre_dat_annual, strain_isolation_times_annual,pad_chain=FALSE,plot_den=FALSE,n_alive = NULL,
                          true_ar=data.frame(group=1,j=strain_isolation_times_annual,AR=n_infections_annual/n_alive_annual),
                          prior_pars = c("prior_version"=2,
                                         "alpha"=1,
                                         "beta"=1)) +
    theme(legend.position="bottom") +
    labs(tag="A") + ggtitle("Estimated attack rates ignoring strain-specific measurement offsets")

p_inf_hists2 <- generate_cumulative_inf_plots(inf_chain, indivs=1:5,real_inf_hist=annual_inf_hist,strain_isolation_times=strain_isolation_times_annual,pad_chain=TRUE)


############################################################
## 3. Annual estimates WITH offset terms
############################################################
chains <- load_mcmc_chains("~/Documents/GitHub/fluscape_serosolver/chains/sim_fluscape_annual/",convert_mcmc=FALSE,burnin = 200000,unfixed = FALSE,thin=1)
chain <- as.data.frame(chains$theta_chain)
inf_chain <- chains$inf_chain

chain <- left_join(chain %>% pivot_longer(-c(sampno,chain_no)) %>% rename(names=name,est_value=value), par_tab_plot)
chain$names <- factor(chain$names, levels=par_tab_plot$names)
chain <- chain %>% filter(fixed == 0) %>% filter(!(names %like% "rho"))
chain$names <- as.character(chain$names)
chain$names <- par_key[chain$names]
chain$names <- factor(chain$names, levels=par_key)
chain_annual <- chain
p_par_estimates3 <- ggplot(chain) + 
    #geom_rect(aes(xmin=lower_bound,xmax=upper_bound,ymin=0,ymax=0)) +
    geom_density(aes(x=est_value,fill="Posterior"),alpha=0.25) + 
    geom_vline(aes(xintercept=values,col="True value"))+ 
    facet_wrap(~names,scales="free") +
    ylab("Density") +
    xlab("Value") +
    scale_color_manual(name="Quantity",values=c("True value"="blue")) +
    scale_fill_manual(name="Quantity",values=c("Posterior"="green")) +
    theme_main+
    labs(tag="B") + ggtitle("Parameter estimates using estimated strain-specific measurement offsets")


chains_offsets_combined <- bind_rows(chain_annual %>% mutate(Model="Using measurement offsets"),
                                     chain_no_offsets%>% mutate(Model="No measurement offsets"))


p_par_estimates_offset_compare <- ggplot(chains_offsets_combined) + 
    geom_density(aes(x=est_value,fill=Model),alpha=0.25) + 
    geom_vline(aes(xintercept=values,col="True value"),linetype="dashed")+ 
    scale_color_manual("",values=c("True value" = "black")) +
    facet_wrap(~names,scales="free") +
    xlab("Value") +
    ylab("Posterior density") +
    theme_main +
    theme(legend.position="bottom")
if(save_plots){
jahR::save_plots(p_par_estimates_offset_compare, save_wd, "compare_offset_par_estimates",width=8,height=6)
}
p_ar3 <- plot_attack_rates(inf_chain, titre_dat_annual, strain_isolation_times_annual,pad_chain=FALSE,plot_den=FALSE,n_alive = NULL,
                           true_ar=data.frame(group=1,j=strain_isolation_times_annual,AR=n_infections_annual/n_alive_annual),
                           prior_pars = c("prior_version"=2,
                                          "alpha"=1,
                                          "beta"=1)) +
    theme(legend.position="bottom") +
    labs(tag="B") + ggtitle("Estimated attack rates using estimated strain-specific measurement offsets")

p_inf_hists3 <- generate_cumulative_inf_plots(inf_chain, indivs=1:5,real_inf_hist=annual_inf_hist,strain_isolation_times=strain_isolation_times_annual,pad_chain=TRUE)

if(save_plots){
jahR::save_plots(p_ar2/p_ar3, save_wd, "no_offsets_ar",width=7,height=8)
jahR::save_plots(p_par_estimates2/p_par_estimates3+ plot_layout(guides="collect"), save_wd, "no_offsets_pars",width=8,height=8)
}

############################################################
## 4. Annual estimates estimating offsets jointly
############################################################
chains <- load_mcmc_chains("~/Documents/GitHub/fluscape_serosolver/chains/sim_fluscape_annual_full/",convert_mcmc=FALSE,burnin = 200000,unfixed = FALSE,thin=1)
chain <- as.data.frame(chains$theta_chain)
inf_chain <- chains$inf_chain

chain <- left_join(chain %>% pivot_longer(-c(sampno,chain_no)) %>% rename(names=name,est_value=value), par_tab_plot)
chain$names <- factor(chain$names, levels=par_tab_plot$names)
chain <- chain %>% filter(fixed == 0)

p_par_estimates4 <- ggplot(chain) + 
    geom_density(aes(x=est_value),fill="blue",alpha=0.25) + 
    geom_vline(aes(xintercept=values,col="True value"))+ 
    facet_wrap(~names,scales="free") +
    theme_main

p_ar4 <- plot_attack_rates(inf_chain, titre_dat_annual, strain_isolation_times_annual,pad_chain=FALSE,plot_den=FALSE,n_alive = NULL,
                           true_ar=data.frame(group=1,j=strain_isolation_times_annual,AR=n_infections_annual/n_alive_annual),
                           prior_pars = c("prior_version"=2,
                                          "alpha"=1,
                                          "beta"=1))

p_inf_hists4 <- generate_cumulative_inf_plots(inf_chain, indivs=1:5,real_inf_hist=annual_inf_hist,strain_isolation_times=strain_isolation_times_annual,pad_chain=TRUE)

if(save_plots){
    jahR::save_plots(p_ar4, save_wd, "joint_ar",width=7,height=4)
    jahR::save_plots(p_par_estimates4, save_wd, "joint_ar_pars",width=8,height=8)
    jahR::save_plots(p_inf_hists4[[1]], save_wd, "joint_ar_inf_hists",width=8,height=8)
}

############################################################
## 5. Quarterly estimates with offsets
############################################################
chains <- load_mcmc_chains("~/Documents/GitHub/fluscape_serosolver/chains/sim_fluscape_main/",convert_mcmc=FALSE,burnin = 200000,unfixed = FALSE,thin=1)
chain <- as.data.frame(chains$theta_chain)
inf_chain <- chains$inf_chain


chain <- left_join(chain %>% pivot_longer(-c(sampno,chain_no)) %>% rename(names=name,est_value=value), par_tab_plot)
chain$names <- factor(chain$names, levels=par_tab_plot$names)
chain <- chain %>% filter(fixed == 0) %>% filter(!(names %like% "rho")) %>% mutate(values = if_else(names == "wane", 0.2, values))
chain_quarter <- chain


chain_quarter$names <- as.character(chain_quarter$names)
chain_quarter$names <- par_key[chain_quarter$names]
chain_quarter$names <- factor(chain_quarter$names, levels=par_key)

p_par_estimates5 <- ggplot(chain) + 
    geom_density(aes(x=est_value),fill="blue",alpha=0.25) + 
    geom_vline(aes(xintercept=values,col="True value"))+ 
    facet_wrap(~names,scales="free") +
    theme_main

p_ar5 <- plot_attack_rates_monthly(inf_chain, titre_dat, strain_isolation_times,pad_chain=FALSE,n_alive = NULL,
                           true_ar=data.frame(group=1,j=strain_isolation_times,AR=n_infections/n_alive))



p_ar5_alt <- plot_attack_rates(inf_chain, titre_dat, strain_isolation_times,pad_chain=FALSE,plot_den=FALSE,n_alive = NULL,
                           true_ar=data.frame(group=1,j=strain_isolation_times,AR=n_infections/n_alive),
                           prior_pars = c("prior_version"=2,
                                          "alpha"=1,
                                          "beta"=1)) +
    theme(legend.position="bottom") +
    labs(tag="B")


p_inf_hists5 <- generate_cumulative_inf_plots(inf_chain, indivs=1:25,number_col=5,real_inf_hist=true_inf_hist,strain_isolation_times=strain_isolation_times,pad_chain=FALSE)
p_inf_hists5_main <- p_inf_hists5[[1]] +
    scale_fill_manual(name="Chain",
                      values=c("1"="green","2"="green",
                               "3"="green","4"="green","5"="green")) + 
    scale_color_manual(name="Chain",values=c("1"="green","2"="green",
                                             "3"="green","4"="green","5"="green"))
############################################################
## 6. Mixed time estimates with offsets
############################################################
chains <- load_mcmc_chains("~/Documents/GitHub/fluscape_serosolver/chains/sim_fluscape_mixed_time2/",convert_mcmc=FALSE,burnin = 200000,unfixed = FALSE,thin=1)
chain <- as.data.frame(chains$theta_chain)
inf_chain <- chains$inf_chain

chain <- left_join(chain %>% pivot_longer(-c(sampno,chain_no)) %>% rename(names=name,est_value=value), par_tab_plot)
chain$names <- factor(chain$names, levels=par_tab_plot$names)
chain <- chain %>% filter(fixed == 0) %>% filter(!(names %like% "rho")) %>% mutate(values = if_else(names == "wane", 0.2, values))


chain_mixed <- chain
chain_mixed$names <- as.character(chain_mixed$names)
chain_mixed$names <- par_key[chain_mixed$names]
chain_mixed$names <- factor(chain_mixed$names, levels=par_key)

p_par_estimates6 <- ggplot(chain) + 
    geom_density(aes(x=est_value),fill="blue",alpha=0.25) + 
    geom_vline(aes(xintercept=values,col="True value"))+ 
    facet_wrap(~names,scales="free") +
    theme_main

strain_isolation_times_mixed <- c(seq(1968*4, 2005*4, by=4),(2005*4+1):max(antigenic_map$inf_times))


p_ar6 <- plot_attack_rates_monthly(inf_chain, titre_dat, strain_isolation_times_mixed,pad_chain=FALSE,n_alive = NULL)
p_ar6_alt <- plot_attack_rates(inf_chain, titre_dat, strain_isolation_times_mixed,pad_chain=FALSE,
                               plot_den=FALSE,n_alive = NULL,
                               true_ar=true_ar_combined,
                               prior_pars = c("prior_version"=2,
                                              "alpha"=1,
                                              "beta"=1)) +
    theme(legend.position="bottom") +
    labs(tag="C")

p_inf_hists6 <- generate_cumulative_inf_plots(inf_chain, indivs=1:25,real_inf_hist=true_inf_hist_mixed,number_col=5,
                                              strain_isolation_times=strain_isolation_times_mixed,pad_chain=FALSE)
p_inf_hists6_main <- p_inf_hists6[[1]] +
    scale_fill_manual(name="Chain",
                      values=c("1"="green","2"="green",
                               "3"="green","4"="green","5"="green")) + 
    scale_color_manual(name="Chain",values=c("1"="green","2"="green",
                                             "3"="green","4"="green","5"="green"))

quarter_x_axis <- scale_x_continuous(limits=c(1968,2016)*4,breaks=c(seq(1970,2016,by=5))*4,
                                     labels=paste0("Q1-",c(seq(1970,2015,by=5))))

quarter_x_axis2 <- scale_x_continuous(limits=c(1968,2016)*4,breaks=c(seq(1975,2016,by=10))*4,
                                     labels=c(seq(1975,2015,by=10)))


p_ar_compare <- (p_ar3 + theme(legend.position="none") + ggtitle("") + labs(tag="A")+ xlab("Time (year)")) + 
    (p_ar5_alt + theme(legend.position="none") + xlab("Time (3-month period)") + quarter_x_axis) +
    (p_ar6_alt  + theme(legend.position="bottom")+ xlab("Time (per year up to 2005, per 3-months thereafter)") + quarter_x_axis) +
    plot_layout(ncol=1,guides="keep")

chains_par_combined <- bind_rows(chain_annual %>% mutate(Model="Annual"),
                                 chain_quarter%>% mutate(Model="Quarter"),
                                 chain_mixed%>% mutate(Model="Mixed"))


p_par_estimates_combined <- ggplot(chains_par_combined %>% 
                               mutate(est_value=if_else(Model =="Annual" & names=="ω",est_value/4,est_value),
                                      values=if_else(Model =="Annual" & names=="ω",values/4,values)
                                      )) + 
    geom_density(aes(x=est_value,fill=Model),alpha=0.25) + 
    geom_vline(aes(xintercept=values,col="True value"),linetype="dashed")+ 
    scale_color_manual("",values=c("True value" = "black")) +
    facet_wrap(~names,scales="free") +
    xlab("Value") +
    ylab("Posterior density") +
    theme_main +
    theme(legend.position="bottom")
if(save_plots){
jahR::save_plots(p_ar_compare, save_wd, "compare_ar",width=8,height=8)
jahR::save_plots(p_par_estimates_combined, save_wd, "compare_par_estimates",width=8,height=6)
jahR::save_plots((p_inf_hists5_main + theme(legend.position="none",axis.text=element_text(size=7)) + 
                      quarter_x_axis2 + labs(tag="A") + scale_y_continuous(breaks=seq(0,12,by=2)))/
                     (p_inf_hists6_main+ theme(legend.position="none",axis.text=element_text(size=7)) + 
                          quarter_x_axis2+ labs(tag="B")+ scale_y_continuous(breaks=seq(0,12,by=2))), save_wd, "compare_inf_hists",width=8,height=10)


}

############################################################
## 7. Mixed time estimates with offsets, Bedford map
############################################################
chains <- load_mcmc_chains("~/Documents/GitHub/fluscape_serosolver/chains/sim_fluscape_wrong_map/",convert_mcmc=FALSE,burnin = 200000,unfixed = FALSE,thin=1)
chain <- as.data.frame(chains$theta_chain)
inf_chain <- chains$inf_chain

chain <- left_join(chain %>% pivot_longer(-c(sampno,chain_no)) %>% rename(names=name,est_value=value), par_tab_plot)
chain$names <- factor(chain$names, levels=par_tab_plot$names)
chain <- chain %>% filter(fixed == 0) %>% filter(!(names %like% "rho")) %>% mutate(values = if_else(names == "wane", 0.2, values))


chain_wrong_map <- chain
chain_wrong_map$names <- as.character(chain_wrong_map$names)
chain_wrong_map$names <- par_key[chain_wrong_map$names]
chain_wrong_map$names <- factor(chain_wrong_map$names, levels=par_key)

p_ar7_alt <- plot_attack_rates(inf_chain, titre_dat, strain_isolation_times_mixed,pad_chain=FALSE,
                               plot_den=FALSE,n_alive = NULL,
                               true_ar=true_ar_combined,
                               prior_pars = c("prior_version"=2,
                                              "alpha"=1,
                                              "beta"=1)) +
    theme(legend.position="none") +
    labs(tag="B")

p_inf_hists7 <- generate_cumulative_inf_plots(inf_chain, indivs=1:5,real_inf_hist=true_inf_hist,strain_isolation_times=strain_isolation_times,pad_chain=FALSE)


############################################################
## 8. Mixed time estimates with offsets, smooth map
############################################################
chains <- load_mcmc_chains("~/Documents/GitHub/fluscape_serosolver/chains/sim_fluscape_smooth_map/",convert_mcmc=FALSE,burnin = 200000,unfixed = FALSE,thin=1)
chain <- as.data.frame(chains$theta_chain)
inf_chain <- chains$inf_chain

chain <- left_join(chain %>% pivot_longer(-c(sampno,chain_no)) %>% rename(names=name,est_value=value), par_tab_plot)
chain$names <- factor(chain$names, levels=par_tab_plot$names)
chain <- chain %>% filter(fixed == 0) %>% filter(!(names %like% "rho")) %>% mutate(values = if_else(names == "wane", 0.2, values))


chain_smooth_map <- chain
chain_smooth_map$names <- as.character(chain_smooth_map$names)
chain_smooth_map$names <- par_key[chain_smooth_map$names]
chain_smooth_map$names <- factor(chain_smooth_map$names, levels=par_key)

p_ar8_alt <- plot_attack_rates(inf_chain, titre_dat, strain_isolation_times_mixed,pad_chain=FALSE,
                               plot_den=FALSE,n_alive = NULL,
                               true_ar=true_ar_combined,
                               prior_pars = c("prior_version"=2,
                                              "alpha"=1,
                                              "beta"=1)) +
    theme(legend.position="none") +
    labs(tag="B")

p_inf_hists8 <- generate_cumulative_inf_plots(inf_chain, indivs=1:5,real_inf_hist=true_inf_hist,strain_isolation_times=strain_isolation_times,pad_chain=FALSE)



############################################################
## Compare estimates using different maps
############################################################

p_ar_compare_maps <- (p_ar6_alt  + labs(tag="A") + theme(legend.position="none")+ xlab("Time (per year up to 2005, per 3-months thereafter)") + quarter_x_axis) +
    (p_ar7_alt  + labs(tag="B") + theme(legend.position="none")+ xlab("Time (per year up to 2005, per 3-months thereafter)") + quarter_x_axis) +
    (p_ar8_alt  + labs(tag="C") + theme(legend.position="bottom")+ xlab("Time (per year up to 2005, per 3-months thereafter)") + quarter_x_axis) +
    plot_layout(ncol=1,guides="keep")


chains_par_combined_maps <- bind_rows(chain_mixed%>% mutate(Model="Correct map"),
                                      chain_wrong_map %>% mutate(Model = "Bedford map"),
                                      chain_smooth_map %>% mutate(Model = "Smoothed map"))


p_par_estimates_combined_map <- ggplot(chains_par_combined_maps) + 
    geom_density(aes(x=est_value,fill=Model),alpha=0.25) + 
    geom_vline(aes(xintercept=values,col="True value"),linetype="dashed")+ 
    scale_color_manual("",values=c("True value" = "black")) +
    facet_wrap(~names,scales="free") +
    xlab("Value") +
    ylab("Posterior density") +
    theme_main +
    theme(legend.position="bottom")
if(save_plots){
jahR::save_plots(p_ar_compare_maps, save_wd, "compare_ar_maps",width=8,height=8)
jahR::save_plots(p_par_estimates_combined_map, save_wd, "compare_par_estimates_maps",width=8,height=6)
}