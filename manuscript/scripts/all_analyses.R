########################################
## RUN THROUGH ALL EXPLORATORY ANALYSES 
## Purpose: source all scripts and generate plots required for raw data analyses
## Date: 15 March 2024
## Author: James Hay
########################################
library(coda)
library(data.table)
library(plyr)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(gtable)
library(mgcv)
library(itsadug)
library(tidyverse)
library(patchwork)
library(viridis)
library(ggpubr)
library(NatParksPalettes)
library(lattice)
library(sp)
library(raster)
#library(maptools)
library(ncf)
library(Hmisc)
library(corrplot)
library(gganimate)
library(cowplot)
library(itsadug)
select <- dplyr::select
mutate <- dplyr::mutate
rename <- dplyr::rename

theme_main <- theme_bw()

ggsave_jah <- function(p, wd, save_name, width, height){
    message("Saving to ", paste0(wd, save_name,".pdf/png"))
    ggsave(p, filename=paste0(wd, save_name,".pdf"),width=width,height=height)
    ggsave(p, filename=paste0(wd, save_name,".png"),width=width,height=height,units='in',dpi=300)
    message("Success")
}

convert_years_to_quarters <- function(times){
    times <- sapply(times, function(x){
        if(x %% 1 == 0){
            paste0("Q1-",x)
        } else if(x %% 1 == 0.25){
            paste0("Q2-",x-0.25)
        } else if(x%%1 == 0.5){
            paste0("Q3-",x-0.5)
        } else {
            paste0("Q4-",x-0.75)
        }
    })
    times
}



devtools::load_all("~/Documents/GitHub/serosolver")
#library(serosolver)

main_wd <- "~/Documents/GitHub/fluscape_infection_histories/manuscript/"
figure_wd <- paste0(main_wd, "/figures/")
setwd(main_wd)
source(paste0(main_wd,"scripts/aux/function_plot.R"))

## Where are the MCMC chains saved?
chain_wd<- "~/Documents/GitHub/fluscape_serosolver/chains/fluscape_main/"
chain_wd_alternative <- "~/Documents/GitHub/fluscape_serosolver/chains/fluscape_offset_estimation2/"
use_alternative_prior <- FALSE
use_base_results <- TRUE

sim_chain_wd <- "~/Documents/GitHub/fluscape_serosolver/chains/sim_fluscape_main_supplement_1/"
sim_chain_wd_misspecified <- "~/Documents/GitHub/fluscape_serosolver/chains/sim_fluscape_no_offsets_supplement_1/"
sim_chain_wd_wrong_map <- "~/Documents/GitHub/fluscape_serosolver/chains/sim_fluscape_wrong_map_supplement_1//"
sim_chain_wd_smooth_map <- "~/Documents/GitHub/fluscape_serosolver/chains/sim_quarterly_fluscape_supplement_2_smooth_tmp//"
sim_chain_wd_smooth_map_correct <-  "~/Documents/GitHub/fluscape_serosolver/chains/sim_quarterly_fluscape_supplement_2_main/"
base_chain_wd <- "~/Documents/GitHub/fluscape_serosolver/chains/fluscape_base/"
vietnam_chain_wd <- "~/Documents/GitHub/fluscape_serosolver/chains/ha_nam/"

## Where is the fluscape Git repo saved?
fluscape_wd <- "~/Documents/GitHub/Fluscape/"
#data_file_loc <- "~/Google Drive/Influenza/serosolver/all_data/fluscape_data_4.csv"
data_file_loc <- "~/Documents/GitHub/fluscape_infection_histories/data/fluscape_data_4_resolution.csv"
buckets <- 4 ## Set to 1 for annual version, 4 for quarterly version
burnin_use <- 20000000
burnin_use_alternative <- 2000000
remove_duplicate_infections <- FALSE
rerun_regressions <- FALSE
#################################
##  Load in ALL of the fluscape data
##     - diseases
##     - fluscape_dat
##     - household_dat_all
##     - ili_last_year
##     - ili_last_visit
##     - part_info_bmi
##     - smoking_all
##     - vaccine_last_year
##     - vaccine_last_visit
if(!file.exists("r_data/fluscape_dat.RData")){
    source("scripts/aux/extract_fluscape_dat.R")
  ## Run the age randomisation, Lat/Long jittering code now
  ## Density information of locations also removed
  save(fluscape_dat, file="r_data/fluscape_dat.RData")
} else {
    load("r_data/fluscape_dat.RData")
}
needed_fluscape_dat <- unique(fluscape_dat[,c("individual","Participant_ID","visit", "dens.1","dens.9","DOB","order","age_group","age","birth_cohort")])
#################################
## Run GAM models on titres
## This won't work without the main fluscape git repo
if(rerun_regressions){
    source("scripts/aux/data_regressions.R")
}
break
#################################
## Load in MCMC chains - only do this once, as slow
## Objects that come from this:
##     - inf_chain1/inf_chain: the full infection history chain
##     - theta_chain: the full mcmc chain for antibody kinetics parameters
if(!file.exists("r_data/measurement_quarter.RData")){
    source("scripts/aux/extract_infection_histories.R")
  
    if(remove_duplicate_infections){
      ## Find consecutive infections
      double_infs <- inf_chain %>% 
        arrange(chain_no, sampno, i, j) %>% 
        group_by(chain_no, sampno, i) %>% 
        mutate(infection_new = ifelse(!is.na(x) & x == 1 & !is.na(lag(x,1)) & lag(x, 1) == 1,0, x))
      
      double_infs <- double_infs %>% ungroup() %>% 
        group_by(chain_no, sampno, i) %>% 
        mutate(start_inf = infection_new == 1 & x == 1 & (!is.na(lead(x,1)) & lead(x,1) == 1))
      
      repeat_infs_summary <- double_infs %>% ungroup() %>% filter(x == 1) %>% group_by(i, sampno, chain_no) %>% dplyr::summarize(n_infs = sum(x),n_distinct_infs=sum(infection_new), n_run_infections = sum(start_inf))
      repeat_infs_summary <- repeat_infs_summary %>% mutate(prop_repeat = n_run_infections/n_infs)
     
      ## Remove runs of consecutive infections
      removed_double_infs <- inf_chain %>% 
        arrange(chain_no, sampno, i, j) %>% 
        group_by(chain_no, sampno, i) %>% 
        mutate(infection_new = ifelse(!is.na(x) & x == 1 & !is.na(lag(x,1)) & lag(x, 1) == 1, 
                                      0, x)) %>%
        dplyr::select(-x) %>% rename(x=infection_new)
      inf_chain_old <- inf_chain
      inf_chain <- as.data.table(removed_double_infs %>% ungroup())
      inf_chain <- inf_chain[,c("i","j","x","sampno","chain_no")]
      rm(removed_double_infs)
      rm(inf_chain_old)
    }
    
    save(inf_chain, theta_chain, par_tab, antigenic_map, titre_dat, strain_isolation_times,
         file=paste0(main_wd,"/r_data/measurement_quarter.RData"))
} else {
    load("r_data/measurement_quarter.RData")
    colnames(antigenic_map)[3] <- "inf_times"
}

## Load MCMC chain for base results
if(use_base_results){
  chain_wd <- base_chain_wd
  
  if(!file.exists("r_data/base_quarter.RData")){
  source("scripts/aux/extract_infection_histories.R")
  
  if(remove_duplicate_infections ){
    ## Remove runs of consecutive infections
    removed_double_infs <- inf_chain %>% 
      arrange(chain_no, sampno, i, j) %>% 
      group_by(chain_no, sampno, i) %>% 
      mutate(infection_new = ifelse(!is.na(x) & x == 1 & !is.na(lag(x,1)) & lag(x, 1) == 1, 
                                    0, x)) %>%
      dplyr::select(-x) %>% rename(x=infection_new)
    inf_chain_old <- inf_chain
    inf_chain <- as.data.table(removed_double_infs %>% ungroup())
    inf_chain <- inf_chain[,c("i","j","x","sampno","chain_no")]
    rm(removed_double_infs)
    rm(inf_chain_old)
  }
  
  save(inf_chain, theta_chain, par_tab, antigenic_map, titre_dat, strain_isolation_times,
       file=paste0(main_wd,"/r_data/base_quarter.RData"))
} else {
  load("r_data/base_quarter.RData")
  colnames(antigenic_map)[3] <- "inf_times"
}
}


#################################
## Load in MCMC chains for alternative prior
if(use_alternative_prior){
  if(!file.exists("r_data/alternative_prior.RData")){
    burnin_use <- burnin_use_alternative
    chain_wd <- chain_wd_alternative
    
    source("scripts/aux/extract_infection_histories.R")
   
    save(inf_chain, theta_chain, par_tab, antigenic_map, titre_dat, strain_isolation_times,
         file=paste0(main_wd,"/r_data/alternative_prior.RData"))
  } else {
    load("r_data/alternative_prior.RData")
    colnames(antigenic_map)[3] <- "inf_times"
  }
}

## Supplementary figure with individual fits, and also residuals and RMSE
source("scripts/aux/plot_individual_titre_fits.R")
source("scripts/aux/plot_individual_titre_fits_meas.R")

## Look at the distribution of infection "runs"
source("scripts/aux/infection_run_summaries.R")



## Plot some MCMC convergence statistics
## Assess convergence
theta_chains <- load_theta_chains(location=chain_wd,thin=1,burnin=20000000,par_tab=par_tab,unfixed=TRUE,convert_mcmc=TRUE)
effectiveSize(theta_chains$list)
gelman.diag(lapply(theta_chains$list, function(x) x[,2:(ncol(x)-4)]))
plot_posteriors_theta(as.data.frame(theta_chains$chain),par_tab,samples=100)
p_inf_chain_mcmc <- plot_infection_history_chains_time(inf_chain,years=1:25)
p_inf_chain_mcmc <- plot_infection_history_chains_time(inf_chain,years=61:80,pad_chain=FALSE)
p_inf_chain_mcmc[[1]]

## Plot attack rates and individual infection probabilities
source("scripts/aux/plot_individual_infections.R")

## Proportion seroconverted
seroconv <- cbind(year=seroconv_by_year[,1],signif(Hmisc::binconf(seroconv_by_year[,2],seroconv_by_year[,3])*100,3))
seroconv[seroconv[,1] %in% c(2010,2012,2014),] %>% as_tibble() %>% mutate(res=paste0(PointEst, "% (",Lower,"%-",Upper,"%)"))


## Attack rates per year and reinfected
source("scripts/aux/get_attack_rates.R")
res_attack_rates %>% mutate(year=c(2010,2011,2012,2013,2014))
res_reinfected
res_total_infections %>% mutate(infections=infections/100)

############################
## Convert infection histories to annual buckets


## Spatial analysis of attack rates
## Plot attack rates ordered by location
## Spline correlograms for how attack rates are correlated in space
source("scripts/aux/plot_attack_rates_loc.R")

## Create gif
source("scripts/aux/spatial_attack_rates_gif.R")

## plot_attack_rates_age is by-age stats
source("scripts/aux/plot_attack_rates_age.R")
source("scripts/aux/plot_infection_age_trends.R")
source("scripts/aux/plot_infection_age_trends_alternative_prior.R")

## Titre reconstruction
source("scripts/aux/reconstruct_titers.R")

## Get antibody kinetics parameters
## note that ESS is smaller than it should be because only a subset of posterior samples are kept
## from loading the chains
## Convert waning rate to years
theta_chain %>% mutate(wane=1/(wane*4)) %>% pull(wane) %>% mean
theta_chain %>% mutate(wane=1/(wane*4)) %>% pull(wane) %>% quantile(0.025)
theta_chain %>% mutate(wane=1/(wane*4)) %>% pull(wane) %>% quantile(0.975)
par_estimates <- plot_posteriors_theta(theta_chain %>% mutate(wane = wane*mu_short*4), par_tab,plot_corr=FALSE,plot_mcmc=FALSE)$results
par_key <- c("mu"="µl","mu_short"="µs","sigma1"="σl","sigma2"="σs","wane"="ω","tau"="τ","error"="ε")
par_estimates <- par_estimates %>% filter(names %in% c("mu","mu_short","sigma1","sigma2","wane","tau","error")) %>% select(names, estimate) %>% rename("Parameter"="names","Estimate"="estimate")
par_estimates$Parameter <- par_key[par_estimates$Parameter]
write_csv(par_estimates,file="results/fluscape_parameter_estimates.csv")


## Plot antigenic map used
spar_fon <- 0.3
antigenic_coords_fonville <- read_csv("~/Documents/GitHub/fluscape_infection_histories/data/antigenic_map_coords/fonville_antigenic_coordinates.csv")

antigenic_coords_fonville$virus <- factor(antigenic_coords_fonville$virus,levels=antigenic_coords_fonville$virus)
fit_dat_fonville <- generate_antigenic_map_flexible(antigenic_coords_fonville,buckets=4,clusters=NULL,use_clusters=FALSE,spar_fon)
p_antigenic_map <- ggplot() + 
    geom_line(data=fit_dat_fonville,aes(x=x_coord,y=y_coord)) +
    geom_point(data=fit_dat_fonville,aes(x=x_coord,y=y_coord),size=0.5) +
    geom_text(data=fit_dat_fonville%>%filter(inf_times %in% range(fit_dat_fonville$inf_times)),
              aes(x=x_coord-2,y=y_coord,label=floor(inf_times/4))) +
    geom_point(data=antigenic_coords_fonville,aes(x=X,y=Y,col=virus),size=1.5) +
    geom_label(data=antigenic_coords_fonville,aes(x=X-2.5,y=Y,label=virus,col=virus),size=3) +
    scale_color_viridis_d(name="Strain") +
    xlab("Antigenic distance (x-coordinate)") +
    ylab("Antigenic distance (y-coordinate)") +
    theme_main
ggsave_jah(p_antigenic_map, figure_wd, "antigenic_map",width=7,height=4)

############################
## Vietnam fitting results
## Load in MCMC chains - only do this once, as slow
## Objects that come from this:
##     - inf_chain1/inf_chain: the full infection history chain
##     - theta_chain: the full mcmc chain for antibody kinetics parameters
buckets <- 1 ## Set to 1 for annual version, 4 for quarterly version
if(!file.exists("r_data/vietnam_fits.RData")){
    chain_wd <- vietnam_chain_wd
    burnin_use <- 1500000
    source("scripts/aux/extract_infection_histories.R")
    save(inf_chain, theta_chain, par_tab, antigenic_map, titre_dat, strain_isolation_times,
         file=paste0(main_wd,"r_data/vietnam_fits.RData"))
} else {
    load("r_data/vietnam_fits.RData")
    colnames(antigenic_map)[3] <- "inf_times"
}
## Get antibody kinetics parameters
## note that ESS is smaller than it should be because only a subset of posterior samples are kept
## from loading the chains
## How long for short-term response to subside?
theta_chain %>% mutate(wane=1/wane) %>% pull(wane) %>% mean
theta_chain %>% mutate(wane=1/wane) %>% pull(wane) %>% quantile(0.025)
theta_chain %>% mutate(wane=1/wane) %>% pull(wane) %>% quantile(0.975)

par_estimates_vietnam <- plot_posteriors_theta(theta_chain %>% mutate(wane=mu_short*wane), par_tab,plot_corr=FALSE,plot_mcmc=FALSE)$results
par_key <- c("mu"="µl","mu_short"="µs","sigma1"="σl","sigma2"="σs","wane"="ω","tau"="τ","error"="ε")
par_estimates_vietnam <- par_estimates_vietnam %>% filter(names %in% c("mu","mu_short","sigma1","sigma2","wane","tau","error")) %>% select(names, estimate) %>% rename("Parameter"="names","Estimate"="estimate")
par_estimates_vietnam$Parameter <- par_key[par_estimates_vietnam$Parameter]
write_csv(par_estimates_vietnam,file="results/vietnam_parameter_estimates.csv")

vietnam_dob <- read.csv(paste0(main_wd,"../data/vietnam_ages.csv"))%>% mutate(individual = 1:n()) %>% 
  left_join(titre_dat %>% select(-DOB))
n_alive_group_vietnam <- get_n_alive_group(vietnam_dob, strain_isolation_times)
ars_vietnam <- get_attack_rates(inf_chain, strain_isolation_times, titre_dat, n_alive_group_vietnam,FALSE)

## Get comparable AR for fluscape:
inf_chain %>% mutate(j = j + 1967) %>% group_by(j,sampno, chain_no) %>% dplyr::summarize(infected=sum(x)) %>% left_join(data.frame(j=1968:2012, n_alive=n_alive_group_vietnam[1,])) %>% mutate(prop = infected/n_alive) %>% group_by(sampno, chain_no) %>% dplyr::summarize(median_ar = median(prop)) %>% ungroup() %>% dplyr::summarize(med_med=median(median_ar),lower=quantile(median_ar,0.025), upper=quantile(median_ar,0.975))

ggplot(ars_vietnam) + geom_pointrange(aes(x=j,ymin=lower,ymax=upper,y=median))
ar_estimates_annual

p_ar_comparison <- ar_estimates_annual %>%  rename(fluscape_ar = med_ar,fluscape_lower=lower,fluscape_upper=upper) %>% 
  left_join(ars_vietnam %>% rename(vietnam_ar=median,vietnam_lower=lower,vietnam_upper=upper)) %>%
  ggplot() + 
  geom_rect(data=data.frame(xmin=1967.5,xmax=1968.5,ymin=-0.2,ymax=1.1),fill="green",alpha=0.25,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)) +
  geom_rect(data=data.frame(xmin=1984.5,xmax=1985.5,ymin=-0.2,ymax=1.1),fill="yellow",alpha=0.25,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)) +
  geom_rect(data=data.frame(xmin=1991.5,xmax=1992.5,ymin=-0.2,ymax=1.1),fill="yellow",alpha=0.25,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)) +
  geom_rect(data=data.frame(xmin=1994.5,xmax=1997.5,ymin=-0.2,ymax=1.1),fill="yellow",alpha=0.25,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)) +
  geom_rect(data=data.frame(xmin=1999.5,xmax=2004.5,ymin=-0.2,ymax=1.1),fill="yellow",alpha=0.25,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)) +
  geom_rect(data=data.frame(xmin=1988.5,xmax=1989.5,ymin=-0.2,ymax=1.1),fill="green",alpha=0.25,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)) +
  geom_rect(data=data.frame(xmin=2008.5,xmax=2009.5,ymin=-0.2,ymax=1.1),fill="green",alpha=0.25,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)) +
  geom_rect(data=data.frame(xmin=2009.5,xmax=2012.5,ymin=-0.2,ymax=1.1),fill="yellow",alpha=0.25,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)) +
  geom_pointrange(aes(x=j-0.1,ymin=fluscape_lower,ymax=fluscape_upper,y=fluscape_ar,col="Guangzhou, China"),size=0.2,fatten=0.2) +
  geom_pointrange(aes(x=j+0.1,ymin=vietnam_lower,ymax=vietnam_upper,y=vietnam_ar,col="Ha Nam, Vietnam"),size=0.2,fatten=0.2) +
  theme_classic() +
  coord_cartesian(ylim=c(0,1))+
  scale_color_manual(name="Location",values=c("Guangzhou, China"="red","Ha Nam, Vietnam"="grey30")) +
  scale_y_continuous(expand=c(.01,.01),breaks=seq(0,1,by=0.2)) +
  scale_x_continuous(breaks=seq(1970,2015,by=5))+
  ylab("Estimated proportion infected\n at least once per year") +
  xlab("Date (years)") +
  theme(legend.position="bottom")
ggsave_jah(p_ar_comparison,wd = "figures/",save_name = "compare_ars",width=8,height=3)
## Plot attack rates and individual infection probabilities
source("scripts/aux/plot_individual_titre_fits_vietnam.R")

############################
## Simulation recovery results
## Load in MCMC chains - only do this once, as slow
## Objects that come from this:
##     - inf_chain1/inf_chain: the full infection history chain
##     - theta_chain: the full mcmc chain for antibody kinetics parameters
if(!file.exists("r_data/sim_measurement_quarter.RData")){
    chain_wd <- sim_chain_wd
    burnin_use <- 500000
    source("scripts/aux/extract_infection_histories.R")
    save(inf_chain, theta_chain, par_tab, antigenic_map, titre_dat, strain_isolation_times,
         file=paste0(main_wd,"r_data/sim_measurement_quarter.RData"))
} else {
    load("r_data/sim_measurement_quarter.RData")
    colnames(antigenic_map)[3] <- "inf_times"
}

par_estimates_sim <- plot_posteriors_theta(theta_chain %>% mutate(wane = wane*mu_short), par_tab,plot_corr=FALSE,plot_mcmc=FALSE)$results
par_key <- c("mu"="µl","mu_short"="µs","sigma1"="σl","sigma2"="σs","wane"="ω","tau"="τ","error"="ε")
par_estimates_sim <- par_estimates_sim %>% filter(names %in% c("mu","mu_short","sigma1","sigma2","wane","tau","error")) %>% select(names, estimate) %>% rename("Parameter"="names","Estimate"="estimate")
par_estimates_sim$Parameter <- par_key[par_estimates_sim$Parameter]
write_csv(par_estimates_sim,file="results/sim_parameter_estimates.csv")


## Plot recovered infection histories
## Plot simulated data and fits
source("scripts/aux/plot_individual_titre_fits_sim.R")

if(!file.exists("r_data/sim_measurement_quarter_misspecified.RData")){
  chain_wd <- sim_chain_wd_misspecified
  burnin_use <- 500000
  source("scripts/aux/extract_infection_histories.R")
  save(inf_chain, theta_chain, par_tab, antigenic_map, titre_dat, strain_isolation_times,
       file=paste0(main_wd,"r_data/sim_measurement_quarter_misspecified.RData"))
} else {
  load("r_data/sim_measurement_quarter_misspecified.RData")
  colnames(antigenic_map)[3] <- "inf_times"
}
sim_chain_wd <- sim_chain_wd_misspecified
source("scripts/aux/plot_individual_titre_fits_sim_misspecified.R")


if(!file.exists("r_data/sim_measurement_quarter_wrong_map.RData")){
  chain_wd <- sim_chain_wd_wrong_map
  burnin_use <- 500000
  source("scripts/aux/extract_infection_histories.R")
  save(inf_chain, theta_chain, par_tab, antigenic_map, titre_dat, strain_isolation_times,
       file=paste0(main_wd,"r_data/sim_measurement_quarter_wrong_map.RData"))
} else {
  load("r_data/sim_measurement_quarter_wrong_map.RData")
  colnames(antigenic_map)[3] <- "inf_times"
}
sim_chain_wd <- sim_chain_wd_wrong_map
source("scripts/aux/plot_individual_titre_fits_sim_wrong_map.R")


compare_ar_annual_sim <- bind_rows(compare_ar_annual_base, compare_ar_annual_wrong_map, compare_ar_annual_no_offsets)
compare_ar_annual_sim$Model <- factor(compare_ar_annual_sim$Model,levels=c("Main","No offsets","Bedford map"))
p_sim_ar_annual <-  ggplot(compare_ar_annual_sim) + 
  geom_vline(xintercept=1:48 - 0.5,col="grey70",linewidth=0.1) +
  geom_pointrange(aes(x = time, y=med_y,ymin = lower2, ymax = upper2,col=Model), alpha = 1,position=position_dodge(width=0.75),fatten=0.75) +
  geom_point(aes(x=time,y=y_true,col="True value"),size=0.75) +
  ylab("Estimated annual per\n capita incidence") +
  xlab("Date") +
  scale_x_continuous(breaks=seq(1970,2015,by=5)-1968,labels=seq(1970,2015,by=5),expand=c(0.01,0.01)) + 
  theme_pubr()+
  scale_y_continuous(limits=c(0,1),expand=c(0,0)) +
  scale_color_manual(values=c("Main"=brewer.pal(n = 8, name = "Dark2")[1], 
                              "No offsets"=brewer.pal(n = 8, name = "Dark2")[2],
                              "Bedford map"=brewer.pal(n = 8, name = "Dark2")[3],
                              "True value"="black")) +
  theme(legend.title=element_text(size=8),
        legend.text=element_text(size=8),
        legend.position="bottom",
        legend.margin = margin(-1,-1,-3,-1),
        axis.title=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=8),
        axis.text.y=element_text(size=8),
        plot.tag = element_text(face="bold"),
        strip.text=element_blank(),strip.background = element_blank(),
        plot.margin=margin(r=15,t=5,l=5))

ggsave_jah(p_sim_ar_annual,figure_wd,"p_sim_recovery_ar_annual_compare",7,3)


if(!file.exists("r_data/sim_measurement_smooth_map.RData")){
  chain_wd <- sim_chain_wd_smooth_map
  burnin_use <- 500000
  source("scripts/aux/extract_infection_histories.R")
  save(inf_chain, theta_chain, par_tab, antigenic_map, titre_dat, strain_isolation_times,
       file=paste0(main_wd,"r_data/sim_measurement_smooth_map.RData"))
} else {
  load("r_data/sim_measurement_smooth_map.RData")
  colnames(antigenic_map)[3] <- "inf_times"
}
sim_chain_wd <- sim_chain_wd_smooth_map
source("scripts/aux/plot_individual_titre_fits_sim_smooth_map.R")



if(!file.exists("r_data/sim_measurement_smooth_map_correct.RData")){
  chain_wd <- sim_chain_wd_smooth_map_correct
  burnin_use <- 500000
  source("scripts/aux/extract_infection_histories.R")
  save(inf_chain, theta_chain, par_tab, antigenic_map, titre_dat, strain_isolation_times,
       file=paste0(main_wd,"r_data/sim_measurement_smooth_map_correct.RData"))
} else {
  load("r_data/sim_measurement_smooth_map_correct.RData")
  colnames(antigenic_map)[3] <- "inf_times"
}
sim_chain_wd <- sim_chain_wd_smooth_map_correct
source("scripts/aux/plot_individual_titre_fits_sim_smooth_map_correct.R")
