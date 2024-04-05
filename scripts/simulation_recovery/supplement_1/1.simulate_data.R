######################################################
## SIMULATE DATASET FOR SIM-RECOVERY
## Author: James Hay
## Date: 05 Jan 2024
## Summary: simulates serosurvey data matching the fluscape data with smooth antigenic map, random measurement offsets

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

serosolver_wd <- "~/Documents/GitHub/serosolver/"
devtools::load_all(serosolver_wd)
#library(serosolver)

buckets <- 4

setwd(main_wd)
print(paste0("In directory: ", main_wd))
print(paste0("Saving to: ", save_wd))

set.seed(1)

if(!dir.exists(save_wd)) dir.create(save_wd,recursive = TRUE)


## Simulation parameters
n_indivs <- 1000
n_groups <- 1
n_samps <- 2
repeats <- 2
samp_min <- 2009*buckets
samp_max <- 2015*buckets
year_min <- 1968*buckets
year_max <- 2015*buckets
age_min <- 5*buckets
age_max <- 75*buckets

sampled_viruses <- seq(year_min, year_max, by=2*buckets)
sampling_times <- seq(samp_min, samp_max, by=1)

## Antigenic map -- for simulation use all time periods, for estimation we will merge the historical times
antigenic_map <- read.csv(antigenic_map_file)
antigenic_map <- antigenic_map[antigenic_map$inf_times >= year_min & antigenic_map$inf_times <= year_max,]
strain_isolation_times_simulation <- c(seq(1968*4, max(antigenic_map$inf_times)))
antigenic_map_simulation <- antigenic_map[antigenic_map$inf_times %in% strain_isolation_times_simulation,]


## Set up parameter table
par_tab <- read.csv("inputs/par_tab_cr.csv",stringsAsFactors=FALSE)
par_tab <- par_tab[par_tab$names != "phi",]
par_tab[par_tab$names %in% c("alpha","beta"),c("values")] <- c(1,1)
par_tab[par_tab$names == "wane",c("values")] <- 0.2

par_tab_rhos <- read.csv("inputs/par_tab_offset_estimates.csv",stringsAsFactors = FALSE)
par_tab_rhos$values <- rnorm(nrow(par_tab_rhos),0, 0.5)

par_tab <- bind_rows(par_tab, par_tab_rhos)
annual_part <- length(seq(1968*4, 2005*4, by=4))
measurement_indices <- rep(1:48,each=buckets)
measurement_indices <- measurement_indices[1:length(strain_isolation_times_simulation)]


## Simulate realistic-ish attack rates
sim_inf_pars=c("mean"=0.15/buckets,"sd"=1.5,large_first_year=TRUE,"bigMean"=0.6/buckets)
attack_rates <- simulate_attack_rates(strain_isolation_times_simulation, sim_inf_pars["mean"],sim_inf_pars["sd"],TRUE,sim_inf_pars["bigMean"])
attack_rates[1] <- 0.6
attack_rates[attack_rates > 0.6] <- 0.6
plot(attack_rates)

sim_data <- simulate_data(par_tab=par_tab, group=1, n_indiv=n_indivs, 
                          buckets=1,
                          strain_isolation_times=strain_isolation_times_simulation,
                          measured_strains=sampled_viruses,
                          sampling_times=sampling_times, nsamps=n_samps, 
                          antigenic_map=antigenic_map_simulation, 
                          titre_sensoring=0, ## Randomly censor 0% of measurements
                          age_min=age_min,age_max=age_max,
                          attack_rates=attack_rates, repeats=repeats,
                          measurement_indices = measurement_indices,
                          add_noise=TRUE)


plot_data(sim_data$data,sim_data$infection_histories,strain_isolation_times_simulation,n_indivs=10)

titre_dat <- sim_data$data
titre_dat <- titre_dat %>% left_join(sim_data$ages)
titre_dat <- titre_dat %>% arrange(individual, samples, virus, run)

## Subset to have fewer repeats in older samples, only repeats in second sample
titre_dat <- titre_dat %>% group_by(individual) %>% 
    mutate(samp_no = if_else(samples == max(samples), 2, 1)) %>% 
    filter(samp_no == 2 | (samp_no == 1 & run == 1))

## Create antigenic map for fitting real data
strain_isolation_times <- c(seq(year_min, 2005*4, by=4),8021:year_max)
antigenic_map_real <- antigenic_map[antigenic_map$inf_times %in% strain_isolation_times,]

## Create measurement bias indices, this bit is tricky
measurement_indices_fit <- measurement_indices[match(strain_isolation_times, strain_isolation_times_simulation)]

## Create annual versions of data
titre_dat_annual <- titre_dat
titre_dat_annual$virus <- floor(titre_dat_annual$virus/4)
titre_dat_annual$samples <- floor(titre_dat_annual$samples/4)
titre_dat_annual$DOB <- floor(titre_dat_annual$DOB/4)

## Save titre data
write_csv(titre_dat, file=titre_dat_filename)
write_csv(titre_dat_annual, file=titre_dat_filename_annual)

## Save parameter table
write_csv(par_tab, file=par_tab_filename)
## Save attack rates
write_csv(sim_data$attack_rates, file=attack_rates_filename)
## Save infection histories
write_csv(as.data.frame(sim_data$infection_histories), file=infection_histories_filename)