################################################################
## CALCULATE ATTACK RATE STATISTICS FOR PAPER
################################################################
save_name <- if_else(remove_duplicate_infections, "Table2","TableS2")
## --------- Sample sizes over time ------------
## Need to manage who is "alive" in each time period
## Find if individual was available to be infected at each timepoint
alive_effective <- titre_dat %>% 
  select(individual, DOB,samples) %>% 
  distinct() %>% 
  group_by(individual) %>% 
  filter(samples == max(samples)) %>%
  expand_grid(time=strain_isolation_times)%>% 
  group_by(individual) %>% 
  mutate(alive=as.numeric(time >= DOB & time <= samples)) %>% 
  mutate(impute = as.numeric(time > samples)) %>%
  mutate(born = as.numeric(time >= DOB)) %>%
  select(individual,time,alive,born,impute) %>% 
  rename(i=individual,j=time) %>%
  arrange(i,j)

alive_effective %>% ggplot() + geom_tile(aes(x=j,y=i,fill=alive))
alive_effective %>% ggplot() + geom_tile(aes(x=j,y=i,fill=born))
alive_effective %>% ggplot() + geom_tile(aes(x=j,y=i,fill=impute))

## Denominator for reinfections is individuals who were "alive" for all 4 quarters in the year
alive_effective_reinfections <- alive_effective %>% 
  filter(alive==TRUE) %>% 
  mutate(year=floor(j/4)) %>% 
  group_by(i, year) %>% 
  dplyr::summarize(n_contribution = n()) %>%
  mutate(include=if_else(n_contribution==4, TRUE, FALSE))


## --------- Calculate reinfections ------------
inf_chain$year <- floor(1968 + inf_chain$j/4 - 0.25) ## Convert time to years
## Group by years and get total number of infections per year
setkey(inf_chain, "sampno","i","chain_no","year")
all_reinfections <-  inf_chain[,list(total_infs=sum(x)),by=key(inf_chain)]
## Only include those who were alive for all 4 quarters of the year
all_reinfections <- left_join(all_reinfections, alive_effective_reinfections) %>% filter(include==TRUE)

## Find those who were reinfected within a year and get proportion reinfected
all_reinfections_summaries <- all_reinfections %>% 
  mutate(reinfection = as.numeric(total_infs > 1)) %>% 
  group_by(sampno, chain_no, year) %>% 
  dplyr::summarize(N=n(), N_reinfection=sum(reinfection)) %>%
  mutate(prop_reinfected=N_reinfection/N)

## Get posterior median and 95% CrI for reinfections for each year
all_reinfections_stats <- all_reinfections_summaries %>% 
  group_by(year) %>% 
  dplyr::summarize(med=median(prop_reinfected*100),lower=quantile(prop_reinfected*100,0.025),upper=quantile(prop_reinfected*100,0.975))

##
all_reinfections_stats %>%
  mutate(label = paste0(signif(med,3),"% (", signif(lower,3), "%-",signif(upper,3),"%)")) %>%
  filter(year %in% 2010:2014) %>%
  pull(label) %>%
  print()

## Save values used for table
all_reinfections_summaries %>% 
  ungroup() %>% 
  select(-chain_no) %>% 
  rename(`Posterior draw`=sampno,Year=year,N=N,`Number of reinfections`=N_reinfection,`Proportion reinfected`=prop_reinfected) %>%
write.csv(file=paste0("~/Documents/GitHub/fluscape_infection_histories/data/figure_data/",save_name,"ii.csv"),row.names=FALSE)

## What proportion of reinfections occurred since 2008?
print("What proportion of reinfections occurred since 2008?")
all_reinfections %>% mutate(reinf_recent = total_infs > 1 & year >2008, reinf=total_infs > 1) %>%
  group_by(sampno) %>% dplyr::summarize(reinfs_recent=sum(reinf_recent),reinfs=sum(reinf)) %>%
  mutate(prop=reinfs_recent/reinfs) %>% 
  dplyr::summarize(med=median(prop),lower=quantile(prop,0.025),upper=quantile(prop,0.975)) %>%
  print()

## --------- Calculate attack rates with imputed infections ------------
## Impute infection states for censored individuals (after they left study)
## Find individuals who were born and prior to final sample time -- these are the individuals we estimated infection states for
## Get number alive and number of infections for each MCMC sample
n_alive_time <- alive_effective %>% filter(alive == 1) %>% group_by(j) %>% tally()
n_alive_time <- n_alive_time %>% mutate( j = j - min(j) + 1)
n_infections_time <- inf_chain %>% group_by(j,sampno,chain_no) %>% dplyr::summarize(n=sum(x))
n_infections_time <- n_infections_time %>% rename(n_infections=n)

## Merge with infection state posteriors and flag which time periods need to have new states imputed for (individuals who were alive but after their final sample collection)
inf_chain <- inf_chain %>% left_join(n_alive_time) %>% left_join(n_infections_time)
inf_chain <- inf_chain %>% left_join(alive_effective %>% mutate(j = j - min(j) + 1))
inf_chain_impute <- inf_chain %>% filter(impute == 1)

## Draw new infection states for those flagged for imputation
inf_chain_impute$x <- rbinom(nrow(inf_chain_impute),size = 1, 
                             prob=(inf_chain_impute$n_infections + unique(theta_chain$alpha))/
                               (inf_chain_impute$n + unique(theta_chain$alpha) + unique(theta_chain$beta)))

## Move back into inf_chain
inf_chain <- bind_rows(inf_chain %>% filter(impute != 1), inf_chain_impute)

## Find proportion who were alive (born) and number of infections for each draw and time
attack_rates <- inf_chain %>% filter(impute == 0) %>% 
  group_by(j, sampno, chain_no) %>% dplyr::summarize(n_alive=sum(alive),n_infections=sum(x)) %>%
  mutate(ar=n_infections/n_alive)

## Calculate posterior median and 95% CrI for each time period
attack_rates_summaries <- attack_rates %>% 
  mutate(ar=n_infections/n_alive) %>% 
  group_by(j) %>% 
  dplyr::summarize(med=median(ar),low=quantile(ar,0.025),high=quantile(ar,0.975))

ggplot(attack_rates_summaries) + geom_pointrange(aes(x=j,y=med,ymin=low,ymax=high))


## Compare attack rates from 2010 onwards to all times
attack_rates_summaries %>% mutate(year =  (1968*4 + (j-1))/4) %>% 
  mutate(recent=if_else(year >= 2009.75,TRUE,FALSE)) %>% 
  group_by(recent) %>%
  dplyr::summarize(med_ar=median(med),lower=quantile(med,0.025),upper=quantile(med,0.975)) %>%
  print()


## Get overall median attack rate (reported in paper)
print("Median overall attack rate")
attack_rates %>% group_by(sampno) %>% dplyr::summarize(med_ar=median(ar)) %>% 
  ungroup() %>% dplyr::summarize(med=median(med_ar),low=quantile(med_ar,0.025),high=quantile(med_ar,0.975)) %>%print()

print("Lowest AR")
attack_rates_summaries %>% filter(med == min(med)) %>% print()
print("Highest AR")
attack_rates_summaries %>% filter(med == max(med)) %>% print()
print("1985 AR")
attack_rates_summaries %>% mutate(year = (1968*4 + (j-1))/4) %>% filter(year == 1985) %>% print()

## Get number infected at least once per year
setkey(inf_chain, "sampno","i","year","chain_no")
inf_chain_annual <- inf_chain[,list(infected=as.numeric(any(x==1))),by=key(inf_chain)]

## Find number infected at least once per year
setkey(inf_chain_annual, "sampno","year","chain_no")
ars_annual <- inf_chain_annual[,list(total_infs=sum(infected)),by=key(inf_chain_annual)]

## Get maximum alive in each year
n_alive_annual <- alive_effective %>% filter(born==TRUE) %>% mutate(year=floor(j/4)) %>% group_by(i, year) %>% dplyr::summarize(n_contribution = n()) %>% mutate(alive=as.numeric(n_contribution > 0)) %>% group_by(year) %>% tally()

ars_annual <- ars_annual %>% left_join(n_alive_annual) %>% mutate(ar = total_infs/n)
ars_annual %>% group_by(year) %>% dplyr::summarize(med=median(ar),low=quantile(ar,0.025),high=quantile(ar,0.975)) %>% tail()

ars_annual %>% 
  group_by(sampno,chain_no) %>% 
  dplyr::summarize(med=median(ar),low=quantile(ar,0.025),high=quantile(ar,0.975)) %>%
  ungroup() %>% dplyr::summarize(med1=median(med),low=quantile(med,0.025),high=quantile(med,0.975))

ar_estimates_annual_summary <- ars_annual %>% group_by(year) %>% dplyr::summarize(med_ar = median(ar),lower=quantile(ar, 0.025),upper=quantile(ar,0.975))

if(remove_duplicate_infections){
  save(ar_estimates_annual_summary,file=paste0(main_wd,"r_data/fluscape_ar_estimates.RData"))
  save(ars_annual,file=paste0(main_wd,"r_data/fluscape_ar_estimates_draws.RData"))
}

ars_annual  %>% select(year,total_infs,n,ar,sampno) %>% 
  rename(`Proportion infected at least once`=ar,`N alive`=n,`Posterior draw`=sampno,Year=year) %>% arrange(Year, `Posterior draw`) %>% 
  write.csv(file=paste0("~/Documents/GitHub/fluscape_infection_histories/data/figure_data/",save_name,"i.csv"),row.names=FALSE)

ggplot(ar_estimates_annual_summary) + geom_pointrange(aes(x=year,y=med_ar,ymin=lower,ymax=upper))

## Get total infections between 2010-2014 inclusive from the per-quarter chain, excluding imputed infections
print("Total infections between 2010-2014 inclusive from the per-quarter chain, excluding imputed infections")
total_infections <- inf_chain %>% 
  filter(impute == 0) %>%
  filter(year %in% 2010:2014) %>% 
  #filter(impute != 1) %>%
  group_by(sampno,i) %>% 
  dplyr::summarize(N_infections=sum(x))

total_infections_dists <- total_infections %>% group_by(sampno) %>% dplyr::summarize(N=n(), 
                                                           `No infections`=100*sum(N_infections==0)/N,
                                                           `1 infection`=100*sum(N_infections==1)/N,
                                                           `2 infections`=100*sum(N_infections==2)/N,
                                                           `3 infections`=100*sum(N_infections==3)/N,
                                                           `4 infections`=100*sum(N_infections==4)/N,
                                                           `5+ infections`=100*sum(N_infections >= 5)/N)
print(apply(total_infections_dists[,3:ncol(total_infections_dists)],2,function(x) quantile(x,c(0.5,0.025,0.975))))

## Get total infections between 2010-2014 inclusive from the annualized chain, including imputed infections
print("Total infections between 2010-2014 inclusive from the annualized chain, including imputed infections")
total_infections <- inf_chain_annual %>% 
  filter(year %in% 2010:2014) %>% 
  #filter(impute != 1) %>%
  group_by(sampno,i) %>% 
  dplyr::summarize(N_infections=sum(infected))

total_infections_dists <- total_infections %>% group_by(sampno) %>% dplyr::summarize(N=n(), 
                                                                                     `No infections`=100*sum(N_infections==0)/N,
                                                                                     `1 infection`=100*sum(N_infections==1)/N,
                                                                                     `2 infections`=100*sum(N_infections==2)/N,
                                                                                     `3 infections`=100*sum(N_infections==3)/N,
                                                                                     `4 infections`=100*sum(N_infections==4)/N,
                                                                                     `5+ infections`=100*sum(N_infections >= 5)/N)
print(apply(total_infections_dists[,3:ncol(total_infections_dists)],2,function(x) quantile(x,c(0.5,0.025,0.975))))

## Get total infections between 2010-2014 inclusive from the annualized chain, excluding imputed infections
print("Total infections between 2010-2014 inclusive from the annualized chain, excluding imputed infections")
setkey(inf_chain, "sampno","i","chain_no","year")
inf_chain_annual_alt <-  inf_chain[impute==0,list(infected=as.numeric(any(x==1))),by=key(inf_chain)]
setkey(inf_chain_annual_alt, "sampno","i")
total_infections <- inf_chain_annual_alt[year %in% 2010:2014, list(N_infections=sum(infected)),by=key(inf_chain_annual_alt)]

total_infections  %>% rename(`Number of years infected at least once`=N_infections,`Individual`=i,`Posterior draw`=sampno) %>% 
  write.csv(file=paste0("~/Documents/GitHub/fluscape_infection_histories/data/figure_data/",save_name,"iii.csv"),row.names=FALSE)

total_infections_dists <- total_infections %>% group_by(sampno) %>% dplyr::summarize(N=n(), 
                                                                                     `No infections`=100*sum(N_infections==0)/N,
                                                                                     `1 infection`=100*sum(N_infections==1)/N,
                                                                                     `2 infections`=100*sum(N_infections==2)/N,
                                                                                     `3 infections`=100*sum(N_infections==3)/N,
                                                                                     `4 infections`=100*sum(N_infections==4)/N,
                                                                                     `5+ infections`=100*sum(N_infections >= 5)/N)
print(apply(total_infections_dists[,3:ncol(total_infections_dists)],2,function(x) quantile(x,c(0.5,0.025,0.975))))




