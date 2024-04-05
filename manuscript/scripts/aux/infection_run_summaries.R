#inf_chain_new <- inf_chain
#inf_chain <- inf_chain_old
library(dplyr)
library(purrr)
library(lubridate)

infection_runs <- identify_run_lengths(inf_chain)

vacc_states_bounds <- read_csv("~/Documents/GitHub/fluscape_serosolver/manuscript/results/vacc_state_bounds.csv") %>%
  rename(Participant_ID = Full.ID)

## Merge with Fluscape data to check what the predictors are
infection_runs_tmp <- infection_runs %>% mutate(is_run = if_else(run_length == 1, 0, 1)) %>% select(i, sampno, start_time, is_run,infection_index) %>%
  mutate(time = 1968*4 + start_time - 1) %>%
  rename(individual = i) %>%
  left_join(fluscape_dat %>% select(individual, DOB, Participant_ID) %>% distinct()) %>%
  mutate(age_at_infection = (time - DOB)/4)
infection_runs_tmp$age_group <- cut(infection_runs_tmp$age_at_infection,breaks=c(0,10,20,30,40,50, 60,100),include.lowest=TRUE)

infection_runs_tmp <-  infection_runs_tmp %>%
  mutate(time = 1968 + (start_time-1)/4) %>%
  mutate(year = lubridate::ymd(floor(time), truncated = 4L)) %>%
  mutate(month = time - floor(time)) %>%
  mutate(date = case_when(month == 0 ~ year,
                          month == 0.25 ~ year + months(3),
                          month == 0.5 ~ year + months(6),
                          month == 0.75 ~ year + months(9))) %>%
  select(-c(month,year))


## Some individuals have multiple vaccine windows
multiple_windows_ids <- vacc_states_bounds %>% group_by(Participant_ID) %>% tally() %>% filter(n > 1) %>% pull(Participant_ID)
single_windows_ids <- vacc_states_bounds %>% group_by(Participant_ID) %>% tally() %>% filter(n == 1) %>% pull(Participant_ID)


infection_runs_tmp_unvacc <- infection_runs_tmp %>% filter(!(Participant_ID %in% unique(vacc_states_bounds$Participant_ID)))
infection_runs_tmp_unvacc$within_1_year_flag <- FALSE

infection_runs_tmp_vacc_single <- infection_runs_tmp %>% 
  filter(Participant_ID %in% single_windows_ids) %>%
  left_join(vacc_states_bounds %>% filter(Participant_ID %in% single_windows_ids) %>% select(individual,upper)) %>%
  mutate(within_1_year_flag = if_else(date>upper & date < (upper + years(1)),TRUE,FALSE))

## Trickier for multiple infection windows
infection_runs_tmp_vacc_multiple <- infection_runs_tmp %>% 
  filter(Participant_ID %in% multiple_windows_ids) %>%
  left_join(
vacc_states_bounds %>% filter(Participant_ID %in% multiple_windows_ids) %>% 
  select(Participant_ID,individual,change,upper) %>%
  group_by(Participant_ID,individual) %>% mutate(change = 1:n()) %>%
  pivot_wider(names_from=change,values_from=upper)
) %>%
  mutate(within_1_year_flag1 = if_else(date>`1` & date < (`1` + years(1)),TRUE,FALSE),
         within_1_year_flag2 = if_else(!is.na(`2`) & date>`2` & date < (`2` + years(1)),TRUE,FALSE),
         within_1_year_flag3 = if_else(!is.na(`3`) & date>`3` & date < (`3` + years(1)),TRUE,FALSE)) %>%
  mutate(within_1_year_flag = within_1_year_flag1 | within_1_year_flag2 | within_1_year_flag3) %>%
  select(-c(within_1_year_flag1,within_1_year_flag2,within_1_year_flag3,`1`,`2`,`3`))

infection_runs_tmp_new <- bind_rows(infection_runs_tmp_unvacc, infection_runs_tmp_vacc_single, infection_runs_tmp_vacc_multiple) %>%
  rename(recently_vaccinated = within_1_year_flag)

## Raw proportions
infection_runs_summary_age_tmp <- infection_runs_tmp_new %>% 
  #filter(date <= "2005-01-01") %>%
  group_by(age_group, sampno) %>% 
  dplyr::summarize(prop=sum(is_run)/n(),N=n()) %>%  
  group_by(age_group) %>% 
  dplyr::summarize(mean_prop=mean(prop),mean_N=mean(N),
                   lower=quantile(prop,0.025),upper=quantile(prop,0.975)) %>% 
  arrange(age_group) %>%
  mutate(Infections="All")
infection_runs_summary_age_tmp1 <- infection_runs_tmp_new %>% 
  filter(date <= "2009-01-01") %>%
  group_by(age_group, sampno) %>% 
  dplyr::summarize(prop=sum(is_run)/n(),N=n()) %>%  
  group_by(age_group) %>% 
  dplyr::summarize(mean_prop=mean(prop),mean_N=mean(N),
                   lower=quantile(prop,0.025),upper=quantile(prop,0.975)) %>% 
  arrange(age_group)%>%
  mutate(Infections="Pre 2009-01-01")

p_raw_age <- ggplot(bind_rows(infection_runs_summary_age_tmp,infection_runs_summary_age_tmp1)) +
  geom_pointrange(aes(x=age_group,y=mean_prop,ymin=lower,ymax=upper,col=Infections),position=position_dodge(width=0.5)) +
  scale_y_continuous(limits=c(0,0.2)) +
  theme(legend.position='bottom') +
  theme_minimal() +
  scale_color_viridis_d() +
  xlab("Age group at time of infection") +
  ylab("Proportion of infections which are runs")
ggsave("~/Documents/GitHub/fluscape_serosolver/manuscript/figures/runs_prop_age.png",height=4,width=7)
ggsave("~/Documents/GitHub/fluscape_serosolver/manuscript/figures/runs_prop_age.pdf",height=4,width=7)


samps <- unique(infection_runs_tmp$sampno)[1:100]
all_ests <- NULL
newdata1 <- expand_grid(age_group=unique(infection_runs_tmp_new$age_group),time=2009,recently_vaccinated=unique(infection_runs_tmp_new$recently_vaccinated))
for(i in seq_along(samps)){
  samp <- samps[i]
  print(i)
  ## What is the probability of it being an infection run, given infection
  fit <- gam(data=infection_runs_tmp_new %>% filter(sampno ==samp), 
             is_run ~ as.factor(age_group) + recently_vaccinated + s(time),family="binomial")
 
  all_ests[[i]] <- cbind(newdata1,prob=predict(fit, newdata=newdata1,type="link",se=TRUE)) %>%
    mutate(predicted_prob=plogis(prob.fit),
           lower=plogis(prob.fit - (1.96 * prob.se.fit)),
           upper=plogis(prob.fit + (1.96 * prob.se.fit))) %>%
    mutate(sampno = samp)
           
}
all_ests <- do.call("bind_rows",all_ests)

p_all_ests <- all_ests %>% group_by(age_group,recently_vaccinated) %>%
  dplyr::summarize(upper1 = quantile(upper,0.975),lower1=quantile(lower,0.025),mean=mean(predicted_prob)) %>%
  mutate(recently_vaccinated = if_else(recently_vaccinated,"Vaccinated within the past year","Not vaccinated within the past year")) %>%
  ggplot() + geom_pointrange(aes(x=age_group,ymin=lower1,ymax=upper1,y=mean)) +
  facet_wrap(~recently_vaccinated,ncol=1) +
  theme_minimal() +
  xlab("Age group") +
  ylab("Probability of infection being a run")

ggsave("~/Documents/GitHub/fluscape_serosolver/manuscript/figures/run_predictors.png",height=6,width=7)
ggsave("~/Documents/GitHub/fluscape_serosolver/manuscript/figures/run_predictors.pdf",height=6,width=7)

infection_runs_tmp_summary <-  infection_runs_tmp_new %>% 
  group_by(individual,infection_index) %>% 
  dplyr::summarize(prob_run = sum(is_run)/n(), 
                   recently_vaccinated = sum(recently_vaccinated)/n(),
                   median_time=median(time),median_age=median(age_at_infection))%>%
  mutate(is_run = prob_run > 0.5, is_vaccinated=recently_vaccinated>0.5)
infection_runs_tmp_summary$age_group <- cut(infection_runs_tmp_summary$median_age,breaks=c(0,10,20,30,40,50, 60,100),include.lowest=TRUE)

fit <- gam(data=infection_runs_tmp_summary %>% filter(median_time > 2009), is_run ~ as.factor(age_group) + is_vaccinated,family="binomial")
summary(fit)

cbind(expand_grid(age_group = unique(infection_runs_tmp_summary$age_group),median_time=2009,
                  is_vaccinated=c(FALSE,TRUE)),prob=predict(fit, 
                            newdata=expand_grid(age_group = unique(infection_runs_tmp_summary$age_group),median_time=2009,
                                                is_vaccinated=c(FALSE,TRUE)),type="link",se=TRUE)) %>%
  mutate(predicted_prob=plogis(prob.fit),
         lower=plogis(prob.fit - (1.96 * prob.se.fit)),
         upper=plogis(prob.fit + (1.96 * prob.se.fit))) %>%
  ggplot() + geom_pointrange(aes(x=age_group,col=is_vaccinated,y=predicted_prob,ymin=lower,ymax=upper)) +
  facet_wrap(~is_vaccinated,scales="free")

infection_runs_tmp %>% filter(sampno ==2) %>% group_by(is_run) %>% dplyr::summarize(mean_age = mean(age_at_infection))
infection_runs_tmp %>% filter(sampno ==2) %>% group_by(is_run) %>% dplyr::summarize(median_index = median(infection_index))


summary_inf_runs <- summarize_run_lengths(infection_runs)


# Overall distribution of runs --------------------------------------------

p_distribution <- infection_runs %>% group_by(sampno, run_length) %>% tally() %>% 
  group_by(run_length) %>% dplyr::summarize(med=median(n),upper=quantile(n,0.975),lower=quantile(n,0.025)) %>%
  ggplot() + geom_pointrange(aes(x=run_length,y=med,ymin=lower,ymax=upper),size=0.25,fatten=0.25) +
  scale_x_continuous(breaks=seq(0,12,by=1)) +
  scale_y_continuous(breaks=seq(0,10000,by=1000)) +
  ylab("Count") + xlab("Length of infection run") + theme_bw()

p_distribution_sub <- infection_runs %>% filter(run_length > 1) %>% group_by(sampno, run_length) %>% tally() %>% 
  group_by(run_length) %>% dplyr::summarize(med=median(n),upper=quantile(n,0.975),lower=quantile(n,0.025)) %>%
  ggplot() + geom_pointrange(aes(x=run_length,y=med,ymin=lower,ymax=upper),size=0.25,fatten=0.25) +
  scale_x_continuous(breaks=seq(0,12,by=1)) +
  #scale_y_continuous(breaks=seq(0,10000,by=1000)) +
  ylab("Count") + xlab("Length of infection run") + theme_bw()

ggsave_jah(p_distribution/p_distribution_sub,"figures/infection_runs/","infection_run_dist",width = 7,height=7)



# Distribution in time ----------------------------------------------------
## When in time do these infections occur?
## Well, they match the attack rates probably...
## So find the proportion of infections which are starts of runs?
summaries <- infection_runs %>% group_by(sampno, start_time) %>% tally() %>% 
  left_join(infection_runs %>% filter(run_length > 1) %>% group_by(sampno, start_time) %>% tally() %>% rename(n_run=n))
summaries <- summaries %>% mutate(n_run = if_else(is.na(n_run),0,n_run))
summaries <- summaries %>% mutate(prop = n_run/n) %>% group_by(start_time) %>% dplyr::summarize(mean_y=mean(prop))


n_alive <- get_n_alive(titre_dat, strain_isolation_times)

total_infs <- inf_chain %>% filter(x == 1) %>% group_by(sampno, j) %>% tally()
total_infs <- total_infs %>% left_join(tibble(j=seq_along(n_alive),n_alive=n_alive))
total_infs <- total_infs %>% mutate(ar = n/n_alive)
total_infs_summary <- total_infs %>% group_by(j) %>% dplyr::summarize(mean_ar=mean(ar),lower=quantile(ar,0.025),upper=quantile(ar,0.975))
p_timing <- ggplot(total_infs_summary) + 
  geom_ribbon(aes(x=j,ymin=lower,ymax=upper),fill="red",alpha=0.25) +
  geom_line(aes(x=j,y=mean_ar,col="Attack rate"))+
  geom_line(data=summaries,aes(x=start_time,y=mean_y,col="Start of infection run")) + 
  scale_y_continuous(limits=c(0,1)) + 
  scale_x_continuous(labels=c(seq(1970,2015,by=5)),breaks=c(seq(1970,2015,by=5)*4 - 1968*4)) +
  geom_hline(yintercept=0.5) +
  ylab("Attack rate/proportion of\n infections which are runs") +
  xlab("Time since start") + theme_bw() +
  scale_color_manual(name="",values=c("Attack rate"="red","Start of infection run"="blue")) +
  theme(legend.position=c(0.8,0.8))

ggsave_jah(p_timing,"figures/infection_runs/","infection_run_timing",width = 8,height=4)

## Not particularly interesting -- it's exactly correlated with the ARs because we have to re-sample infections around time periods of high incidence. More interesting I think is where in individual's lives there are more runs.


# Distribution in age -----------------------------------------------------
infection_runs_age <- infection_runs %>% 
  select(i, sampno, run_length,start_time,infection_index) %>%
  rename(j = start_time) %>% 
  left_join(titre_dat %>% select(individual, DOB) %>% distinct() %>% rename(i=individual))
infection_runs_age$j <- strain_isolation_times[infection_runs_age$j]
infection_runs_age$age <- infection_runs_age$j - infection_runs_age$DOB

## Get distribution of run lengths by age, compared to distribution of all infections
infection_runs_age_summary <- infection_runs_age %>% 
  filter(run_length > 1) %>%
  mutate(age=age/4) %>% mutate(age=floor(age)) %>% group_by(sampno, age) %>% tally() %>% mutate(n_reinfs=sum(n))  %>% 
  mutate(prop=n/n_reinfs) %>%
  group_by(age) %>% 
  dplyr::summarize(median_y=median(prop),lower=quantile(prop,0.025),upper=quantile(prop,0.975)) 

infection_all_age_summary <- infection_runs_age %>% 
  mutate(age=age/4) %>% mutate(age=floor(age)) %>% group_by(sampno, age) %>% tally() %>% mutate(n_reinfs=sum(n))  %>% 
  mutate(prop=n/n_reinfs) %>%
  group_by(age) %>% 
  dplyr::summarize(median_y=median(prop),lower=quantile(prop,0.025),upper=quantile(prop,0.975)) 

## Suggests to me that the infection runs are most common at younger ages to account for age-specific boosting.  
p_age_distribution <- ggplot(infection_runs_age_summary) + 
  geom_ribbon(aes(x=age,ymin=lower,ymax=upper,fill="Runs"),alpha=0.25) + 
  geom_line(aes(x=age,y=median_y,col="Runs")) +
  geom_ribbon(data=infection_all_age_summary,aes(x=age,ymin=lower,ymax=upper,fill="All"),alpha=0.25) + 
  geom_line(data=infection_all_age_summary,aes(x=age,y=median_y,col="All")) +
  ylab("Proportion of infections") + xlab("Age at start of infection (year)") +theme_bw()

## Is this true even for those born after 1968? Yes -- even if restricting to recent years
## Get distribution of run lengths by age, compared to distribution of all infections
infection_runs_age_summary_subset <- infection_runs_age %>% filter(DOB > (1985*4)) %>% 
  filter(run_length > 1) %>%
  mutate(age=age/4) %>% mutate(age=floor(age)) %>% group_by(sampno, age) %>% tally() %>% mutate(n_reinfs=sum(n))  %>% 
  mutate(prop=n/n_reinfs) %>%
  group_by(age) %>% 
  dplyr::summarize(median_y=median(prop),lower=quantile(prop,0.025),upper=quantile(prop,0.975)) 

infection_all_age_summary_subset <- infection_runs_age %>% filter(DOB > (1985*4)) %>% 
  mutate(age=age/4) %>% mutate(age=floor(age)) %>% group_by(sampno, age) %>% tally() %>% mutate(n_reinfs=sum(n))  %>% 
  mutate(prop=n/n_reinfs) %>%
  group_by(age) %>% 
  dplyr::summarize(median_y=median(prop),lower=quantile(prop,0.025),upper=quantile(prop,0.975)) 

## Suggests to me that the infection runs are most common at younger ages to account for age-specific boosting.  
p_age_distribution_recent <- ggplot(infection_runs_age_summary_subset) + 
  geom_ribbon(aes(x=age,ymin=lower,ymax=upper,fill="Runs"),alpha=0.25) + 
  geom_line(aes(x=age,y=median_y,col="Runs")) +
  geom_ribbon(data=infection_all_age_summary_subset,aes(x=age,ymin=lower,ymax=upper,fill="All"),alpha=0.25) + 
  geom_line(data=infection_all_age_summary_subset,aes(x=age,y=median_y,col="All")) +
  ylab("Proportion of infections") + xlab("Age at start of infection (year)") +theme_bw()+
  ggtitle("Restricting to individuals born after 1985")


p_comb <- ggplot(bind_rows(
  infection_runs_age_summary %>% mutate(Subset="All individuals"),
  infection_runs_age_summary_subset %>% mutate(Subset="Born after 1985")
  )
  )+ geom_ribbon(aes(x=age,ymin=lower,ymax=upper,fill=Subset),alpha=0.25) + geom_line(aes(x=age,y=median_y,col=Subset)) +
  ylab("Proportion of infection runs") + xlab("Age at start of infection run (year)")+theme_classic() +
  theme(legend.position=c(0.8,0.8))



## Which n'th infection run are most of the runs? Vast majority are first infections
p_nth_distribution <- summary_inf_runs %>% filter(median_run_length > 1) %>% ggplot() + geom_histogram(aes(x=infection_index),binwidth=1,fill="grey70",col="black") +
  scale_y_continuous(breaks=seq(0,60,by=5)) +xlab("Infection number") + ylab("Count") + theme_classic()
  #ggtitle("Infection number of infection history runs of 1s (using median run length)") 

ggsave_jah((p_comb + labs(tag="A"))/(p_nth_distribution + labs(tag="B")),"figures/infection_runs/","infection_run_summary",width = 8,height=8)


ggsave_jah(p_age_distribution/p_nth_distribution,"figures/infection_runs/","infection_run_age_dist",width = 8,height=8)
ggsave_jah(p_age_distribution_recent,"figures/infection_runs/","infection_run_age_dist_recent",width = 8,height=4)

# Individual plots --------------------------------------------------------
inf_chain <- as.data.table(inf_chain)
n_strain <- max(inf_chain$j)
data.table::setkey(inf_chain, "i", "sampno", "chain_no")
## For each individual, how many infections did they have in each sample in total?
n_inf_chain <- inf_chain[, list(V1 = sum(x)), by = key(inf_chain)]
## Get quantiles on total number of infections per indiv across all samples
indiv_hist <- plyr::ddply(n_inf_chain, .(i), function(x) quantile(x$V1, c(0.025, 0.5, 0.975)))
colnames(indiv_hist) <- c("i", "lower", "median", "upper")

## Plot all individual antibody profiles
unique_indivs <- unique(inf_chain$i)
subvector_length <- 25
# Split the vector into contiguous sub-vectors
unique_indivs_broken <- split(unique_indivs, ceiling(seq_along(unique_indivs)/subvector_length))
pdf(paste0(main_wd, "/figures/inf_run_summaries.pdf"))
for(index in seq_along(unique_indivs_broken)){
  use_ids <- unique_indivs_broken[[index]]
  p <- ggplot(summary_inf_runs %>% filter(i %in% use_ids)) +
    geom_pointrange(aes(x=infection_index,y=median_run_length,ymin=lower95_run_length,ymax=upper95_run_length), col="black",size=0.1,fatten=0.1,shape=21) +
    geom_pointrange(data=indiv_hist %>% filter(i %in% use_ids),aes(y=0, x=median,xmin=lower,xmax=upper), col="darkgreen",size=1,fatten=1,shape=21) +
    scale_y_continuous(expand=c(0,0),limits=c(-1,10),labels=c("Total", seq(2,10,by=2)),breaks=seq(0,10,by=2)) +
    scale_x_continuous(breaks=seq(0,25,by=5)) +
    facet_wrap(~i) +
    theme_bw() +
    xlab("Infection run number") +
    ylab("Infection run length")
  print(p)
}
dev.off()

## Subset of individuals who are predicted to have definite runs of >1 infections
use_i <- summary_inf_runs %>% filter(median_run_length > 1) %>% pull(i) %>% unique()
subvector_length <- 25
# Split the vector into contiguous sub-vectors
unique_indivs_broken <- split(use_i, ceiling(seq_along(use_i)/subvector_length))
pdf(paste0(main_wd, "/figures/inf_run_summaries_high.pdf"))
for(index in seq_along(unique_indivs_broken)){
  use_ids <- unique_indivs_broken[[index]]
  p <- ggplot(summary_inf_runs %>% filter(i %in% use_ids)) +
    geom_pointrange(aes(x=infection_index,y=median_run_length,ymin=lower95_run_length,ymax=upper95_run_length), col="black",size=0.1,fatten=0.1,shape=21) +
    geom_pointrange(data=indiv_hist %>% filter(i %in% use_ids),aes(y=0, x=median,xmin=lower,xmax=upper), col="darkgreen",size=1,fatten=1,shape=21) +
    scale_y_continuous(expand=c(0,0),limits=c(-1,10),labels=c("Total", seq(2,10,by=2)),breaks=seq(0,10,by=2)) +
    scale_x_continuous(breaks=seq(0,25,by=5)) +
    facet_wrap(~i) +
    theme_bw() +
    xlab("Infection run number") +
    ylab("Infection run length")
  print(p)
}
dev.off()

ps <- plot_infection_history_chains_indiv(inf_chain, indivs=use_i, pad_chain=FALSE)

## Create key to make better sample time labels
time_key <- strain_isolation_times/4#titre_dat %>% select(samples) %>% distinct() %>% arrange(samples) %>% mutate(samples=samples/4) %>% pull(samples)
time_key <- convert_years_to_quarters(time_key)
names(time_key) <- strain_isolation_times
## Create key to make better virus labels
virus_key1 <- fluscape_dat %>% select(virus,Virus) %>% distinct()
virus_key <- virus_key1 %>% pull(Virus)
names(virus_key) <- virus_key1 %>% pull(virus)


subvector_length <- 5
# Split the vector into contiguous sub-vectors
unique_indivs_broken <- split(use_i, ceiling(seq_along(use_i)/subvector_length))
pdf(paste0(main_wd, "/figures/titre_fits_high_inf_runs.pdf"))
for(index in seq_along(unique_indivs_broken)){
  use_indivs <- unique_indivs_broken[[index]]
  DOBs <- titre_dat %>% filter(individual %in% use_indivs) %>% select(individual,DOB) %>% distinct()
  
  p_titre_fits <- plot_infection_histories_long_mod(chain=theta_chain,infection_histories=inf_chain,
                                                    titre_dat=titre_dat,individuals=use_indivs,
                                                    antigenic_map=antigenic_map,
                                                    strain_isolation_times = strain_isolation_times,
                                                    mu_indices=rep(1:48,each=4),
                                                    measurement_indices_by_time = rep(1:48,each=4),
                                                    par_tab=par_tab,time_key=time_key,virus_key=virus_key)
  colnames(antigenic_map)[3] <- "inf_times"
  p_inf_hists <- generate_cumulative_inf_plots(inf_chain, 0, use_indivs, ages = DOBs,strain_isolation_times = antigenic_map$inf_times)
  p_cumu_infhist <- p_inf_hists[[1]] + scale_x_continuous(breaks=seq(1970,2015,by=5)*4,labels=seq(1970,2015,by=5)) + 
    scale_fill_manual(name="",values=c(`1`="#D55E00")) + 
    scale_color_manual(name="",values=c(`1`="#D55E00"))+
    theme_pubr()+
    theme(legend.title=element_text(size=7),
          legend.text=element_text(size=7),
          legend.position="none",
          legend.margin = margin(-1,-1,-3,-1),
          axis.title=element_text(size=10),
          strip.text=element_text(size=7),
          axis.text.x=element_text(angle=45,hjust=1,size=7),
          axis.text.y=element_text(size=7),
          plot.tag = element_text(face="bold"),
          plot.margin=margin(r=15,t=5,l=5)) +labs(tag="B")
  
  p_titre_fit_comb <- (p_titre_fits+labs(tag="A") + theme(plot.tag = element_text(face="bold"))) + p_cumu_infhist + plot_layout(ncol=2,widths=c(2.5,1))
  print(p_titre_fit_comb)
}
dev.off() 

## Particularly notable individuals
pdf(paste0(main_wd, "/figures/titre_fits_967.pdf"),height=4,width=8)
DOBs <- titre_dat %>% filter(individual ==967) %>% select(individual,DOB) %>% distinct()

p_titre_fits <- plot_infection_histories_long_mod(chain=theta_chain,infection_histories=inf_chain,
                                                  titre_dat=titre_dat,individuals=967,
                                                  antigenic_map=antigenic_map,
                                                  strain_isolation_times = strain_isolation_times,
                                                  mu_indices=rep(1:48,each=4),
                                                  measurement_indices_by_time = rep(1:48,each=4),
                                                  par_tab=par_tab,time_key=time_key,virus_key=virus_key)
colnames(antigenic_map)[3] <- "inf_times"
p_inf_hists <- generate_cumulative_inf_plots(inf_chain, 0, 967, ages = DOBs,strain_isolation_times = antigenic_map$inf_times)
p_cumu_infhist <- p_inf_hists[[1]] + scale_x_continuous(breaks=seq(1970,2015,by=5)*4,labels=seq(1970,2015,by=5)) + 
  scale_fill_manual(name="",values=c(`1`="#D55E00")) + 
  scale_color_manual(name="",values=c(`1`="#D55E00"))+
  theme_pubr()+
  theme(legend.title=element_text(size=7),
        legend.text=element_text(size=7),
        legend.position="none",
        legend.margin = margin(-1,-1,-3,-1),
        axis.title=element_text(size=10),
        strip.text=element_text(size=7),
        axis.text.x=element_text(angle=45,hjust=1,size=7),
        axis.text.y=element_text(size=7),
        plot.tag = element_text(face="bold"),
        plot.margin=margin(r=15,t=5,l=5)) +labs(tag="B")

p_titre_fit_comb <- (p_titre_fits+labs(tag="A") + theme(plot.tag = element_text(face="bold"))) + p_cumu_infhist + plot_layout(ncol=2,widths=c(2.5,1))
print(p_titre_fit_comb)
dev.off() 

## Massive boosts which subside to very low levels
massive_recent_boosts_all <- c(1,54, 147,512,578,684, 967)
massive_recent_boosts <- c(1,54, 147,512,578)
use_indivs <- massive_recent_boosts

pdf(paste0(main_wd, "/figures/inf_run_titre_fits_recent_boosts.pdf"))
p <- ggplot(summary_inf_runs %>% filter(i %in% use_indivs)) +
    geom_pointrange(aes(x=infection_index,y=median_run_length,ymin=lower95_run_length,ymax=upper95_run_length), col="black",size=0.1,fatten=0.1,shape=21) +
    geom_pointrange(data=indiv_hist %>% filter(i %in% use_indivs),aes(y=0, x=median,xmin=lower,xmax=upper), col="darkgreen",size=1,fatten=1,shape=21) +
    scale_y_continuous(expand=c(0,0),limits=c(-1,10),labels=c("Total", seq(2,10,by=2)),breaks=seq(0,10,by=2)) +
    scale_x_continuous(breaks=seq(0,25,by=5)) +
    facet_wrap(~i) +
    theme_bw() +
    xlab("Infection run number") +
    ylab("Infection run length")
  print(p)
dev.off()


pdf(paste0(main_wd, "/figures/titre_fits_recent_boosts.pdf"))
DOBs <- titre_dat %>% filter(individual %in% use_indivs) %>% select(individual,DOB) %>% distinct()
  
p_titre_fits <- plot_infection_histories_long_mod(chain=theta_chain,infection_histories=inf_chain,
                                                    titre_dat=titre_dat,individuals=use_indivs,
                                                    antigenic_map=antigenic_map,
                                                    strain_isolation_times = strain_isolation_times,
                                                    mu_indices=rep(1:48,each=4),
                                                    measurement_indices_by_time = rep(1:48,each=4),
                                                    par_tab=par_tab,time_key=time_key,virus_key=virus_key)
colnames(antigenic_map)[3] <- "inf_times"
p_inf_hists <- generate_cumulative_inf_plots(inf_chain, 0, use_indivs, ages = DOBs,strain_isolation_times = antigenic_map$inf_times)
p_cumu_infhist <- p_inf_hists[[1]] + scale_x_continuous(breaks=seq(1970,2015,by=5)*4,labels=seq(1970,2015,by=5)) + 
  scale_fill_manual(name="",values=c(`1`="#D55E00")) + 
  scale_color_manual(name="",values=c(`1`="#D55E00"))+
  theme_pubr()+
  theme(legend.title=element_text(size=7),
        legend.text=element_text(size=7),
        legend.position="none",
        legend.margin = margin(-1,-1,-3,-1),
        axis.title=element_text(size=10),
        strip.text=element_text(size=7),
        axis.text.x=element_text(angle=45,hjust=1,size=7),
        axis.text.y=element_text(size=7),
        plot.tag = element_text(face="bold"),
        plot.margin=margin(r=15,t=5,l=5)) +labs(tag="B")

p_titre_fit_comb <- (p_titre_fits+labs(tag="A") + theme(plot.tag = element_text(face="bold"))) + p_cumu_infhist + plot_layout(ncol=2,widths=c(2.5,1))
print(p_titre_fit_comb)
dev.off() 

## Early life boosts
early_life_boosts <- c(590, 610, 654, 863, 891, 1086)
use_indivs <- early_life_boosts

pdf(paste0(main_wd, "/figures/inf_run_titre_fits_early_life.pdf"))
p <- ggplot(summary_inf_runs %>% filter(i %in% use_indivs)) +
  geom_pointrange(aes(x=infection_index,y=median_run_length,ymin=lower95_run_length,ymax=upper95_run_length), col="black",size=0.1,fatten=0.1,shape=21) +
  geom_pointrange(data=indiv_hist %>% filter(i %in% use_indivs),aes(y=0, x=median,xmin=lower,xmax=upper), col="darkgreen",size=1,fatten=1,shape=21) +
  scale_y_continuous(expand=c(0,0),limits=c(-1,10),labels=c("Total", seq(2,10,by=2)),breaks=seq(0,10,by=2)) +
  scale_x_continuous(breaks=seq(0,25,by=5)) +
  facet_wrap(~i) +
  theme_bw() +
  xlab("Infection run number") +
  ylab("Infection run length")
print(p)
dev.off()


pdf(paste0(main_wd, "/figures/titre_fits_early_life_boosts.pdf"))
DOBs <- titre_dat %>% filter(individual %in% use_indivs) %>% select(individual,DOB) %>% distinct()

p_titre_fits <- plot_infection_histories_long_mod(chain=theta_chain,infection_histories=inf_chain,
                                                  titre_dat=titre_dat,individuals=use_indivs,
                                                  antigenic_map=antigenic_map,
                                                  strain_isolation_times = strain_isolation_times,
                                                  mu_indices=rep(1:48,each=4),
                                                  measurement_indices_by_time = rep(1:48,each=4),
                                                  par_tab=par_tab,time_key=time_key,virus_key=virus_key)
colnames(antigenic_map)[3] <- "inf_times"
p_inf_hists <- generate_cumulative_inf_plots(inf_chain, 0, use_indivs, ages = DOBs,strain_isolation_times = antigenic_map$inf_times)
p_cumu_infhist <- p_inf_hists[[1]] + scale_x_continuous(breaks=seq(1970,2015,by=5)*4,labels=seq(1970,2015,by=5)) + 
  scale_fill_manual(name="",values=c(`1`="#D55E00")) + 
  scale_color_manual(name="",values=c(`1`="#D55E00"))+
  theme_pubr()+
  theme(legend.title=element_text(size=7),
        legend.text=element_text(size=7),
        legend.position="none",
        legend.margin = margin(-1,-1,-3,-1),
        axis.title=element_text(size=10),
        strip.text=element_text(size=7),
        axis.text.x=element_text(angle=45,hjust=1,size=7),
        axis.text.y=element_text(size=7),
        plot.tag = element_text(face="bold"),
        plot.margin=margin(r=15,t=5,l=5)) +labs(tag="B")

p_titre_fit_comb <- (p_titre_fits+labs(tag="A") + theme(plot.tag = element_text(face="bold"))) + p_cumu_infhist + plot_layout(ncol=2,widths=c(2.5,1))
print(p_titre_fit_comb)
dev.off() 

## Plot a few individuals for supplement
use_i_births <- c(863)
p_titre_fits_births <- plot_infection_histories_long_mod(chain=theta_chain,infection_histories=inf_chain%>% mutate(j = j - (1968*4) + 1),
                                                  titre_dat=titre_dat,individuals=use_i_births,
                                                  antigenic_map=antigenic_map,
                                                  nsamp = 1000,
                                                  strain_isolation_times = strain_isolation_times,
                                                  mu_indices=rep(1:48,each=4),
                                                  measurement_indices_by_time = rep(1:48,each=4),
                                                  par_tab=par_tab,time_key=time_key,virus_key=virus_key)

use_i_boosts <- c(147)
p_titre_fits_boosts <- plot_infection_histories_long_mod(chain=theta_chain,infection_histories=inf_chain%>% mutate(j = j - (1968*4) + 1),
                                                         titre_dat=titre_dat,individuals=use_i_boosts,
                                                         antigenic_map=antigenic_map,
                                                         nsamp = 1000,
                                                         strain_isolation_times = strain_isolation_times,
                                                         mu_indices=rep(1:48,each=4),
                                                         measurement_indices_by_time = rep(1:48,each=4),
                                                         par_tab=par_tab,time_key=time_key,virus_key=virus_key)


p_main1 <- (p_comb + labs(tag="A") + theme(legend.position=c(0.6,0.8))) + (p_nth_distribution + labs(tag="B"))
p_main2 <- (p_titre_fits_births+ labs(tag="A") + theme(plot.tag = element_text(face="bold"))) + (p_titre_fits_boosts+ labs(tag="B")+ theme(plot.tag = element_text(face="bold"))) + plot_layout(ncol=1)
ggsave_jah(p_main1,"figures/infection_runs/","infection_run_summary1",width = 7,height=3)
ggsave_jah(p_main2,"figures/infection_runs/","infection_run_summary2",width = 7,height=8)

