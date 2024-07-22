load("r_data/measurement_quarter_runs.RData")
colnames(antigenic_map)[3] <- "inf_times"

## Create key to make better sample time labels
time_key <- titre_dat %>% select(samples) %>% distinct() %>% arrange(samples) %>% mutate(samples=samples) %>% pull(samples)
time_key <- convert_years_to_quarters(times=(time_key/4))
names(time_key) <- unique(titre_dat$samples)
## Create key to make better virus labels
virus_key1 <- fluscape_dat %>% select(virus,Virus) %>% distinct() %>% mutate(virus = virus)
virus_key <- virus_key1 %>% pull(Virus)
names(virus_key) <- virus_key1 %>% pull(virus)

use_indivs <- sample(unique(titre_dat$individual),3)
use_indivs <- c(709, 712, 947)#c(604, 674, 774)
DOBs <- titre_dat %>% filter(individual %in% use_indivs) %>% select(individual,DOB) %>% distinct()

DOBs <- DOBs %>% mutate(DOB = floor(DOB/4)*4)
titre_dat <- titre_dat %>% mutate(DOB = floor(DOB/4)*4)

plot_indiv_key <- data.frame(i = use_indivs, i_new = seq_along(use_indivs))


p_titre_fits <- plot_infection_histories_long_mod(chain=theta_chain,inf_chain %>% filter(i %in% use_indivs) %>%
                                                    left_join(plot_indiv_key) %>% select(-i) %>% 
                                                    rename(i = i_new),
                                                  titre_dat=titre_dat %>% filter(individual %in% use_indivs) %>%
                                                    left_join(plot_indiv_key %>% rename(individual = i)) %>% select(-individual) %>% rename(individual = i_new),
                                                  individuals=seq_along(use_indivs),
                                                  antigenic_map=antigenic_map,
                                                  strain_isolation_times = strain_isolation_times,
                                                  mu_indices=rep(1:48,each=4),
                                                  measurement_indices_by_time = rep(1:48,each=4),
                                                  par_tab=par_tab,time_key=time_key,virus_key=virus_key)
p_titre_fits_data <- p_titre_fits[[2]]
p_titre_fits <- p_titre_fits[[1]]
p_titre_fits <- p_titre_fits + labs(tag="B") + theme(plot.tag=element_text(face="bold")) #+ facet_wrap(individual~samples_label,nrow=1)


## Extract data for plot
plot_data1 <- p_titre_fits_data[[1]] %>% select(individual, virus, titre,lower,run, median,upper,samples,samples_label,virus_label) %>% rename(lower_obs=lower,upper_obs=upper,median_obs=median) 
plot_data2 <-p_titre_fits_data[[2]]%>% mutate(infection_time=variable/4)
plot_data3 <- p_titre_fits_data[[3]] %>% select(individual, virus, lower, median,upper,samples) %>% distinct()
plot_data4 <- left_join(plot_data1,plot_data3, by=c("individual","virus","samples")) %>% select(-c(samples,virus)) %>%
  rename(`Lower 95% prediction interval`=lower_obs, `Posterior median observation`=median_obs,
         `Upper 95% prediction interval`=upper_obs, `Lower 95% CrI`=lower, `Upper 95% CrI`=upper, `Posterior median`=median,`Virus`=virus_label,`Sample time`=`samples_label`,`Repeat number`=run)
plot_data2 <- plot_data2 %>% rename(`Sample time`=`samples_label`) %>% select(-samples) %>% rename(`Posterior probability of infection`=value)

write.csv(plot_data4,"~/Documents/GitHub/fluscape_infection_histories//data/figure_data/FigS24Ai.csv",row.names=FALSE)
write.csv(plot_data2,"~/Documents/GitHub/fluscape_infection_histories//data/figure_data/FigS24Aii.csv",row.names=FALSE)

load("r_data/base_quarter_runs.RData")
colnames(antigenic_map)[3] <- "inf_times"

DOBs <- DOBs %>% mutate(DOB = floor(DOB/4)*4)
titre_dat <- titre_dat %>% mutate(DOB = floor(DOB/4)*4)

load("r_data/base_quarter_runs.RData")
colnames(antigenic_map)[3] <- "inf_times"

p_titre_fits1 <- plot_infection_histories_long_mod(chain=theta_chain,inf_chain %>% filter(i %in% use_indivs) %>%
                                                     left_join(plot_indiv_key) %>% select(-i) %>% 
                                                     rename(i = i_new),
                                                   titre_dat=titre_dat %>% filter(individual %in% use_indivs) %>%
                                                     left_join(plot_indiv_key %>% rename(individual = i)) %>% select(-individual) %>% rename(individual = i_new),
                                                   individuals=seq_along(use_indivs),
                                                  antigenic_map=antigenic_map,
                                                  strain_isolation_times = strain_isolation_times,
                                                  mu_indices=NULL,
                                                  measurement_indices_by_time = NULL,
                                                  par_tab=par_tab,time_key=time_key,virus_key=virus_key,
                                                  fill_color = "cornflowerblue")
p_titre_fits1_data <- p_titre_fits1[[2]]
p_titre_fits1 <- p_titre_fits1[[1]]
p_titre_fits1 <- p_titre_fits1 + labs(tag="A") + theme(plot.tag=element_text(face="bold"))#+ facet_wrap(individual~samples_label,nrow=1)


## Extract data for plot
plot_data1 <- p_titre_fits1_data[[1]] %>% select(individual, virus, titre,lower,run, median,upper,samples,samples_label,virus_label) %>% rename(lower_obs=lower,upper_obs=upper,median_obs=median) 
plot_data2 <- p_titre_fits1_data[[2]] %>% mutate(infection_time=variable/4)
plot_data3 <- p_titre_fits1_data[[3]] %>% select(individual, virus, lower, median,upper,samples) %>% distinct()
plot_data4 <- left_join(plot_data1,plot_data3, by=c("individual","virus","samples")) %>% select(-c(samples,virus)) %>%
  rename(`Lower 95% prediction interval`=lower_obs, `Posterior median observation`=median_obs,
         `Upper 95% prediction interval`=upper_obs, `Lower 95% CrI`=lower, `Upper 95% CrI`=upper, `Posterior median`=median,`Virus`=virus_label,`Sample time`=`samples_label`,`Repeat number`=run)
plot_data2 <- plot_data2 %>% rename(`Sample time`=`samples_label`) %>% select(-samples) %>% rename(`Posterior probability of infection`=value)

write.csv(plot_data4,"~/Documents/GitHub/fluscape_infection_histories//data/figure_data/FigS22Bi.csv",row.names=FALSE)
write.csv(plot_data2,"~/Documents/GitHub/fluscape_infection_histories//data/figure_data/FigS22Bii.csv",row.names=FALSE)


p_meas_fit <- p_titre_fits1 + p_titre_fits + plot_layout(nrow=1)
ggsave_jah(p_meas_fit,figure_wd,"titre_fits_meas",width=8,height=6)
