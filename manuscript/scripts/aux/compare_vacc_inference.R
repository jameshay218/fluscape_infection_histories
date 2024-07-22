######################################################
## PURPOSE: uses the tidied vaccine status data and compares to inferred infection histories using serosolver
## Author: James Hay
## Last edited: 22 July 2024

library(tidyverse)
library(patchwork)
summarize <- dplyr::summarize
summarise <- dplyr::summarise

## Function to plot individual-level observed antibody profiles
plot_real_data <- function(titre_dat, indivs){
  unique_viruses <- unique(titre_dat$virus)
  
  titre_dat1 <- titre_dat %>% filter(individual %in% indivs) %>% 
    group_by(individual,samples,virus,DOB) %>% 
    summarize(titre=median(titre))
  plot_dat2 <- titre_dat1 %>% 
    filter(individual %in% indivs) %>%
    group_by(individual) %>% 
    filter(samples == max(samples)) %>%
    ungroup() %>%
    mutate(label=paste0("i=", individual, "; t=",samples))
  
  plot_dat1 <- titre_dat1 %>% 
    filter(individual %in% indivs) %>%
    group_by(individual) %>% 
    filter(samples == min(samples)) %>%
    ungroup() %>%
    mutate(label=paste0("i=", individual, "; t=",samples))
  
  birth_dats <- titre_dat1 %>% ungroup() %>% 
    select(individual, DOB) %>% distinct() %>% 
    mutate(xmin=min(titre_dat$virus),ymin=0,ymax=8) %>%
    mutate(DOB = pmax(1968, DOB))
  
  
  sample_dats <- titre_dat1 %>% ungroup() %>% 
    select(individual, samples) %>% distinct() %>%
    arrange(individual, samples) %>%
    group_by(individual) %>% mutate(samp_no=1:n()) %>%
    mutate(samp_no = if_else(samp_no == 1, "First sample","Second sample"))
  
  titre_differences <- plot_dat1 %>% 
    select(individual, virus, titre, DOB) %>% mutate(samp=1) %>% 
    bind_rows(plot_dat2 %>% select(individual, virus, titre, DOB) %>% mutate(samp=2)) %>% 
    tidyr::pivot_wider(names_from=samp,values_from=titre) %>%
    mutate(change=`2`-`1`)
  
  p1 <- ggplot() + 
    geom_rect(data=birth_dats,aes(xmin=xmin,xmax=DOB,ymin=ymin,ymax=ymax),fill="grey70") +
    geom_ribbon(data=plot_dat1,aes(x=virus,ymax=titre,ymin=0,fill="First sample"),col="black",alpha=0.5) + 
    geom_ribbon(data=plot_dat2,aes(x=virus,ymax=titre,ymin=0,fill="Second sample"),col="black",alpha=0.5) + 
    geom_vline(data=sample_dats,aes(xintercept=samples,col=samp_no),linetype="dashed",linewidth=0.5) +
    facet_wrap(~individual)+
    scale_y_continuous(limits=c(0,8),breaks=seq(0,8,by=2),expand=c(0,0)) +
    scale_x_continuous(limits=c(1967,2016),breaks=unique(titre_dat$virus),expand=c(0,0))+
    scale_fill_manual(name="Sample",values=c("Second sample"="red","First sample"="blue")) +
    scale_color_manual(name="Sample",values=c("Second sample"="red","First sample"="blue")) +
    theme_minimal() +
    theme(axis.text.x=element_text(size=5,angle=45,hjust=1),
          legend.position="bottom",
          strip.text=element_text(size=6),
          panel.grid=element_blank()) +
    ylab("log titre") +
    xlab("Strain")
  titre_differences <- titre_differences %>% mutate(virus = factor(virus,levels=unique_viruses))
  titre_differences_up <- titre_differences %>% filter(change >= 0)
  titre_differences_down <- titre_differences %>% filter(change <= 0)
  birth_dats1 <- birth_dats %>% mutate(DOB = as.factor(DOB))
  p2 <- ggplot() + 
    geom_vline(data=birth_dats1,aes(xintercept=DOB,col="Birth")) +
    geom_bar(data=titre_differences_up,aes(x=virus,y=change),stat="identity",alpha=1,fill="orange") + 
    geom_bar(data=titre_differences_down,aes(x=virus,y=change),stat="identity",alpha=1,fill="darkgreen") +
    geom_vline(data=sample_dats %>% mutate(samples = as.factor(samples)),aes(xintercept=samples,col=samp_no),linetype="dashed",linewidth=0.25) +
    geom_hline(yintercept=0) +
    geom_hline(yintercept=c(-2,2),linetype="dashed") +
    facet_wrap(~individual)+
    scale_y_continuous(limits=c(-8,8),breaks=seq(-8,8,by=4),expand=c(0,0)) +
    scale_x_discrete(limits=c(unique(titre_dat1$virus),"2015"))+
    scale_color_manual(name="Sample",values=c("Second sample"="red","First sample"="blue","Birth"="grey70")) +
    theme_minimal() +
    theme(axis.text.x=element_text(size=5,angle=45,hjust=1),
          legend.position="bottom",
          strip.text=element_text(size=6),
          panel.grid=element_blank()) +
    ylab("Fold change in titre") +
    xlab("Strain")
  return(list(p1,p2))
}
save_wd <- "~/Documents/GitHub/fluscape_serosolver/manuscript/figures/"

## Load in inferred infection states -- use the versions prior to thinning
inf_states <- read_csv("~/Google Drive/My Drive/Fluscape Infection Histories/fluscape_inf_history_posteriors_all.csv")

## Get fluscape IDs of fitted individuals
fitted_ids <- unique(inf_states$Participant_ID)

## Read in fluscape data
load("~/Documents/GitHub/fluscape_infection_histories/manuscript/r_data/fluscape_dat.RData")

## Load in vaccination states
vacc_states_easy <- read_csv("~/Documents/GitHub/fluscape_serosolver/manuscript/results/vaccination_statuses_easy.csv")
vacc_states_cleaned <- read_csv("~/Documents/GitHub/fluscape_serosolver/manuscript/results/vaccination_statuses.csv")

## Sort out individual ID keys
inf_state_key <- inf_states %>% select(individual,Participant_ID) %>% distinct()
vacc_state_key <- vacc_states_cleaned %>% select(Full.ID,id) %>% rename(Participant_ID=Full.ID,individual_vacc=id) %>% distinct()
all_ids_key <- left_join(inf_state_key,vacc_state_key) %>% rename(id=individual_vacc)
vacc_states_easy <- vacc_states_easy %>% left_join(all_ids_key) %>% filter(Participant_ID %in% fitted_ids) %>% 
  select(-c(Participant_ID, id))
vacc_states_cleaned <- vacc_states_cleaned %>% left_join(all_ids_key) %>% filter(Participant_ID %in% fitted_ids) %>% 
  select(-c(Participant_ID, id))
vacc_states <- vacc_states_cleaned

## Get sample times
time_key_vacc <- vacc_states %>% select(individual,PART_SAMPLE_TIME) %>% distinct()
time_key_titres <- fluscape_dat %>% select(individual, samples) %>% distinct() %>%
  mutate(year = lubridate::ymd(floor(samples/4), truncated = 4L)) %>%
  mutate(month = samples/4 - floor(samples/4)) %>%
  mutate(date = case_when(month == 0 ~ year,
                          month == 0.25 ~ year + months(3),
                          month == 0.5 ~ year + months(6),
                          month == 0.75 ~ year + months(9))) 

## Individuals with vaccination states  
use_ids <- unique(vacc_states$individual)

## Filter inferred infection states to those with confirmed vaccination
inf_states <- inf_states %>% filter(individual %in% use_ids)
inf_states_expanded <- expand_grid(individual=unique(inf_states$individual),
                                   time=unique(inf_states$time),
                                   sampno=unique(inf_states$sampno),chain_no=1)%>% 
  left_join(inf_states) %>% 
  mutate(infection_state = if_else(is.na(infection_state),0,infection_state)) %>%
  mutate(year = lubridate::ymd(floor(time), truncated = 4L)) %>%
  mutate(month = time - floor(time)) %>%
  mutate(date = case_when(month == 0 ~ year,
                          month == 0.25 ~ year + months(3),
                          month == 0.5 ~ year + months(6),
                          month == 0.75 ~ year + months(9)))
## Expand all infection states
tmp <- expand_grid(individual=unique(inf_states$individual),
                   time=unique(inf_states$time),
                   sampno=unique(inf_states$sampno),chain_no=1)%>% 
  left_join(inf_states) %>% 
  mutate(infection_state = if_else(is.na(infection_state),0,infection_state)) %>% 
  group_by(individual) %>% 
  arrange(individual,time) %>% 
  group_by(individual, sampno) %>% 
  mutate(cumu_inf_hist=cumsum(infection_state)) %>%
  group_by(individual,time) %>%
  dplyr::summarize(mean_inf = median(cumu_inf_hist),
                   low=quantile(cumu_inf_hist,0.025),
                   high=quantile(cumu_inf_hist,0.975)) %>%
  mutate(year = lubridate::ymd(floor(time), truncated = 4L)) %>%
  mutate(month = time - floor(time)) %>%
  mutate(date = case_when(month == 0 ~ year,
                          month == 0.25 ~ year + months(3),
                          month == 0.5 ~ year + months(6),
                          month == 0.75 ~ year + months(9)))


## Calculate posterior probability of infection in each 3-month window
prob_infection <- expand_grid(individual=unique(inf_states$individual),
              time=unique(inf_states$time),
              sampno=unique(inf_states$sampno),chain_no=1)%>% 
  left_join(inf_states) %>% 
  mutate(infection_state = if_else(is.na(infection_state),0,infection_state)) %>%
  group_by(individual,time) %>% dplyr::summarize(prob_inf = sum(infection_state)/1000) %>%
  ungroup() %>%
  mutate(year = lubridate::ymd(floor(time), truncated = 4L)) %>%
  mutate(month = time - floor(time)) %>%
  mutate(date = case_when(month == 0 ~ year,
                          month == 0.25 ~ year + months(3),
                          month == 0.5 ~ year + months(6),
                          month == 0.75 ~ year + months(9))) 

## Posterior probability of infection in each time period
p_inf_hists <- ggplot(prob_infection %>% filter(individual %in% use_ids) %>% filter(date > "2009-01-01")) + 
  geom_ribbon(aes(x=date,ymax=prob_inf,ymin=0),fill="grey70",col="black") + 
  geom_vline(data=time_key_titres%>% filter(individual %in% use_ids),aes(xintercept=date,linetype = "Serum sample")) +
  facet_wrap(~individual) + 
  ylab("Posterior probability of infection") + xlab("Date") + theme_minimal() +
  scale_linetype_manual(name="",values=c("Serum sample"="dashed")) +
  theme(legend.position=c(0.8,0.05),
        axis.text.x=element_text(angle=45,hjust=1,size=7),
        plot.background = element_rect(fill='white'),
        axis.text.y=element_text(size=7))
ggsave(filename=paste0(save_wd, "/posterior_prob_infection_for_vacc.png"),height=7,width=8)
ggsave(filename=paste0(save_wd, "/posterior_prob_infection_for_vacc.pdf"),height=7,width=8)

## Cumulative infection history
ggplot(tmp %>%  filter(individual %in% use_ids)) + 
  geom_ribbon(aes(x=date,ymin=low,ymax=high),alpha=0.25) + 
  geom_line(aes(x=date,y=mean_inf)) + 
  facet_wrap(~individual)

## Tidy data to allow comparison to raw titre data
id_key <- fluscape_dat %>% select(individual) %>% distinct()
prob_infection_plot <- prob_infection %>% 
  filter(individual %in% use_ids) %>% 
  filter(date > "2005-01-01") 
time_key_titres_plot <- time_key_titres %>% filter(individual %in% use_ids)
time_key_vacc_plot <- time_key_vacc%>% filter(individual %in% use_ids)

## Plot inferred infections vs. vaccine windows
vacc_states_plot <- vacc_states%>% filter(date > "2005-01-01",date<=max(tmp$date)) %>% left_join(id_key) %>% 
  left_join(time_key_titres_plot %>% group_by(individual) %>% filter(samples == max(samples)) %>% rename(sample = date)) %>% filter(date <= sample)

use_ids_updated <- vacc_states_plot %>% group_by(individual) %>% summarize(any_vaccinated = any(value == "Vaccinated")) %>% filter(any_vaccinated) %>% pull(individual)

## Next, see if infections land within a vaccination window
## Convert vaccine status dates into quarters
## Get boundaries on vaccine windows
vacc_states_plot1 <- vacc_states_plot %>%
  arrange(individual, date) %>%
  filter(individual %in% use_ids_updated) %>%
  group_by(individual) %>% 
  mutate(change=0) %>%
  mutate(change = if_else(value != lag(value,1) | is.na(lag(value,1)), change+1, change)) %>%
  mutate(change = cumsum(change)) %>%
  mutate(boundary_lower = value != lag(value,1) & value == "Vaccinated") %>% 
  mutate(boundary_lower = if_else(is.na(boundary_lower),TRUE,boundary_lower)) %>%
  mutate(boundary_upper = value != lead(value,1) & value == "Vaccinated") %>% 
  mutate(boundary_upper = if_else(is.na(boundary_upper),TRUE,boundary_upper)) %>%
  arrange(individual,change,date)

vacc_states_bounds <- vacc_states_plot1 %>% filter(boundary_lower == TRUE | boundary_upper==TRUE)  %>% 
  filter(value =="Vaccinated") %>%
  select(date, individual, boundary_lower,boundary_upper, change) %>% 
  arrange(individual, change) %>%
  mutate(boundary=case_when(boundary_lower == TRUE~"lower",
                            boundary_upper==TRUE~"upper"))%>%
  select(-c(boundary_lower,boundary_upper)) %>%
  pivot_wider(values_from=date,names_from=boundary,id_cols=c(individual,change,individual))
  
vacc_states_bounds <- vacc_states_bounds %>% ungroup() %>% mutate(lower = floor_date(lower,"quarter"),upper=ceiling_date(upper,"quarter"))

write_csv(vacc_states_bounds,"~/Documents/GitHub/fluscape_serosolver/manuscript/results/vacc_state_bounds.csv")


## Merge with infection states
inf_states_expanded_summary <- inf_states_expanded %>% 
  filter(individual %in% unique(vacc_states_bounds$individual)) %>% 
  select(individual, infection_state, date, sampno,chain_no) %>%
  filter(date >= "2005-01-01") %>% 
  left_join(vacc_states_bounds) %>% 
  filter(date >= lower - months(3), date <= upper + months(3)) %>%
  group_by(individual, change, sampno) %>% 
  summarize(y=any(infection_state == 1)) %>%
  group_by(individual, change) %>% 
  summarize(prop=sum(y)/n(),N=n())

inf_states_expanded_summary %>% mutate(confident = prop > 0.25) %>% group_by(confident) %>% tally()

ids_has_alignment <- inf_states_expanded_summary %>% filter(prop > 0.25) %>% arrange(individual) %>% pull(individual)
ids_has_misalignment<- inf_states_expanded_summary %>% filter(prop < 0.25) %>%  arrange(individual) %>% pull(individual)



unique_indivs <- ids_has_alignment
subvector_length <- 5
# Split the vector into contiguous sub-vectors
unique_indivs_broken <- split(unique_indivs, ceiling(seq_along(unique_indivs)/subvector_length))

## Plot inferred infection histories against vaccine windows next to each individual's antibody data
fluscape_dat_plot <- read.csv("~/Documents/GitHub/fluscape_infection_histories/data/fluscape_data_1_resolution.csv")

pdf("~/Documents/GitHub/fluscape_serosolver/manuscript/figures/vaccination_status_compare_aligned.pdf")
for(index in seq_along(unique_indivs_broken)){
  use_ids1 <- unique_indivs_broken[[index]]
  p <- plot_real_data(fluscape_dat_plot, use_ids1)[[1]] + facet_wrap(~individual,ncol=1) + theme(legend.title=element_text(size=6),legend.text=element_text(size=6)) + ylab("log HI titre")
  p_rhs <- ggplot(vacc_states_plot %>% filter(individual %in% use_ids1) %>%
                    mutate(value = if_else(value=="Vaccinated","Vaccination reported","No vaccination reported"))) + 
    geom_tile(aes(x=as.Date(date),fill=value,y=0.5),height=1) +
    geom_ribbon(data=prob_infection_plot%>% filter(individual %in% use_ids1),aes(x=date,ymax=prob_inf,ymin=0),fill="grey70",alpha=1,col="black") +
    geom_vline(data=time_key_titres_plot%>% filter(individual %in% use_ids1),aes(xintercept=as.Date(date),linetype="Serum sample")) +
    
    facet_wrap(~individual,ncol=1) + theme_minimal() +
    scale_linetype_manual(name="Data update",values=c("Serum sample"="dashed","Vaccination status update"="dotted")) +
    scale_fill_manual(name="Vaccination status",values=c("No vaccination reported"="purple","Vaccination reported"="orange")) +
    theme(axis.text.x=element_text(size=5,angle=45,hjust=1),
          legend.position="bottom",
          strip.text=element_text(size=6),
          panel.grid=element_blank(),
          legend.title=element_text(size=6),legend.text=element_text(size=6),
          legend.direction="vertical") +
    ylab("Posterior probability of infection") +
    xlab("Date")
  print(p | p_rhs)
  ggsave(filename=paste0(save_wd, "/vacc_compare/vacc_compare_set_aligned_",index,".png"),p | p_rhs,height=8,width=8)
}
dev.off()

unique_indivs <- ids_has_misalignment
subvector_length <- 5
# Split the vector into contiguous sub-vectors
unique_indivs_broken <- split(unique_indivs, ceiling(seq_along(unique_indivs)/subvector_length))
## Plot inferred infection histories against vaccine windows next to each individual's antibody data
pdf("~/Documents/GitHub/fluscape_serosolver/manuscript/figures/vaccination_status_compare_misaligned.pdf")
for(index in seq_along(unique_indivs_broken)){
  use_ids1 <- unique_indivs_broken[[index]]
  p <- plot_real_data(fluscape_dat_plot, use_ids1)[[1]] + facet_wrap(~individual,ncol=1) + theme(legend.title=element_text(size=6),legend.text=element_text(size=6)) + ylab("log HI titre")
  p_rhs <- ggplot(vacc_states_plot %>% filter(individual %in% use_ids1) %>%
                    mutate(value = if_else(value=="Vaccinated","Vaccination reported","No vaccination reported"))) + 
    geom_tile(aes(x=as.Date(date),fill=value,y=0.5),height=1) +
    geom_ribbon(data=prob_infection_plot%>% filter(individual %in% use_ids1),aes(x=date,ymax=prob_inf,ymin=0),fill="grey70",alpha=1,col="black") +
    geom_vline(data=time_key_titres_plot%>% filter(individual %in% use_ids1),aes(xintercept=as.Date(date),linetype="Serum sample")) +
    
    facet_wrap(~individual,ncol=1) + theme_minimal() +
    scale_linetype_manual(name="Data update",values=c("Serum sample"="dashed","Vaccination status update"="dotted")) +
    scale_fill_manual(name="Vaccination status",values=c("No vaccination reported"="purple","Vaccination reported"="orange")) +
    theme(axis.text.x=element_text(size=5,angle=45,hjust=1),
          legend.position="bottom",
          strip.text=element_text(size=6),
          panel.grid=element_blank(),
          legend.title=element_text(size=6),legend.text=element_text(size=6),
          legend.direction="vertical") +
    ylab("Posterior probability of infection") +
    xlab("Date")
  print(p | p_rhs)
  ggsave(filename=paste0(save_wd, "/vacc_compare/vacc_compare_set_misaligned_",index,".png"),p | p_rhs,height=8,width=8)
}
dev.off()

vacc_data_to_save1 <- vacc_states_plot %>% mutate(value = if_else(value=="Vaccinated","Vaccination reported","No vaccination reported")) %>% select(individual, date, value) %>% rename("Vaccination status"=value) %>% filter(individual %in% c(ids_has_alignment,ids_has_misalignment))
vacc_data_to_save2 <- prob_infection_plot %>% select(individual,date, prob_inf) %>% rename(`Posterior probability of infection`=prob_inf)
vacc_data_to_save3 <- time_key_titres_plot %>% select(individual,date) %>% rename(`Date of serum sample collection`=date)

save(vacc_data_to_save1,vacc_data_to_save2,vacc_data_to_save3,file="~/Documents/GitHub/fluscape_infection_histories/data/figure_data/vacc_figure.RData")

## Get random windows matching the size of these windows -- do the same analyses and see what proportion of windows contain infections
## i.e., we expect much less than 77%

## Choose 90 random individuals
## Choose random time points after birth
break
inf_states_all <- read_csv("~/Google Drive/My Drive/Fluscape Infection Histories/fluscape_inf_history_posteriors_all.csv") %>%
  mutate(year = lubridate::ymd(floor(time), truncated = 4L)) %>%
  mutate(month = time - floor(time)) %>%
  mutate(date = case_when(month == 0 ~ year,
                          month == 0.25 ~ year + months(3),
                          month == 0.5 ~ year + months(6),
                          month == 0.75 ~ year + months(9)))

sensitivities_random <- numeric(100)
for(i in 1:100){
  print(i)
  random_windows <- fluscape_dat %>% select(Participant_ID, individual,DOB,samples) %>%
    distinct() %>%
    mutate(year = lubridate::ymd(floor(samples/4), truncated = 4L)) %>%
    mutate(month = (samples/4 - floor(samples/4))) %>%
    mutate(samples = case_when(month == 0 ~ year,
                                month == 0.25 ~ year + months(3),
                                month == 0.5 ~ year + months(6),
                                month == 0.75 ~ year + months(9))) %>%
    group_by(Participant_ID, individual,DOB) %>% 
    filter(samples == max(samples)) %>%
    ungroup() %>%
    mutate(year = lubridate::ymd(floor(DOB/4), truncated = 4L)) %>%
    mutate(month = DOB/4 - floor(DOB/4)) %>%
    mutate(DOB = case_when(month == 0 ~ year,
                            month == 0.25 ~ year + months(3),
                            month == 0.5 ~ year + months(6),
                            month == 0.75 ~ year + months(9))) %>%
    #mutate(DOB = pmax(DOB,"2009-01-01")) %>%
    select(-c(year,month)) %>%
    sample_n(90) %>% 
    group_by(Participant_ID, individual) %>%
    mutate(lower = list(sample(seq(DOB,samples,by="3 months"),1))) %>%
    unnest(lower) %>%
    bind_cols(vacc_states_bounds %>% mutate(width = upper - lower) %>% select(width)) %>%
    mutate(upper = lower + width) %>%
    mutate(upper = pmin(upper, samples)) %>%
    select(-width) %>%
    arrange(individual)
  
  random_inf_states_expanded_summary <- inf_states_all %>% 
    filter(individual %in% unique(random_windows$individual)) %>%
    select(individual, infection_state, date, sampno,chain_no) %>%
    left_join(random_windows %>% ungroup() %>% select(individual, lower, upper)) %>% 
    filter(date >= lower - months(3), date <= upper + months(3)) %>%
    group_by(individual, sampno) %>% 
    summarize(y=any(infection_state == 1)) %>%
    group_by(individual) %>% 
    summarize(prop=sum(y)/1000,N=1000)
  
  sensitivities_random[i] <- random_inf_states_expanded_summary %>% 
    full_join(expand_grid(individual=unique(random_windows$individual))) %>%
    mutate(prop=if_else(is.na(prop),0,prop)) %>%
    mutate(confident = prop > 0.25) %>% group_by(confident) %>% tally() %>%
    pivot_wider(values_from=n,names_from=confident) %>%
    mutate(sens = `TRUE`/(`TRUE` + `FALSE`)) %>%
    pull(sens)
}
mean(sensitivities_random)
quantile(sensitivities_random,c(0.025,0.975))
