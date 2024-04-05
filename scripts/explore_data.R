######################################################
## EXPLORE RAW DATA AND ANTIBODY KINETICS MODEL
## Author: James Hay
## Date: 12 July 2023
## Summary: demonstrates how the antibody kinetics model builds antibody profiles over the life course and creates all individual-level antibody profiles

#library(serosolver)
devtools::load_all("~/Documents/GitHub/serosolver")
library(ggplot2)
library(dplyr)

main_wd <- "~/Documents/GitHub/fluscape_infection_histories//"

## Function to plot example infection history and antibody landscapes over time
plot_model <- function(theta,
                       infection_history,
                       antigenic_map,
                       strain_isolation_times,
                       measured_strains_timeseries=strain_isolation_times,
                       add_noise = TRUE,
                       DOB = NULL,
                       buckets=1){
  infection_history[1:(which(strain_isolation_times == DOB)-1)] <- 0
  sampling_times <- strain_isolation_times[seq(1, length(strain_isolation_times),by=buckets)]
  y <- simulate_individual(theta,infection_history,antigenic_map,
                           sampling_times=sampling_times,
                           strain_isolation_times,
                           strain_isolation_times,
                           measurement_bias=NULL, measurement_indices=NULL,
                           add_noise=add_noise,repeats=1,
                           DOB)
  
  colnames(y) <- c("t","virus","titre","repeat")
  
  inf_times <- strain_isolation_times[which(infection_history==1)]
  infection_history_melted <- tidyr::expand_grid(inf_times=inf_times,t=sampling_times) %>%
    dplyr::filter(inf_times <= t)
  
  
  p1 <- ggplot(as.data.frame(y) %>% filter(virus %in% measured_strains_timeseries)) + 
    geom_rect(xmin=min(strain_isolation_times),xmax=DOB,ymin=0,ymax=8,fill="grey70") +
    geom_line(aes(x=t,y=titre,col=virus,group=virus)) +
    geom_vline(data=data.frame(inf_times=strain_isolation_times[which(infection_history==1)]),aes(xintercept=inf_times),linetype='dashed',col="black") +
    scale_y_continuous(limits=c(0,8),breaks=seq(0,8,by=2),expand=c(0,0)) +
      scale_x_continuous(breaks=seq(1970, 2015, by=10*buckets),expand=c(0,0))+
      xlab("Time of sample") +
    ylab("log titre") +
    theme_bw() +
    scale_color_viridis_c()
  y <- as.data.frame(y)
  y_shifted <- as.data.frame(y)
  y_shifted$t <- y_shifted$t + buckets
  y_shifted <- y_shifted[y_shifted$t <= max(sampling_times),]
  
  colnames(y)[3] <- "titre"
  colnames(y_shifted)[3] <- "previous_titre"
  
  y_comb <- y %>% left_join(y_shifted)
  
  y_boosts <- y_comb %>% filter(t %in% inf_times)
  y_wanes <- y_comb %>% filter(!(t %in% inf_times))
  
  
  
  p2 <- ggplot() + 
    geom_rect(data=data.frame(t=sampling_times,xmin=min(strain_isolation_times),xmax=DOB,ymin=0,ymax=8),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="grey70") +
      geom_rect(data=data.frame(t=sampling_times),aes(xmin=t,xmax=max(strain_isolation_times),ymin=0,ymax=8),fill='grey70') +
    geom_ribbon(data=y_boosts,aes(x=virus,ymax=titre,ymin=0,fill="Current landscape"),col="black") + 
    geom_ribbon(data=y_boosts,aes(x=virus,ymax=previous_titre,ymin=0,fill="Previous landscape"),col="black") + 
    geom_ribbon(data=y_wanes,aes(x=virus,ymax=previous_titre,ymin=0,fill="Previous landscape"),col="black") + 
    geom_ribbon(data=y_wanes,aes(x=virus,ymax=titre,ymin=0,fill="Current landscape"),col="black") + 
    geom_vline(data=infection_history_melted,aes(xintercept=inf_times), linetype='dashed',col="black") +
    facet_wrap(~t)+
    scale_y_continuous(limits=c(0,8),breaks=seq(0,8,by=2),expand=c(0,0)) +
      scale_x_continuous(breaks=seq(1970, 2015, by=10*buckets),expand=c(0,0))+
      theme_bw() +
    theme(axis.text.x=element_text(size=7,angle=45,hjust=1)) +
    scale_fill_manual(name="Landscape",values=c("Previous landscape"="red","Current landscape"="blue")) +
    ylab("log titre") +
    xlab("Strain")
  p2
  return(list(p1, p2))
  
}
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

## Demonstrate how model predicts changes in antibody landscapes over time
par_tab <- read.csv(paste0(main_wd,"inputs/par_tab_base.csv"))

theta <- par_tab$values
names(theta) <- par_tab$names
theta["wane"] <- 0.25
antigenic_map <- read.csv(paste0(main_wd, "/data/antigenic_maps/antigenic_map_fonville_annual_continuous.csv"))
buckets <- 1
strain_isolation_times <- antigenic_map$inf_times
infection_history <- rep(0, length(strain_isolation_times))
infection_history[seq(1,length(strain_isolation_times),by=10*buckets)] <- 1
p_model <- plot_model(theta,infection_history,antigenic_map,
           strain_isolation_times,
           measured_strains_timeseries=strain_isolation_times[which(infection_history == 1)],
           add_noise=FALSE,DOB=1980*buckets,
           buckets=buckets)

jahR::save_plots(p_model[[2]],paste0(main_wd, "/figures/model/"),"example_antibody_profile",10,8)
jahR::save_plots(p_model[[1]],paste0(main_wd, "/figures/model/"),"example_antibody_kinetics",8,4)

## Plot all individual antibody profiles
fluscape_dat <- read.csv(paste0(main_wd,"data/fluscape_data_1_resolution.csv"))


fluscape_dat_latest <- fluscape_dat %>%group_by(individual) %>% filter(samples == max(samples))

## Plot titers against strains that circulated prior to birth as function of time between circulation and birth
fluscape_dat_latest %>% filter(virus <= DOB) %>% mutate(t_since = virus - DOB) %>%
    ggplot() + geom_jitter(aes(x=t_since,y=titre),height=0.25,width=0.25) +
    geom_smooth(aes(x=t_since,y=titre))

## Titres are pretty high to stuff that circulated up to 10 years prior to birth. Many individuals have elevated titers to even older strains.
## By the time we're at 30-40 years prior to birth, titres are rare, though there is a little blip at ~37 years prior to birth.
fluscape_dat_latest %>% filter(virus <= DOB) %>% mutate(t_since = virus - DOB) %>%
    group_by(t_since) %>% #summarize(med_titre=median(titre),upper=quantile(titre,0.75),lower=quantile(titre,0.25)) %>%
    ggplot() + geom_boxplot(aes(x=t_since,y=titre,group=t_since)) 

## Look at stuff that circulated after the sample was taken
## Titers are pretty elevated against titers up to 5 years into the future.
fluscape_dat %>% group_by(individual) %>% filter(samples == min(samples)) %>% filter(virus > samples) %>%
    mutate(t_future = virus - samples) %>%
    ggplot() + geom_boxplot(aes(x=t_future,y=titre,group=t_future))


## Find people that had a particularly large increase in titer across the board, and a particularly large decrease
titre_changes <- fluscape_dat %>% 
    left_join(fluscape_dat %>% group_by(individual) %>% select(individual, samples) %>% distinct() %>% mutate(samp_no = 1:n())) %>%
    filter(virus > 2005) %>%
    group_by(individual,samp_no, virus) %>%
    summarize(titre=mean(titre)) %>%
    pivot_wider(values_from=titre,names_from=samp_no) %>%
    mutate(titre_diff = `2`-`1`) %>%
    group_by(individual) %>%
    summarize(mean_diff = mean(titre_diff)) 
titre_changes %>% filter(mean_diff > 2) %>% pull(individual) -> large_changes
titre_changes %>% filter(mean_diff <= 2) %>% pull(individual) -> small_changes


## Of people who had a boost, elevated titers for the strain 3 years into the future
## Look at stuff that circulated after the sample was taken
fluscape_dat %>% filter(individual %in% large_changes) %>% 
    group_by(individual) %>% filter(samples == max(samples)) %>% filter(virus > samples) %>%
    mutate(t_future = virus - samples) %>%
    ggplot() + geom_boxplot(aes(x=t_future,y=titre,group=t_future))+
    scale_y_continuous(limits=c(0,8))

## Of people who had a fairly small change, still pretty high against future strains
## Look at stuff that circulated after the sample was taken
fluscape_dat %>% filter(individual %in% small_changes) %>% 
    group_by(individual) %>% filter(samples == max(samples)) %>% filter(virus > samples) %>%
    mutate(t_future = virus - samples) %>%
    ggplot() + geom_boxplot(aes(x=t_future,y=titre,group=t_future)) +
    scale_y_continuous(limits=c(0,8))
titre_dat <- fluscape_dat
## Plot all individual antibody profiles
unique_indivs <- unique(titre_dat$individual)
subvector_length <- 16
# Split the vector into contiguous sub-vectors
unique_indivs_broken <- split(unique_indivs, ceiling(seq_along(unique_indivs)/subvector_length))
pdf(paste0(main_wd, "/figures/raw_data/comb.pdf"))
for(index in seq_along(unique_indivs_broken)){
  use_ids <- unique_indivs_broken[[index]]
  p <- plot_real_data(titre_dat, use_ids)[[1]]
  print(p)
  #gsave(filename=paste0("~/Documents/GitHub/fluscape_infection_histories//figures/raw_data/raw_data_",min(use_ids),"-",max(use_ids),".png"),
  #      p, width=10,height=8,units="in",dpi=300)
}
dev.off()

pdf(paste0(main_wd, "/figures/raw_data/comb_change.pdf"))
for(index in seq_along(unique_indivs_broken)){
  use_ids <- unique_indivs_broken[[index]]
  p <- plot_real_data(titre_dat, use_ids)[[2]]
  print(p)
  #gsave(filename=paste0("~/Documents/GitHub/fluscape_infection_histories//figures/raw_data/raw_data_",min(use_ids),"-",max(use_ids),".png"),
  #      p, width=10,height=8,units="in",dpi=300)
}
dev.off()

## Plotting subset of younger-individuals for whom DOB is after 1990
young_indivs <- titre_dat %>% filter(DOB > 1990) %>% select(individual) %>% distinct() %>% pull(individual)
pdf(paste0(main_wd, "/figures/raw_data/young_comb.pdf"))
plot_real_data(titre_dat,young_indivs[1:25])[[1]]
plot_real_data(titre_dat,young_indivs[26:50])[[1]]
plot_real_data(titre_dat,young_indivs[51:75])[[1]]
plot_real_data(titre_dat,young_indivs[76:length(young_indivs)])[[1]]
dev.off()

pdf(paste0(main_wd, "/figures/raw_data/young_comb_change.pdf"))
plot_real_data(titre_dat,young_indivs[1:25])[[2]]
plot_real_data(titre_dat,young_indivs[26:50])[[2]]
plot_real_data(titre_dat,young_indivs[51:75])[[2]]
plot_real_data(titre_dat,young_indivs[76:length(young_indivs)])[[2]]
dev.off()
