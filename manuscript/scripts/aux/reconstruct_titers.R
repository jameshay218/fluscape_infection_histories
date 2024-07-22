#inf_chain <- inf_chain_old

n_samps <- 1000
use_indivs <- c(77, 232, 328, 398, 726, 753, 795, 1000)
min_time <- 2005*buckets
max_time <- max(titre_dat$samples)
use_indivs_all <- unique(titre_dat$individual)
## Generate titre data matrix to use to solve titres against circulating strain at time of possible infection (before infection occurs)
expanded_titre_dat <- expand.grid("individual"=use_indivs_all,"group"=1,"virus"=strain_isolation_times,"run"=1)
expanded_titre_dat <- merge(expanded_titre_dat, unique(titre_dat[,c("individual","DOB","Participant_ID")]))
expanded_titre_dat$samples <- expanded_titre_dat$virus
expanded_titre_dat$titre <- 0
expanded_titre_dat <- expanded_titre_dat[order(expanded_titre_dat$individual, expanded_titre_dat$samples),]

## Comment out to solve all individuals
inf_chain <- inf_chain %>% filter(i %in% use_indivs)
expanded_titre_dat <- expanded_titre_dat %>% filter(individual %in% use_indivs)
use_indivs_all <- use_indivs

## This gives median and 95% CI on model predicted titres, probability of infection in each quarter, 
## MLE titre predictions and MLE infection histories
## Not that strain-specific measurement bias is not used, as this is an observation process
## This is another bottleneck
tmp_preboost <- get_titre_predictions(theta_chain, inf_chain, expanded_titre_dat, 
                                      individuals=use_indivs_all, antigenic_map, par_tab=par_tab,
                                      nsamp=n_samps,
                                      measurement_indices_by_time = NULL,
                                      titre_before_infection = TRUE,titres_for_regression=TRUE)

## Extract titre predictions and get age at sample
summary_titre_pred <- tmp_preboost$summary_titres
summary_titre_pred$age_at_infection <- summary_titre_pred$samples - summary_titre_pred$DOB
#summary_titre_pred <- summary_titre_pred[summary_titre_pred$age_at_infection >= 0,]
summary_titre_pred$age_at_inf_year <- summary_titre_pred$age_at_infection/buckets
summary_titre_pred$age_group_at_inf <- cut(summary_titre_pred$age_at_inf_year, c(0,10,20,30,40,50,60,100),include.lowest = TRUE)
summary_titre_pred$year <- floor(summary_titre_pred$samples/buckets)
summary_titre_pred <- data.table(summary_titre_pred)


## Get titre predictions for each MCMC samp and melt
titre_preds <- tmp_preboost$all_predictions

summary_titre_pred <- cbind(summary_titre_pred, data.table(titre_preds))
summary_titre_pred <- melt(summary_titre_pred, id.vars=c("individual", "group", "virus", "run", "DOB", "Participant_ID", 
                                                         "samples", "titre", "lower", "lower_50", "median", "upper_50", 
                                                         "upper", "max", "age_at_infection", "age_at_inf_year", "age_group_at_inf", 
                                                         "year"))
colnames(summary_titre_pred)[c(ncol(summary_titre_pred)-1, ncol(summary_titre_pred))] <- c("sampno","log_titre")

## Convert to integers and truncate at titre of 8
summary_titre_pred$floor_log_titre <- floor(summary_titre_pred$log_titre)
summary_titre_pred$floor_log_titre[summary_titre_pred$floor_log_titre > 8] <- 8

## Remove titres before birth
summary_titre_pred <- summary_titre_pred[summary_titre_pred$age_at_infection >= 0,]
summary_titre_pred$sampno <- as.numeric(summary_titre_pred$sampno)


## Now extract corresponding infection histories
inf_hists <- tmp_preboost$all_inf_hist
combined_inf_hist <- do.call("rbind", inf_hists)
combined_inf_hist <- data.table(combined_inf_hist)
combined_inf_hist$individual <- rep(use_indivs_all, n_samps)
combined_inf_hist$sampno <- rep(1:n_samps, each=length(use_indivs_all))
combined_inf_hist_melted <- melt(combined_inf_hist, id.vars=c("individual","sampno"))

combined_inf_hist_melted$samples <- as.numeric(combined_inf_hist_melted$variable)
combined_inf_hist_melted$samples <- combined_inf_hist_melted$samples + 1968*buckets - 1
colnames(combined_inf_hist_melted) <- c("individual","sampno","var","infection","samples")


all_combined_data <- merge(summary_titre_pred, combined_inf_hist_melted[,c("individual","infection","samples","sampno")],
                           by=c("individual","samples","sampno"))
all_combined_data <- all_combined_data[all_combined_data$samples <= max_time & all_combined_data$samples >= min_time,]
positive_infections <- all_combined_data[all_combined_data$infection == 1,]

## Find average age at first infection
combined_inf_hist_melted <- merge(combined_inf_hist_melted, 
                                  unique(summary_titre_pred[,c("individual","sampno","DOB")]),
                                  by=c("individual","sampno"))
combined_inf_hist_melted$age_at_inf_year <- (combined_inf_hist_melted$samples - combined_inf_hist_melted$DOB)/buckets
## Can only look at this for individuals who were born since 1968
first_infection <- ddply(combined_inf_hist_melted[combined_inf_hist_melted$DOB >= 1968*buckets,], .(sampno, individual), function(x){
    min(x[x$infection == 1,"age_at_inf_year"])
})
print("Age at first infection (years): ")
print(c(mean(first_infection$V1), quantile(first_infection$V1, c(0.5,0.025,0.975))))



second_infection <- ddply(combined_inf_hist_melted[combined_inf_hist_melted$DOB >= 1968*buckets,], .(sampno, individual), function(x){
    infs <- x[x$infection == 1,"age_at_inf_year"]
    infs[2]
    #min(x[x$infection == 1,"age_at_inf_year"])
})
print("Age at second infection (years): ")
print(c(mean(second_infection$V1), quantile(second_infection$V1, c(0.5,0.025,0.975))))

# Plot titres vs. circulating strains and infection -----------------------
## Get predicted titre ranges
tmp <- as.data.table(all_combined_data[,c("individual","samples","sampno","log_titre","DOB")])
setkey(tmp, individual,samples,DOB)
titre_preds_tmp <- tmp[,list(median=median(log_titre),
                             lower=quantile(log_titre,c(0.025)),
                             upper=quantile(log_titre,c(0.975))),by=key(tmp)]

## Get distribution of times of infection
tmp <- as.data.table(all_combined_data[,c("individual","samples","infection","DOB")])
setkey(tmp, individual,samples,DOB)
infection_predictions <- tmp[,list(prob_inf=sum(infection)/n_samps),by=key(tmp)]

used_i <- use_indivs
x_breaks <- seq(min_time, max_time+buckets, by=buckets)
x_labels <- paste0("Q1-",floor(min_time/buckets):(floor((max_time+buckets)/buckets)))
seq_along_labels <- seq(1,length(x_breaks),by=4)
x_breaks <- x_breaks[seq_along_labels]
x_labels <- x_labels[seq_along_labels]

strain_isolation_times_tmp <- strain_isolation_times[strain_isolation_times <= max_time & strain_isolation_times >= min_time]

## Get masks for each individual
DOBs <- get_DOBs(titre_dat)[,2]
age_mask <- create_age_mask(DOBs,strain_isolation_times_tmp)
strain_mask <- create_strain_mask(titre_dat, strain_isolation_times_tmp)
masks <- cbind(age_mask, strain_mask)
masks <- data.frame(masks)
masks$individual <- 1:nrow(masks)
masks$birth <- masks$age_mask + 1968*buckets-1
masks$sample <- masks$strain_mask + 1968*buckets -1


titre_preds_tmp[titre_preds_tmp$individual %in% used_i,] %>% mutate(DOB = floor(DOB/4)*4) %>%
  mutate(`Year of birth`=DOB/4,`Observation time`=samples/4) %>%
  rename(`Individual`=individual,`Posterior median predicted titre`=median,`Lower 95% CrI`=lower,`Upper 95% CrI`=upper) %>%
  select(-c(samples,DOB)) %>%
  write.csv(file="~/Documents/GitHub/fluscape_infection_histories/data/figure_data/Fig4.csv",row.names = FALSE)

hi_titre_labels <- 5*2^seq(0,8,by=1)
hi_titre_labels <- paste0("1:",hi_titre_labels)
hi_titre_labels[1] <- "<1:10"
hi_titre_labels[9] <- "≥1:1280"

p2 <- ggplot(titre_preds_tmp[titre_preds_tmp$individual %in% used_i,]) +  
  geom_rect(data=masks[masks$individual %in% used_i,],aes(xmin=1968*buckets,xmax=birth,ymin=0,ymax=9),
            fill="grey70")+ 
  geom_rect(data=masks[masks$individual %in% used_i,],aes(xmax=2020*buckets,xmin=sample,ymin=0,ymax=9),
            fill="grey70")+
  geom_vline(data=infection_predictions[infection_predictions$individual %in% used_i,],
             aes(xintercept=samples,col=prob_inf))+
  geom_label(data=titre_preds_tmp[titre_preds_tmp$individual %in% used_i,] %>% select(individual,DOB) %>%
               mutate(DOB=paste0("Year of birth: ", floor(DOB/4)))%>%distinct(),aes(x=1980*4,y=7,label=DOB),size=3)  +
  geom_ribbon(aes(x=samples,ymin=lower,ymax=upper),fill="lightblue")+
  geom_line(aes(x=samples,y=median),col="blue") + 
  theme_classic() +
  scale_y_continuous(expand=c(0,0),breaks=seq(0,8,by=1),labels=hi_titre_labels)+
  scale_color_gradient2(low="white",mid="darkorange",high="red",limits=c(0,1),midpoint=0.5)+
  guides(col=guide_colorbar(title="Predicted probability that infection occurred",
                            title.position = "top",
                            label.position = "bottom",
                            direction="horizontal"))+
  coord_cartesian(ylim=c(0,8))+
  scale_x_continuous(breaks=x_breaks,labels=x_labels,expand=c(0.01,0.01))+
  facet_wrap(~individual,ncol=3) +
  theme(axis.text.x=element_text(angle=90,vjust=0.5, size=8),
        legend.key.width = unit(0.8,"cm"),
        legend.position=c(0.85,0.15),
        strip.text=element_blank(),
        panel.spacing = unit(1, "lines"),
        axis.text.y=element_text(family="sans",size=8,colour="black"),
        axis.title.y=element_text(family="sans",size=12,colour="black"),
        axis.title.x=element_text(family="sans",size=12,colour="black"),
        legend.text=element_text(size=7,family="sans"),
        legend.title=element_text(size=7,family="sans"),
        legend.key=element_rect(color="black",fill="none")) +
  xlab("Circulation time") + ylab("HI titre (pre infection)")
ggsave_jah(p2, figure_wd, "reconstructed_titers",width=8,height=6)


## Remove runs of consecutive infections
removed_double_infs <- all_combined_data %>% 
    arrange(individual, sampno, samples) %>% 
    group_by(sampno, individual) %>% 
    mutate(infection_new = ifelse(!is.na(infection) & infection == 1 & !is.na(lag(infection,1)) & lag(infection, 1) == 1, 
                                  0, infection)) 

removed_double_infs <- removed_double_infs %>% rename(infection_old=infection,infection=infection_new)

## Total number of infections
all_combined_data %>% filter(infection == 1) %>% ungroup() %>% group_by(sampno) %>% dplyr::summarize(n_inf = sum(infection)) %>% ungroup() %>% dplyr::summarize(y=median(n_inf),lower=quantile(n_inf, 0.025),upper=quantile(n_inf,0.975))
removed_double_infs %>% filter(infection == 1) %>% ungroup() %>% group_by(sampno) %>% dplyr::summarize(n_inf = sum(infection)) %>% ungroup() %>% dplyr::summarize(y=median(n_inf),lower=quantile(n_inf, 0.025),upper=quantile(n_inf,0.975))

## Subset to desired time range
removed_double_infs <- removed_double_infs[removed_double_infs$samples <= max_time & removed_double_infs$samples >= min_time,]

## Proportion of states with an infection in them
overall_comparison_nodoubles <- removed_double_infs %>% group_by(sampno, floor_log_titre) %>% 
    dplyr::summarize(n_infected=sum(infection),n_tot=n()) %>%
    mutate(V1=n_infected/n_tot)
tmp <- overall_comparison_nodoubles %>% filter(floor_log_titre == 0) %>% select(sampno, V1) %>% rename(zero_titre=V1)

## Compare infection probs across all times and ages NO REPEATS
overall_comparison_nodoubles1 <- merge(overall_comparison_nodoubles, tmp)
overall_comparison_nodoubles1 <- data.table(overall_comparison_nodoubles1)
overall_comparison_nodoubles1$V1 <- overall_comparison_nodoubles1$V1/overall_comparison_nodoubles1$zero_titre
setkey(overall_comparison_nodoubles1, floor_log_titre)
overall_comparison_nodoubles1 <- overall_comparison_nodoubles1[,list(median=median(V1),
                                                                     lower=quantile(V1, c(0.025)),
                                                                     upper=quantile(V1,c(0.975))),
                                                               by=key(overall_comparison_nodoubles1)]
colnames(overall_comparison_nodoubles1) <- c("log HI titre","median","lower","upper")
overall_comparison_nodoubles1$`Age at time of infection` <- "All"

#################
## AGE
## Compare infection probs across all times by age NO REPEATS
age_comparison_nodoubles <- removed_double_infs %>% 
    group_by(sampno, age_group_at_inf, floor_log_titre) %>% 
    dplyr::summarize(n_infected=sum(infection),n_tot=n()) %>%
    mutate(V1=n_infected/n_tot)
tmp_age <- age_comparison_nodoubles %>% filter(floor_log_titre == 0) %>% 
    select(sampno, age_group_at_inf,V1) %>% rename(zero_titre=V1)

## Compare infection probs across all times and ages NO REPEATS
overall_comparison_nodoubles_age <- merge(age_comparison_nodoubles, tmp_age)
overall_comparison_nodoubles_age <- data.table(overall_comparison_nodoubles_age)
overall_comparison_nodoubles_age$V1 <- overall_comparison_nodoubles_age$V1/overall_comparison_nodoubles_age$zero_titre
setkey(overall_comparison_nodoubles_age, floor_log_titre, age_group_at_inf)
overall_comparison_nodoubles_age <- overall_comparison_nodoubles_age[,list(median=median(V1),
                                                                     lower=quantile(V1, c(0.025)),
                                                                     upper=quantile(V1,c(0.975))),
                                                               by=key(overall_comparison_nodoubles_age)]
colnames(overall_comparison_nodoubles_age) <- c("log HI titre","Age at time of infection","median","lower","upper")


###########
## PLOT NO RUNS
noruns_comparison_all <- rbind(overall_comparison_nodoubles1, overall_comparison_nodoubles_age)

cols <- c("black",colorRampPalette(c("royalblue", "orange", "red"))(length(unique(noruns_comparison_all$`Age at time of infection`))))
noruns_comparison_all$`Age at time of infection` <- factor(noruns_comparison_all$`Age at time of infection`, 
                                                           levels=c("All","[0,10]", "(10,20]", "(20,30]", "(30,40]", "(40,50]", "(50,60]","(60,100]"))

p_raw_noruns <- ggplot(noruns_comparison_all) + 
    geom_hline(yintercept=0.5,linetype="dashed",col="grey")+
    geom_vline(xintercept=3,linetype="dashed",col="grey") +
    geom_hline(yintercept=1,linetype="dashed",col="black")+
    geom_ribbon(aes(x=`log HI titre`,ymin=lower,ymax=upper,fill=`Age at time of infection`),alpha=0.2) +
    geom_line(aes(x=`log HI titre`,y=median,col=`Age at time of infection`)) + 
    theme_classic() +   
    ylab("Relative risk of infection") +
    xlab("log HI titre pre infection") +
    scale_color_manual(values=cols)+
    scale_fill_manual(values=cols)+
    scale_x_continuous(expand=c(0,0), limits=c(0,8.1),breaks=seq(0,8,by=1),labels=c(seq(0,7,by=1),"≥8"))+
    facet_wrap(~`Age at time of infection`) + 
    guides(col = guide_legend(title.position = "top",
                              label.position = "bottom",
                              nrow = 2,byrow=TRUE)) +
    theme(legend.direction = "horizontal",
          legend.position=c(0.85,0.15),
          axis.text.x=element_text(family="sans",size=8),
          axis.text.y=element_text(family="sans",size=8),
          legend.text=element_text(family="sans",size=8),
          legend.title=element_text(family="sans",size=8)
    )+
    coord_cartesian(ylim=c(0,1.5))

ggsave_jah(p_raw_noruns, figure_wd, "age_titre_relationship_noruns",height=6,width=7.2)

noruns_comparison_all %>% filter(`log HI titre` == 3)
