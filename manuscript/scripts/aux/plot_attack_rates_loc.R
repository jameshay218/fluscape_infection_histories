## Note -- this needs some functions from the fluscape GeneralUtility.R file

run_name <- "annual"
buckets <- 1
cutoff_time <- 1968*buckets
strain_isolation_times_quarter <- strain_isolation_times
strain_isolation_times <- unique(floor(strain_isolation_times/(4/buckets)))
## Bucket into years rather than per-quarter for spatial analyses

times <- seq_along(strain_isolation_times)
times <- times[which(strain_isolation_times >= cutoff_time)]
strain_isolation_times_tmp <- strain_isolation_times[strain_isolation_times >= cutoff_time]

###############################################################
######### Using most up-do-date latitude/longitude by location. Prior to V1 it's all V1 data, then if it is updated in a subsequent visit,
######### uses the new information.
###############################################################
## Pull out information on locations
locs <- fluscape_dat[,c("LOC_ID", "LOC_Lat", "LOC_Long", "URBAN", "DIST_FRM_GZ", "dens.1", "dens.9")] %>% distinct()
locs <- locs %>% dplyr::arrange(DIST_FRM_GZ) %>% dplyr::mutate(rank.DIST_FRM_GZ=1:n())

loc_dat_complete <- load.and.merge.locs.V1.V2.V3(topdir=fluscape_wd,trim=FALSE)
loc_dat_complete <- apply(loc_dat_complete, 2, as.numeric) %>% as_tibble()
loc_dat_complete <- loc_dat_complete %>% pivot_longer(-c(LOC_ID))
loc_dat_complete <- loc_dat_complete %>% 
    mutate(Visit = if_else(grepl("\\.V[23]$", name), sub(".*\\.([V][23])$", "\\1", name), "V1"),
                                                name = sub("\\.V[23]$", "", name))


## NOTE - updated to be per year
## Find first time period of each location ID
earliest_sample_per_visit_per_loc <- fluscape_dat %>% 
    select(LOC_ID, raw_visit, samples) %>% 
    distinct() %>% 
    mutate(samples = floor(samples*(buckets/4))) %>%
    group_by(LOC_ID, raw_visit) %>% 
    filter(samples==min(samples)) %>% ## Find first sample per location per visit
    arrange(LOC_ID, raw_visit) %>% rename(Visit=raw_visit)

loc_dat_complete <- loc_dat_complete %>% pivot_wider(id_cols=c(LOC_ID, Visit), names_from="name",values_from="value")  %>% left_join(earliest_sample_per_visit_per_loc)
loc_dat_complete <- loc_dat_complete %>% fill(URBAN, .direction="down")
#loc_dat_complete <- loc_dat_complete %>% fill(LOC_Label, .direction="down")

all_samples <- expand_grid(LOC_ID = unique(loc_dat_complete$LOC_ID), samples=strain_isolation_times)
loc_dat_complete <- all_samples %>% full_join(loc_dat_complete)
loc_dat_complete <- loc_dat_complete %>% select(-Visit) %>% fill(c(LOC_Lat, LOC_Long, URBAN, DIST_FRM_GZ, dens.1, dens.9),
                                                       .direction="updown")    
loc_dat_complete <- loc_dat_complete %>% left_join(locs %>% select(LOC_ID, rank.DIST_FRM_GZ))

loc_dat_complete <- loc_dat_complete %>% filter(!is.na(samples)) %>% 
    mutate(samples = samples - min(samples) + 1) %>% rename(time=samples)
###############################################################

## Pull out individual-level data we need
needed_fluscape_dat <- unique(fluscape_dat[,c("individual","Participant_ID","URBAN","DIST_FRM_GZ","LOC_ID",
                                              "dens.1","dens.9","DOB","age_group","age","birth_cohort")])
needed_fluscape_dat <- data.table(needed_fluscape_dat)

## Get masks for each individual
DOBs <- floor(get_DOBs(titre_dat)[,2]/(4/buckets))
age_mask <- create_age_mask(DOBs,strain_isolation_times_tmp)
strain_mask <- create_strain_mask(titre_dat%>%mutate(samples=floor(samples/(4/buckets))),  strain_isolation_times_tmp)
masks <- cbind(age_mask, strain_mask)
masks <- data.frame(masks)
masks$individual <- 1:nrow(masks)
masks <- merge(masks,needed_fluscape_dat[,c("individual","LOC_ID")],by="individual")

## Calculate infection histories and attack rates from the chain
infection_histories <- inf_chain
rm(inf_chain)
colnames(infection_histories) <- c("individual","time","x","sampno","chain_no")
## Convert time incidences to annual if needed, leave as quarterly otherwise. Convert to correct scale
infection_histories <- infection_histories %>% mutate(time = floor((time-1)*(buckets/4)) + 1)
setkey(infection_histories, "sampno","individual","time","chain_no")
infection_histories <- infection_histories[,list(x=as.numeric(any(x==1))),by=key(infection_histories)]

#infection_histories <- infection_histories[infection_histories$time >= min(strain_isolation_times),]
#infection_histories$time <- infection_histories$time - min(strain_isolation_times) + 1
#times <- times - min(times) + 1


infection_histories$individual <- as.numeric(infection_histories$individual)
masks$individual <- as.numeric(masks$individual)
masks <- data.table(masks)
infection_histories <- merge(infection_histories, masks[,c("individual","LOC_ID")],by="individual")


# AR stats by location ----------------------------------------------------
## Look at stats by location
n_alive_loc <- ddply(masks, ~LOC_ID, function(x) {
  sapply(times, function(y){
    nrow(x[x$age_mask <= y & x$strain_mask >= y,])
  })
})
n_alive_loc <- melt(n_alive_loc,id.vars = "LOC_ID")
colnames(n_alive_loc) <- c("LOC_ID","time","n")
n_alive_loc$time <- as.integer(n_alive_loc$time)

## Overall AR
n_alive <- get_n_alive(titre_dat %>% mutate(DOB=floor(DOB/(4/buckets)),samples=floor(samples/(4/buckets))),strain_isolation_times)
n_alive <- data.frame("time"=seq_along(strain_isolation_times),"n"=n_alive)

## Stratifying by location ID
data.table::setkey(infection_histories, "time", "sampno")

## Number of samples with a 1 divided by total samples to get AR
total_infs <- infection_histories[, list(V1 = sum(x)), by = key(infection_histories)]
total_infs <- merge(total_infs, n_alive, by=c("time"))
total_infs$ar <- total_infs$V1/total_infs$n
total_infs[is.nan(total_infs$ar),"ar"] <- NA
total_infs <- total_infs[,c("time","sampno","ar")]
data.table::setkey(total_infs, "time")
total_infs$sampno <- as.numeric(total_infs$sampno)
total_infs <- total_infs[complete.cases(total_infs),]
ar_estimates <- total_infs[,list(lower_quantile=quantile(ar,c(0.025)),
                                 median=quantile(ar,c(0.5)),
                                 upper_quantile=quantile(ar,c(0.975)),
                                 precision=1/var(ar),
                                 var=sd(ar)/mean(ar)
),by=key(total_infs)]


## Stratifying by location ID
data.table::setkey(infection_histories, "time", "LOC_ID","sampno")
## Number of samples with a 1 divided by total samples to get AR
total_infs <- infection_histories[, list(V1 = sum(x)), by = key(infection_histories)]
total_infs <- merge(total_infs, n_alive_loc, by=c("time","LOC_ID"))
total_infs$ar <- total_infs$V1/total_infs$n
total_infs[is.nan(total_infs$ar),"ar"] <- NA
total_infs <- total_infs[,c("time","LOC_ID","sampno","ar","n","V1")]
total_infs <- total_infs %>% filter(n > 0)
data.table::setkey(total_infs, "time", "LOC_ID","n")
total_infs$sampno <- as.numeric(total_infs$sampno)
total_infs <- total_infs[complete.cases(total_infs),]
ar_estimates_by_loc <- total_infs[,list(lower_quantile=quantile(ar,c(0.025)),
                                 median=quantile(ar,c(0.5)),
                                 upper_quantile=quantile(ar,c(0.975)),
                                 precision=1/var(ar),
                                 var=sd(ar)/mean(ar)
                                 ),by=key(total_infs)]
## Get coefficient of variation of attack rates for each time point, within each posterior draw
total_infs <- total_infs %>% filter(is.finite(ar))
setkey(total_infs, "time","sampno")
coef_var <- total_infs[,list(mean=mean(ar),sd=sd(ar),coef1=sd(ar)/mean(ar)),by=key(total_infs)]
coef_var[is.nan(coef_var$coef1),"coef1"] <- NA
coef_var_overall <- ddply(coef_var, ~sampno, function(x) quantile(x$coef1,c(0.5,0.025,0.975),na.rm=TRUE))
coef_var_overall_res <- quantile(coef_var_overall$`50%`,c(0.5,0.025,0.975))

coef_var <- melt(coef_var, id.vars=c("time","sampno"))
coef_var_by_time_res <- ddply(coef_var, .(variable,time), function(x) quantile(x$value, c(0.5,0.025,0.975),na.rm=TRUE))

var_key <- c("coef1"="Coefficient of variation","sd"="Standard deviation of location-specific attack rates","mean"="Mean of location-specific attack rates")
coef_var_by_time_res$variable <- var_key[as.character(coef_var_by_time_res$variable)]

coef_var_by_time_res %>% select(-c(`2.5%`,`97.5%`)) %>% pivot_wider(names_from= variable,values_from=`50%`) %>% select(c(`Mean of location-specific attack rates`,`Coefficient of variation`)) %>% as.matrix() %>% cor(method="pearson")
coef_var_by_time_res %>% select(-c(`2.5%`,`97.5%`)) %>% pivot_wider(names_from= variable,values_from=`50%`) %>% ggplot() + geom_point(aes(x=`Mean of location-specific attack rates`,y=`Coefficient of variation`))
## Create a simulation for how much variation would be expected if each location's attack rate was drawn from
## the same distribution
## So we have attack rates for 40 locations, each of population size 25, and we want the proportion infected
## in each location
n <- 10000
coef_sims <- numeric(n)
sds <- numeric(n)
for(i in 1:n){
  ar <- runif(1)
  ars <- rbinom(n=40,size=25,prob=ar)/25
  coef_sims[i] <- sd(ars)/mean(ars)
  sds[i] <- sd(ars)
}
coef_sims_null <- melt(data.frame(c(mean(coef_sims,na.rm=TRUE), quantile(coef_sims, c(0.025,0.975),na.rm=TRUE))))
coef_sims_null$variable <- "Coefficient of variation"
overall_var_dat <- data.frame(variable="Coefficient of variation",y=coef_var_overall_res[1])

x_breaks <- seq(cutoff_time, 2016*buckets,by=4*buckets) - cutoff_time + 1
x_labels <- paste0("Q1-",seq(floor(cutoff_time/buckets),2016, by=4)) 

ar_var_p <- ggplot(coef_var_by_time_res) + 
  geom_ribbon(aes(x=time,ymin=`2.5%`,ymax=`97.5%`,fill=variable),alpha=0.2) + 
  geom_line(aes(x=time,y=`50%`,col=variable))  + 
  geom_hline(data=coef_sims_null[c(1,3),],aes(yintercept=value),linetype="dashed",col="grey40") +
  scale_x_continuous(expand=c(0.01,0.01),breaks=x_breaks,labels=x_labels) + 
  theme_classic()+
  geom_hline(data=overall_var_dat,aes(yintercept=y,col=variable),linewidth=1)+
  ylab("Posterior attack rate estimate") +
  xlab("Circulation time")+
  theme(legend.position="none",
        axis.text.x=element_text(family="sans",size=7,colour="black",angle=90,vjust=0.5),
        axis.title.y=element_text(family="sans",size=12,colour="black"),
        axis.text.y=element_text(family="sans",size=7,colour="black"),
        strip.background = element_blank())+
  facet_wrap(~variable,ncol=1,scales="free_y") +
    scale_fill_manual(values=natparks.pals("SouthDowns",3)) +
    scale_color_manual(values=natparks.pals("SouthDowns",3)) 
ar_var_p

## Get publication plot of CoV, mean etc

ar_var_pA <- ggplot(coef_var_by_time_res %>% filter(variable == "Coefficient of variation")) + 
    geom_ribbon(aes(x=time,ymin=`2.5%`,ymax=`97.5%`),fill=natparks.pals("SouthDowns",3)[1],alpha=0.2) + 
    geom_line(aes(x=time,y=`50%`),col=natparks.pals("SouthDowns",3)[1])  + 
    geom_hline(data=coef_sims_null[c(1,3),],aes(yintercept=value),linetype="dashed",col="grey40") +
    scale_x_continuous(expand=c(0.01,0.01),breaks=x_breaks,labels=x_labels) + 
    theme_classic()+
    geom_hline(data=overall_var_dat,aes(yintercept=y),col=natparks.pals("SouthDowns",3)[1],linewidth=1)+
    ylab("Coefficient of variation") +
    xlab("Circulation time")+
    theme(legend.position="none",
          axis.text.x=element_text(family="sans",size=7,colour="black",angle=90,vjust=0.5),
          axis.title.y=element_text(family="sans",size=12,colour="black"),
          axis.text.y=element_text(family="sans",size=7,colour="black"),
          plot.tag = element_text(face="bold"),
          strip.background = element_blank())+
    labs(tag="A") +
    jahR::theme_no_x_axis()

ar_var_pB <- ggplot(coef_var_by_time_res %>% filter(variable == "Mean of location-specific attack rates")) + 
    geom_ribbon(aes(x=time,ymin=`2.5%`,ymax=`97.5%`),fill=natparks.pals("SouthDowns",3)[2],alpha=0.2) + 
    geom_line(aes(x=time,y=`50%`),col=natparks.pals("SouthDowns",3)[2])  + 
    scale_x_continuous(expand=c(0.01,0.01),breaks=x_breaks,labels=x_labels) + 
    theme_classic()+
    ylab("Posterior mean annual\n attack rate") +
    xlab("Circulation time")+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25)) +
    theme(legend.position="none",
          axis.text.x=element_text(family="sans",size=7,colour="black",angle=90,vjust=0.5),
          axis.title.y=element_text(family="sans",size=12,colour="black"),
          axis.text.y=element_text(family="sans",size=7,colour="black"),
          plot.tag = element_text(face="bold"),
          strip.background = element_blank())+
    scale_fill_manual(values=natparks.pals("SouthDowns",3)) +
    scale_color_manual(values=natparks.pals("SouthDowns",3)) +
    labs(tag="B") +
    jahR::theme_no_x_axis()

ar_var_pC <- ggplot(coef_var_by_time_res %>% filter(variable == "Standard deviation of location-specific attack rates")) + 
    geom_ribbon(aes(x=time,ymin=`2.5%`,ymax=`97.5%`),fill=natparks.pals("SouthDowns",3)[3],alpha=0.2) + 
    geom_line(aes(x=time,y=`50%`),col=natparks.pals("SouthDowns",3)[3])  + 
    scale_x_continuous(expand=c(0.01,0.01),breaks=x_breaks,labels=x_labels) + 
    theme_classic()+
    ylab("Standard deviation") +
    xlab("Circulation time")+
    theme(legend.position="none",
          axis.text.x=element_text(family="sans",size=7,colour="black",angle=90,vjust=0.5),
          axis.title.y=element_text(family="sans",size=12,colour="black"),
          axis.text.y=element_text(family="sans",size=7,colour="black"),
          plot.tag = element_text(face="bold"),
          strip.background = element_blank())+
    scale_fill_manual(values=natparks.pals("SouthDowns",3)) +
    scale_color_manual(values=natparks.pals("SouthDowns",3)) +
    labs(tag="C") 


ar_loc_p <- ggplot(ar_estimates_by_loc) + 
  geom_ribbon(aes(x=time,ymin=lower_quantile,ymax=upper_quantile),fill="red",alpha=0.2) + 
  geom_line(aes(x=time,y=median),col="red") +
  facet_wrap(~LOC_ID,ncol=8) + 
  theme_classic() + 
  scale_x_continuous(expand=c(0.01,0.01),breaks=x_breaks,labels=x_labels) + 
  scale_y_continuous(expand=c(0.01,0.01),limits=c(0,1)) +
  ylab("Estimated quarterly \nincidence per capita") +
  xlab("") +
  theme( 
    axis.text.x=element_text(family="sans",size=7,colour="black",angle=90,vjust=0.5),
    #axis.text.x=element_blank(),
    axis.title.y=element_text(family="sans",size=10,colour="black"),
    #panel.grid.major=element_line(colour="gray40"),
    #axis.ticks.x=element_blank(),
    #plot.margin=margin(c(10,0,-10,0)),
    axis.text.y=element_text(family="sans",size=7,colour="black")
    #strip.background = element_blank(),
    #strip.text.x = element_blank()
  )
ar_loc_p
#fit_quantiles_comb %>% rename(`Distance (km)`=x,`Aggregation of infection histories/metric`=ver2,`Included time period`=ver) %>%
#  write.csv(file="~/Documents/GitHub/fluscape_infection_histories/data/figure_data/FigS9D.csv",row.names=FALSE)
## Get average attack rate over time
data.table::setkey(total_infs, "sampno","LOC_ID")
loc_quarterly_ar <- total_infs[,list(median(ar)),by=key(total_infs)]
data.table::setkey(loc_quarterly_ar, "LOC_ID")
average_quarterly_ar <- loc_quarterly_ar[, list(lower_quantile=quantile(V1,0.025),
                                                   median=quantile(V1,0.5),
                                                   upper_quantile=quantile(V1,0.975),
                                                   sd=sd(V1)
                                                   ),by=key(loc_quarterly_ar)]

locs <- locs[order(locs$rank.DIST_FRM_GZ),]
locs$rank.trav.min <- order(locs$trav.min)
locs <- locs[order(locs$rank.trav.min),]
loc_quarterly_ar$LOC_ID <- factor(loc_quarterly_ar$LOC_ID,levels=locs$LOC_ID,ordered=TRUE)
p_average_ar <- ggplot(loc_quarterly_ar) + 
  geom_violin(aes(x=as.factor(LOC_ID),y=V1),fill="grey70",alpha=0.5,
              trim=FALSE,scale="width",
              draw_quantiles = c(0.025,0.5,0.975)) + 
  #scale_y_continuous(limits=c(0.05,0.2), breaks=seq(0.05,0.2,by=0.005),labels=seq(0.05,0.2,by=0.005),expand=c(0,0)) +
  #scale_y_continuous(limits=c(0.04,0.1), breaks=seq(0.04,0.1,by=0.005),labels=seq(0.04,0.1,by=0.005),expand=c(0,0))+
  #facet_wrap(~LOC_ID,ncol=8) +
  theme_classic() + 
  theme(#plot.margin = unit(c(0,0.5,0,0),"cm"),
    axis.text.x=element_text(size=7,family="sans",color="black"),
    axis.text.y=element_text(size=7,family="sans",color="black"),
    axis.title.x=element_text(size=10,family="sans",color="black"),
    axis.title.y=element_text(size=10,family="sans",color="black")) +
  xlab("Location (increasing distance from Guangzhou)") +
  ylab("Median annual attack rate")
p_average_ar

p_average_ar_hist <- ggplot(average_quarterly_ar) + 
  geom_histogram(aes(x=median),binwidth=0.005,fill="grey70",col="black") +
  #scale_x_continuous(limits=c(0.05,0.2), breaks=seq(0.05,0.2,by=0.005),labels=seq(0.05,0.2,by=0.005),expand=c(0,0)) +
  #scale_x_continuous(limits=c(0.04,0.1), breaks=seq(0.04,0.1,by=0.005),labels=seq(0.04,0.1,by=0.005),expand=c(0,0))+
  theme_classic() +
  coord_flip() +
  ylab("Count") +
  xlab("Median of mean posterior\n quarterly attack rates") +
  theme(#plot.margin = unit(c(0,0.5,0,0),"cm"),
        axis.text.x=element_text(size=7,family="sans",color="black"),
        axis.text.y=element_text(size=7,family="sans",color="black"),
        axis.title.x=element_text(size=10,family="sans",color="black"),
        axis.title.y=element_text(size=10,family="sans",color="black"))
p_average_ar_hist


p_ar_loc_densities <- ggplot(loc_quarterly_ar) + 
  geom_density(aes(x=V1),fill="blueviolet") + 
  facet_wrap(~LOC_ID,ncol=8) +
  theme_classic() +
  xlab("Location (increasing distance from Guangzhou)") +
  ylab("Mean quarterly attack rate")

ar_estimates <- merge(ar_estimates_by_loc, locs, by="LOC_ID")
locs <- locs[order(locs$rank.DIST_FRM_GZ),]
ar_estimates$LOC_ID <- factor(ar_estimates$LOC_ID,levels=locs$LOC_ID,ordered=TRUE)
x_breaks <- seq(cutoff_time, 2015*buckets,by=4*buckets)
x_labels <- paste0("Q1-",seq(ceiling(cutoff_time*(1/buckets)),2014, by=4))
p_infection <- ggplot(ar_estimates) + 
  geom_tile(aes(y=LOC_ID,x=time,fill=median)) + 
  theme_bw() + 
  scale_fill_gradient2(low="blue",mid="yellow",high="red",midpoint=0.5,na.value="gray40",limits=c(0,1),
                       guide=guide_colorbar(title="Attack rate (AR, posterior median)", 
                                            direction="vertical",title.position="right",
                                            label.position="left")) +
  scale_x_continuous(expand=c(0.01,0.01),breaks=x_breaks-cutoff_time*buckets + 1,labels=x_labels) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color="black"),
        panel.background = element_rect(fill="gray40"),
        axis.text.x=element_text(family="sans",size=6,colour="black",angle=45,hjust=1),
        axis.text.y=element_text(family="sans",size=6,colour="black"),
        axis.title.y=element_text(family="sans",size=8,colour="black"),
        legend.text=element_text(size=6,family="sans"),
        legend.title=element_text(size=6,family="sans",angle=-90,hjust=0.5),
        legend.key.height=unit(0.5,"in"),
        legend.key.width=unit(0.1,"in"),
        legend.position="right",
        legend.box.background = element_rect(color="black"),
       plot.margin=margin(c(10,10,-10,10)),
        legend.box.margin=margin(c(1,1,1,1)),
        legend.key=element_rect(color="black"))+
  ylab("Location ID (closest to Guangzhou at bottom)") +
  xlab("")
p_infection


p_var <- ggplot(ar_estimates) + 
  geom_tile(aes(y=LOC_ID,x=time,fill=var)) + 
  theme_bw() + 
  scale_fill_gradient2(low="springgreen4",mid="cornsilk1",high="purple4",na.value="gray40",midpoint=max(ar_estimates$var,na.rm=TRUE)/2,
                       guide=guide_colorbar(title="Coefficient of variation of posterior AR estimate", direction="vertical",
                                            title.position="right",
                                            label.position="left")) +
    scale_x_continuous(expand=c(0.01,0.01),breaks=x_breaks-cutoff_time*buckets + 1,labels=x_labels) +
    #scale_y_continuous(expand=c(0.01,0.01),breaks=c(age_group_counts[2:length(age_group_counts)], max(wow$individual)), 
  #                   labels=c(10,20,30,40,50,60,100)) + #+
  theme(panel.grid = element_blank(),
        axis.line = element_line(color="black"),
        panel.background = element_rect(fill="gray40"),
        axis.text.x=element_text(family="sans",size=6,colour="black",angle=45,hjust=1),
        axis.text.y=element_text(family="sans",size=6,colour="black"),
        axis.title.y=element_text(family="sans",size=8,colour="black"),
        legend.text=element_text(size=6,family="sans"),
        legend.title=element_text(size=6,family="sans",angle=-90,hjust=0.5),
        legend.key.height=unit(0.5,"in"),
        legend.key.width=unit(0.1,"in"),
        legend.position="right",
        legend.box.background = element_rect(color="black"),
        plot.margin=margin(c(0,10,-10,10)),
        legend.box.margin=margin(c(1,1,1,1)),
        legend.key=element_rect(color="black"))+
  ylab("Location ID (closest to Guangzhou at bottom)") +
  xlab("")
p_var


print("Coefficient of variation simulation: ")
print(coef_sims_null)

print("Calculated coefficient of variation: ")
print(coef_var_overall_res)

p1 <- ggplot(loc_quarterly_ar) + geom_density(aes(x=V1))+facet_wrap(~LOC_ID)
p2 <- ggplot(loc_quarterly_ar) + geom_line(aes(x=sampno, y=V1))+facet_wrap(~LOC_ID)

cairo_pdf(paste0(figure_wd,run_name,"_trace.pdf"),width=6,height=7.5,family="sans")
print(plot_grid(p1,p2,nrow=2, labels=c("A","B")),align="hv")
dev.off()

p_supp1 <- plot_grid(p_infection,p_var,nrow=2, labels=c("A","B"),align="hv")
p_supp2 <- plot_grid(
    plot_grid(p_average_ar, p_average_ar_hist,rel_widths=c(4,1),align="hv",labels=c("A","B")),
    ar_loc_p, ncol=1, rel_heights = c(1,2),labels=c("","C")
)

save(ar_var_p, file=paste0(figure_wd, "/",run_name,"_coef_var_space.RData"))
ggsave_jah(p_supp1, figure_wd,paste0(run_name, "_loc_ar"),width=6,height=7.5)
ggsave_jah(p_supp2, figure_wd, paste0(run_name,"_loc_ar_averages"),width=7.5,height=8.75)

if(FALSE){
  ar_estimates %>% select(LOC_ID,time,n,lower_quantile,median,upper_quantile,precision,var) %>%
    mutate(time = (time - 1 + 1968*4)/4) %>%
    rename(`Lower 95% CrI`=lower_quantile,`Posterior median`=median,`Upper 95% CrI`=upper_quantile,Precision=precision,`Coefficient of variation`=var) %>%
    rename(Date=time) %>%
    write.csv(file="~/Documents/GitHub/fluscape_infection_histories/data/figure_data/FigS7.csv",row.names=FALSE)
}

###############################################################
## Spatial correlation -- using posterior medians
###############################################################
## Spline correlograms/spatial correlation for seroconversion

seroconv_loc <- fluscape_dat %>% select(individual, Participant_ID, LOC_ID, visit, change, Year) %>% filter(visit == "Second visit") %>% group_by(Participant_ID,individual,  LOC_ID,Year) %>% dplyr::summarize(mean_change = mean(change)) %>% mutate(seroconv = mean_change >= 2) %>% group_by(LOC_ID, Year) %>% dplyr::summarize(seroconv=sum(seroconv)/n())

seroconv_loc <- seroconv_loc %>% left_join(locs) %>% select(Year, seroconv,LOC_Lat,LOC_Long) %>% pivot_wider(names_from=Year,values_from=seroconv)

res_seroconv <- Sncf(y=seroconv_loc$LOC_Lat, x=seroconv_loc$LOC_Long,
                     z=seroconv_loc[,4:ncol(seroconv_loc)],na.rm=TRUE,latlon=TRUE)
plot(res_seroconv)

res_seroconv_recent <- Sncf(y=seroconv_loc$LOC_Lat, x=seroconv_loc$LOC_Long,
                     z=seroconv_loc[,(ncol(seroconv_loc)-3):ncol(seroconv_loc)],na.rm=TRUE,latlon=TRUE)
plot(res_seroconv_recent)


break
ar_estimates_by_loc_old <- ar_estimates_by_loc
ar_estimates_by_loc <- ar_estimates_by_loc %>% left_join(locs)
ars_long <- as.data.frame(dcast(ar_estimates_by_loc[,c("LOC_ID","time","median")], LOC_ID~time))
ars_long <- ars_long[,2:ncol(ars_long)]

res_all <- Sncf(y=ar_estimates_by_loc$LOC_Lat, x=ar_estimates_by_loc$LOC_Long,
     z=ars_long,na.rm=TRUE,latlon=TRUE)

## Quick look at relationship between AR and population density of location, urban/rural and travel time to GZ
tmp <- total_infs %>% group_by(LOC_ID,time,n) %>% dplyr::summarize(y=round(median(V1))) %>% left_join(locs) %>% as.data.frame()%>% filter(n > 0)

## First just GLM assuming linear relationship
fit <- glm(y/n ~ as.factor(time) + trav.min + URBAN + log10(dens.1),data=tmp,#tmp %>% filter(sampno == 4),
           family=binomial(link="logit"),weights=n)
summary(fit)
exp(fit$coefficients)
gam.check(fit)
#plot(fit)

newdata <- expand_grid(time=1,dens.1=log10(seq(min(tmp$dens.1),max(tmp$dens.1),length.out=100)),URBAN=0,trav.min=seq(min(tmp$trav.min),max(tmp$trav.min),length.out=100))
newdata$predictions <- predict(fit,newdata=newdata,type="response")
ggplot(newdata %>% filter(trav.min==min(trav.min))) + geom_line(aes(x=dens.1,y=predictions)) + scale_y_continuous(limits=c(0,1))
ggplot(newdata %>% filter(dens.1==min(dens.1))) + geom_line(aes(x=trav.min,y=predictions)) + scale_y_continuous(limits=c(0,1))
## Maybe a very slight relationship with travel distance, but a tiny effect

## For recent times only
## First just GLM assuming linear relationship
fit <- glm(y/n ~ as.factor(time) + trav.min + URBAN + log10(dens.1),data=tmp %>% filter(time >= 41),#tmp %>% filter(sampno == 4),
           family=binomial(link="logit"),weights=n)
summary(fit)
exp(fit$coefficients)
gam.check(fit)
#plot(fit)
## Definitely no relationship

## Second, using GAM to test non-linear relationship with these factors
fit_gam <- gam(y/n ~ as.factor(time) + s(trav.min) + URBAN + s(log10(dens.1)),data=tmp,
           family=binomial(link="logit"),weights=n,method="REML")

## Again, small relationship with travel distance
summary(fit_gam)
exp(fit_gam$coefficients)
gam.check(fit_gam)
plot(fit_gam)

newdata <- expand_grid(time=1,dens.1=log10(seq(min(tmp$dens.1),max(tmp$dens.1),length.out=100)),URBAN=0,trav.min=seq(min(tmp$trav.min),max(tmp$trav.min),length.out=100))
newdata$predictions <- predict(fit_gam,newdata=newdata,type="response")
ggplot(newdata %>% filter(trav.min==min(trav.min))) + geom_line(aes(x=dens.1,y=predictions)) + scale_y_continuous(limits=c(0,1))
ggplot(newdata %>% filter(dens.1==min(dens.1))) + geom_line(aes(x=trav.min,y=predictions)) + scale_y_continuous(limits=c(0,1))

#### Look at spline correlograms
ar_ests_recent <- ar_estimates_by_loc %>% filter(time >= (2009 - 1968))
ars_long_recent <- as.data.frame(dcast(ar_ests_recent[,c("LOC_ID","time","median")], LOC_ID~time))
ars_long_recent <- ars_long_recent[,2:ncol(ars_long_recent)]

res_recent <- Sncf(y=ar_ests_recent$LOC_Lat, x=ar_ests_recent$LOC_Long,
                z=ars_long_recent,na.rm=TRUE,latlon=TRUE,resamp=1000)

ar_ests_old <- ar_estimates_by_loc %>% filter(time < (2009 - 1968))
ars_long_old <- as.data.frame(dcast(ar_ests_old[,c("LOC_ID","time","median")], LOC_ID~time))
ars_long_old <- ars_long_old[,2:ncol(ars_long_old)]

res_old <- Sncf(y=ar_ests_old$LOC_Lat, x=ar_ests_old$LOC_Long,
                   z=ars_long_old,na.rm=TRUE,latlon=TRUE,resamp=1000)


## Spatial correlation -- using posterior samples
total_infs <- total_infs %>% left_join(locs)
n_samps <- 100
use_samps <- sample(unique(total_infs$sampno),n_samps)
all_res <- vector("list",n_samps)
all_res_recent <- vector("list",n_samps)
all_res_old <- vector("list",n_samps)

for(i in seq_along(use_samps)){
    ## All attack rates
    tmp <- total_infs %>% filter(sampno == use_samps[i])
    tmp_ars <- as.data.frame(dcast(tmp[,c("LOC_ID","time","ar")], LOC_ID~time,value.var="ar"))
    tmp_ars <- tmp_ars[,2:ncol(tmp_ars)]
    tmp <- tmp %>% select(LOC_ID, LOC_Lat,LOC_Long) %>% distinct()
    tmp_res <- Sncf(y=tmp$LOC_Lat, x=tmp$LOC_Long,z=tmp_ars,na.rm=TRUE,latlon=TRUE,resamp=100)
    all_res[[i]] <- data.frame(x=tmp_res$real$predicted$x[1,], y=tmp_res$real$predicted$y[1,], samp=use_samps[i])
}

for(i in seq_along(use_samps)){
    ## Recent times
    tmp_recent <- total_infs %>% filter(sampno == use_samps[i]) %>% filter(time >= (2009*buckets - 1968*buckets))
    ars_long_recent <- as.data.frame(dcast(tmp_recent[,c("LOC_ID","time","ar")], LOC_ID~time,value.var="ar"))
    ars_long_recent <- ars_long_recent[,2:ncol(ars_long_recent)]
    tmp_recent <- tmp_recent %>% select(LOC_ID, LOC_Lat,LOC_Long) %>% distinct()
    
    res_recent <- Sncf(y=tmp_recent$LOC_Lat, x=tmp_recent$LOC_Long,z=ars_long_recent,na.rm=TRUE,latlon=TRUE,resamp=100)
    all_res_recent[[i]] <- data.frame(x=res_recent$real$predicted$x[1,], y=res_recent$real$predicted$y[1,], samp=use_samps[i])
}

for(i in seq_along(use_samps)){
    ## Older times
    tmp_old <- total_infs %>% filter(sampno == use_samps[i]) %>% filter(time < (2009*buckets - 1968*buckets))
    ars_long_old <- as.data.frame(dcast(tmp_old[,c("LOC_ID","time","ar")], value.var="ar",formula=LOC_ID~time))
    ars_long_old <- ars_long_old[,2:ncol(ars_long_old)]
    tmp_old <- tmp_old %>% select(LOC_ID, LOC_Lat,LOC_Long) %>% distinct()
    
    res_old <- Sncf(y=tmp_old$LOC_Lat, x=tmp_old$LOC_Long,z=ars_long_old,na.rm=TRUE,latlon=TRUE,resamp=100)
    all_res_old[[i]] <- data.frame(x=res_old$real$predicted$x[1,], y=res_old$real$predicted$y[1,], samp=use_samps[i])
}

fits_all <- do.call("bind_rows",all_res)
fit_quantiles <- ddply(fits_all, ~x, function(tmp) quantile(tmp$y, c(0.5,0.025,0.975)))
fits_recent <- do.call("bind_rows",all_res_recent)
fit_quantiles_recent <- ddply(fits_recent, ~x, function(tmp) quantile(tmp$y, c(0.5,0.025,0.975)))
fits_old <- do.call("bind_rows",all_res_old)
fit_quantiles_old <- ddply(fits_old, ~x, function(tmp) quantile(tmp$y, c(0.5,0.025,0.975)))

fit_quantiles_all <- bind_rows(fit_quantiles %>% mutate(ver="All times"),
                               fit_quantiles_recent%>% mutate(ver="Recent time (2009 onwards)"),
                               fit_quantiles_old%>% mutate(ver="Historic (pre-2009)"))

seroconv_dat <- data.frame(x=res_seroconv$real$predicted$x[1,], "50%" =res_seroconv$real$predicted$y[1,],
                           "2.5%"=res_seroconv$boot$boot.summary$predicted$y[2,],
                           "97.5%"=res_seroconv$boot$boot.summary$predicted$y[10,],ver2="All strains",ver="Seroconversion")
seroconv_recent_dat <- data.frame(x=res_seroconv_recent$real$predicted$x[1,], 
                                  "50%"=res_seroconv_recent$real$predicted$y[1,],
                                  "2.5%"=res_seroconv_recent$boot$boot.summary$predicted$y[2,],
                                  "97.5%"=res_seroconv_recent$boot$boot.summary$predicted$y[10,],ver2="Recent strains\n( A/Victoria/2009 onwards)",ver="Seroconversion")

colnames(seroconv_dat)[2:4] <- c("50%","2.5%","97.5%")
colnames(seroconv_recent_dat)[2:4] <- c("50%","2.5%","97.5%")
break
save(fit_quantiles_all,
     file=paste0(figure_wd,"/corr_plot_",run_name,".RData"))
combine<-TRUE
if(combine){
  fit_quantiles_all_annual <- fit_quantiles_all
  fit_quantiles_all_annual$ver2 <- "annual"
  load(paste0(figure_wd,"/corr_plot_","quarter",".RData"))
  fit_quantiles_comb <- bind_rows(fit_quantiles_all_annual,fit_quantiles_all)
} else {
  fit_quantiles_comb <- fit_quantiles_all
  fit_quantiles_comb$ver2 <- run_name
}
run_name_key <- c("per-quarter"="Per-quarter","annual"="Annual")
fit_quantiles_comb$ver2 <- run_name_key[fit_quantiles_comb$ver2]

fit_quantiles_comb <- bind_rows(fit_quantiles_comb, seroconv_dat,seroconv_recent_dat)

p_corr <- ggplot(fit_quantiles_comb) + 
    geom_hline(yintercept=0) +
    geom_ribbon(aes(x=x,ymin=`2.5%`,ymax=`97.5%`,fill=ver2),alpha=0.25) + 
    geom_line(aes(x=x,y=`50%`,col=ver2)) + 
    scale_y_continuous(limits=c(-1,1),breaks=seq(-1,1,by=0.25)) + 
    scale_x_continuous(expand=c(0,0)) +
    facet_wrap(~ver, ncol=1) +
    theme_classic() + 
    ylab("Correlation coefficient") +
    xlab("Distance (km)") +
  scale_fill_manual(name="Aggregation of infection\n histories/metric",values=c("Annual"="red","Per-quarter"="blue",
                                                                       "All strains"="purple","Recent strains\n( A/Victoria/2009 onwards)"="orange")) +
  scale_color_manual(name="Aggregation of infection\n histories/metric",values=c("Annual"="red","Per-quarter"="blue",
                                                                               "All strains"="purple","Recent strains\n( A/Victoria/2009 onwards)"="orange")) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
    theme(legend.position="bottom",
          legend.text = element_text(size=7),
          legend.title=element_text(size=7),
          panel.grid.major = element_line(linewidth=0.1,color="grey70"),
          axis.text.x=element_text(family="sans",size=7,colour="black"),
          axis.title.y=element_text(family="sans",size=12,colour="black"),
          axis.text.y=element_text(family="sans",size=7,colour="black"),
          plot.tag=element_text(face='bold'),
          strip.background = element_blank()) +
    labs(tag="D")


p_lhs <- ar_var_pA / ar_var_pB / ar_var_pC
fig_ar_corr <- p_lhs | p_corr
#fit_quantiles_comb %>% rename(`Distance (km)`=x,`Aggregation of infection histories/metric`=ver2,`Included time period`=ver) %>%
#  write.csv(file="~/Documents/GitHub/fluscape_infection_histories/data/figure_data/FigS9D.csv",row.names=FALSE)

ggsave_jah(fig_ar_corr, figure_wd,paste0("location_ar_corr_",run_name),width=8,height=8)
