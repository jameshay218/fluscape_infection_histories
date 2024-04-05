library(serosolver)
times <- seq_along(strain_isolation_times)
needed_fluscape_dat <- unique(fluscape_dat[,c("individual","Participant_ID",
                                              "dens.1","dens.9","DOB","order","age_group","age","birth_cohort")])
needed_fluscape_dat <- data.table(needed_fluscape_dat)
DOBs <- get_DOBs(titre_dat)[,2]
age_mask <- create_age_mask(DOBs,strain_isolation_times)
strain_mask <- create_strain_mask(titre_dat, strain_isolation_times)
masks <- cbind(age_mask, strain_mask)
masks <- data.frame(masks)
masks$individual <- 1:nrow(masks)
masks <- merge(masks,needed_fluscape_dat[,c("individual","birth_cohort")],by="individual")


n_alive_group <- ddply(masks, ~birth_cohort, function(x) {
  sapply(times, function(y){
    nrow(x[x$age_mask <= y & x$strain_mask >= y,])
  })
})

n_alive_group <- melt(n_alive_group)
colnames(n_alive_group) <- c("birth_cohort","time","n")
n_alive_group$time <- as.integer(n_alive_group$time)

infection_histories <- inf_chain
colnames(infection_histories) <- c("individual","time","x","sampno","chain_no")
infection_histories$individual <- as.numeric(infection_histories$individual)
masks$individual <- as.numeric(masks$individual)
masks <- data.table(masks)
infection_histories <- merge(infection_histories, masks[,c("individual","birth_cohort")],by="individual")
data.table::setkey(infection_histories, "time", "birth_cohort","sampno","chain_no")

## Number of samples with a 1 divided by total samples
total_infs <- infection_histories[, list(V1 = sum(x)), by = key(infection_histories)]
total_infs <- merge(total_infs, n_alive_group, by=c("time","birth_cohort"))
total_infs$ar <- total_infs$V1/total_infs$n
total_infs[is.nan(total_infs$ar),"ar"] <- NA
total_infs <- total_infs[,c("time","birth_cohort","sampno","ar")]
data.table::setkey(total_infs, "time", "birth_cohort")
total_infs$sampno <- as.numeric(total_infs$sampno)
total_infs <- total_infs[complete.cases(total_infs),]
ar_estimates <- total_infs[,list(lower_quantile=quantile(ar,c(0.025)),
                                 median=quantile(ar,c(0.5)),
                                 upper_quantile=quantile(ar,c(0.975)),
                                 precision=1/var(ar),
                                 var=sd(ar)
                                 ),by=key(total_infs)]


x_breaks <- seq(1968*4, 2015*4,by=8)
x_labels <- paste0("Q1-",seq(1968,2015, by=2))
p_infection <- ggplot(ar_estimates) + 
  #geom_tile(aes(x=order,y=virus,fill=`log titre`)) + 
  geom_tile(aes(y=birth_cohort,x=time,fill=median)) + 
  theme_bw() + 
  scale_fill_gradient2(low="blue",mid="yellow",high="red",midpoint=0.5,na.value="gray40",limits=c(0,1),
                       guide=guide_colorbar(title="Attack rate (AR, posterior median)", direction="horizontal",title.position="top",
                                            label.position="bottom")) +
  scale_x_continuous(expand=c(0.01,0.01),breaks=x_breaks-1968*4,labels=x_labels) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color="black"),
        panel.background = element_rect(fill="gray40"),
        axis.text.x=element_text(family="sans",size=7,colour="black",angle=45,hjust=1),
        axis.text.y=element_text(family="sans",size=8,colour="black"),
        axis.title.y=element_text(family="sans",size=8,colour="black"),
        legend.text=element_text(size=6,family="sans"),
        legend.title=element_text(size=6,family="sans",angle=0,hjust=0.5),
        legend.key.height=unit(0.1,"in"),
        legend.key.width=unit(0.35,"in"),
        legend.justification = c(0,0),
        legend.position=c(0, 0),
        legend.box.background = element_rect(color="black"),
        plot.margin=margin(c(0,10,-10,10)),
        plot.tag=element_text(face="bold") ,
        legend.box.margin=margin(c(1,1,1,1)),
        legend.key=element_rect(color="black"))+
  ylab("Birth cohort") +
  xlab("") +
labs(tag="B") 
p_infection


p_var <- ggplot(ar_estimates) + 
  #geom_tile(aes(x=order,y=virus,fill=`log titre`)) + 
  geom_tile(aes(y=birth_cohort,x=time,fill=var)) + 
  theme_bw() + 
  scale_fill_gradient2(low="springgreen4",mid="cornsilk1",high="purple4",na.value="gray40",midpoint=max(ar_estimates$var)/2,
                       guide=guide_colorbar(title="Standard deviation of posterior AR estimate", direction="horizontal",
                                            title.position="top",title.hjust = 0.5,
                                            label.position="bottom")) +
  scale_x_continuous(expand=c(0.01,0.01),breaks=x_breaks-1968*4,labels=x_labels) +
  #scale_y_continuous(expand=c(0.01,0.01),breaks=c(age_group_counts[2:length(age_group_counts)], max(wow$individual)), 
  #                   labels=c(10,20,30,40,50,60,100)) + #+
  theme(panel.grid = element_blank(),
        axis.line = element_line(color="black"),
        panel.background = element_rect(fill="gray40"),
        axis.text.x=element_text(family="sans",size=7,colour="black",angle=45,hjust=1),
        axis.text.y=element_text(family="sans",size=8,colour="black"),
        axis.title.y=element_text(family="sans",size=8,colour="black"),
        legend.text=element_text(size=6,family="sans"),
        legend.title=element_text(size=6,family="sans",angle=0,hjust=0),
        legend.key.height=unit(0.1,"in"),
        legend.key.width=unit(0.35,"in"),
        legend.justification = c(0,0),
        legend.position=c(0, 0),
        plot.tag=element_text(face="bold") ,
        legend.box.background = element_rect(color="black"),
        plot.margin=margin(c(0,10,-10,10)),
        legend.box.margin=margin(c(1,1,1,1)),
        legend.key=element_rect(color="black"))+
    ylab("Birth cohort") +
    labs(tag="C") +
    xlab("")
p_var


p_ar_age <- ggplot(ar_estimates) + 
  #geom_tile(aes(x=order,y=virus,fill=`log titre`)) +
  geom_rect(xmin=min(titre_dat$samples) - 1968*4, xmax=max(ar_estimates$time),ymin=-1,ymax=2,fill="gray90",alpha=0.1) +
  geom_ribbon(aes(fill=birth_cohort,x=time,ymin=lower_quantile,ymax=upper_quantile),alpha=0.1) +
  geom_line(aes(col=birth_cohort,x=time,y=median)) + 
  theme_bw() + 
  #scale_fill_gradient2(low="blue",mid="yellow",high="red",midpoint=0.5,na.value="gray40",limits=c(0,1),
  #                     guide=guide_colorbar(title="Attack rate (AR, posterior median)", direction="horizontal",title.position="top",
  #                                          label.position="bottom")) +
  scale_x_continuous(expand=c(0.01,0.01),breaks=x_breaks-1968*4,labels=x_labels) +
  scale_y_continuous(expand=c(0,0)) +
  scale_color_viridis_d(name="Birth cohort", option="turbo") +
  scale_fill_viridis_d(name="Birth cohort", option="turbo") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color="black"),
        #panel.background = element_rect(fill="gray40"),
        axis.text.x=element_text(family="sans",size=7,colour="black",angle=45,hjust=1),
        axis.text.y=element_text(family="sans",size=8,colour="black"),
        axis.title.y=element_text(family="sans",size=8,colour="black"),
        legend.text=element_text(size=6,family="sans"),
        legend.title=element_text(size=6,family="sans",angle=0,hjust=0.5),
        legend.key.height=unit(0.1,"in"),
        legend.key.width=unit(0.35,"in"),
        legend.justification = c(1,1),
        legend.position=c(1, 1),
        legend.box.background = element_rect(color="black"),
        #plot.margin=margin(c(0,10,-10,10)),
        plot.tag=element_text(face="bold") ,
        legend.box.margin=margin(c(1,1,1,1)),
        legend.key=element_rect(color="black"))+
  ylab("Per capita incidence per quarter") +
  xlab("")
p_ar_age

ggsave_jah(p_infection/p_var, figure_wd, "age_ar",width=7,height=8)
ggsave_jah(p_ar_age, figure_wd, "age_ar_simple",width=7,height=4)




