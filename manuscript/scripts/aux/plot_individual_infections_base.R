needed_fluscape_dat <- unique(fluscape_dat[,c("individual","Participant_ID","dens.1","dens.9","DOB","order","age_group","age","birth_cohort")])
needed_fluscape_dat <- data.table(needed_fluscape_dat)
DOBs <- get_DOBs(titre_dat)[,2]
age_mask <- create_age_mask(DOBs,strain_isolation_times)
strain_mask <- create_strain_mask(titre_dat, strain_isolation_times)
masks <- cbind(age_mask, strain_mask)
masks <- data.frame(masks)
masks$individual <- 1:nrow(masks)

#########################################
## POSTERIOR DENSITIES
## Extract posterior probability of infection for each individual for each time
data.table::setkey(inf_chain, "i", "j")
max_sampno <- length(unique(inf_chain$samp))
## Number of samples with a 1 divided by total samples
densities <- inf_chain[, list(V1 = sum(x) / max_sampno), by = key(inf_chain)]

colnames(densities) <- c("individual","time","Posterior probability of infection")

wow <- merge(needed_fluscape_dat, densities,by="individual")
wow <- merge(wow, data.table(masks),by="individual")

wow[wow$time < wow$age_mask | wow$time > wow$strain_mask,"Posterior probability of infection"] <- NA


len.cols <- length(0:8)
cols <- colorRampPalette(c("royalblue", "khaki1", "red"))(len.cols)
## Generate age group breaks
tmp <- unique(wow[,c("individual","age_group")])
age_group_counts <- c(0,cumsum(ddply(tmp, ~age_group, nrow)$V1))
age_group_counts <- age_group_counts[1:(length(age_group_counts)-1)]

x_breaks <- seq(1968*4, 2015*4,by=4)
x_labels <- paste0("Q1-",seq(1968,2015, by=1))
wow$time <- wow$time + 7871

p_infection <- ggplot(wow) + 
  #geom_tile(aes(x=order,y=virus,fill=`log titre`)) + 
  geom_tile(aes(y=order,x=time,fill=`Posterior probability of infection`)) + 
  theme_bw() + 
  scale_fill_gradient2(low="blue",mid="yellow",high="red",midpoint=0.5,na.value="gray40",
                       guide=guide_colorbar(direction="horizontal",title.position="top",
                                            label.position="bottom")) +
  scale_x_continuous(expand=c(0.01,0.01),breaks=x_breaks,labels=x_labels) +
  scale_y_continuous(expand=c(0.01,0.01),breaks=c(age_group_counts[2:length(age_group_counts)], max(wow$individual)), 
                     labels=c(10,20,30,40,50,60,100)) + #+
  theme(panel.grid = element_blank(),
        axis.line = element_line(color="black"),
        panel.background = element_rect(fill="gray40"),
        axis.text.x=element_text(family="sans",size=8,colour="black",angle=45,hjust=1),
        axis.text.y=element_text(family="sans",size=8,colour="black"),
        axis.title.y=element_text(family="sans",size=10,colour="black"),
        legend.text=element_text(size=6,family="sans"),
        legend.title=element_text(size=6,family="sans",angle=0,hjust=0.5),
        legend.key.height=unit(0.1,"in"),
        legend.key.width=unit(0.3,"in"),
        legend.justification = c(0,0),
        legend.position=c(0, 0),
        legend.box.background = element_rect(color="black"),
        plot.margin=margin(c(0,10,-10,10)),
        legend.box.margin=margin(c(1,1,1,1)),
        legend.key=element_rect(color="black"))+
  ylab("Individual") +
  xlab("")
p_infection

colnames(inf_chain)[4] <- "sampno"

first_sample <- min(titre_dat$samples)
first_sample <- 8020
inf_chain_tmp <- inf_chain[inf_chain$j >= (first_sample - 7872),]
inf_chain_tmp$j <- inf_chain_tmp$j - min(inf_chain_tmp$j) + 1
strain_isolation_times_tmp <- strain_isolation_times[strain_isolation_times >= first_sample]

virus_stuff <- unique(fluscape_dat[,c("virus","Virus")])


ar_p1 <- plot_attack_rates_monthly(inf_chain,titre_dat, strain_isolation_times,pad_chain=FALSE,add_box=TRUE)
#,fill_col="red",add_box=TRUE,x_box_min=min(titre_dat$samples)
ar_p <- ar_p1 + theme_classic() + 
  scale_x_continuous(expand=c(0.01,0.01),breaks=x_breaks,labels=x_labels) + 
  scale_y_continuous(expand=c(0.01,0.01),limits=c(0,1)) +
  geom_point(data=virus_stuff,aes(x=virus,y=1),shape=8,size=1) +
 # geom_text(data=virus_stuff,aes(x=virus,y=0.82,label=Virus),angle=45,size=3.2,hjust=0,vjust=1.2)+
  ylab("Estimated quarterly \nincidence per capita") +
  xlab("") +
  theme( 
    axis.text.x=element_text(family="sans",size=8,colour="black",angle=45,hjust=1),
    #axis.text.x=element_blank(),
    axis.title.y=element_text(family="sans",size=10,colour="black"),
    #panel.grid.major=element_line(colour="gray40"),
    #axis.ticks.x=element_blank(),
    plot.margin=margin(c(10,0,-10,0)),
    axis.text.y=element_text(family="sans",size=8,colour="black"), 
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
ar_p

first_sample <- min(titre_dat$samples)
first_sample <- 8020
inf_chain_tmp <- inf_chain[inf_chain$j >= (first_sample - 7872),]
inf_chain_tmp$j <- inf_chain_tmp$j - min(inf_chain_tmp$j) + 1
strain_isolation_times_tmp <- strain_isolation_times[strain_isolation_times >= first_sample]

x_breaks2 <- strain_isolation_times_tmp
#x_breaks2 <- seq(1968*4, 2015*4,by=1)
#x_breaks2 <- c(x_breaks2,8061,8062,8063)
x_labels2_tmp <- seq(2005, 2015,by=1)
x_labels2_tmp <- rep(x_labels2_tmp, each=4)
x_labels2 <- paste0(c("Q1-","Q2-","Q3-","Q4-"),x_labels2_tmp)
x_labels2 <- x_labels2[seq_along(x_breaks2)]

p_fig2 <- cowplot::plot_grid(ar_p, p_infection, nrow=2,align="hv",axis="lr",rel_heights = c(1,1.8),labels=c("A","B"))

ggsave_jah(p_fig2, figure_wd,"base_combined_infections_fig2",width=7.2,height=6)

infection_histories <- inf_chain
colnames(infection_histories) <- c("individual","time","x","sampno","chain_no")
infection_histories$individual <- as.numeric(infection_histories$individual)
masks$individual <- as.numeric(masks$individual)
masks <- data.table(masks)
infection_histories <- merge(infection_histories, masks[,c("individual")],by="individual")
data.table::setkey(infection_histories, "time", "sampno")
n_alive <- get_n_alive(titre_dat, strain_isolation_times)
n_alive <- data.frame("time"=unique(infection_histories$time),"n"=n_alive)
## Number of samples with a 1 divided by total samples
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


data.table::setkey(total_infs, "sampno")
ar_estimates_overall <- total_infs[,list(lower_quantile=quantile(ar,c(0.025)),
                                         median=quantile(ar,c(0.5)),
                                         upper_quantile=quantile(ar,c(0.975)),
                                         precision=1/var(ar),
                                         var=sd(ar)/mean(ar)
),by=key(total_infs)]

ar_estimates$time <- ar_estimates$time + min(strain_isolation_times) -1
ar_estimates_recent <- ar_estimates[ar_estimates$time >= min(titre_dat$samples),]

total_infs$time <- total_infs$time + min(strain_isolation_times) -1
total_infs_recent <- total_infs[total_infs$time >= min(titre_dat$samples),]

data.table::setkey(total_infs_recent, "sampno")
ar_estimates_overall_recent <- total_infs_recent[,list(lower_quantile=quantile(ar,c(0.025)),
                                         median=quantile(ar,c(0.5)),
                                         upper_quantile=quantile(ar,c(0.975)),
                                         precision=1/var(ar),
                                         var=sd(ar)/mean(ar)
),by=key(total_infs_recent)]
