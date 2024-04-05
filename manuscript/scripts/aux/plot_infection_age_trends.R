## Where to read MCMC chains from?
inf_chain1 <- inf_chain
inf_chain1 <- inf_chain1[inf_chain1$x == 1,]
colnames(inf_chain1)[colnames(inf_chain1) == "samp"] <- "sampno"
inf_chain1$chain_no <- 1
buckets <- 4

p1 <- plot_number_infections(inf_chain,FALSE)  +  
  geom_point(aes(x=individual + 1, y=median),col="red",size=0.1) +
  theme_classic() + 
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(size=8,family="sans",color="black"),
        axis.text.y=element_text(size=8,family="sans",color="black"),
        axis.title.x=element_text(size=8,family="sans",color="black"),
        axis.title.y=element_text(size=8,family="sans",color="black"),
                                  plot.tag=element_text(face="bold"))+
            labs(tag="A") +
  scale_y_continuous(limits=c(0,25),breaks=seq(0,25,by=5),expand=c(0,0))
p1

## Fitted data
## For each individual, how many infections did they have in each sample in total?
n_strain <- max(inf_chain1$j)
data.table::setkey(inf_chain1, "i", "sampno", "chain_no")
n_inf_chain <- inf_chain1[, list(V1 = sum(x)), by = key(inf_chain1)]

## Get quantiles on total number of infections per indiv across all samples
indiv_hist <- plyr::ddply(n_inf_chain, .(i), function(x) quantile(x$V1, c(0.025, 0.5, 0.975)))
colnames(indiv_hist) <- c("individual", "lower", "median", "upper")
indiv_hist <- indiv_hist[order(indiv_hist$median), ]
indiv_hist$individual <- 1:nrow(indiv_hist)
indiv_hist$ver <- "data"

## Prior
## For each individual, how many infections did they have in each sample in total?
# n_strain_prior <- max(inf_chain_prior1$j)
# data.table::setkey(inf_chain_prior1, "i", "sampno", "chain_no")
# n_inf_chain_prior <- inf_chain_prior1[, list(V1 = sum(x)), by = key(inf_chain_prior1)]
# 
# ## Get quantiles on total number of infections per indiv across all samples
# indiv_hist_prior <- plyr::ddply(n_inf_chain_prior, .(i), function(x) quantile(x$V1, c(0.025, 0.5, 0.975)))
# colnames(indiv_hist_prior) <- c("individual", "lower", "median", "upper")
# indiv_hist_prior <- indiv_hist_prior[order(indiv_hist_prior$median), ]
# indiv_hist_prior$individual <- 1:nrow(indiv_hist_prior)
# indiv_hist_prior$ver <- "prior"
# 
# indiv_hist <- rbind(indiv_hist, indiv_hist_prior)

DOBs <- get_DOBs(titre_dat)[,2]
age_mask <- create_age_mask(DOBs,strain_isolation_times)
strain_mask <- create_strain_mask(titre_dat, strain_isolation_times)
masks <- cbind(age_mask, strain_mask,DOBs)
masks <- data.frame(masks)
masks$individual <- 1:nrow(masks)
masks <- merge(masks,needed_fluscape_dat[,c("individual","age_group")] %>% distinct(),by="individual")
masks$n_exposure_years <- (masks$strain_mask - masks$age_mask + 1)/buckets


indiv_hist <- merge(indiv_hist, masks[,c("individual","n_exposure_years")],by="individual")
indiv_hist$median_density <- indiv_hist$median/indiv_hist$n_exposure_years

## Histogram of total number of infections
## Can add prior to this
p2 <- gghistogram(indiv_hist,x="median",y="..density..",fill="grey70",
                  binwidth=1,
                  xlab="Total number of lifetime infections",
                  ylab="Density") +  
  coord_flip()+
  #scale_x_continuous(limits=c(0,25)) +
  #scale_y_continuous(limits=c(0,0.15),breaks=seq(0,0.15,by=0.05),expand=c(0,0)) +
  theme(plot.margin = unit(c(0,0.5,0,0),"cm"),
        axis.text.x=element_text(size=8,family="sans",color="black"),
        plot.tag=element_text(face="bold"),
        axis.text.y=element_text(size=8,family="sans",color="black"),
        axis.title.x=element_text(size=8,family="sans",color="black"),
        axis.title.y=element_text(size=8,family="sans",color="black")) +
    labs(tag="B")
p2

## Histogram of infection frequency per year
p4 <- ggplot(indiv_hist %>%
         mutate(median_density = pmin(median_density,1)) %>%
         mutate(median_density=median_density*10)) +
  geom_density(aes(x=median_density),fill="grey70") +
  coord_cartesian(xlim=c(0,10)) +
  scale_x_continuous(breaks=seq(0,10,by=2),limits=c(0,10),labels=c(0,2,4,6,8,"10+")) +
  coord_flip()+
  #scale_x_continuous(labels=seq(0,3,by=0.25)) +
  scale_y_continuous(expand=c(0,0)) +
  theme_classic() +
  #scale_y_continuous(limits=c(0,0.15),breaks=seq(0,0.15,by=0.05),expand=c(0,0)) +
  theme(plot.margin = unit(c(0,0.5,0,0),"cm"),
        axis.text.x=element_text(size=8,family="sans",color="black"),
        axis.text.y=element_text(size=8,family="sans",color="black"),
        axis.title.x=element_text(size=8,family="sans",color="black"),
        axis.title.y=element_text(size=8,family="sans",color="black"),
        plot.tag=element_text(face="bold")) +
  labs(tag="D") +
  xlab("Number of infections per decade alive") +
  ylab("Density")

if(FALSE){
p4 <- gghistogram(indiv_hist %>% mutate(median_density=median_density*10),x="median_density",fill="grey70", y="..density..",
                  binwidth=0.25,
                  xlab="Number of infections per decade alive",
                  ylab="Density") +  
  coord_cartesian(xlim=c(0,5)) +
  scale_x_continuous(breaks=seq(0,5,by=1),limits=c(0,5)) +
  coord_flip()+
  #scale_x_continuous(labels=seq(0,3,by=0.25)) +
  scale_y_continuous(expand=c(0,0)) +
  #scale_y_continuous(limits=c(0,0.15),breaks=seq(0,0.15,by=0.05),expand=c(0,0)) +
  theme(plot.margin = unit(c(0,0.5,0,0),"cm"),
        axis.text.x=element_text(size=8,family="sans",color="black"),
        axis.text.y=element_text(size=8,family="sans",color="black"),
        axis.title.x=element_text(size=8,family="sans",color="black"),
        axis.title.y=element_text(size=8,family="sans",color="black"),
        plot.tag=element_text(face="bold")) +
    labs(tag="D")
}
p4

print("Number of infections per decade alive")
print(quantile(indiv_hist$median_density*10,c(0.025,0.5,0.975)))

## Put age at time of infection into inf_chain
colnames(masks)[1] <- "i"
masks$i <- as.numeric(masks$i)
inf_chain_age <- merge(inf_chain1,masks %>% select(i, DOBs, strain_mask) %>% distinct(),by="i")
inf_chain_age$j <- inf_chain_age$j - 1 + 1968*buckets
inf_chain_age$strain_mask <- inf_chain_age$strain_mask - 1 + 1968*buckets

inf_chain_age$age_at_infection <- (inf_chain_age$j - inf_chain_age$DOBs)/buckets

## Find duration each individual was in each age group i.e., time at risk in each age group
ranges <- titre_dat %>% select(individual, DOB,samples) %>% distinct() %>% group_by(individual,DOB) %>% filter(samples == max(samples))
ranges <- expand_grid(ranges, j=strain_isolation_times)
ranges <- ranges %>% group_by(individual) %>% filter(DOB <= j, j <= samples)
ranges <- ranges %>% mutate(age_at_infection = (j - DOB)/buckets)

#age_breaks <- c(0,10,20,40,60,100)
age_breaks <- c(0,5,10,20,30,40,50,60,70,80,100)
ranges$age_group_at_infection <- cut(ranges$age_at_infection,
                                            breaks=age_breaks,
                                     include.lowest=TRUE)
widths <- ranges %>% group_by(individual,age_group_at_infection) %>% dplyr::summarize(widths=n()) %>% distinct()

inf_chain_age$age_group_at_infection <- cut(inf_chain_age$age_at_infection,
                                            breaks=age_breaks,include.lowest=TRUE)

#age_widths <- data.frame(age_group_at_infection=unique(inf_chain_age$age_group_at_infection),widths=c(20,20,40,10,10))
inf_chain_age <- merge(inf_chain_age, widths %>% rename(i = individual),by=c("age_group_at_infection","i")) %>% arrange(i, j) %>% filter(widths >= buckets*2)


## Get number of infections per age group
setkey(inf_chain_age, "i","age_group_at_infection","sampno","chain_no","widths")
inf_chain_age1 <- inf_chain_age# %>% filter(j < min(titre_dat$samples))
setkey(inf_chain_age1, "i","age_group_at_infection","sampno","chain_no","widths")
no_infs_age <- inf_chain_age1[,list(no_infs=sum(x)),by=key(inf_chain_age1)]
no_infs_age$no_infs_per_year <- no_infs_age$no_infs/no_infs_age$widths
no_infs_age_summary <- no_infs_age %>% group_by(i, age_group_at_infection, widths) %>% dplyr::summarize(y=median(no_infs_per_year))

age_labels <- levels(no_infs_age_summary$age_group_at_infection)
age_labels[1] <- "<5"
age_labels[length(age_labels)] <- c("80+")

samp_sizes <- no_infs_age_summary %>% group_by(age_group_at_infection) %>% tally() %>% mutate("Time period"="All")
no_infs_age_summary <- no_infs_age_summary %>% mutate("Time period"="All")


## Get number of infections per age group for pre Fluscape time periods
inf_chain_age_early <- inf_chain_age %>% filter(j < min(titre_dat$samples))
setkey(inf_chain_age_early, "i","age_group_at_infection","sampno","chain_no","widths")
no_infs_age_early <- inf_chain_age_early[,list(no_infs=sum(x)),by=key(inf_chain_age_early)]
no_infs_age_early$no_infs_per_year <- no_infs_age_early$no_infs/no_infs_age_early$widths
no_infs_age_summary_early <- no_infs_age_early %>% group_by(i, age_group_at_infection, widths) %>% dplyr::summarize(y=median(no_infs_per_year))

age_labels <- levels(no_infs_age_summary_early$age_group_at_infection)
age_labels[1] <- "<5"
age_labels[length(age_labels)] <- c("80+")

samp_sizes_early <- no_infs_age_summary_early %>% group_by(age_group_at_infection) %>% tally() %>% mutate("Time period"="<2009")
no_infs_age_summary_early <- no_infs_age_summary_early %>% mutate("Time period"="<2009")

p3 <- ggplot(no_infs_age_summary_early) + 
  geom_boxplot(aes(y=y*buckets*10,x=age_group_at_infection),outlier.size=0.1,fill="grey70") + 
  geom_text(data=samp_sizes_early,y=10.25, aes(x=age_group_at_infection,label=paste0("N=",n)),size=3) +
  theme_classic() +
  xlab("Age group at time of infection (years)") +
  ylab("Infections per decade spent in age group") +
  coord_cartesian(ylim=c(0,10.25)) +
  scale_y_continuous(breaks=seq(0,10,by=2),labels=seq(0,10,by=2))+
  scale_x_discrete(labels=age_labels) +
  theme(legend.position="none",
        axis.text.x=element_text(size=8,family="sans",color="black"),
        axis.text.y=element_text(size=8,family="sans",color="black"),
        plot.tag=element_text(face="bold"),
        axis.title.x=element_text(size=8,family="sans",color="black"),
        axis.title.y=element_text(size=8,family="sans",color="black")) +
  labs(tag="C")
p3


## Get number of infections per age group for during Fluscape time periods
inf_chain_age_late <- inf_chain_age %>% filter(j >= min(titre_dat$samples))
setkey(inf_chain_age_late, "i","age_group_at_infection","sampno","chain_no","widths")
no_infs_age_late <- inf_chain_age_late[,list(no_infs=sum(x)),by=key(inf_chain_age_late)]
no_infs_age_late$no_infs_per_year <- no_infs_age_late$no_infs/no_infs_age_late$widths
no_infs_age_summary_late <- no_infs_age_late %>% group_by(i, age_group_at_infection, widths) %>% dplyr::summarize(y=median(no_infs_per_year))

age_labels <- levels(no_infs_age_summary_late$age_group_at_infection)
age_labels[1] <- "<5"
age_labels[length(age_labels)] <- c("80+")

samp_sizes_late <- no_infs_age_summary_late %>% group_by(age_group_at_infection) %>% tally() %>% mutate("Time period"="≥2009")
no_infs_age_summary_late <- no_infs_age_summary_late %>% mutate("Time period"="≥2009")


## Combined
no_infs_age_summary_comb <- bind_rows(no_infs_age_summary, no_infs_age_summary_early,no_infs_age_summary_late)
no_infs_age_summary_comb$`Time period` <- factor(no_infs_age_summary_comb$`Time period`,levels=c("<2009","≥2009","All"))
samp_sizes_comb <- bind_rows(samp_sizes, samp_sizes_late, samp_sizes_early)
samp_sizes_comb$`Time period` <- factor(samp_sizes_comb$`Time period`,levels=c("<2009","≥2009","All"))
p3_alt <- ggplot(no_infs_age_summary_comb) + 
  geom_boxplot(aes(y=y*buckets*10,x=age_group_at_infection,fill=`Time period`),outlier.size=0.1,position="dodge",alpha=0.25) + 
  #geom_text(data=samp_sizes_comb,y=5.25, aes(x=age_group_at_infection,label=paste0("N=",n),col=`Time period`),size=3) +
  theme_classic() +
  xlab("Age group at time of infection (years)") +
  ylab("Infections per decade spent in age group") +
  coord_cartesian(ylim=c(0,10.25)) +
  scale_y_continuous(breaks=seq(0,10,by=2),labels=seq(0,10,by=2))+
  scale_x_discrete(labels=age_labels) +
  scale_fill_manual(name="Time period",values=c("All"="blue","<2009"="grey70","≥2009"="orange")) +
  theme(legend.position=c(0.9,0.8),
        axis.text.x=element_text(size=8,family="sans",color="black"),
        axis.text.y=element_text(size=8,family="sans",color="black"),
        plot.tag=element_text(face="bold"),
        axis.title.x=element_text(size=8,family="sans",color="black"),
        axis.title.y=element_text(size=8,family="sans",color="black"))
p3_alt




p_lhs <- p1 / p3
p_rhs <- p2 / p4


fig_age <-  (p_lhs | p_rhs) + plot_layout(widths=c(3,1),ncol=2)

ggsave_jah(fig_age,figure_wd, "age_trends",width=8,height=7)
ggsave_jah(p3_alt,figure_wd, "age_trends_comparison",width=7,height=4)
if(FALSE){
## Quick check to see if any participant variables affect infection frequency
library(brms)
tmp <- no_infs_age %>% as.data.table()
setkey(tmp, "i","age_group_at_infection","widths")
tmp1 <- tmp[,list(y=floor(median(no_infs))),by=key(tmp)]

tmp1 <- tmp1 %>% rename(individual=i) %>% 
  left_join(fluscape_dat[,c("individual")] %>% distinct()) %>%
  mutate(TRUE_OCC_STATUS = if_else(age_group_at_infection %in% c("[0,5]","(5,10]","(10,20]"),5,TRUE_OCC_STATUS))%>%#,
         #TRUE_OCC_STATUS = if_else(TRUE_OCC_STATUS == 4 & !(age_group_at_infection %in% c("(70,80]","(80,100]")),999,TRUE_OCC_STATUS)) %>%
  mutate(y = pmin(y, widths),
         y = floor(y)) %>%
  mutate(y=as.integer(y),
         widths=as.integer(widths))

## Gender
## Age at infection
## Distance from GZ (not super reliable)
## SES
## Occupation (if recent)
## Household members (if recent)
## Animals in household
## Running water

fit1 <- glm(y/widths ~ age_group_at_infection + as.factor(PART_GENDER) + DIST_FRM_GZ,
           family=binomial(),
           weights=widths,
           data=tmp1 %>%  drop_na())
summary(fit1)

fit <- brm(y | trials(widths) ~ age_group_at_infection + as.factor(PART_GENDER) + as.factor(URBAN) + DIST_FRM_GZ + as.factor(TRUE_OCC_STATUS),
           family=binomial(),
           data=tmp1 %>%  drop_na())
summary(fit)
exp(coef(fit))
}