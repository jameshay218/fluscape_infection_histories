library(ncf)
library(Hmisc)
library(corrplot)
library(gganimate)
library(data.table)
library(coda)
library(data.table)
library(plyr)
library(bayesplot)
library(ggplot2)

library("lattice")
library("sp")
library("raster")
#library("maptools")

#devtools::load_all("~/Documents/Fluscape/serosolver")
#setwd("~/Drive/Influenza/serosolver/FluScape_epidemiology/")
buckets <- 4
source(paste0(fluscape_wd,"manuscripts/paired_titer/rCodes/GeneralUtility.r"))
source(paste0(fluscape_wd,"manuscripts/paired_titer/rCodes/GeneralUtility_BY.r"))

#locs <- read.csv("data/loc_data_clean.csv")
locs <- load.and.merge.locs.V1.V2.V3(topdir=fluscape_wd)
load("r_data/fluscape_dat.RData")

## Need to plot with jittered locations
locs_use <- read.csv("~/Documents/GitHub/fluscape_serosolver/data/loc_data_rnd.csv")

locs <- locs %>% select(-c(LOC_Lat,LOC_Long)) %>% left_join(locs_use %>% select(-c(LOC_Lat_original,LOC_Long_original))) %>% rename(LOC_Lat=LOC_Lat_rnd,
                                                                                                                                    LOC_Long=LOC_Long_rnd)


## Some potentially useful functions
euc_dist <- function(x1,x2,y1,y2){
  sqrt((x2 - x1)^2 + (y2-y1)^2)
}

unique_locs <- unique(locs$LOC_ID)
distances <- matrix(ncol=length(unique_locs),nrow=length(unique_locs))
for(i in seq_along(unique_locs)){
  x1 <- locs[locs$LOC_ID == unique_locs[i],"LOC_Long"]
  y1 <- locs[locs$LOC_ID == unique_locs[i],"LOC_Lat"]
  
  for(j in seq_along(unique_locs)){
    x2 <- locs[locs$LOC_ID == unique_locs[j],"LOC_Long"]
    y2 <- locs[locs$LOC_ID == unique_locs[j],"LOC_Lat"]
    distances[i,j] <- euc_dist(x1, x2, y1, y2)
  }
}

load("r_data/measurement_quarter_runs.RData")
cutoff_time <- 1968*buckets
max_time <- max(titre_dat$samples)
#load("r_data/annual_meas_bb.RData")
#load("r_data/base_quarter_estimates.RData")

times <- seq_along(strain_isolation_times)
times <- times[which(strain_isolation_times >= cutoff_time & strain_isolation_times <= max_time)]
strain_isolation_times_tmp <- strain_isolation_times[strain_isolation_times >= cutoff_time & strain_isolation_times <= max_time]

needed_fluscape_dat <- unique(fluscape_dat[,c("individual","Participant_ID","URBAN","DIST_FRM_GZ","LOC_ID",
                                              "dens.1","dens.9","DOB","age_group","age","birth_cohort")])
needed_fluscape_dat <- data.table(needed_fluscape_dat)

devtools::load_all("~/Documents/GitHub//serosolver")
## Get masks for each individual
DOBs <- get_DOBs(titre_dat)[,2]
age_mask <- create_age_mask(DOBs,strain_isolation_times_tmp)
strain_mask <- create_strain_mask(titre_dat, strain_isolation_times_tmp)
masks <- cbind(age_mask, strain_mask)
masks <- data.frame(masks)
masks$individual <- 1:nrow(masks)
masks <- merge(masks,needed_fluscape_dat[,c("individual","LOC_ID")],by="individual")
## Create a new object in case we make a mistake
infection_histories <- inf_chain
colnames(infection_histories) <- c("individual","time","x","sampno","chain_no")
infection_histories <- infection_histories[infection_histories$time >= min(times) & infection_histories$time <= max(times),]
## Re-center so that first time is 1
infection_histories$time <- infection_histories$time - min(times) + 1
times <- times - min(times) + 1

## Merge with age masks
infection_histories$individual <- as.numeric(infection_histories$individual)
masks$individual <- as.numeric(masks$individual)
masks <- data.table(masks)
infection_histories <- merge(infection_histories, masks[,c("individual","LOC_ID")],by="individual")


# AR stats by location ----------------------------------------------------
## Get number alive in each location
n_alive_loc <- ddply(masks, ~LOC_ID, function(x) {
  sapply(times, function(y){
    nrow(x[x$age_mask <= y & x$strain_mask >= y,])
  })
})
n_alive_loc <- melt(n_alive_loc,id.vars = "LOC_ID")
colnames(n_alive_loc) <- c("LOC_ID","time","n")
n_alive_loc$time <- as.integer(n_alive_loc$time)

## Stratify infection histories by location ID
data.table::setkey(infection_histories, "time", "LOC_ID","sampno")

## For each MCMC sample, number of individuals with a 1 for 
## infection state divided by total number alive gives AR
total_infs <- infection_histories[, list(V1 = sum(x)), by = key(infection_histories)]
total_infs <- merge(total_infs, n_alive_loc, by=c("time","LOC_ID"))
total_infs$ar <- total_infs$V1/total_infs$n
total_infs[is.nan(total_infs$ar),"ar"] <- NA
total_infs <- total_infs[,c("time","LOC_ID","sampno","ar")]
data.table::setkey(total_infs, "time", "LOC_ID")
total_infs$sampno <- as.numeric(total_infs$sampno)
total_infs <- total_infs[complete.cases(total_infs),]
ar_estimates_by_loc <- total_infs[,list(lower_quantile=quantile(ar,c(0.025)),
                                        median=quantile(ar,c(0.5)),
                                        upper_quantile=quantile(ar,c(0.975)),
                                        precision=1/var(ar),
                                        var=sd(ar)/mean(ar)
),by=key(total_infs)]
colnames(ar_estimates_by_loc)[4] <- "median"

ars <- ar_estimates_by_loc[,c("time","LOC_ID","median")]
ars <- merge(ars, locs, by="LOC_ID")         


# Population density data
pop <- read.asciigrid(paste0(fluscape_wd,"data/landscan/prd_sim_ascii_n_44384655.txt" ))
x1 = raster(paste0(fluscape_wd,"data/landscan/fluscape_wide.gri"))
margin <- 0.05 # how much to extend the density area by?
minlocx <- min(locs$LOC_Long)
maxlocx <- max(locs$LOC_Long)
minlocy <- min(locs$LOC_Lat)
maxlocy <- max(locs$LOC_Lat)
ext <- extent(minlocx-margin,maxlocx+margin,minlocy-margin,maxlocy+margin) 
x2 <- crop(x1,ext)

raster_for_plot <- as(x2, "SpatialPixelsDataFrame")
raster_for_plot <- as.data.frame(raster_for_plot)
colnames(raster_for_plot) <- c("pop","LOC_Long","LOC_Lat")
raster_for_plot[raster_for_plot$pop == 0,"pop"] <- 0.001
head(ars)
xbreaks <- seq(1968, 2015, by=1)
xbreaks <- rep(xbreaks, each=4)
xlabels <- paste0(c("Q1-","Q2-","Q3-","Q4-"),xbreaks)
ars$time1 <- xlabels[ars$time]
times <- unique(ars$time)
use_times <- times[seq(1,length(times),by=33)]
ars$xlab <- xlabels[ars$time]
ars$time1 <- factor(ars$time1, levels=unique(ars$time1))

p_static <- ggplot(ars[ars$time %in% use_times,]) + 
  geom_tile(data=raster_for_plot,aes(x=LOC_Long,y=LOC_Lat,fill=log10(pop)),alpha=0.5)+
  geom_point(aes(y=LOC_Lat,x=LOC_Long,size=median, col=median)) +
 # transition_time(time) + 
  scale_color_gradient2(low="blue",mid="yellow",high="red",midpoint=0.5,na.value="gray40",limits=c(0,1),
                        breaks=c(0,0.5,1),
                        guide=guide_colorbar(title="AR (posterior median)",
                                             direction="horizontal",title.position="top",
                                             label.position="bottom")) +
  scale_fill_gradient2(low="grey",mid="cornsilk1",high="darkgreen",na.value="gray40",midpoint=2,limits=c(0,6),
                       breaks=seq(0,12,by=2),
                       guide=guide_colorbar(title="Population density (log10, km^-1)",
                                            direction="horizontal",title.position="top",
                                            label.position="bottom")) +
  scale_size_continuous(limits=c(0,1),range=c(2,20),guide="none")+
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))+
  theme_bw() + 
  facet_wrap(~time1,ncol=2) +
  ylab("Latitude") +
  xlab("Longitude") +
  theme(legend.position="bottom",
        axis.title.x=element_text(size=12,family="sans"),
        axis.title.y=element_text(size=12,family="sans"),
        axis.text=element_blank(),axis.ticks = element_blank(),
        legend.text=element_text(size=8,family="sans"),
        legend.title=element_text(size=10,family="sans",hjust=0),
        legend.justification = "center",
        legend.key.width=unit(1.2,"cm"),
        legend.box.background = element_rect(colour = "black",fill="white"),
        plot.title=element_text(size=14)) 
p_static
#ars$time <- as.numeric(ars$time)
#ars$time1 <- factor(ars$time1,levels=unique(ars$time1[ars$time]))
#x_labels <- as.character(ars$time1)

p <- ggplot(ars) + 
  geom_tile(data=raster_for_plot,aes(x=LOC_Long,y=LOC_Lat,fill=log10(pop)),alpha=0.5)+ 
  geom_point(aes(y=LOC_Lat,x=LOC_Long,size=median, col=median)) +
  scale_color_gradient2(low="blue",mid="yellow",high="red",midpoint=0.5,na.value="gray40",limits=c(0,1),
                        breaks=c(0,0.5,1),
                        guide=guide_colorbar(title="AR (posterior median)",
                                                           direction="horizontal",title.position="top",
                                                           label.position="bottom")) +
  scale_fill_gradient2(low="grey",mid="cornsilk1",high="darkgreen",na.value="gray40",midpoint=2,limits=c(0,6),
                       breaks=seq(0,12,by=2),
                       guide=guide_colorbar(title="Population density (log10, km^-1)",
                                            direction="horizontal",title.position="top",
                                            label.position="bottom")) +
  scale_size_continuous(limits=c(0,1),range=c(2,20),guide="none")+
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))+
  theme_bw() + 
  ylab("Latitude") +
  xlab("Longitude") +
  theme(legend.position=c(0.2,0.8),
        axis.title.x=element_text(size=12,family="sans"),
        axis.title.y=element_text(size=12,family="sans"),
        axis.text=element_blank(),axis.ticks = element_blank(),
        legend.text=element_text(size=8,family="sans"),
        legend.title=element_text(size=10,family="sans",hjust=0),
        legend.justification = "center",
        legend.key.width=unit(1.2,"cm"),
        legend.box.background = element_rect(colour = "black",fill="white"),
        plot.title=element_text(size=14)) +
    transition_manual(time1) +
  labs(title="Time: {current_frame}") 
#p
library(gifski)
library(gganimate)
anim_p <- gganimate::animate(p,nframes=length(unique(ars$time)),duration=60,height=500,width=600,
                             renderer=gifski_renderer())#length(unique(ars$time)))#,height = 8, width =10,units="in",res=300)
#animate(p,duration=30,nframes=length(unique(ars$time))*2)
anim_save(paste0(figure_wd,"/ar_time_space_meas.gif"), anim_p)
ggsave_jah(p_static, figure_wd, "static_snapshot",6,7)
# 
# ars <- dcast(ars, LOC_ID~time)
# ars <- ars[,2:ncol(ars)]
# fit <- Sncf(x=locs$LOC_Long,y=locs$LOC_Lat,z=ars,na.rm=TRUE,latlon=TRUE)
# plot(fit)
# summary(fit)
# use_ars <- seq_along(ars)[seq(1,ncol(ars),by=25)]
# par(mfrow=c(ceiling(length(use_ars)/2),mfcol=2))
# for(i in seq_along(use_ars)){
#   spatial.plot(x=locs$LOC_Long,y=locs$LOC_Lat,z=ars[,i],main=strain_isolation_times[use_ars[i]]/4)
# }
# par(mfrow=1,mfcol=1)
# 
