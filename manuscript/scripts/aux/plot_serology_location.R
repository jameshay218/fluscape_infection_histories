setwd("~/Documents/GitHub/fluscape")
#source("~/Google Drive/Drive/Influenza/serosolver/FluScape_epidemiology/scripts/extract_fluscape_dat.R")
load("~/Documents/GitHub/fluscape_infection_histories/manuscript/r_data/fluscape_dat.RData")

fluscape_dat[fluscape_dat$change > 3 & !is.na(fluscape_dat$change), "change"] <- 4
fluscape_dat[fluscape_dat$change < -3 & !is.na(fluscape_dat$change), "change"] <- -4
fluscape_dat[fluscape_dat$change == 4 & !is.na(fluscape_dat$change),"change"] <- "≥4"
fluscape_dat[fluscape_dat$change == "-4" & !is.na(fluscape_dat$change),"change"] <- "≤-4"
change_levels <- c("≤-4",as.character(seq(-3,3,by=1)),"≥4")
fluscape_dat$change <- factor(fluscape_dat$change, levels=change_levels,ordered=TRUE)
colnames(fluscape_dat)[colnames(fluscape_dat) == "change"] <- "Change in titre"

library("lattice")
library("sp")
library("raster")
#library("maptools")
library("RColorBrewer")
# read in china census data
path = "data/China census 2008/"
census.hh.size = read.csv( paste(path,"china2008.household.size.csv",sep=""), header=TRUE )
census.age.demogr = read.csv( paste(path,"china2008.demography.csv",sep=""), header=TRUE )

# read in population density data
path = "data/landscan/"
require("sp")
pop <- read.asciigrid( paste(path,"prd_sim_ascii_n_44384655.txt",sep="") )

# read in boundary & coastline map data
path = "data/Noaa maps/"
boundaries <- read.table(paste(path,"boundaries.dat",sep=""), header=FALSE, 
                         sep="", stringsAsFactors=FALSE, na.strings = c("NA",">"), fill=TRUE )
coast <- read.table(paste(path,"coast.dat",sep=""), header=FALSE, 
                    sep="", stringsAsFactors=FALSE, na.strings = c("NA",">"), fill=TRUE )
rivers <- read.table(paste(path,"rivers.dat",sep=""), header=FALSE, 
                     sep="", stringsAsFactors=FALSE, na.strings = c("NA",">"), fill=TRUE )

locs <- read.csv("~/Google Drive/My Drive/Influenza/serosolver/FluScape_epidemiology/data/loc_data_clean.csv",stringsAsFactors=FALSE)

fluscape_dat <- merge(fluscape_dat, locs[,c("LOC_ID","LOC_Lat_rnd","LOC_Long_rnd","n.samples.V1",
                                            "badge.long","badge.lat","rank.DIST_FRM_GZ","quintile")],by="LOC_ID")

## Need to plot with jittered locations
locs_use <- read.csv("~/Documents/GitHub/fluscape_serosolver/data/loc_data_rnd.csv")
locs <- locs %>% select(-c(LOC_Lat,LOC_Long)) %>% left_join(locs_use %>% select(-c(LOC_Lat_original,LOC_Long_original))) %>% rename(LOC_Lat=LOC_Lat_rnd,
                                                                                                                                    LOC_Long=LOC_Long_rnd)

fluscape_dat <- fluscape_dat %>% select(-c(LOC_Lat,LOC_Long)) %>% left_join(locs_use %>% select(-c(LOC_Lat_original,LOC_Long_original))) %>% rename(LOC_Lat=LOC_Lat_rnd,
                                                                                                                                    LOC_Long=LOC_Long_rnd)



## Order by age by group
tmp <- unique(fluscape_dat[,c("individual","DOB","quintile")])
tmp <- as.data.frame(tmp %>% arrange(quintile, desc(DOB)) %>% group_by(quintile) %>% mutate(counter=row_number(desc(DOB))))
fluscape_dat <- merge(fluscape_dat, tmp)



len.cols <- length(0:8)
cols <- colorRampPalette(c("royalblue", "khaki1", "red"))(len.cols)
age_labels <- c(0,10,20,30,40,50,60)
age_labels_text <- levels(fluscape_dat$age_group)
ps <- list()
tmp <- list()
age_group_counts_group <- list()
age_group_labels_group <- list()
n_indivs <- list()
age_label_rect_dats <- list()
age_label_ps <- list()
loc_plots <- list()

#y <- age_group_counts1
visit <- "First visit"

## Group 3
for(k in unique(fluscape_dat$quintile)){
  print(k)
  tmp <- unique(fluscape_dat[,c("individual","age_group","quintile")])
  tmp <- tmp[tmp$quintile == k,]
  n_indivs[[k]] <- nrow(tmp)
  age_group_counts_bottom <- c(0,cumsum(ddply(tmp, ~age_group, nrow)$V1))
  age_group_counts_bottom <- age_group_counts_bottom[1:(length(age_group_counts_bottom)-1)]
  age_group_counts_group[[k]] <- age_group_counts_bottom
  age_group_counts_top <- cumsum(ddply(tmp, ~age_group, nrow)$V1)
  
  age_group_rect <- data.frame(ymin=age_group_counts_bottom, ymax=age_group_counts_top)
  age_group_rect$age_group <- age_labels_text[which(unique(tmp$age_group) %in% age_labels_text)]
  age_group_rect$age_group <- factor(age_group_rect$age_group, levels=levels(fluscape_dat$age_group),ordered=TRUE)
  age_label_rect_dats[[k]] <- age_group_rect
  
  loc_plots[[k]] <- 
    ggplot(locs) + geom_point(aes(y=LOC_Lat, x=LOC_Long),col="gray70",size=1) + 
    geom_point(data=locs[locs$quintile == k,],aes(y=LOC_Lat, x=LOC_Long),col="gray5",size=1) +
    theme_void() +
    theme(panel.border=element_rect(color="gray5",fill=NA,size=1),
          plot.margin=unit(c(0.1,0.1,0.1,0.11),"cm"))
  
  
  age_label_ps[[k]] <- ggplot(age_group_rect) + 
    geom_rect(aes(xmin=0,xmax=1, ymin=ymin,ymax=ymax, fill=age_group)) + 
    scale_y_continuous(expand=c(0,0), breaks=age_group_counts_group[[k]]) +
    theme_void() + 
    scale_fill_manual(values=brewer.pal(8,name="Set2"),guide=guide_legend(title="Age group",
                                         reverse=TRUE, direction="vertical",title.position="left",
                                                      label.position="right"))+
    #scale_fill_brewer(palette="Set2") + 
    theme(legend.text=element_text(size=6,family="sans"),
          legend.title=element_text(size=8,family="sans",angle=90,hjust=0.5),
          legend.key.size=unit(0.15,"in"),
          legend.justification = "center",
          legend.key=element_rect(color="black")) 
    
    
  leg_age <- get_legend(age_label_ps[[k]])
  age_label_ps[[k]] <- age_label_ps[[k]] + theme(legend.position = "none",
                                                 plot.margin=unit(c(1,0,1,0.01),"cm"))
  
  age_group_labels_group[[k]] <- age_labels[levels(tmp$age_group) %in% unique(tmp$age_group)]

  ps[[k]] <- ggplot(fluscape_dat[fluscape_dat$visit == visit & fluscape_dat$quintile == k,]) + 
  #geom_tile(aes(x=order,y=virus,fill=`log titre`)) + 
  geom_tile(aes(y=counter,x=Full_name,fill=`log titre`)) + 
  theme_bw() + 
    scale_fill_manual(values=cols, guide=guide_legend(title="Log titre",
                                                      reverse=TRUE, direction="vertical",title.position="left",
                                                      label.position="right"))+
  scale_x_discrete(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0), breaks=age_group_counts_group[[k]], labels=age_group_labels_group[[k]]) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill="grey40"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(family="sans",size=6,colour="black",angle=90,hjust=1,vjust=0.5),
        legend.text=element_text(size=6,family="sans"),
        legend.title=element_text(size=8,family="sans",angle=90,hjust=0.5),
        legend.key.size=unit(0.15,"in"),
        legend.justification = "center",
        legend.key=element_rect(color="black")) +
  xlab("") +
  ylab("")
  leg <- get_legend(ps[[k]])
  plot(leg)
  ps[[k]] <- ps[[k]] + theme(legend.position="none",
                             plot.margin=unit(c(1,0.1,0,0),"cm"))
  
  ps[[k]] <- plot_grid(age_label_ps[[k]], ps[[k]] ,nrow=1,rel_widths=c(1,10),axis="bt",align="h")
  #ps[[k]] <- plot_grid(ps[[k]], loc_plots[[k]], nrow=2,rel_heights=c(6,1),align="hv")
}

blank_p <- ggplot(data.frame()) + geom_blank() + theme_void()
max_height <- max(unlist(n_indivs))
for(k in seq_along(ps)){
  rel_height1 <- n_indivs[[k]]
  rel_height2 <- max_height - rel_height1
  print(rel_height1)
  if(n_indivs[[k]] == min(unlist(n_indivs))){
    ps[[k]] <- plot_grid(plot_grid(leg_age,leg,nrow=1), ps[[k]], nrow=2,rel_heights=c(rel_height2, rel_height1))
  } else {
    if(rel_height1 != max_height){
      ps[[k]] <- plot_grid(blank_p, ps[[k]], nrow=2,rel_heights=c(rel_height2, rel_height1))
    }
  }
}
pdf("~/Documents/GitHub/fluscape_infection_histories/figures/serology_by_loc_v1.pdf",height=8,width=7.2)
plot_grid(plot_grid(plotlist=ps,nrow=1), plot_grid(plotlist=loc_plots,nrow=1),nrow=2,rel_heights=c(6,1))
dev.off()
png("~/Documents/GitHub/fluscape_infection_histories/figures/serology_by_loc_v1.png",height=8,width=7.2,res=300,units="in")
plot_grid(plot_grid(plotlist=ps,nrow=1), plot_grid(plotlist=loc_plots,nrow=1),nrow=2,rel_heights=c(6,1))
dev.off()
