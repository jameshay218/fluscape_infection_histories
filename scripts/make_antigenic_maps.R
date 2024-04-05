######################################################
## ANTIGENIC MAPS
## Author: James Hay
## Date: 12 July 2023
## Summary: plots various forms of the H3N2 antigenic map and creates versions used in model fitting

library(ggplot2)
library(cowplot)
library(ggrepel)
library(plyr)
library(RColorBrewer)
library(tidyverse)
library(patchwork)
devtools::load_all("~/Documents/GitHub/serosolver") ## library(serosolver)

dataWD <- "~/Documents/GitHub/fluscape_infection_histories//data/"
saveWD <- "~/Documents/GitHub/fluscape_infection_histories//figures/antigenic_maps/"

## Read in the various coordinates
virus_key <- c("HK68"=1968, "EN72"=1972, "VI75"=1975, "TX77"=1977, "BK79"=1979, "SI87"=1987, 
               "BE89"=1989, "BJ89"=1989,"BE92"=1992, "WU95"=1995, "SY97"=1997, "FU02"=2002, 
               "CA04"=2004, "WI05"=2005, "PE06"=2009)

fluscape_summary_coords <- read.csv(paste0(dataWD, "/antigenic_map_coords/fonville_antigenic_coordinates.csv")) %>% 
    rename(x_coord=X,y_coord=Y, year=Strain) %>%
    mutate(dataset="Fluscape strains") %>%
    arrange(year) %>%
    select(year, dataset,x_coord,y_coord)
fluscape_summary_coords$x_coord <- fluscape_summary_coords$x_coord - mean(fluscape_summary_coords$x_coord)
fluscape_summary_coords$y_coord <- fluscape_summary_coords$y_coord - mean(fluscape_summary_coords$y_coord)

vietnam_coords <- read.csv(paste0(dataWD, "/antigenic_map_coords/vietnam_antigenic_coords.csv")) %>%
    separate(viruses, into = c("prefix", "year"), sep = "/(?=[^/]+$)", extra = "drop", remove = FALSE) %>%
    mutate(year = sub("_.*", "", year)) %>%
    mutate(year = as.numeric(year)) %>%
    mutate(year = if_else(year < 2000 & year > 50, year + 1900, year)) %>%
    mutate(year = if_else(year < 50, year + 2000, year)) %>%
    rename(x_coord=AG_y, y_coord=AG_x)%>%
    mutate(dataset="Fonville et al. 2014")%>%
    arrange(year) %>%
    select(year, dataset,x_coord,y_coord)

vietnam_coords$x_coord <- vietnam_coords$x_coord - mean(vietnam_coords$x_coord)
vietnam_coords$y_coord <- vietnam_coords$y_coord - mean(vietnam_coords$y_coord)

bedford_coords <- read_tsv(paste0(dataWD, "/antigenic_map_coords/elife-01914-fig2-data1-v1.tsv")) %>% 
    filter(lineage == "H3N2") %>%
    rename(x_coord=ag1, y_coord=ag2) %>%
    mutate(dataset="Bedford et al. 2014")%>%
    arrange(year) %>%
    select(year, dataset,x_coord,y_coord)

bedford_coords$x_coord <- bedford_coords$x_coord - mean(bedford_coords$x_coord)
bedford_coords$y_coord <- bedford_coords$y_coord - mean(bedford_coords$y_coord)

## Combine datasets and re-center to mean in each antigenic dimension
comb_coords <- bind_rows(fluscape_summary_coords, vietnam_coords, bedford_coords) %>% 
    mutate(year=as.factor(year)) %>%
    group_by(dataset) 

summary_coords <- comb_coords %>% group_by(dataset,year) %>% summarize(x_mean=mean(x_coord),y_mean=mean(y_coord))

## Look at raw data
p_raw <- ggplot(comb_coords) +
    geom_point(aes(x=x_coord,y=y_coord,col=year)) + 
    geom_line(data=summary_coords,aes(x=x_mean,y=y_mean),col="black",linewidth=0.5) +
    geom_point(data=summary_coords,aes(x=x_mean,y=y_mean,fill=year),col="black",size=3,pch=21) +
    facet_wrap(~dataset,ncol=1) +
    theme_minimal() +
    xlab("Antigenic dimension 1") + ylab("Antigenic dimension 2")
jahR::save_plots(p_raw, saveWD,"raw_data",width=10,height=10)    

tmp <- summary_coords %>% group_by(dataset) %>% arrange(year) %>% 
    mutate(year = as.numeric(year)) %>%
    do(as.data.frame(as.matrix(dist(.[,c("x_mean","y_mean")],method="euclidean")))) %>% 
    mutate(index = 1:n()) %>% pivot_longer(cols=-c("dataset","index")) %>%
    mutate(name = as.numeric(name))

year_key <- comb_coords %>% select(dataset, year) %>% distinct() %>% group_by(dataset) %>% arrange(year) %>% mutate(index = 1:n())

tmp <- tmp %>% drop_na() %>% left_join(year_key) %>% left_join(year_key %>% rename(year_name = year, name = index))

p_pairwise <- ggplot(tmp) + 
    geom_tile(aes(x=year,y=year_name,fill=value)) +
    facet_wrap(~dataset,ncol=1,scales="free") +
    scale_fill_viridis_c(name="Euclidean distance") +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=45,hjust=1,size=6),
          axis.text.y=element_text(size=6)) +
    xlab("Virus isolation year") + ylab("Virus isolation year") 
jahR::save_plots(p_pairwise, saveWD,"pairwise_distances",width=8,height=10)    


summary_coords <- summary_coords %>% mutate(year = as.numeric(as.character(year))) %>% 
    group_by(dataset) %>% mutate(dist = sqrt((x_mean - lag(x_mean,1))^2 + (y_mean - lag(y_mean,1))^2)) %>%
    mutate(dist = if_else(is.na(dist),0,dist)) %>%
    mutate(cumu_dist = cumsum(dist)) %>%
    mutate(cumu_dist = cumu_dist/max(cumu_dist))

p_cumu_dist <- ggplot(summary_coords) + geom_line(aes(x=year,y=cumu_dist,col=dataset)) +
    xlab("Virus isolation year") +
    ylab("Relative antigenic distance from root") +
    scale_color_manual(name="Dataset",values=c("Bedford et al. 2014"="blue","Fonville et al. 2014"="red","Fluscape strains"="darkgreen")) +
    theme_minimal() + 
    theme(legend.position="bottom",
          panel.grid.major=element_line(linewidth=0.1,color="grey50"),
          panel.grid.minor=element_line(linewidth=0.1,color="grey50")) +
    #scale_y_continuous(breaks=seq(-10,10,by=2)) + scale_x_continuous(breaks=seq(-50,50,by=5)) +
    ylab("Cumulative antigenic distance") + xlab("Virus isolation year")
jahR::save_plots(p_cumu_dist, saveWD,"cumulative_distance",width=8,height=5)    

## Function to plot map and coordinates
plot_fit <- function(fit_dat=NULL, coord_data=NULL,buckets, clusters=NULL){
    if(!is.null(clusters)){
        clusters1 <- clusters[rep(seq_len(nrow(clusters)),each=buckets),]
        clusters1$year <- seq(1968*buckets, length.out=nrow(clusters1))
        clusters1$first_year <- clusters1$first_year*buckets
        clusters1 <- clusters1[,c("cluster_used","year")]
        colnames(clusters1)[2] <- "inf_times"
        fit_dat <- merge(fit_dat, clusters1, by="inf_times")
        if(!is.null(coord_data)) coord_data <- merge(coord_data, clusters, by="inf_times")
    } else {
        fit_dat$cluster_used <- 1
        if(!is.null(coord_data)) coord_data$cluster_used <- 1
    }
    p <- ggplot(fit_dat) + 
        geom_point(aes(x=x_coord,y=y_coord)) + 
        geom_line(aes(x=x_coord,y=y_coord))
    if(!is.null(coord_data)){
        p <- p +
            geom_point(data=coord_data,aes(x=X,y=Y, col=as.factor(cluster_used)),size=0.8)
    }
    yearMin <- min(fit_dat$inf_times)
    yearMax <- max(fit_dat$inf_times)
    yearLabels <- seq(yearMin, yearMax, by=buckets)
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))
    colourCount = length(unique(fit_dat$cluster_used))
    p <- p + geom_label_repel(data=fit_dat[fit_dat$inf_times %in% yearLabels, ],
                              aes(x=x_coord,y=y_coord,label=inf_times/buckets,col=as.factor(cluster_used)),
                              label.size=NA,size=4,fill=NA) +
        scale_colour_manual(values=getPalette(colourCount))+
        theme(legend.position="none")
    p
}

## Function to create duplicate rows within a year if desired. 
## Use this to assume no antigenic change within a year (ie. same for 4 consecutive quarters)
replicate_years_f <- function(fit_dat,replicate_buckets){
    fit_dat <- fit_dat[rep(row.names(fit_dat),each=replicate_buckets),]
    yearMin <- min(fit_dat$inf_times)*replicate_buckets
    yearMax <- (max(fit_dat$inf_times)+1)*replicate_buckets
    years <- seq(yearMin, yearMax-1,by=1)
    fit_dat$inf_times <- years
    fit_dat
    
}

# FOR QUARTERLY CLUSTERED WITHIN YEAR
# buckets <- 1; replicate_buckets <- 4; use_clusters <- FALSE; replicate_years <- TRUE

# FOR QUARTERLY CLUSTERED BY CLUSTER
# buckets <- 4; replicate_buckets <- 4; use_clusters <- TRUE; replicate_years <- TRUE

# FOR QUARTERLY CONTINUOUS
# buckets <- 4; replicate_buckets <- 4; use_clusters <- FALSE; replicate_years <- FALSE

# FOR ANNUAL CONTINUOUS
# buckets <- 1; replicate_buckets <- 1; use_clusters <- FALSE; replicate_years <- FALSE

# FOR ANNUAL CLUSTERED
# buckets <- 1; replicate_buckets <- 1; use_clusters <- TRUE, replicate_years <- TRUE

spar_bed <- 0.8 ## Smoothing parameter for the Bedford antigenic map
spar_fon <- 0.3 ## Smoothing parameter for the Fluscape strains map
spar_viet <- 0.8 ## Smoothing parameter for the Fluscape strains map

## Generate maps for each set of assumptions
combs <- list(
    c(buckets=1,replicate_buckets=4, use_clusters=FALSE,replicate_years=TRUE),
    c(buckets=4,replicate_buckets=4, use_clusters=TRUE,replicate_years=TRUE),
    c(buckets=4,replicate_buckets=4, use_clusters=FALSE,replicate_years=FALSE),
    c(buckets=1,replicate_buckets=1, use_clusters=FALSE,replicate_years=FALSE),
    c(buckets=1,replicate_buckets=1, use_clusters=TRUE,replicate_years=TRUE)
    
)
for(index in seq_along(combs)){
    comb <- combs[[index]]
    print(comb)
    buckets <- comb[1]
    replicate_buckets <- comb[2]
    use_clusters <- comb[3]
    replicate_years <- comb[4]
    if(buckets == 1 & replicate_buckets == 4){
        filename <- "quarterly_annual_clusters"
    } else if(buckets == 4 & replicate_buckets == 4 & use_clusters == TRUE){
        filename <- "quarterly_clustered"
    } else if(buckets == 4 & replicate_buckets == 4 & use_clusters == FALSE){
        filename <- "quarterly_continuous"
    } else if(buckets == 1 & replicate_buckets == 1 & use_clusters == FALSE){
        filename <- "annual_continuous"
    } else {
        filename <- "annual_clustered"
    }
    
    if(use_clusters){
        filename <- paste0(filename,"_clustered")
    }
    print(filename)
    save_maps <- FALSE
    
    ## Clusters assigned based on Xiangjun et al. Nat Comms 2012 
    clusters <- read.csv(paste0(dataWD, "antigenic_map_coords/xiangjun2012_clusters.csv"))
    
    ## Antigenic maps
    ## Aim is to generate coordinates for all viruses that will be used, either annual or quarterly
    ## Coordinates of the Fluscape strains as per Fonville et al. Science 2014
    
    ## Create the smooth antigenic maps from both sets of coordinates
    colnames(bedford_coords) <- c("Strain","Dataset","X","Y")
    fit_dat_bedford <- generate_antigenic_map_flexible(bedford_coords,buckets,clusters,use_clusters,spar_bed)
    colnames(fluscape_summary_coords) <- c("Strain","Dataset","X","Y")
    fit_dat_fluscape <- generate_antigenic_map_flexible(fluscape_summary_coords,buckets,clusters,use_clusters,spar_fon)
    colnames(vietnam_coords) <- c("Strain","Dataset","X","Y")
    fit_dat_vietnam <- generate_antigenic_map_flexible(vietnam_coords,buckets,clusters,use_clusters,spar_viet)
    
    ## If specified, duplicate rows within a year
    if(replicate_years & !use_clusters){
      fit_dat_bedford <- replicate_years_f(fit_dat_bedford,replicate_buckets)
      fit_dat_fluscape <- replicate_years_f(fit_dat_fluscape,replicate_buckets)
      fit_dat_vietnam <- replicate_years_f(fit_dat_vietnam,replicate_buckets)
    }
    
    if(!is.null(clusters)) colnames(clusters)[2] <- "inf_times"
    colnames(bedford_coords)[1] <- colnames(fluscape_summary_coords)[1] <- colnames(vietnam_coords)[1] <- "inf_times"
    
    ## Plot maps
    p1 <- plot_fit(fit_dat_bedford,bedford_coords,replicate_buckets, clusters) + 
        theme_minimal() + theme(legend.position="none") + ggtitle("Bedford et al. 2014") +
        ylab("Antigenic dimension 1") + xlab("Antigenic dimension 2")
    p2 <- plot_fit(fit_dat_fluscape,fluscape_summary_coords,replicate_buckets, clusters)+ 
        theme_minimal() + theme(legend.position="none") + ggtitle("Fluscape strains")+
        ylab("Antigenic dimension 1") + xlab("Antigenic dimension 2")
    p3 <- plot_fit(fit_dat_vietnam,vietnam_coords,replicate_buckets, clusters)+ 
        theme_minimal() + theme(legend.position="none") + ggtitle("Fonville et al. 2014")+
        ylab("Antigenic dimension 1") + xlab("Antigenic dimension 2")
    p_main1 <- p1/p2/p3
    
    
    if(save_maps){
        write.table(fit_dat_bedford, paste0(dataWD,"/antigenic_maps/antigenic_map_bedford_",filename,".csv"),sep=",",row.names=FALSE)
        write.table(fit_dat_fluscape, paste0(dataWD,"/antigenic_maps/antigenic_map_fonville_",filename,".csv"),sep=",",row.names=FALSE)
        write.table(fit_dat_vietnam, paste0(dataWD,"/antigenic_maps/antigenic_map_vietnam_",filename,".csv"),sep=",",row.names=FALSE)
        jahR::save_plots(p_main1, saveWD, paste0(filename,"_maps"),width=8,height=10)
    }
    
    comb_fits <- bind_rows(fit_dat_bedford %>% mutate(Dataset="Bedford et al. 2014"),
                           fit_dat_vietnam %>% mutate(Dataset = "Fonville et al. 2014"),
                           fit_dat_fluscape %>% mutate(Dataset="Fluscape strains"))
    
    p_compare <- ggplot(comb_fits) + 
        geom_line(aes(x=x_coord,y=y_coord,col=Dataset)) + 
        scale_color_manual(name="Dataset",values=c("Bedford et al. 2014"="blue","Fonville et al. 2014"="red","Fluscape strains"="darkgreen")) +
        theme_minimal() + 
        theme(legend.position="bottom") +
        #scale_y_continuous(breaks=seq(-10,10,by=2)) + scale_x_continuous(breaks=seq(-50,50,by=5)) +
        ylab("Antigenic dimension 1") + xlab("Antigenic dimension 2") +
        ggtitle("Comparison of antigenic map")
    
    comb_fits <- comb_fits %>% group_by(Dataset) %>%
        mutate(dist = sqrt((x_coord-lag(x_coord,1))^2 + (y_coord-lag(y_coord,1))^2))%>% 
        mutate(dist=ifelse(is.na(dist),0,dist)) %>%
        mutate(cdist = cumsum(dist))
    
    p_compare_dist <- ggplot(comb_fits) + 
        geom_line(aes(x=inf_times,y=cdist,col=Dataset)) + 
        scale_color_manual(name="Dataset",values=c("Bedford et al. 2014"="blue","Fonville et al. 2014"="red","Fluscape strains"="darkgreen")) +
        theme_minimal() + 
        theme(legend.position="bottom") +
        ylab("Antigenic distance from root") + xlab("Time") +
        ggtitle("Cumulative antigenic distance")
    p_main2 <- p_compare/p_compare_dist
    if(filename=="quarterly_clustered_clustered"){
        comb_fits_punctuated <- bind_rows(fit_dat_bedford %>% mutate(Dataset="Bedford et al. 2014"),
                                          fit_dat_fluscape %>% mutate(Dataset="Fluscape strains"))
        
        comb_fits_punctuated2 <- comb_fits %>% group_by(Dataset) %>%
            mutate(dist = sqrt((x_coord-lag(x_coord,1))^2 + (y_coord-lag(y_coord,1))^2))%>% 
            mutate(dist=ifelse(is.na(dist),0,dist)) %>%
            mutate(cdist = cumsum(dist))
        
        pA <- p2
        #pB <- p1
    }
    
    if(filename=="quarterly_continuous"){
        comb_fits_continuous <- bind_rows(fit_dat_bedford %>% mutate(Dataset="Bedford et al. 2014"),
                                          fit_dat_fluscape %>% mutate(Dataset="Fluscape strains"))
        
        comb_fits_continuous2 <- comb_fits %>% group_by(Dataset) %>%
            mutate(dist = sqrt((x_coord-lag(x_coord,1))^2 + (y_coord-lag(y_coord,1))^2))%>% 
            mutate(dist=ifelse(is.na(dist),0,dist)) %>%
            mutate(cdist = cumsum(dist))
        pB <- p1
        pC <- p2
    }
    if(save_maps){
        jahR::save_plots(p_main2, saveWD, paste0(filename,"_compare"),width=10,height=10)
    }
}

p_appendix_compare <- (pC + labs(tag="A") + ggtitle("Map used for simulation, smoothed") + theme(plot.tag=element_text(face="bold"))) + 
    (pB+ labs(tag="B") +ggtitle("Map using Bedford et al. 2014 coordinates")+ theme(plot.tag=element_text(face="bold"))) + 
    (pA+ labs(tag="C")+ ggtitle("Map used for simulation, clustered")+ theme(plot.tag=element_text(face="bold"))) +
    plot_layout(ncol=1)
if(save_maps){
    jahR::save_plots(p_appendix_compare, saveWD, paste0(filename,"_compare_appendix"),width=10,height=10)
}
