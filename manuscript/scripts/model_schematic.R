library(tidyverse)
library(serosolver)
library(patchwork)
setwd("~/Documents/GitHub/fluscape_serosolver/")
main_theme <- theme_classic
## Create model schematic
## 1. Example attack rates
set.seed(1)
ts <- 1:10
ARs <- data.frame(t=ts, y=simulate_attack_rates(ts))
p_ar <- ggplot(ARs) + geom_line(aes(x=t,y=y)) + 
    scale_x_continuous(breaks=ts) +
    ylab("Probability of infection") +
    labs(tag="A") +
    main_theme() +
    xlab("Time") +
    scale_y_continuous(expand=c(0,0),limits=c(0,0.35),breaks=seq(0,0.35,by=0.1))

## Simulate example infection history
inf_hist <- rep(0, 10)
inf_hist[c(2,8)] <- 1

## Simulate kinetics
par_tab <- read.csv("inputs/par_tab_base.csv")
theta <- par_tab$values
names(theta) <- par_tab$names
theta["wane"] <- 0.5
theta["sigma2"] <- 0.01
antigenic_map <- read.csv("data/antigenic_maps/antigenic_map_fonville_annual_continuous.csv")
antigenic_map$inf_times <- antigenic_map$inf_times - 1967
antigenic_map <- antigenic_map %>% dplyr::filter(inf_times %in% (ts*4))
antigenic_map$inf_times <- ts
titers <- as.data.frame(simulate_individual(theta, inf_hist, antigenic_map, 1:10,1:10,c(2,8),add_noise=FALSE))
colnames(titers) <- c("samples","virus","titer","repeat")
titers <- titers %>% mutate(virus = if_else(virus == 2, "B","H"))
titers$virus <- as.factor(titers$virus)

p_titer1 <- ggplot(titers %>% filter(virus == "B")) + 
    geom_line(aes(x=samples,y=titer,group=virus,col=virus)) +
    labs(tag="D") + 
    geom_hline(yintercept=0,linetype="dotted",col="grey60") +
    ylab('HI titer') +
    scale_x_continuous(breaks=ts) +
    scale_y_continuous(limits=c(0,5),breaks=seq(0,5,by=1)) +
    scale_color_manual(name="Strain",values=c("B"="blue","H"="orange")) +
    main_theme() +
    theme(axis.text.x=element_blank(), axis.title.x=element_blank(),axis.ticks.x=element_blank(),
          legend.position=c(0.9,0.9))

p_titer2 <- ggplot(titers) + geom_line(aes(x=samples,y=titer,group=virus,col=virus,
                                           alpha=virus,linetype=virus)) +
    scale_linetype_manual(name="Strain",values=c("B"="dashed","H"="solid")) +
    scale_alpha_manual(name="Strain",values=c("B"=0.5,"H"=1))+
    geom_hline(yintercept=0,linetype="dotted",col="grey60") +
    ylab('HI titer') +
    labs(tag="E") +
    scale_color_manual(name="Strain",values=c("B"="blue","H"="orange")) +
    scale_x_continuous(breaks=ts) +
    scale_y_continuous(limits=c(0,5),breaks=seq(0,5,by=1)) +
    xlab("Time") +
    main_theme() +
    theme(legend.position=c(0.9,0.9))


antigenic_map$inf_times <- LETTERS[antigenic_map$inf_times]
p_map <- ggplot(antigenic_map) + 
    geom_point(aes(x=x_coord,y=y_coord,col=inf_times),size=2) +
    scale_color_manual(name="Strain",values=c("A"="grey40", "B"="blue", "C"="grey40", "D"="grey40", 
                                               "E"="grey40", "F"="grey40", "G"="grey40", "H"="orange", "I"="grey40", "J"="grey40")) +
    main_theme() +
    xlab("Antigenic dimension 1") + ylab("Antigenic dimension 2") +
    scale_x_continuous(breaks=seq(-15,15,by=5)) +
    scale_y_continuous(breaks=seq(-4,6,by=2)) +
    labs(tag="B") +
    theme(legend.position="none")

p_serum1 <- ggplot(titers %>% filter(samples == 6)) + 
    geom_bar(aes(x=virus,y=titer,group=virus,fill=virus),stat="identity",col="black") +
    geom_hline(yintercept=0,linetype="dotted",col="grey60") +
    ylab('HI titer') +
    labs(tag="F") +
    scale_fill_manual(name="Strain",values=c("B"="blue","H"="orange")) +
    scale_y_continuous(limits=c(0,5),breaks=seq(0,5,by=1)) +
    main_theme() +
    xlab("Strain") +
    theme(legend.position="none")
p_serum2 <- ggplot(titers %>% filter(samples == 9)) + 
    geom_bar(aes(x=virus,y=titer,group=virus,fill=virus),stat="identity",col="black") +
    geom_hline(yintercept=0,linetype="dotted",col="grey60") +
    ylab('HI titer') +
    labs(tag="G") +
    scale_fill_manual(name="Strain",values=c("B"="blue","H"="orange")) +
    scale_y_continuous(limits=c(0,5),breaks=seq(0,5,by=1)) +
    main_theme() +
    xlab("Strain") +
    theme(legend.position="none")

p_main <- ({{p_ar + p_titer1 + p_titer2 + plot_layout(guides="keep", ncol=1, heights=c(1,2,2))}  | 
    {p_map / {(p_serum1 | p_serum2)}}} + plot_layout(widths=c(1.5,1)))

jahR::save_plots(p_main,"figures/model/","schematic",10,8)

