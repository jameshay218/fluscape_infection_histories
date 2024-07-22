
#################################
## Plot some example titre fits
plot_infection_histories_long_mod <- function(chain, infection_histories, titre_dat,
                                              individuals, antigenic_map=NULL, 
                                              strain_isolation_times=NULL, par_tab,
                                              nsamp = 1000,
                                              mu_indices = NULL,expand_titre_dat=FALSE,
                                              measurement_indices_by_time = NULL,
                                              time_key=NULL,virus_key=NULL) {
    individuals <- individuals[order(individuals)]
    ## Generate titre predictions
    titre_preds <- get_titre_predictions(
        chain, infection_histories, titre_dat, individuals,
        antigenic_map, strain_isolation_times, 
        par_tab, nsamp, FALSE, mu_indices,
        measurement_indices_by_time,
        expand_titredat=expand_titre_dat
    )
    
    ## Use these titre predictions and summary statistics on infection histories
    to_use <- titre_preds$predicted_observations
    model_preds <- titre_preds$predictions
    to_use$individual <- individuals[to_use$individual]
    
    inf_hist_densities <- titre_preds$histories
    inf_hist_densities$xmin <- inf_hist_densities$variable-0.5
    inf_hist_densities$xmax <- inf_hist_densities$variable+0.5
    
    max_titre <- max(titre_dat$titre)
    min_titre <- min(titre_dat$titre)
    
    max_x <- max(inf_hist_densities$variable) + 5
    time_range <- range(inf_hist_densities$variable)
    
    if(!is.null(time_key)){
        to_use$samples_label <- time_key[as.character(to_use$samples)]
        to_use$samples_label <- factor(to_use$samples_label, levels=time_key)
        model_preds$samples_label <- time_key[as.character(model_preds$samples)]
        model_preds$samples_label <- factor(model_preds$samples_label, levels=time_key)
        titre_dat$samples_label <- time_key[as.character(titre_dat$samples)]
        titre_dat$samples_label <- factor(titre_dat$samples_label, levels=time_key)
    } else{
        to_use$samples_label <- to_use$samples
        model_preds$samples_label <- model_preds$samples
        titre_dat$samples_label <- titre_dat$samples
    }
    if(!is.null(virus_key)){
        to_use$virus_label <- virus_key[as.character(to_use$virus)]
        to_use$virus_label <- factor(to_use$virus_label, levels=virus_key)
        model_preds$virus_label <- time_key[as.character(model_preds$virus)]
        model_preds$virus_label <- factor(model_preds$virus_label, levels=virus_key)
        titre_dat$virus_label <- time_key[as.character(titre_dat$virus)]
        titre_dat$virus_label <- factor(titre_dat$virus_label, levels=virus_key)
    } else {
        to_use$virus_label <- to_use$virus
        model_preds$virus_label <- model_preds$virus
        titre_dat$virus_label <- titre_dat$virus
    }
    
    inf_hist_densities <- full_join(titre_dat %>% filter(individual %in% individuals) %>% 
                                        select(individual, samples_label, samples) %>% distinct(),inf_hist_densities)%>% 
      filter(variable <= samples)
    
    all_plot_data <- list(to_use,inf_hist_densities, model_preds[model_preds$individual %in% individuals,],
                          titre_dat[titre_dat$individual %in% individuals,])
    
    titre_pred_p <- ggplot(to_use) +
        geom_rect(data=inf_hist_densities,
                  aes(xmin=xmin,xmax=xmax,fill=value),ymin=min_titre-1,ymax=max_titre+2)+
        geom_ribbon(aes(x=virus,ymin=lower, ymax=upper),alpha=0.4, fill="#009E73",size=0.2)+
        geom_ribbon(data=model_preds[model_preds$individual %in% individuals,], 
                    aes(x=virus,ymin=lower,ymax=upper),alpha=0.7,fill="#009E73",size=0.2) + 
        geom_line(data=model_preds, aes(x=virus, y=median),linetype="dotted",color="grey10")+
        geom_rect(ymin=max_titre,ymax=max_titre+2,xmin=0,xmax=max_x,fill="grey70")+
        geom_rect(ymin=min_titre-2,ymax=min_titre,xmin=0,xmax=max_x,fill="grey70")+
        scale_x_continuous(expand=c(0,0),breaks=as.numeric(names(virus_key)),labels=virus_key) +
        scale_fill_gradient(low="white",high="#D55E00",limits=c(0,1),name="Posterior probability of infection")+
        guides(fill=guide_colourbar(title.position="top",title.hjust=0.5,label.position = "bottom",
                                    barwidth=10,barheight = 0.5, frame.colour="black",ticks=FALSE)) +
        geom_point(data=titre_dat[titre_dat$individual %in% individuals,], aes(x=virus, y=titre),shape=23, 
                   col="black",size=1)+
        ylab("log HI titre") +
        xlab("Strain") +
        theme_pubr()+
        theme(legend.title=element_text(size=7),
              legend.text=element_text(size=7),
              legend.margin = margin(-1,-1,-3,-1),
              axis.title=element_text(size=10),
              strip.text=element_text(size=7),
              axis.text.x=element_text(angle=45,hjust=1,size=5),
              axis.text.y=element_text(size=7),
              plot.margin=margin(r=15,t=5,l=5))+
        coord_cartesian(ylim=c(min_titre,max_titre+1),xlim=time_range) +
        scale_y_continuous(breaks=seq(min_titre,max_titre+2,by=2)) +
        facet_grid(individual~samples_label)
    list(titre_pred_p,all_plot_data)
}
## Create key to make better sample time labels
time_key <- strain_isolation_times#titre_dat %>% select(samples) %>% distinct() %>% arrange(samples) %>% mutate(samples=samples) %>% pull(samples)
time_key <- convert_years_to_quarters(time_key)
names(time_key) <- strain_isolation_times
## Create key to make better virus labels
virus_key <- seq(1970, 2015,by=5)
names(virus_key) <- virus_key

use_indivs <- c(10, 15, 17, 32, 56)
DOBs <- titre_dat %>% filter(individual %in% use_indivs) %>% select(individual,DOB) %>% distinct()

p_titre_fits <- plot_infection_histories_long_mod(chain=theta_chain,infection_histories=inf_chain,
                                                  titre_dat=titre_dat,individuals=use_indivs,
                                                  antigenic_map=antigenic_map,
                                                  strain_isolation_times = strain_isolation_times,
                                                  mu_indices=rep(1:48,each=4),
                                                  measurement_indices_by_time = rep(1:48,each=4),
                                                  par_tab=par_tab,time_key=time_key,virus_key=virus_key)


plot_data <- p_titre_fits[[2]]
p_titre_fits <- p_titre_fits[[1]]

## Extract data for plot
plot_data1 <- plot_data[[1]] %>% select(individual, virus, titre,lower,run, median,upper,samples,samples_label) %>% rename(lower_obs=lower,upper_obs=upper,median_obs=median) 
plot_data2 <-plot_data[[2]]%>% mutate(infection_time=variable)
plot_data3 <- plot_data[[3]] %>% select(individual, virus, lower, median,upper,samples) %>% distinct()
plot_data4 <- left_join(plot_data1,plot_data3, by=c("individual","virus","samples")) %>% select(-c(samples)) %>%
  rename(`Lower 95% prediction interval`=lower_obs, `Posterior median observation`=median_obs,
         `Upper 95% prediction interval`=upper_obs, `Lower 95% CrI`=lower, `Upper 95% CrI`=upper, `Posterior median`=median,`Virus`=virus,`Sample time`=`samples_label`,`Repeat number`=run)
plot_data2 <- plot_data2 %>% rename(`Sample time`=`samples_label`) %>% select(-samples) %>% rename(`Posterior probability of infection`=value)

write.csv(plot_data4,"~/Documents/GitHub/fluscape_infection_histories//data/figure_data/FigS5Ai.csv",row.names=FALSE)
write.csv(plot_data2,"~/Documents/GitHub/fluscape_infection_histories//data/figure_data/FigS5Aii.csv",row.names=FALSE)

colnames(antigenic_map)[3] <- "inf_times"

p_titre_fits$data %>% select(individual, samples_label, virus_label, titre, lower,lower_50,median,upper_50,upper)

colnames(antigenic_map)[3] <- "inf_times"
p_inf_hists <- generate_cumulative_inf_plots(inf_chain, 0, use_indivs, ages = DOBs,strain_isolation_times = antigenic_map$inf_times)
p_cumu_infhist <- p_inf_hists[[1]] + scale_x_continuous(breaks=seq(1970,2015,by=5)*buckets,labels=seq(1970,2015,by=5)) + 
    scale_fill_manual(name="",values=c(`1`="#D55E00")) + 
    scale_color_manual(name="",values=c(`1`="#D55E00"))+
    theme_pubr()+
    theme(legend.title=element_text(size=7),
          legend.text=element_text(size=7),
          legend.position="none",
          legend.margin = margin(-1,-1,-3,-1),
          axis.title=element_text(size=10),
          strip.text=element_text(size=7),
          axis.text.x=element_text(angle=45,hjust=1,size=7),
          axis.text.y=element_text(size=7),
          plot.tag = element_text(face="bold"),
          plot.margin=margin(r=15,t=5,l=5)) +labs(tag="B")
p_titre_fits

p_cumu_infhist_data <- p_cumu_infhist$data %>% dplyr::select(-chain_no) %>% dplyr::mutate(`Circulation time`=(as.numeric(as.character(variable)))) %>% rename(`Lower 95% CrI`=lower,
                                                                                                                                                              `Upper 95% CrI`=upper,
                                                                                                                                                              `Posterior median cumulative infections`=`median`)
write.csv(p_cumu_infhist_data,"~/Documents/GitHub/fluscape_infection_histories//data/figure_data/FigS5B.csv",row.names=FALSE)

## Find RMSE
if(FALSE){
model_preds <- get_titre_predictions(chain=theta_chain,infection_histories=inf_chain,
                                     titre_dat=titre_dat,individuals=unique(titre_dat$individual),
                                     antigenic_map=antigenic_map,
                                     strain_isolation_times = strain_isolation_times,
                                     add_residuals=TRUE,
                                     for_res_plot=FALSE,
                                     expand_titredat=FALSE,
                                     mu_indices=rep(1:48,each=4),
                                     measurement_indices_by_time = rep(1:48,each=4),
                                     par_tab=par_tab)
residuals <- model_preds$residuals$`50%`
RMSE <- sqrt(sum((residuals^2))/length(residuals))
print(paste0("RMSE: ",signif(RMSE,3)))


vietnam_dob <- read.csv(paste0(main_wd,"data/vietnam_ages.csv"))%>% mutate(individual = 1:n()) %>% 
    left_join(titre_dat %>% select(-DOB))
n_alive <- matrix(get_n_alive(vietnam_dob, antigenic_map$inf_times),nrow=1)
colnames(n_alive) <- antigenic_map$inf_times
n_alive <- as.data.frame(n_alive)
n_alive$group <- 1

p_vietnam_ar <- plot_attack_rates(inf_chain %>% mutate(group=1),titre_dat %>% mutate(group=1),
                                      antigenic_map$inf_times,pad_chain=TRUE,n_alive=n_alive)

## Plot estimated attack rate vs. real
p_vietnam_ar <- p_vietnam_ar +  ylab("Per capita incidence per year") +
    scale_x_continuous(breaks=seq(1970,2015,by=5)*4,labels=seq(1970,2015,by=5)) + theme_pubr()+
    theme(legend.title=element_text(size=7),
          legend.text=element_text(size=7),
          legend.position="none",
          legend.margin = margin(-1,-1,-3,-1),
          axis.title=element_text(size=10),
          axis.text.x=element_text(angle=45,hjust=1,size=7),
          axis.text.y=element_text(size=7),
          plot.tag = element_text(face="bold"),
          strip.text=element_blank(),strip.background = element_blank(),
          plot.margin=margin(r=15,t=5,l=5)) +labs(tag="C")
}

p_titre_fit_comb <- (p_titre_fits+labs(tag="A") + theme(plot.tag = element_text(face="bold"))) + 
    (p_cumu_infhist + facet_wrap(~individual,nrow=1)) + 
    plot_layout(ncol=1,heights=c(4,1))
ggsave_jah(p_titre_fit_comb,figure_wd,"vietnam_titre_fits_example",width=7,height=8)
