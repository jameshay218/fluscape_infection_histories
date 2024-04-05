#################################
## Plot some example titre fits
plot_infection_histories_long_mod <- function(chain, infection_histories, titre_dat,
                                              individuals, antigenic_map=NULL, 
                                              strain_isolation_times=NULL, par_tab,
                                              nsamp = 1000,
                                              mu_indices = NULL,expand_titre_dat=FALSE,
                                              measurement_indices_by_time = NULL,
                                              time_key=NULL,virus_key=NULL,
                                              fill_color="#009E73") {
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
                                        select(individual, samples_label, samples) %>% distinct(),inf_hist_densities) %>% 
      filter(variable <= samples)
    titre_pred_p <- ggplot(to_use) +
          geom_rect(data=inf_hist_densities,
                  aes(xmin=xmin,xmax=xmax,fill=value),ymin=min_titre-1,ymax=max_titre+2)+
          geom_rect(data=to_use %>% ## Get DOB to plot pre-birth times
                      select(individual,DOB,samples_label) %>% 
                      distinct() %>% 
                      mutate(DOB_min=min(strain_isolation_times)-4) %>%
                      group_by(individual,samples_label) %>% mutate(min=pmin(DOB_min,DOB),max=pmax(DOB_min,DOB)),
                    aes(xmin=min,xmax=max,ymin=0,ymax=8),fill="purple",alpha=0.25) +
        geom_ribbon(aes(x=virus,ymin=lower, ymax=upper),alpha=0.4, fill=fill_color,size=0.2)+
        geom_ribbon(data=model_preds[model_preds$individual %in% individuals,], 
                    aes(x=virus,ymin=lower,ymax=upper),alpha=0.7,fill=fill_color,size=0.2) + 
        geom_line(data=model_preds, aes(x=virus, y=median),linetype="dotted",color="grey10")+
        geom_rect(ymin=max_titre,ymax=max_titre+2,xmin=0,xmax=max_x,fill="grey70")+
        geom_rect(ymin=min_titre-2,ymax=min_titre,xmin=0,xmax=max_x,fill="grey70")+
        scale_x_continuous(expand=c(0.01,0.01),breaks=as.numeric(names(virus_key)),labels=virus_key) +
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
        facet_wrap(individual~samples_label,ncol=2)
    titre_pred_p
}

