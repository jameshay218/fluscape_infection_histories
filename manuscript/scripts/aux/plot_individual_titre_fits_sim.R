## Read in par tab, data, infection history etc
par_tab_sim <- load_start_tab(sim_chain_wd)
sim_data <- load_titre_dat(sim_chain_wd)
#sim_inf_hist <- read_csv(paste0(main_wd,"../data/simulated/sim_quarterly_fluscape_infection_histories.csv"))

sim_inf_hist <- read_csv("~/Documents/GitHub/fluscape_serosolver/data/simulated/sim_quarterly_fluscape_supplement_1_infection_histories.csv")

#real_par_tab <- read_csv(paste0(main_wd,"../data/simulated/sim_quarterly_fluscape_clusters_par_tab.csv"))
real_par_tab <- read_csv("~/Documents/GitHub/fluscape_serosolver/data/simulated/sim_quarterly_fluscape_supplement_1_par_tab.csv")
sim_map <- load_antigenic_map_file(sim_chain_wd)
colnames(antigenic_map)[3] <- "inf_times"


## Dimensions of simulated data
length(unique(sim_data$individual))
length(unique(sim_data$virus))
diff(range(sim_map$inf_times))/4
diff(range(sim_data$samples))/4

## Read in true attack rates
#sim_true_ar <- read_csv(paste0(main_wd,"../data/simulated/sim_quarterly_fluscape_attack_rates.csv"))
sim_true_ar <- read_csv("~/Documents/GitHub/fluscape_serosolver/data/simulated/sim_quarterly_fluscape_supplement_1_attack_rates.csv")

#################################
## Plot some example titre fits
plot_infection_histories_long_mod <- function(chain, infection_histories, titre_dat,
                                              individuals, antigenic_map=NULL, 
                                              strain_isolation_times=NULL, par_tab,
                                              nsamp = 100,
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
                                        select(individual, samples_label,samples) %>% distinct(),inf_hist_densities)%>% 
      filter(variable <= samples)
    
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
        facet_wrap(individual~samples_label,ncol=2)
    titre_pred_p
}
## Create key to make better sample time labels
time_key <- titre_dat %>% select(samples) %>% distinct() %>% arrange(samples) %>% mutate(samples=samples/4) %>% pull(samples)
time_key <- convert_years_to_quarters(time_key)
names(time_key) <- unique(titre_dat$samples)

use_indivs <- sample(unique(titre_dat$individual),3)
use_indivs <- c(22,75,747)
DOBs <- titre_dat %>% filter(individual %in% use_indivs) %>% select(individual,DOB) %>% distinct()

p_titre_fits <- plot_infection_histories_long_mod(chain=theta_chain,infection_histories=inf_chain,
                                                  titre_dat=titre_dat,individuals=use_indivs,
                                                  antigenic_map=antigenic_map,
                                                  nsamp=500,
                                                  strain_isolation_times = strain_isolation_times,
                                                  mu_indices=rep(1:48,each=4),
                                                  measurement_indices_by_time = rep(1:48,each=4),
                                                  par_tab=par_tab,time_key=time_key,virus_key=NULL)
colnames(antigenic_map)[3] <- "inf_times"
p_inf_hists <- generate_cumulative_inf_plots(inf_chain, 0, use_indivs,real_inf_hist=sim_inf_hist,
                                             ages = DOBs,strain_isolation_times = antigenic_map$inf_times)
p_cumu_infhist <- p_inf_hists[[1]] + scale_x_continuous(breaks=seq(1970,2015,by=5)*4,labels=seq(1970,2015,by=5)) + 
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

p_titre_fit_comb <- (p_titre_fits+labs(tag="A") + theme(plot.tag = element_text(face="bold"))) + p_cumu_infhist + plot_layout(ncol=2,widths=c(2.5,1))
ggsave_jah(p_titre_fit_comb,figure_wd,"sim_titre_fits_example",width=7,height=8)

## Find RMSE tajes a very long time
if(FALSE){
model_preds <- get_titre_predictions(chain=theta_chain,infection_histories=inf_chain,
                                     titre_dat=titre_dat,individuals=unique(titre_dat$individual),
                                     antigenic_map=antigenic_map,
                                     strain_isolation_times = strain_isolation_times,
                                     add_residuals=TRUE,
                                     for_res_plot=FALSE,
                                     nsamp=100,
                                     expand_titredat=FALSE,
                                     mu_indices=rep(1:48,each=4),
                                     measurement_indices_by_time = rep(1:48,each=4),
                                     par_tab=par_tab)
residuals <- model_preds$residuals$`50%`
RMSE <- sqrt(sum((residuals^2))/length(residuals))
print(paste0("RMSE: ",signif(RMSE,3)))
}

p_sim_ar <- plot_attack_rates_monthly(inf_chain %>% mutate(group=1),sim_data %>% mutate(group=1),
                                      antigenic_map$inf_times,buckets=4,pad_chain=TRUE)

n_alive <- get_n_alive(sim_data, strain_isolation_times)
sim_true_ar$AR_sim <- sim_true_ar$AR
sim_true_ar$AR[1:189] <- colSums(sim_inf_hist[,1:189])/n_alive



## Calculate infection histories and attack rates from the chain
infection_histories <- inf_chain
colnames(infection_histories) <- c("individual","time","x","sampno","chain_no")
## Convert time incidences to annual
infection_histories <- infection_histories %>% mutate(time = floor((time-1)*(1/4)) + 1)
setkey(infection_histories, "sampno","individual","time","chain_no")
infection_histories <- infection_histories[,list(x=as.numeric(any(x==1))),by=key(infection_histories)]

sim_ar_true_annual <- sum_buckets(colSums(sim_inf_hist),rep(4,48))

est_ars_annual <- get_attack_rates(inf_chain %>% mutate(group=1),unique(floor(antigenic_map$inf_times/4)),
                            titre_dat=sim_data %>% mutate(group=1) %>% mutate(DOB=floor(DOB/4)),
                            n_alive=NULL,by_group = FALSE)
n_alive_annual <- get_n_alive(sim_data %>% mutate(group=1) %>% mutate(DOB=floor(DOB/4)),unique(floor(antigenic_map$inf_times/4)))

compare_ar_annual <- infection_histories %>% group_by(time,sampno, chain_no) %>% 
  dplyr::summarize(y=sum(x))%>% 
  left_join(data.frame(time=1:48,n_alive=n_alive_annual)) %>% mutate(y=y/n_alive) %>% 
  group_by(time) %>% 
  dplyr::summarize(med_y=median(y),lower2=quantile(y,0.025),upper2=quantile(y,0.975),
                   lower=quantile(y,0.25),upper=quantile(y,0.75)
                   ) %>% 
  left_join(data.frame(time=1:48,y_true=sim_ar_true_annual)) %>% 
  left_join(data.frame(time=1:48,n_alive=n_alive_annual)) %>% 
  mutate(y_true=y_true/n_alive)
compare_ar_annual_base <- compare_ar_annual %>% mutate(Model="Main")

p_sim_ar_annual <- 
  ggplot(compare_ar_annual) + 
  #geom_ribbon(aes(x = time, ymin = lower, ymax = upper), fill = "red", alpha = 0.2) +
  geom_pointrange(aes(x = time, y=med_y,ymin = lower2, ymax = upper2), col = "red", alpha = 0.4) +
  
   #geom_point(aes(x=time,y=med_y),col="red") + 
   geom_point(aes(x=time,y=y_true),col="black") +
  ylab("Estimated monthly per capita incidence") +
  xlab("Date") +
  scale_x_continuous(breaks=seq(1970,2015,by=5)-1968,labels=seq(1970,2015,by=5)) + theme_pubr()+
  scale_y_continuous(limits=c(0,1),expand=c(0,0)) +
  theme(legend.title=element_text(size=7),
        legend.text=element_text(size=7),
        legend.position="none",
        legend.margin = margin(-1,-1,-3,-1),
        axis.title=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=7),
        axis.text.y=element_text(size=7),
        plot.tag = element_text(face="bold"),
        strip.text=element_blank(),strip.background = element_blank(),
        plot.margin=margin(r=15,t=5,l=5))


## Plot estimated attack rate vs. real
p_sim_ar <- p_sim_ar + geom_line(data = sim_true_ar, aes(x = year, y = AR),col = "black", size = 0.25) + ylab("Per capita incidence per quarter") +
    scale_x_continuous(breaks=seq(1970,2015,by=5)*4,labels=seq(1970,2015,by=5)) + theme_pubr()+
  scale_y_continuous(limits=c(0,1),expand=c(0,0)) +
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

## Plot density vs. true value
par_key <- c("mu"="µ_l","mu_short"="µ_s","tau"="τ","wane"="ω","sigma1"="σ_l","sigma2"="σ_s","total_infections"="ΣZ","error"="ε")
chains <- load_mcmc_chains(sim_chain_wd, par_tab,burnin=200000,unfixed=TRUE,convert_mcmc=FALSE)
sim_theta_chain <- chains$theta_chain
sim_theta_chain1 <- sim_theta_chain[,c("sampno","chain_no","mu","mu_short","tau","wane","sigma1","sigma2","error","total_infections")]
sim_theta_chain1 <- sim_theta_chain1 %>% pivot_longer(-c(sampno,chain_no))


sim_theta_chain1 %>% group_by(name) %>% dplyr::summarize(med=median(value),lower=quantile(value,0.025),upper=quantile(value,0.975)) %>% distinct()

#sim_theta_chain1 <- sim_theta_chain1 %>% filter(name != "error")

sim_theta_chain1$name <- par_key[sim_theta_chain1$name]
real_par_tab <- read_csv("~/Documents/GitHub/fluscape_serosolver/data/simulated/sim_quarterly_fluscape_supplement_1_par_tab.csv")

real_par_tab <- real_par_tab %>% 
  select(names,values,lower_bound,upper_bound) %>%
  filter(names %in% c("mu","mu_short","tau","wane","sigma1","sigma2","error")) %>% 
    rename(name=names)
real_par_tab$name <- par_key[real_par_tab$name]

real_par_tab <- bind_rows(real_par_tab, tibble(name="ΣZ",values=sum(sim_inf_hist),lower_bound=0,upper_bound=prod(dim(sim_inf_hist))))

equal_breaks <- function(n = 4, s = 0.05, ...){
  function(x){
    # rescaling
    d <- s * diff(range(x)) / (1+2*s)
    round(seq(min(x)+d, max(x)-d, length=n),3)
  }
}

p_sim_recover <- ggplot(sim_theta_chain1) + 
    geom_density(aes(x=value),fill="blue",alpha=0.25) + 
    geom_vline(data=real_par_tab,aes(xintercept=values),linetype="dashed") +
    facet_wrap(~name,scales="free",ncol=2)+ theme_pubr()+
  scale_x_continuous(breaks=equal_breaks(n=4, s=0.05), 
                     expand = c(0.05, 0)) + 
    theme(legend.title=element_text(size=7),
          legend.text=element_text(size=7),
          legend.position="none",
          legend.margin = margin(-1,-1,-3,-1),
          axis.title=element_text(size=10),
          strip.text=element_text(size=7),
          axis.text.x=element_text(size=5),
          axis.text.y=element_text(size=5),
          plot.tag = element_text(face="bold"),
          plot.margin=margin(r=15,t=5,l=5)) +labs(tag="C") + 
    xlab("Value") + ylab("Density") + labs(tag="D")
p_titre_fits <- p_titre_fits+labs(tag="A") + theme(plot.tag = element_text(face="bold"))
p_top <- p_titre_fits + p_cumu_infhist + plot_layout(ncol=2,widths=c(2.5,1))
p_bot <- p_sim_ar + p_sim_recover + plot_layout(ncol=2,widths=c(1.75,1.75))

p_left <- p_titre_fits + p_sim_ar + plot_layout(ncol=1,heights=c(3,1))
p_right <- p_cumu_infhist + p_sim_recover + plot_layout(ncol=1,heights=c(1,1))

p_sim_recovery <- (p_top / p_bot) + plot_layout(heights=c(2.25,1))
p_sim_recovery_alt <- (wrap_elements(p_left) | wrap_elements(p_right))

ggsave_jah(p_sim_recovery_alt,figure_wd,"p_sim_recovery_all",8,8)
ggsave_jah(p_sim_ar_annual,figure_wd,"p_sim_recovery_ar_annual",7,4)
