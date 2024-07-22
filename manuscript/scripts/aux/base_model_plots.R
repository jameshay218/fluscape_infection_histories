#################################
## Load in MCMC chains - only do this once, as slow
## Objects that come from this:
##     - inf_chain1/inf_chain: the full infection history chain
##     - theta_chain: the full mcmc chain for antibody kinetics parameters
if(!file.exists("r_data/base_quarter.RData")){
    chain_wd <- base_chain_wd
    source("scripts/aux/extract_infection_histories.R")
    save(inf_chain, theta_chain, par_tab, antigenic_map, titre_dat, strain_isolation_times,
         file=paste0(main_wd,"/r_data/measurement_quarter.RData"))
} else {
    load("r_data/base_quarter.RData")
    colnames(antigenic_map)[3] <- "inf_times"
}


## Create key to make better sample time labels
time_key <- titre_dat %>% select(samples) %>% distinct() %>% arrange(samples) %>% mutate(samples=samples/4) %>% pull(samples)
time_key <- convert_years_to_quarters(time_key)
names(time_key) <- unique(titre_dat$samples)
## Create key to make better virus labels
virus_key1 <- fluscape_dat %>% select(virus,Virus) %>% distinct()
virus_key <- virus_key1 %>% pull(Virus)
names(virus_key) <- virus_key1 %>% pull(virus)

use_indivs <- sample(unique(titre_dat$individual),5)
use_indivs <- c(372, 736, 801, 929, 1027)
DOBs <- titre_dat %>% filter(individual %in% use_indivs) %>% select(individual,DOB) %>% distinct()

p_titre_fits <- plot_infection_histories_long_mod(chain=theta_chain,infection_histories=inf_chain,
                                                  titre_dat=titre_dat,individuals=use_indivs,
                                                  antigenic_map=antigenic_map,
                                                  nsamp = 500,
                                                  strain_isolation_times = strain_isolation_times,
                                                  mu_indices=NULL,
                                                  measurement_indices_by_time = NULL,
                                                  par_tab=par_tab,time_key=time_key,virus_key=virus_key)[[1]]
colnames(antigenic_map)[3] <- "inf_times"
p_inf_hists <- generate_cumulative_inf_plots(inf_chain, 0, use_indivs, ages = DOBs,strain_isolation_times = antigenic_map$inf_times)
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
ggsave_jah(p_titre_fit_comb,figure_wd,"base_titre_fits_example",width=7,height=8)

## Find RMSE
model_preds <- get_titre_predictions(chain=theta_chain,infection_histories=inf_chain,
                                     titre_dat=titre_dat,individuals=unique(titre_dat$individual),
                                     antigenic_map=antigenic_map,
                                     strain_isolation_times = strain_isolation_times,
                                     add_residuals=TRUE,
                                     for_res_plot=FALSE,
                                     expand_titredat=FALSE,
                                     mu_indices=NULL,
                                     measurement_indices_by_time = NULL,
                                     par_tab=par_tab)
residuals <- model_preds$residuals$`50%`
RMSE <- sqrt(sum((residuals^2))/length(residuals))
print(paste0("RMSE: ",signif(RMSE,3)))

## Residuals for plot
model_preds <- get_titre_predictions(chain=theta_chain,infection_histories=inf_chain,
                                     titre_dat=titre_dat,individuals=unique(titre_dat$individual),
                                     antigenic_map=antigenic_map,
                                     strain_isolation_times = strain_isolation_times,
                                     add_residuals=TRUE,
                                     for_res_plot=TRUE,
                                     expand_titredat=FALSE,
                                     mu_indices=NULL,
                                     measurement_indices_by_time = NULL,
                                     par_tab=par_tab)

res_for_plot <- titre_dat %>% dplyr::select(virus) %>% bind_cols(model_preds[[1]])
res_for_plot <- res_for_plot %>% pivot_longer(-virus) %>% mutate(name=as.numeric(as.factor(name)))
res_for_plot$Strain <- virus_key[as.character(res_for_plot$virus)]
res_for_plot$Strain <- factor(res_for_plot$Strain, levels=virus_key)
p_res_lhs <- ggplot(res_for_plot) + 
    geom_histogram(aes(x=value,y=..ncount..),binwidth=1,fill="grey70",col="black") + 
    ylab("Scaled count") +
    xlab("Titer residual (observed - predicted)") +
    facet_wrap(~Strain,nrow=4) + 
    geom_vline(xintercept=0,col="red",size=0.5) + 
    theme_main +
    theme(plot.tag=element_text(face="bold")) +
    labs(tag="A")
p_res_rhs <- ggplot(res_for_plot) + 
    geom_histogram(aes(x=value,y=..ncount..),binwidth=1,fill="grey70",col="black") + 
    ylab("Scaled count") +
    xlab("Titer residual\n (observed - predicted)") +
    geom_vline(xintercept=0,col="red",size=0.5) + 
    theme_main+
    theme(plot.tag=element_text(face="bold")) +
    labs(tag="B")

p_res_hist <- p_res_lhs + p_res_rhs + plot_layout(widths=c(3,1))
ggsave_jah(p_res_hist,figure_wd,"hist_residuals_base",width=7,height=4)

## Plot attack rates and individual infection probabilities
source("scripts/aux/plot_individual_infections_base.R")

## Get antibody kinetics parameters
## note that ESS is smaller than it should be because only a subset of posterior samples are kept
## from loading the chains
## Convert waning rate to years
par_estimates <- plot_posteriors_theta(theta_chain %>% mutate(wane = wane*4), par_tab,plot_corr=FALSE,plot_mcmc=FALSE)$results
par_key <- c("mu"="µl","mu_short"="µs","sigma1"="σl","sigma2"="σs","wane"="ω","tau"="τ","error"="ε")
par_estimates <- par_estimates %>% filter(names %in% c("mu","mu_short","sigma1","sigma2","wane","tau","error")) %>% select(names, estimate) %>% rename("Parameter"="names","Estimate"="estimate")
par_estimates$Parameter <- par_key[par_estimates$Parameter]
write_csv(par_estimates,file="results/fluscape_parameter_estimates_base.csv")
