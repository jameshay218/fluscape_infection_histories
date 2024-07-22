## Create key to make better sample time labels
time_key <- strain_isolation_times/4#titre_dat %>% select(samples) %>% distinct() %>% arrange(samples) %>% mutate(samples=samples/4) %>% pull(samples)
time_key <- convert_years_to_quarters(time_key)
names(time_key) <- strain_isolation_times#unique(titre_dat$samples)
## Create key to make better virus labels
virus_key1 <- fluscape_dat %>% select(virus,Virus) %>% distinct()
virus_key <- virus_key1 %>% pull(Virus)
names(virus_key) <- virus_key1 %>% pull(virus)

use_indivs <- sample(unique(titre_dat$individual),5)
use_indivs <- c(372, 736, 801, 929, 1027)
DOBs <- titre_dat %>% filter(individual %in% use_indivs) %>% select(individual,DOB) %>% distinct()

DOBs <- DOBs %>% mutate(DOB = floor(DOB/4)*4)
titre_dat <- titre_dat %>% mutate(DOB = floor(DOB/4)*4)

plot_indiv_key <- data.frame(i = use_indivs, i_new = seq_along(use_indivs))


p_titre_fits <- plot_infection_histories_long_mod(chain=theta_chain,
                                                  infection_histories=inf_chain %>% filter(i %in% use_indivs) %>%
                                                    left_join(plot_indiv_key) %>% select(-i) %>% 
                                                    rename(i = i_new),
                                                  titre_dat=titre_dat %>% filter(individual %in% use_indivs) %>%
  left_join(plot_indiv_key %>% rename(individual = i)) %>% select(-individual) %>% rename(individual = i_new),
individuals=seq_along(use_indivs),
                                                  antigenic_map=antigenic_map,
                                                  strain_isolation_times = strain_isolation_times,
                                                  mu_indices=rep(1:48,each=4),
                                                  measurement_indices_by_time = rep(1:48,each=4),
                                                  par_tab=par_tab,time_key=time_key,virus_key=virus_key)


plot_data <- p_titre_fits[[2]]
p_titre_fits <- p_titre_fits[[1]]

## Extract data for plot
plot_data1 <- plot_data[[1]] %>% select(individual, virus, titre,lower,run, median,upper,samples,samples_label,virus_label) %>% rename(lower_obs=lower,upper_obs=upper,median_obs=median) 
plot_data2 <-plot_data[[2]]%>% mutate(infection_time=variable/4)
plot_data3 <- plot_data[[3]] %>% select(individual, virus, lower, median,upper,samples) %>% distinct()
plot_data4 <- left_join(plot_data1,plot_data3, by=c("individual","virus","samples")) %>% select(-c(samples,virus)) %>%
  rename(`Lower 95% prediction interval`=lower_obs, `Posterior median observation`=median_obs,
         `Upper 95% prediction interval`=upper_obs, `Lower 95% CrI`=lower, `Upper 95% CrI`=upper, `Posterior median`=median,`Virus`=virus_label,`Sample time`=`samples_label`,`Repeat number`=run)
plot_data2 <- plot_data2 %>% rename(`Sample time`=`samples_label`) %>% select(-samples) %>% rename(`Posterior probability of infection`=value)

write.csv(plot_data4,"~/Documents/GitHub/fluscape_infection_histories//data/figure_data/FigS4Ai.csv",row.names=FALSE)
write.csv(plot_data2,"~/Documents/GitHub/fluscape_infection_histories//data/figure_data/FigS4Aii.csv",row.names=FALSE)

colnames(antigenic_map)[3] <- "inf_times"

p_titre_fits$data %>% select(individual, samples_label, virus_label, titre, lower,lower_50,median,upper_50,upper)


p_inf_hists <- generate_cumulative_inf_plots(inf_chain%>% filter(i %in% use_indivs) %>%
                                               left_join(plot_indiv_key) %>% select(-i) %>% 
                                               rename(i = i_new), 0, seq_along(use_indivs), ages = DOBs %>% mutate(individual=1:n()),strain_isolation_times = antigenic_map$inf_times)
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
ggsave_jah(p_titre_fit_comb,figure_wd,"titre_fits_example",width=7,height=8)

p_cumu_infhist_data <- p_cumu_infhist$data %>% dplyr::select(-chain_no) %>% dplyr::mutate(`Circulation time`=(as.numeric(as.character(variable))/4)) %>% rename(`Lower 95% CrI`=lower,
                                                                                                                                                                `Upper 95% CrI`=upper,
                                                                                                                                                                `Posterior median cumulative infections`=`median`)
write.csv(p_cumu_infhist_data,"~/Documents/GitHub/fluscape_infection_histories//data/figure_data/FigS4B.csv",row.names=FALSE)

## Find RMSE -- takes a long time
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
}
