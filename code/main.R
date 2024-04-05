
fluscape_serosolver_main <- function(serosolver_wd,mcmc_pars,subset_indivs=NULL,
                                     run_name, 
                                     main_wd,
                                     save_wd, 
                                     chain_wd, 
                                     buckets, prior_version, n_chains,
                                     titre_dat_file, antigenic_map_file,par_tab_file,rho_file,
                                     inf_hist_prior=c(2,5),
                                     solve_likelihood=TRUE,
                                     estimate_measurement_offsets=FALSE,
                                     use_measurement_random_effects=FALSE,append_fixed_rho=FALSE,
                                     create_prior_func=NULL,
                                     proposal_ratios=NULL,
                                     mixed_time_cutoff=NULL,
                                     rerun=TRUE){
    top_wd <- getwd()
    if(!dir.exists(save_wd)) dir.create(save_wd,recursive = TRUE)
    if(!dir.exists(chain_wd)) dir.create(chain_wd,recursive = TRUE)
    
    setwd(main_wd)
    
    print(paste0("Starting in: ", main_wd))
    print(paste0("Saving plots to: ", save_wd))
    print(paste0("Saving chains to: ", chain_wd))
    
    titre_dat <- read.csv(titre_dat_file,stringsAsFactors = FALSE)
    titre_dat <- titre_dat %>% arrange(individual, samples, virus, run)
    if(!("group" %in%colnames(titre_dat))){
        titre_dat$group <- 1
    }
    
    if(!is.null(subset_indivs)){
        message(cat("Subsetting to ", subset_indivs, " individuals"))
        use_indivs <- sample(unique(titre_dat$individual), subset_indivs)
        titre_dat_use <- titre_dat %>% filter(individual %in% use_indivs)
        titre_dat_use_map <- titre_dat_use %>% select(individual) %>% distinct() %>% mutate(individual_new = 1:n())
        titre_dat_use <- titre_dat_use %>% left_join(titre_dat_use_map) %>% rename(individual_old=individual,individual=individual_new) %>% select(-individual_old)
        titre_dat_use_map <- titre_dat_use_map %>% rename(individual_old=individual,individual=individual_new)
        write.csv(titre_dat_use_map, paste0(chain_wd, "/",run_name, "_indiv_id_map.csv"))
        
    } else {
        titre_dat_use <- titre_dat
        titre_dat_use_map <- NULL
    }
    
    ## Antigenic map
    antigenic_map <- read.csv(antigenic_map_file)
    if(!is.null(mixed_time_cutoff)){
        antigenic_map <- antigenic_map[antigenic_map$inf_times >= 1968*buckets & antigenic_map$inf_times <= max(titre_dat_use$samples),]
        strain_isolation_times <- c(seq(1968*buckets, mixed_time_cutoff*buckets, by=buckets),(mixed_time_cutoff*buckets+1):max(antigenic_map$inf_times))
        antigenic_map <- antigenic_map[antigenic_map$inf_times %in% strain_isolation_times,]
    } else {
        antigenic_map <- antigenic_map[antigenic_map$inf_times >= 1968*buckets & antigenic_map$inf_times <= max(titre_dat_use$samples),]
        strain_isolation_times <- antigenic_map$inf_times
    }
    
    ## Set up parameter table
    par_tab <- read.csv(par_tab_file,stringsAsFactors=FALSE)
    if(prior_version != 1){
        par_tab <- par_tab[par_tab$names != "phi",]
    }
    par_tab[par_tab$names %in% c("alpha","beta"),c("values")] <- inf_hist_prior
    
    ## Add entries for each measurement offset
    if(estimate_measurement_offsets){
        measurement_indices <- match(antigenic_map$inf_times, antigenic_map$inf_times)
        for(i in 1:length(measurement_indices)){
            par_tab <- rbind(par_tab, data.frame("names"="rho","values"=rnorm(1,0,1),"fixed"=1,"steps"=0.1,
                                                 "lower_bound"=-3,"upper_bound"=3,"lower_start"=-2,"upper_start"=2,"type"=3))
        }
        fixed <- as.numeric(!(strain_isolation_times %in% unique(titre_dat_use$virus)))
        par_tab[par_tab$names %in% c("rho_mean"),"fixed"] <- 1
        
        if(!use_measurement_random_effects){
            par_tab[par_tab$names == "rho_sd","fixed"] <- 1
        } else {
            par_tab[par_tab$names == "rho_sd","fixed"] <- 0
        }
        par_tab[par_tab$names== "rho","fixed"] <- 1
        par_tab[par_tab$names == "rho",][which(fixed == 0),"fixed"] <- 0
        par_tab[par_tab$names == "rho" & par_tab$fixed == 1, "values"] <- 0
    } else {
        if(append_fixed_rho){
            par_tab_rhos <- read.csv(rho_file,stringsAsFactors = FALSE)
            par_tab <- bind_rows(par_tab, par_tab_rhos)
            measurement_indices <- rep(1:48,each=buckets)
            if(!is.null(mixed_time_cutoff)){
                strain_isolation_times_full <- c(seq(1968*buckets, max(antigenic_map$inf_times)))
                measurement_indices <- measurement_indices[1:length(strain_isolation_times_full)]
                ## Create measurement bias indices, this bit is tricky
                measurement_indices <- measurement_indices[match(strain_isolation_times, strain_isolation_times_full)]
                
            }
        } else {
            measurement_indices <- NULL
            par_tab <- par_tab[par_tab$names != "rho",]
            par_tab <- par_tab[par_tab$names != "rho_mean",]
            par_tab <- par_tab[par_tab$names != "rho_sd",]
        }
    }
    
    f <- create_posterior_func(par_tab,titre_dat_use,antigenic_map=antigenic_map,
                               strain_isolation_times = strain_isolation_times,
                               version=prior_version,
                               solve_likelihood=solve_likelihood,n_alive=NULL,
                               measurement_indices_by_time=measurement_indices)
    
    ## Time runs and use dopar to run multiple chains in parallel
    if(rerun){
        print("Starting fits...")
        t1 <- Sys.time()
        filenames <- paste0(chain_wd, "/",run_name, "_",1:n_chains)
        res <- foreach(x = filenames, .packages = c('data.table','plyr',"dplyr")) %dopar% {
        #for(x in filenames){
            devtools::load_all(serosolver_wd)
            index <- 1
            lik <- -Inf
            inf_hist_correct <- 1
            while((!is.finite(lik) || inf_hist_correct > 0) & index < 100){
                start_tab <- generate_start_tab(par_tab)
                start_inf <- setup_infection_histories_titre(titre_dat_use,strain_isolation_times,
                                                             space=5,titre_cutoff=2,sample_prob=0.9)
                
                inf_hist_correct <- sum(check_inf_hist(titre_dat_use, strain_isolation_times, start_inf))
                y <- f(start_tab$values, start_inf)
                lik <- sum(y[[1]])
                index <- index + 1
                print(index)
            } 
            
            if(index >= 100){
                browser()
                print(lik)
            }
            
            write.csv(start_tab, paste0(x, "_start_tab.csv"))
            write.csv(start_inf, paste0(x, "_start_inf_hist.csv"))
            write.csv(antigenic_map, paste0(x, "_antigenic_map.csv"))
            write.csv(titre_dat_use, paste0(x, "_titre_dat.csv"))
            
            res <- serosolver::run_MCMC(start_tab, titre_dat_use, 
                                        antigenic_map, start_inf_hist=start_inf,filename=x,
                                        CREATE_POSTERIOR_FUNC=create_posterior_func, 
                                        CREATE_PRIOR_FUNC = create_prior_func,
                                        version=prior_version,
                                        measurement_indices=measurement_indices,
                                        measurement_random_effects = use_measurement_random_effects,
                                        mcmc_pars=mcmc_pars,
                                        proposal_ratios=proposal_ratios,
                                        solve_likelihood=solve_likelihood)
        }
        run_time_fast <- Sys.time() - t1
        run_time_fast
    }
    
    ## Read in chains for trace plot
    chains <- load_mcmc_chains(chain_wd,convert_mcmc=TRUE,burnin = mcmc_pars["adaptive_period"],unfixed = TRUE,thin=1)
    pdf(paste0(save_wd,"/",run_name,"_chain.pdf"))
    plot(as.mcmc.list(chains$theta_list_chains))
    dev.off()
    
    ## Read in chains for all other plots
    chains <- load_mcmc_chains(chain_wd,convert_mcmc=FALSE,burnin = mcmc_pars["adaptive_period"],unfixed = FALSE,thin=1)
    chain <- as.data.frame(chains$theta_chain)
    inf_chain <- chains$inf_chain
    
    ## Plot kinetics parameter estimates and number of infections
    p1 <- plot_posteriors_theta(chain, par_tab, burnin=0,samples=100,TRUE,FALSE,FALSE, TRUE,"")
    p_total_inf <- plot_total_number_infections(inf_chain,FALSE)
    ggsave(paste0(save_wd,"/",run_name,"_total_inf.pdf"),p_total_inf[[1]],height=5,width=7,units="in",dpi=300)
    
    ## Plot attack rates
    p_ar <- plot_attack_rates(inf_chain, titre_dat_use, strain_isolation_times,pad_chain=FALSE,plot_den=FALSE,n_alive = NULL,
                              prior_pars = c("prior_version"=prior_version,
                                             "alpha"=inf_hist_prior[1],
                                             "beta"=inf_hist_prior[2]))
    ggsave(paste0(save_wd,"/",run_name,"_ar.pdf"),p_ar,height=7,width=8,units="in",dpi=300)
    
    ## Plot model fits
    titre_preds <- get_titre_predictions(chain = chain[chain$chain_no == 1,], 
                                         infection_histories = inf_chain[inf_chain$chain_no == 1,], 
                                         titre_dat = titre_dat_use, 
                                         measurement_indices_by_time = NULL,
                                         individuals = unique(titre_dat_use$individual),
                                         antigenic_map = antigenic_map, 
                                         par_tab = par_tab,expand_titredat=FALSE)
    
    use_indivs <- sample(unique(titre_dat_use$individual), 9)
    titre_pred_p <- ggplot(titre_preds$predictions[titre_preds$predictions$individual %in% use_indivs,])+
        geom_ribbon(data=titre_preds$predicted_observations[titre_preds$predicted_observations$individual %in% use_indivs,], 
                    aes(x=virus,ymin=lower,ymax=upper),fill="grey90") +
        geom_ribbon(aes(x=virus,ymin=lower, ymax=upper),fill="gray70")+
        geom_line(aes(x=virus, y=median))+
        geom_point(aes(x=virus, y=titre))+
        coord_cartesian(ylim=c(0,8))+
        ylab("log titre") +
        xlab("Time of virus circulation") +
        theme_classic() +
        facet_grid(individual~samples)
    
    ggsave(paste0(save_wd,"/",run_name,"_titre_fits.pdf"),titre_pred_p,height=7,width=8,units="in",dpi=300)
    setwd(top_wd)
}