## Reads in the specified MCMC chains, thins them and pads the infection history chain with missing 0s

## Number of MCMC samples to use
n_samps <- 1000
ignore_chainno <- c()

## Load up data used in fitting
antigenic_map <- load_antigenic_map_file(location=chain_wd)
antigenic_map <- antigenic_map[,colnames(antigenic_map) != "X"]
strain_isolation_times <- antigenic_map$inf_times
titre_dat <- load_titre_dat(location=chain_wd)
titre_dat <- titre_dat[,colnames(titre_dat) != "X"]

## Anonymise titre data
titre_dat$DOB <- floor(titre_dat$DOB/buckets)*buckets
titre_dat$DOB <- pmax(1930*buckets, titre_dat$DOB)
titre_dat$Participant_ID <- titre_dat$individual


par_tab <- load_start_tab(location=chain_wd)
par_tab <- par_tab[,colnames(par_tab) != "X"]

## Read in chains
all_chains <- load_mcmc_chains(location=chain_wd, thin=1,burnin=burnin_use,
                               par_tab=par_tab,unfixed=FALSE,convert_mcmc=TRUE)
## Thin the chains
theta_chain <- as.data.frame(all_chains$theta_chain)
inf_chain <- all_chains$inf_chain

n_samps <- min(1000,nrow(theta_chain))


theta_samps <- unique(theta_chain[,"sampno"])
inf_chain_samps <- unique(inf_chain$sampno)
samps_in_both <- intersect(theta_samps, inf_chain_samps)

max_samps <- min(length(samps_in_both),n_samps)
get_sampnos <- sample(samps_in_both,max_samps)
inf_chain <- inf_chain[inf_chain$sampno %in% get_sampnos,]
theta_chain <- theta_chain[theta_chain$sampno %in% get_sampnos,]
inf_chain <- inf_chain[!(chain_no %in% ignore_chainno),]
theta_chain <- theta_chain[!(theta_chain[,"chain_no"] %in% ignore_chainno),]

list_chains <- all_chains$theta_list_chains
list_chains1 <- lapply(list_chains, function(x) x[,c("mu","mu_short","tau","wane","sigma1",
                                                     "sigma2","error","total_infections","lnlike",
                                                     "likelihood")])
color_scheme_set("viridis")
mcmc_trace(list_chains1)

n_alive <- get_n_alive_group(titre_dat, strain_isolation_times, TRUE)
#plot_infection_history_chains_time(inf_chain,
#                                   0,151:200,n_alive=n_alive,pad_chain=FALSE)[[1]] + 
#  theme(legend.position = "none",axis.text.x=element_text(size=4)) + scale_y_continuous(limits=c(0,1))


samps_all <- unique(inf_chain[,c("sampno","chain_no")])
samps_all$samp <- 1:nrow(samps_all)

samps_all_theta <- unique(theta_chain[,c("sampno","chain_no")])
samps_all_theta$samp <- 1:nrow(samps_all_theta)

theta_chain <- merge(theta_chain, samps_all_theta)
theta_chain$sampno <- theta_chain$samp
#theta_chain$chain_no <- 1

samps <- unique(theta_chain$sampno)
samps <- sample(samps, n_samps)

theta_chain <- theta_chain[theta_chain$sampno %in% samps,]

## Expand inf_chain to have 0 entries as well
inf_chain <- merge(inf_chain, samps_all)
inf_chain <- inf_chain[inf_chain$samp %in% samps,]

inf_chain <- inf_chain[,c("i","j","x","samp")]
colnames(inf_chain)[4] <- "sampno"

is <- unique(inf_chain$i)
js <- unique(inf_chain$j)

samps <- unique(inf_chain$samp)
expanded_values <- data.table::CJ(
  i = is,
  j = js,
  sampno = samps
)
diff_infs <- fsetdiff(expanded_values, inf_chain[, c("i", "j", "sampno")])
diff_infs$x <- 0
inf_chain <- rbind(inf_chain, diff_infs)
inf_chain$chain_no <- 1
## Number of samples with a 1 divided by total samples
colnames(antigenic_map)[3] <- "inf_times"
#save(inf_chain, theta_chain, par_tab, antigenic_map, titre_dat, strain_isolation_times,
#     file="~/Drive/Influenza/serosolver/FluScape_epidemiology/r_data/measurement_annual_bb_new.RData")
