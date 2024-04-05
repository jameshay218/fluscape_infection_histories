by_year <- TRUE

DOBs <- get_DOBs(titre_dat)[,2]
age_mask <- create_age_mask(DOBs,strain_isolation_times)
strain_mask <- create_strain_mask(titre_dat, strain_isolation_times)
masks <- cbind(age_mask, strain_mask)
masks <- data.frame(masks)
masks$i <- 1:nrow(masks)

## Finding cumulative incidence in HK epidemic periods
e1_start <- 2010*4 + 3 - 1
e1_end <- 2010*4 + 4 - 1

e2_start <- 2012*4 + 1 - 1
e2_end <- 2012*4 + 2 - 1

e3_start <- 2013*4 + 3 - 1
e3_end <- 2013*4 + 3 - 1

e4_start <- 2014*4 + 2 - 1
e4_end <- 2014*4 + 3 - 1

e_times <- data.frame("start"=c(e1_start, e2_start, e3_start, e4_start),
                      "end"=c(e1_end,e2_end,e3_end,e4_end))
e_times <- e_times - 1968*4 + 1


e_2011_start <- 2011*4
e_2011_end <- 2011*4 + 3

## Do it by year instead
if(by_year){
  e1_start <- 2010*4 + 1 - 1
  e1_end <- 2010*4 + 4 - 1
  
  e2_start <- 2012*4 + 1 - 1
  e2_end <- 2012*4 + 4 - 1
  
  e3_start <- 2013*4 + 1 - 1
  e3_end <- 2013*4 + 4 - 1
  
  e4_start <- 2014*4 + 1 - 1
  e4_end <- 2014*4 + 4 - 1
}

n_alive_all <- get_n_alive(titre_dat,strain_isolation_times)

## Get number infected at least once per year
inf_chain_annual <- inf_chain
inf_chain_annual <- inf_chain_annual %>% mutate(j = 1968 + floor(j/4))
setkey(inf_chain_annual, "sampno","i","j","chain_no")
inf_chain_annual <- inf_chain_annual[,list(infected=as.numeric(any(x==1))),by=key(inf_chain_annual)]

setkey(inf_chain_annual, "sampno","j","chain_no")
ars_annual <- inf_chain_annual[,list(total_infs=sum(infected)),by=key(inf_chain_annual)]

n_alive_annual <- tibble(
  j =  unique(floor(strain_isolation_times/4)),
  n_alive = get_n_alive(titre_dat %>% mutate(DOB = floor(DOB/4)+1), unique(floor(strain_isolation_times/4))))
ars_annual <- left_join(ars_annual, n_alive_annual)

ars_annual %>% mutate(ar=total_infs/n_alive) %>% group_by(sampno,chain_no) %>% dplyr::summarize(med=median(ar),low=quantile(ar,0.025),high=quantile(ar,0.975)) %>%
  ungroup() %>% dplyr::summarize(med1=median(med),low=quantile(med,0.025),high=quantile(med,0.975))

ar_estimates_annual <- ars_annual %>% mutate(ar = total_infs/n_alive) %>% group_by(j) %>% dplyr::summarize(med_ar = median(ar),lower=quantile(ar, 0.025),upper=quantile(ar,0.975))

save(ar_estimates_annual,file=paste0(main_wd,"r_data/fluscape_ar_estimates.RData"))

## Epidemic 1
e1_infs <- inf_chain[inf_chain$j >= e1_start  - 1968*4 + 1 & inf_chain$j <= e1_end  - 1968*4 + 1,]
e1_times <- e1_start:e1_end - 1968*4 + 1
#n_alive <- get_n_alive(titre_dat, e1_times)
n_alive <- n_alive_all[e1_times]

## How many times were people infected
setkey(e1_infs, "sampno","i")
e1_inf_times <- e1_infs[,list(total_infs=sum(x)),by=key(e1_infs)]

## Some people were infected twice, apparently
table(e1_inf_times$total_infs)
e1_inf_times$any_inf <- e1_inf_times$total_infs
e1_inf_times$any_inf[e1_inf_times$any_inf > 1] <- 1

## Get number of people that were reinfected
e1_reinfs <- e1_inf_times[e1_inf_times$total_infs > 1,]
e1_reinfs <- ddply(e1_reinfs, ~sampno, nrow)
y <- e1_reinfs$V1/n_alive[1] * 100
e1_reinfs_est <- c(mean(y), quantile(y, c(0.025,0.5,0.975)))

## Now get attack rate for this time period
setkey(e1_inf_times,"sampno")
e1_inc <- e1_inf_times[,list(sum(any_inf)),by=key(e1_inf_times)]
e1_ests <- c(mean(e1_inc$V1/max(n_alive)), quantile(e1_inc$V1/n_alive[1], c(0.025,0.5,0.975)))

####################################################
## Epidemic 2
e2_infs <- inf_chain[inf_chain$j >= e2_start  - 1968*4 + 1 & inf_chain$j <= e2_end  - 1968*4 + 1,]
e2_times <- e2_start:e2_end - 1968*4 + 1
#n_alive <- get_n_alive(titre_dat, e2_times)
n_alive <- n_alive_all[e2_times]

## How many times were people infected
setkey(e2_infs, "sampno","i")
e2_inf_times <- e2_infs[,list(total_infs=sum(x)),by=key(e2_infs)]

## Some people were infected twice, apparently
table(e2_inf_times$total_infs)
e2_inf_times$any_inf <- e2_inf_times$total_infs
e2_inf_times$any_inf[e2_inf_times$any_inf > 1] <- 1

## Get number of people that were reinfected
e2_reinfs <- e2_inf_times[e2_inf_times$total_infs > 1,]
e2_reinfs <- ddply(e2_reinfs, ~sampno, nrow)
y <- e2_reinfs$V1/n_alive[1] * 100
e2_reinfs_est <- c(mean(y), quantile(y, c(0.025,0.5,0.975)))

## Now get attack rate for this time period
setkey(e2_inf_times,"sampno")
e2_inc <- e2_inf_times[,list(sum(any_inf)),by=key(e2_inf_times)]
e2_ests <- c(mean(e2_inc$V1/max(n_alive)), quantile(e2_inc$V1/n_alive[1], c(0.025,0.5,0.975)))

####################################################
## Epidemic 3
e3_infs <- inf_chain[inf_chain$j >= e3_start  - 1968*4 + 1 & inf_chain$j <= e3_end  - 1968*4 + 1,]
e3_times <- e3_start:e3_end - 1968*4 + 1
#n_alive <- get_n_alive(titre_dat, e3_times)
n_alive <- n_alive_all[e3_times]

## How many times were people infected
setkey(e3_infs, "sampno","i")
e3_inf_times <- e3_infs[,list(total_infs=sum(x)),by=key(e3_infs)]

## Some people were infected twice, apparently
table(e3_inf_times$total_infs)
e3_inf_times$any_inf <- e3_inf_times$total_infs
e3_inf_times$any_inf[e3_inf_times$any_inf > 1] <- 1

## Get number of people that were reinfected
e3_reinfs <- e3_inf_times[e3_inf_times$total_infs > 1,]
e3_reinfs <- ddply(e3_reinfs, ~sampno, nrow)
y <- (e3_reinfs$V1/max(n_alive)) * 100
e3_reinfs_est <- c(mean(y), quantile(y, c(0.025,0.5,0.975)))

## Now get attack rate for this time period
setkey(e3_inf_times,"sampno")
e3_inc <- e3_inf_times[,list(sum(any_inf)),by=key(e3_inf_times)]
e3_ests <- c(mean(e3_inc$V1/n_alive[1]), quantile(e3_inc$V1/n_alive[1], c(0.025,0.5,0.975)))

####################################################
## Epidemic 4
e4_infs <- inf_chain[inf_chain$j >= e4_start  - 1968*4 + 1 & inf_chain$j <= e4_end  - 1968*4 + 1,]
e4_times <- e4_start:e4_end - 1968*4 + 1
#n_alive <- get_n_alive(titre_dat, e4_times)
n_alive <- n_alive_all[e4_times]

## How many times were people infected
setkey(e4_infs, "sampno","i")
e4_inf_times <- e4_infs[,list(total_infs=sum(x)),by=key(e4_infs)]

## Some people were infected twice, apparently
table(e4_inf_times$total_infs)
e4_inf_times$any_inf <- e4_inf_times$total_infs
e4_inf_times$any_inf[e4_inf_times$any_inf > 1] <- 1


## Get number of people that were reinfected
e4_reinfs <- e4_inf_times[e4_inf_times$total_infs > 1,]
e4_reinfs <- ddply(e4_reinfs, ~sampno, nrow)
y <- (e4_reinfs$V1/max(n_alive)) * 100
e4_reinfs_est <- c(mean(y), quantile(y, c(0.025,0.5,0.975)))

## Now get attack rate for this time period
setkey(e4_inf_times,"sampno")
e4_inc <- e4_inf_times[,list(sum(any_inf)),by=key(e4_inf_times)]
e4_ests <- c(mean(e4_inc$V1/n_alive[1]), quantile(e4_inc$V1/n_alive[1], c(0.025,0.5,0.975)))


e_ests <- rbind(e1_ests, e2_ests, e3_ests, e4_ests)
e_reinf_ests <- rbind(e1_reinfs_est, e2_reinfs_est, e3_reinfs_est, e4_reinfs_est)
colnames(e_ests) <- colnames(e_reinf_ests) <- c("mean","lower_quantile","median","upper_Quantile")

####################################################
## Epidemic 2011
e_2011_infs <- inf_chain[inf_chain$j >= e_2011_start  - 1968*4 + 1 & inf_chain$j <= e_2011_end  - 1968*4 + 1,]
e_2011_times <- e_2011_start:e_2011_end - 1968*4 + 1
#_alive <- get_n_alive(titre_dat, e_2011_times)
n_alive <- n_alive_all[e_2011_times]

## How many times were people infected
setkey(e_2011_infs, "sampno","i")
e_2011_inf_times <- e_2011_infs[,list(total_infs=sum(x)),by=key(e_2011_infs)]

## Some people were infected twice, apparently
table(e_2011_inf_times$total_infs)
e_2011_inf_times$any_inf <- e_2011_inf_times$total_infs
e_2011_inf_times$any_inf[e_2011_inf_times$any_inf > 1] <- 1


## Get number of people that were reinfected
e_2011_reinfs <- e_2011_inf_times[e_2011_inf_times$total_infs > 1,]
e_2011_reinfs <- ddply(e_2011_reinfs, ~sampno, nrow)
y <- (e_2011_reinfs$V1/max(2011)) * 100
e_2011_reinfs_est <- c(mean(y), quantile(y, c(0.025,0.5,0.975)))

## Now get attack rate for this time period
setkey(e_2011_inf_times,"sampno")
e_2011_inc <- e_2011_inf_times[,list(sum(any_inf)),by=key(e_2011_inf_times)]
e_2011_ests <- c(mean(e_2011_inc$V1/n_alive[1]), quantile(e_2011_inc$V1/n_alive[1], c(0.025,0.5,0.975)))

e_ests <- rbind(e1_ests, e2_ests, e3_ests, e4_ests,e_2011_ests)
e_reinf_ests <- rbind(e1_reinfs_est, e2_reinfs_est, e3_reinfs_est, e4_reinfs_est,e_2011_reinfs_est)
colnames(e_ests) <- colnames(e_reinf_ests) <- c("mean","lower_quantile","median","upper_quantile")
e_ests <- data.frame(e_ests)
e_ests$period <- c(2010,2012,2013,2014,2011)  
ggplot(e_ests) + geom_pointrange(aes(x=period,y=median,ymin=lower_quantile,ymax=upper_quantile))
## Join infection status across epidemics
e1_inf_times$epidemic <- 1
e2_inf_times$epidemic <- 2
e3_inf_times$epidemic <- 3
e4_inf_times$epidemic <- 4
e_2011_inf_times$epidemic <- 5

if(by_year){
e_inf_times_all <- rbind(e1_inf_times, e2_inf_times, e3_inf_times, e4_inf_times, e_2011_inf_times)
} else {
  e_inf_times_all <- rbind(e1_inf_times, e2_inf_times, e3_inf_times, e4_inf_times)
}
setkey(e_inf_times_all, "i","sampno")
e_inf_times_all <- e_inf_times_all[,list(total_infs=sum(any_inf)),by=key(e_inf_times_all)]
res <- ddply(e_inf_times_all, ~sampno, function(x) c(sum(x$total_infs==0), sum(x$total_infs==1), 
                                                     sum(x$total_infs==2), sum(x$total_infs==3),
                                                     sum(x$total_infs==4), sum(x$total_infs==5)))
res[,2:ncol(res)] <- res[,2:ncol(res)]/n_alive[1]
inf_number <- t(apply(res[,2:ncol(res)],2, function(x) quantile(x,c(0.025,0.5,0.975))))
colnames(inf_number) <- c("lower","median","upper")
inf_number <- data.frame(inf_number)
inf_number$infections <- 0:5

res_total_infections <- ddply(signif(inf_number*100,3), ~infections, function(x) paste0(x$median,"% (", x$lower, "%-",x$upper,"%)"))

e_ests <- e_ests[order(e_ests$period),]
e_ests$period1 <- 1:nrow(e_ests)
res_attack_rates <- ddply(signif(e_ests*100,3), ~period1, function(x) paste0(x$median,"% (", x$lower_quantile, "%-",x$upper_quantile,"%)"))

n_alive <- get_n_alive(titre_dat,strain_isolation_times)
n_alive <- data.frame("n"=n_alive,"j"=seq_along(strain_isolation_times))

infection_histories <- inf_chain
setkey(infection_histories, j, sampno)
total_infs <- infection_histories[,list(ar=sum(x)),by=key(infection_histories)] 
total_infs <- merge(total_infs, n_alive,by="j")

total_infs$ar <- total_infs$ar/total_infs$n
y <- ddply(total_infs, ~sampno, function(x) median(x$ar))
median_ars <- signif(quantile(y$V1*100, c(0.5,0.025,0.975)), 3)
print("Median overall attack rate")
print(median_ars)

AR_ests <- total_infs %>% dplyr::group_by(j) %>% dplyr::summarize(median=median(ar),lower_quant=quantile(ar,0.025),upper_quant=quantile(ar,0.975)) %>% 
    mutate(year=strain_isolation_times[j]/4)

AR_ests %>% filter(median == min(median))
AR_ests %>% filter(median == max(median))
AR_ests %>% filter(year == 1985)

infection_histories$year <- floor(infection_histories$j/4)
setkey(infection_histories, year,i, sampno)
total_infs_year <- infection_histories[,list(total_infs=sum(x)),by=key(infection_histories)] 
y <- ddply(total_infs_year, ~sampno, function(x){
  n_reinfs_recent <- nrow(x[x$total_infs > 1 & x$year >= 2008-1968,])
  n_reinfs_all <- nrow(x[x$total_infs > 1,])
  n_reinfs_recent/n_reinfs_all
})
## What proportion of reinfections occured since 2008?
print("What proportion of reinfections occured since 2008?")
print(quantile(y$V1*100, c(0.5,0.025,0.975)))

e_reinf_ests <- as.data.frame(e_reinf_ests)
## Reinfected
res_reinfected <- paste0(paste0(signif(e_reinf_ests$median,3),"% (", 
              signif(e_reinf_ests$lower_quantile,3), "%-",
              signif(e_reinf_ests$upper_quantile,3),"%)"))
