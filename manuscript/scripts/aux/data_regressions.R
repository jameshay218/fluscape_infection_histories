## Get seropositivity
fluscape_dat$seropos_v1 <- fluscape_dat$v1_titre >= 3
fluscape_dat$seropos_v2 <- fluscape_dat$v2_titre >= 3

use_virus_year <- c(1968)
use_virus_year <- unique(fluscape_dat$Year)

## Get age at sample and age at circulation
dat_v2 <- fluscape_dat[fluscape_dat$visit == "Second visit" & 
                         fluscape_dat$Year %in% use_virus_year,]
dat_v2$age_raw <- (dat_v2$samples - dat_v2$DOB)/4
dat_v2$age_at_circ <- dat_v2$Year - dat_v2$DOB/4

## Perform some basic regressions
dat_v2$Year <- as.factor(dat_v2$Year)
dat_v2$LOC_ID <- as.factor(dat_v2$LOC_ID)
## Baseline, year intercept
fit_virus <- bam(data=dat_v2, v2_titre ~ Year,family=gaussian())
## Spline on age + year intercept
fit_virus_dob <- bam(data=dat_v2, v2_titre ~  Year + s(age_raw,by=Year),family=gaussian())
## Spline on age, age at circulation + year intercept
fit_virus_dob_age <- bam(data=dat_v2, v2_titre ~  Year + s(age_raw,by = Year) + s(age_at_circ,by = Year),family=gaussian())
## Spline on age, age at circulation + year + location intercept
fit_virus_dob_age_loc <- bam(data=dat_v2, v2_titre ~  Year + LOC_ID + s(age_raw,by=Year) + s(age_at_circ) + s(age_at_circ,by=Year),family=gaussian())
## Spline on age, age at circulation + year + location intercept + visit intercept
fit_virus_dob_age_loc_visit <- bam(data=dat_v2, v2_titre ~  Year + LOC_ID + as.factor(raw_visit) + s(age_raw,by=Year) + s(age_at_circ) + s(age_at_circ,by=Year),family=gaussian())
## Spline on density, age, age at circulation + year + location intercept
fit_virus_dob_age_dens <- bam(data=dat_v2, v2_titre ~  as.factor(Year) + s(age_raw,by=Year) + 
                                s(age_at_circ,by=Year) + s(dens.9,by=Year),family=gaussian())
## Spline on travel time to GZ, age, age at circulation + year + location intercept
fit_virus_dob_age_trav <- bam(data=dat_v2, v2_titre ~  as.factor(Year) + s(age_raw,by=Year) + 
                                s(age_at_circ) + s(trav.min),family=gaussian())


## Compare deviance explained, AIC and BIC
dev_explained <- c(summary(fit_virus)[["dev.expl"]]*100,
                   summary(fit_virus_dob)[["dev.expl"]]*100,
                   summary(fit_virus_dob_age)[["dev.expl"]]*100,
                   summary(fit_virus_dob_age_loc)[["dev.expl"]]*100,
                   summary(fit_virus_dob_age_dens)[["dev.expl"]]*100,
                   summary(fit_virus_dob_age_loc_visit)[["dev.expl"]]*100,
                   summary(fit_virus_dob_age_trav)[["dev.expl"]]*100)
aics <- c(AIC(fit_virus),
          AIC(fit_virus_dob),
          AIC(fit_virus_dob_age),
          AIC(fit_virus_dob_age_loc),
          AIC(fit_virus_dob_age_dens),
          AIC(fit_virus_dob_age_loc_visit),
          AIC(fit_virus_dob_age_trav))
bics <- c(BIC(fit_virus),
          BIC(fit_virus_dob),
          BIC(fit_virus_dob_age),
          BIC(fit_virus_dob_age_loc),
          BIC(fit_virus_dob_age_dens),
          BIC(fit_virus_dob_age_loc_visit),
          BIC(fit_virus_dob_age_trav)
          )
model_names <- as.character(c(fit_virus$formula, fit_virus_dob$formula, fit_virus_dob_age$formula,
                 fit_virus_dob_age_loc$formula,fit_virus_dob_age_dens$formula, fit_virus_dob_age_loc_visit$formula,
                 fit_virus_dob_age_trav$formula))
raw_res <- data.frame("name"=model_names,"deviance"=dev_explained,"aic"=aics,"bic"=bics)

## Regressions on seroconversion
## Get ordered distance from Guangzhou
loc_distances <- unique(fluscape_dat[,c("LOC_ID","DIST_FRM_GZ")])
head(loc_distances)
loc_distances <- loc_distances[order(loc_distances$DIST_FRM_GZ),]
dat <- fluscape_dat[complete.cases(fluscape_dat),]
dat$LOC_ID <- factor(dat$LOC_ID, levels=loc_distances$LOC_ID,ordered=TRUE)

## Use change in titre as outcome
dat <- dat[dat$visit == "Second visit",]

dat$age_raw <- (dat$samples - dat$DOB)/4
dat$age_at_circ <- dat$Year -dat$DOB/4 
dat$Year <- as.factor(dat$Year)

## Same regressions as above
fit_ct_virus <- bam(data=dat, change ~ as.factor(Year),family=gaussian())
fit_ct_virus_dob <- bam(data=dat, change ~  as.factor(Year) + s(age_raw,by=Year),family=gaussian())
fit_ct_virus_dob_age <- bam(data=dat, change ~  as.factor(Year) + s(age_raw,by=Year) + 
                                s(age_at_circ,by=Year),family=gaussian())
fit_ct_virus_dob_age_loc <- bam(data=dat, change ~  as.factor(Year) + s(age_raw,by=Year) + 
                               s(age_at_circ,by=Year) + as.factor(LOC_ID),family=gaussian())
fit_ct_virus_dob_age_loc_visit <- bam(data=dat, change ~  as.factor(Year) + s(age_raw,by=Year) + as.factor(raw_visit) +
                                    s(age_at_circ,by=Year) + as.factor(LOC_ID),family=gaussian())
fit_ct_virus_dob_age_dens <- bam(data=dat, change ~  as.factor(Year) + s(age_raw,by=Year) + 
                                  s(age_at_circ,by=Year) + s(dens.9,by=Year),family=gaussian())
fit_ct_virus_dob_age_trav <- bam(data=dat, change ~  as.factor(Year) + s(age_raw,by=Year) + 
                                   s(age_at_circ,by=Year) + s(trav.min,by=Year),family=gaussian())

dev_explained <- c(summary(fit_ct_virus)[["dev.expl"]]*100,
                   summary(fit_ct_virus_dob)[["dev.expl"]]*100,
                   summary(fit_ct_virus_dob_age)[["dev.expl"]]*100,
                   summary(fit_ct_virus_dob_age_loc)[["dev.expl"]]*100,
                   summary(fit_ct_virus_dob_age_dens)[["dev.expl"]]*100,
                   summary(fit_ct_virus_dob_age_loc_visit)[["dev.expl"]]*100,
                   summary(fit_ct_virus_dob_age_trav)[["dev.expl"]]*100)
aics <- c(AIC(fit_ct_virus),
          AIC(fit_ct_virus_dob),
          AIC(fit_ct_virus_dob_age),
          AIC(fit_ct_virus_dob_age_loc),
          AIC(fit_ct_virus_dob_age_dens),
          AIC(fit_ct_virus_dob_age_loc_visit),
          AIC(fit_ct_virus_dob_age_trav))
bics <- c(BIC(fit_ct_virus),
          BIC(fit_ct_virus_dob),
          BIC(fit_ct_virus_dob_age),
          BIC(fit_ct_virus_dob_age_loc),
          BIC(fit_ct_virus_dob_age_dens),
          BIC(fit_ct_virus_dob_age_loc_visit),
          BIC(fit_ct_virus_dob_age_trav))
model_names <- as.character(c(fit_ct_virus$formula, fit_ct_virus_dob$formula, fit_ct_virus_dob_age$formula,
                              fit_ct_virus_dob_age_loc$formula,
                              fit_ct_virus_dob_age_loc_visit$formula,
                              fit_ct_virus_dob_age_dens$formula, 
                              fit_ct_virus_dob_age_trav$formula))

ct_res <- data.frame("name"=model_names,"deviance"=dev_explained,"aic"=aics,"bic"=bics)

ct_res <- ct_res %>% mutate(delta_aic = aic - min(aic), delta_bic = bic - min(bic)) %>% arrange(delta_aic)
raw_res <- raw_res %>% mutate(delta_aic = aic - min(aic), delta_bic = bic - min(bic))%>% arrange(delta_aic)

write.table(raw_res, paste0(figure_wd, "/data_v2_regression.csv"),sep=",",row.names=FALSE)
write.table(ct_res, paste0(figure_wd, "/data_change_regression.csv"),sep=",",row.names=FALSE)


## Plots
## V2 titres
dat_new_v2_raw <- expand_grid(DOB=seq(2014-80,2014,by=1), 
                              LOC_ID=levels(dat_v2$LOC_ID)[1],
                              Year=unique(dat_v2$Year)) %>% 
    mutate(age_raw=2014-DOB,age_at_circ=as.numeric(as.character(Year))-DOB)
tmp <- predict(fit_virus_dob_age_loc, newdata=dat_new_v2_raw,se.fit=TRUE)
dat_new_v2_raw$prediction <- tmp$fit
dat_new_v2_raw$upper <- tmp$fit + (2*tmp$se.fit)
dat_new_v2_raw$lower <- tmp$fit - (2*tmp$se.fit)
dat_new_v2_raw$model <- "HI titre at second sample"

## Change titres
dat_new_v2_ct <- expand_grid(DOB=seq(2014-80,2014,by=1), 
                             LOC_ID=levels(dat_v2$LOC_ID)[1],
                             Year=unique(dat_v2$Year)) %>% 
    mutate(age_raw=2014-DOB,age_at_circ=as.numeric(as.character(Year))-DOB)

tmp <- predict(fit_ct_virus_dob_age_loc, newdata=dat_new_v2_ct,se.fit=TRUE)
dat_new_v2_ct$prediction <- tmp$fit
dat_new_v2_ct$upper <- tmp$fit + (2*tmp$se.fit)
dat_new_v2_ct$lower <- tmp$fit - (2*tmp$se.fit)
dat_new_v2_raw$model <- "Change in HI titre"

## Combine
dat_all <- bind_rows(dat_new_v2_raw, dat_new_v2_ct)

virus_key <- fluscape_dat %>% select(Full_name, Year) %>% distinct()
virus_key1 <- virus_key$Full_name
names(virus_key1) <- virus_key$Year
dat_new_v2_raw$Year1 <- virus_key1[dat_new_v2_raw$Year]
dat_new_v2_ct$Year1 <- virus_key1[dat_new_v2_ct$Year]

theme_consistent <-  theme(axis.text.x=element_text(colour="black",size=6),
                           axis.text.y=element_text(family="sans",colour="black",size=6),
                           axis.title = element_text(family="sans",size=7),
                           strip.text=element_text(size=5,family="sans"),
                           legend.text=element_text(size=6,family="sans"),
                           panel.spacing=unit(1,"lines"))


## Raw titres
p_gam_raw_age <- ggplot(dat_new_v2_raw) + 
    geom_line(aes(x=age_raw,y=prediction,col=Year1)) +
    guides(color=guide_legend(ncol=1)) +
    scale_color_viridis_d(name="Strain") +
    ylab("Predicted log HI titre") +
    xlab("Age at sampling, years") +
    coord_cartesian(ylim=c(0,8)) +
    scale_y_continuous(breaks=seq(0,8,by=2)) +
    theme_minimal() +
    theme_consistent + 
    theme(legend.position="none") +
    labs(tag="A") +
    guides(color=guide_legend(ncol=1))

p_gam_raw_circ <- ggplot(dat_new_v2_raw) + 
    geom_line(aes(x=age_at_circ,y=prediction,col=Year1))  +
    scale_color_viridis_d(name="Strain") +
    guides(color=guide_legend(ncol=1)) +
    ylab("Predicted log HI titre") +
    xlab("Age at isolation of strain, years") +
    coord_cartesian(ylim=c(0,8)) +
    scale_y_continuous(breaks=seq(0,8,by=2)) +
    theme_minimal() +
    theme_consistent +
    theme(legend.position="none") +
    labs(tag="A") +
    guides(color=guide_legend(ncol=1))

p_gam_raw_circ_indiv <- ggplot(dat_new_v2_raw) + 
    geom_ribbon(aes(x=age_at_circ,ymin=lower,ymax=upper,fill=Year1),alpha=0.25)+ 
    geom_line(aes(x=age_at_circ,y=prediction,col=Year1)) +
    guides(color=guide_legend(ncol=1),fill=guide_legend(ncol=1)) +
    scale_fill_viridis_d(name="Strain") +
    scale_color_viridis_d(name="Strain") +
    ylab("Predicted log HI titre") +
    xlab("Age at isolation of strain, years") +
    coord_cartesian(ylim=c(0,8)) +
    scale_y_continuous(breaks=seq(0,8,by=2)) +
    theme_minimal() +
    theme_consistent +
    theme(legend.position="right") +
    facet_wrap(~Year) +
    labs(tag="B") +
    guides(color=guide_legend(ncol=1),fill=guide_legend(ncol=1))

p_gam_raw_age_indiv <- ggplot(dat_new_v2_raw) + 
    geom_ribbon(aes(x=age_raw,ymin=lower,ymax=upper,fill=Year1),alpha=0.25)+ 
    geom_line(aes(x=age_raw,y=prediction,col=Year1)) +
    guides(color=guide_legend(ncol=1),fill=guide_legend(ncol=1)) +
    scale_fill_viridis_d(name="Strain") +
    scale_color_viridis_d(name="Strain") +
    ylab("Predicted log HI titre") +
    xlab("Age at sampling, years") +
    coord_cartesian(ylim=c(0,8)) +
    scale_y_continuous(breaks=seq(0,8,by=2)) +
    theme_minimal() +
    theme_consistent +
    theme(legend.position="right") +
    facet_wrap(~Year)+
    labs(tag="B")+
    guides(color=guide_legend(ncol=1),fill=guide_legend(ncol=1))

## Change in titres
p_gam_ct_age <- ggplot(dat_new_v2_ct) + 
    geom_line(aes(x=age_raw,y=prediction,col=Year1)) +
    scale_color_viridis_d(name="Strain") +
    guides(color=guide_legend(ncol=1)) +
    ylab("Predicted change in log HI titre") +
    xlab("Age at sampling, years") +
    coord_cartesian(ylim=c(0,8)) +
    scale_y_continuous(breaks=seq(0,8,by=2)) +
    theme_minimal() +
    theme_consistent +
    theme(legend.position="none") +
    labs(tag="A")+
    guides(color=guide_legend(ncol=1))

p_gam_ct_circ <- ggplot(dat_new_v2_ct) + 
    geom_line(aes(x=age_at_circ,y=prediction,col=Year1))  +
    guides(color=guide_legend(ncol=1)) +
    scale_color_viridis_d(name="Strain") +
    ylab("Predicted change in log HI titre") +
    xlab("Age at isolation of strain, years") +
    coord_cartesian(ylim=c(0,8)) +
    scale_y_continuous(breaks=seq(0,8,by=2)) +
    theme_minimal() +
    theme_consistent +
    theme(legend.position="none") +
    labs(tag="A")

p_gam_ct_circ_indiv <- ggplot(dat_new_v2_ct) + 
    geom_ribbon(aes(x=age_at_circ,ymin=lower,ymax=upper,fill=Year1),alpha=0.25)+ 
    geom_line(aes(x=age_at_circ,y=prediction,col=Year1)) +
    guides(color=guide_legend(ncol=1),fill=guide_legend(ncol=1)) +
    scale_fill_viridis_d(name="Strain") +
    scale_color_viridis_d(name="Strain") +
    ylab("Predicted change in log HI titre") +
    xlab("Age at isolation of strain, years") +
    coord_cartesian(ylim=c(0,8)) +
    scale_y_continuous(breaks=seq(0,8,by=2)) +
    theme_minimal() +
    theme_consistent +
    theme(legend.position="right") +
    facet_wrap(~Year) +
    labs(tag="B")

p_gam_ct_age_indiv <- ggplot(dat_new_v2_ct) + 
    geom_ribbon(aes(x=age_raw,ymin=lower,ymax=upper,fill=Year1),alpha=0.25)+ 
    geom_line(aes(x=age_raw,y=prediction,col=Year1)) +
    guides(color=guide_legend(ncol=1),fill=guide_legend(ncol=1)) +
    scale_fill_viridis_d(name="Strain",guide=guide_legend(ncol=1)) +
    scale_color_viridis_d(name="Strain",guide=guide_legend(ncol=1)) +
    ylab("Predicted change in log HI titre") +
    xlab("Age at sampling, years") +
    coord_cartesian(ylim=c(0,8)) +
    scale_y_continuous(breaks=seq(0,8,by=2)) +
    theme_minimal() +
    theme_consistent +
    theme(legend.position="right") +
    facet_wrap(~Year)+
    labs(tag="B")

p1 <- p_gam_raw_age + p_gam_raw_age_indiv + plot_layout(ncol=1,heights=c(1,2.5),guides = "collect")
p2 <- p_gam_raw_circ + p_gam_raw_circ_indiv + plot_layout(ncol=1,heights=c(1,2.5),guides = "collect")
p3 <- p_gam_ct_age + p_gam_ct_age_indiv + plot_layout(ncol=1,heights=c(1,2.5),guides = "collect")
p4 <- p_gam_ct_circ + p_gam_ct_circ_indiv + plot_layout(ncol=1,heights=c(1,2.5),guides = "collect")

p_alt <- (p_gam_raw_age+labs(tag="A") ) + (p_gam_raw_circ+labs(tag="B")+ theme(legend.position="right")) + 
    (p_gam_ct_age+labs(tag="C")) + (p_gam_ct_circ+labs(tag="D")) + plot_layout(ncol=2,guides="collect",byrow=FALSE)

ggsave_jah(p1, wd = figure_wd,save_name = "v2_titre_regression_age",height=7,width=7)
ggsave_jah(p2, wd = figure_wd,save_name = "v2_titre_regression_agecirc",height=7,width=7)
ggsave_jah(p3, wd = figure_wd,save_name = "change_titre_regression_age",height=7,width=7)
ggsave_jah(p4, wd = figure_wd,save_name = "change_titre_regression_agecirc",height=7,width=7)
ggsave_jah(p_alt, wd = figure_wd,save_name = "gams_combined",height=6,width=8)
