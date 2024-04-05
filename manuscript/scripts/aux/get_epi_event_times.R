## Find if individual had a fever/ILI in the past year
tmp_dat <- part_info[part_info$variable %in% c("PART_FEV_ILI_LAST_YR","PART_SYMPT_1YR_MULTIPLE") & !is.na(part_info$value), ]
ili_last_year <- ddply(tmp_dat, .(Participant_ID, visit), function(tmp){ 
  x <-  unique(tmp[,colnames(tmp) == unique(tmp$visit)])
  fever_end <- x
  if(length(x) > 0 && !is.na(fever_end)){
    fever_start <- fever_end - buckets
    fever_presence <- tmp[tmp$variable %in% c("PART_FEV_ILI_LAST_YR","PART_SYMPT_1YR_MULTIPLE"),"value"] 
  } else {
    fever_start <- NA
    fever_presence <- NA
  }
  data.frame(fever_end=fever_end, fever_start=fever_start, fever_presence=fever_presence)
})
 
## Find if participant smoked in the last month ever (ie. ever smoker)
tmp_dat <- part_info[part_info$variable == "PART_SMOKE_LAST_MO" & !is.na(part_info$value), ]
smoking <- ddply(tmp_dat, .(Participant_ID, visit), function(tmp){ 
  x <-  unique(tmp[,colnames(tmp) == unique(tmp$visit)])
  smoke_end <- x
  if(length(x) > 0 && !is.na(smoke_end)){
    smoke_start <- smoke_end - buckets/12
    smoke_presence <- tmp[tmp$variable == "PART_SMOKE_LAST_MO","value"] 
  } else {
    smoke_start <- NA
    smoke_presence <- NA
  }
  data.frame(smoke_end=smoke_end, smoke_start=smoke_start, smoke_presence=smoke_presence)
})
smoking_all <- ddply(smoking, .(Participant_ID), function(x) any(x$smoke_presence == 1))
smoking_all <- smoking_all[complete.cases(smoking_all),]


## Find if individual was vaccinated in the past year
tmp_dat <- part_info[part_info$variable == "PART_VAC_LAST_YR" & !is.na(part_info$value), ]
vaccine_last_year <- ddply(tmp_dat, .(Participant_ID, visit), function(tmp){ 
  x <-  unique(tmp[,colnames(tmp) == unique(tmp$visit)])
  vaccine_end <- x
  if(length(x) > 0 && !is.na(vaccine_end)){
    vaccine_start <- vaccine_end - buckets
    vaccine_presence <- tmp[tmp$variable == "PART_VAC_LAST_YR","value"] 
  } else {
    vaccine_start <- NA
    vaccine_presence <- NA
  }
  data.frame(vaccine_end=vaccine_end, vaccine_start=vaccine_start, vaccine_presence=vaccine_presence)
})

## Find if fev/ili since last visit
tmp_dat <- part_info[part_info$variable %in% c("PART_FEV_ILI_LAST_VISIT","PART_SYMPT_LSV_MULTIPLE") & !is.na(part_info$value), ]
ili_last_visit <- ddply(tmp_dat, .(Participant_ID, visit), function(tmp){ 
  x <-  unique(tmp[,colnames(tmp) == unique(tmp$visit)])
  ili_end <- x
  if(length(x) > 0 && !is.na(ili_end)){
    ili_start <- unique(tmp[,"last_visit_time"])
    ili_presence <- tmp[tmp$variable %in% c("PART_FEV_ILI_LAST_VISIT","PART_SYMPT_LSV_MULTIPLE"),"value"] 
  } else {
    ili_start <- NA
    ili_presence <- NA
  }
  data.frame(ili_end=ili_end, ili_start=ili_start, ili_presence=ili_presence)
})


## Find if vac since last visit
## Find if individual was vaccinated in the past year
tmp_dat <- part_info[part_info$variable == "PART_VAC_LAST_VISIT" & !is.na(part_info$value), ]
vaccine_last_visit <- ddply(tmp_dat, .(Participant_ID, visit), function(tmp){ 
  x <-  unique(tmp[,colnames(tmp) == unique(tmp$visit)])
  vaccine_end <- x
  if(length(x) > 0 && !is.na(vaccine_end)){
    vaccine_start <- unique(tmp[,"last_visit_time"])
    vaccine_presence <- tmp[tmp$variable == "PART_VAC_LAST_VISIT","value"] 
  } else {
    vaccine_start <- NA
    vaccine_presence <- NA
  }
  data.frame(vaccine_end=vaccine_end, vaccine_start=vaccine_start, vaccine_presence=vaccine_presence)
})




## Find if individual has ever been vaccinated
## 1=never 2=last time this year 3=last onelast year 4=last one was 2-5years 5= last one more than 5years ago 9=unsure 888=decline
## 1=never 2=last time this year 3=last onelast year 4=last one was 2-5years 5= last one more than 5years ago 9=unsure 888=decline
tmp_dat <- part_info[part_info$variable == "PART_LAST_VAC" & !is.na(part_info$value), ]

vaccine_ever <- tmp_dat %>% ungroup() %>% select(Participant_ID,visit,value,variable)%>% 
  filter(!(value %in% c(9, 888))) %>%
  group_by(visit,Participant_ID) %>% mutate(any_vacc = as.numeric(value %in% c(2,3,4,5)),                                                    
                                       never_vacc = as.numeric(value == 1),
                                       recent_vacc = as.numeric(value %in% c(2,3))) %>%
  select(Participant_ID, visit, any_vacc,never_vacc,recent_vacc)

## Merge epi events together
all_epi_events <- 
  expand_grid(visit = c("V1","V2","V3","V4"), Participant_ID = unique(part_info$Participant_ID)) %>% 
  left_join(vaccine_ever) %>%
  left_join(vaccine_last_visit %>% 
              rename(vaccine_last_visit = vaccine_presence, vaccine_last_visit_start = vaccine_start, vaccine_last_visit_end = vaccine_end) %>%
              ungroup()) %>%
  left_join(ili_last_visit %>% rename(ili_last_visit = ili_presence, ili_last_visit_start = ili_start, ili_last_visit_end = ili_end)) %>%
  left_join(vaccine_last_year %>%rename(vaccine_last_year = vaccine_presence, vaccine_last_year_start = vaccine_start, vaccine_last_year_end = vaccine_end)) %>%
  left_join(ili_last_year %>%rename(ili_last_year = fever_presence, ili_last_year_start = fever_start, ili_last_year_end = fever_end))%>% 
  arrange(Participant_ID, visit) 

all_epi_events <- all_epi_events %>%
  left_join(part_info %>% select(Participant_ID,V1,V2,V3,V4)  %>% pivot_longer(-Participant_ID) %>% rename(date=value, visit=name) %>% distinct())

## Some participants gave conflicting responses/there was some data input error
all_epi_events %>% filter(never_vacc == 1 & (vaccine_last_year == 1 | vaccine_last_visit == 1))

## Find out if individual has ever had one of these diseases
ever_had_heart_disease <- ddply(part_info, ~Participant_ID, function(x) any(x[x$variable == "PART_DX_HEART_DISEASE","value"] == 1))
colnames(ever_had_heart_disease)[2] <- c("heart_dx")
ever_had_diabetes_disease <- ddply(part_info, ~Participant_ID, function(x) any(x[x$variable == "PART_DX_DIABETES","value"] == 1))
colnames(ever_had_diabetes_disease)[2]  <- c("diabetes_dx")
ever_had_hypertension_disease <- ddply(part_info, ~Participant_ID, function(x) any(x[x$variable == "PART_DX_HYPERTENSION","value"] == 1))
colnames(ever_had_hypertension_disease)[2]  <- c("hypertension_dx")
ever_had_cancer_disease <- ddply(part_info, ~Participant_ID, function(x) any(x[x$variable == "PART_DX_CANCER","value"] == 1))
colnames(ever_had_cancer_disease)[2]  <- c("cancer_dx")
ever_had_lung_disease <- ddply(part_info, ~Participant_ID, function(x) any(x[x$variable == "PART_DX_CHRONIC_LUNG","value"] == 1))
colnames(ever_had_lung_disease)[2]  <- c("lung_dx")
ever_had_allergies_disease <- ddply(part_info, ~Participant_ID, function(x) any(x[x$variable == "PART_DX_ALLERGIES","value"] == 1))
colnames(ever_had_allergies_disease)[2]  <- c("allergies_dx")

diseases <- merge(ever_had_heart_disease, ever_had_diabetes_disease)
diseases <- merge(diseases, ever_had_hypertension_disease)
diseases <- merge(diseases, ever_had_cancer_disease)
diseases <- merge(diseases, ever_had_lung_disease)
diseases <- merge(diseases, ever_had_allergies_disease)
