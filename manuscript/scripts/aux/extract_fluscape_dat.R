get_datetime <- TRUE

## Get general utility functions from fluscape SVN
source(paste0(fluscape_wd,"manuscripts/paired_titer/rCodes/GeneralUtility.r"))
source(paste0(fluscape_wd,"manuscripts/paired_titer/rCodes/GeneralUtility_BY.r"))

## The data used for MCMC fitting
fluscape_dat1 <- read.csv("~/Documents/GitHub/fluscape_serosolver/data/fluscape_data_4_resolution.csv",stringsAsFactors = FALSE)

## Strains used
strains <- read.csv(paste0(main_wd, "data/strains.csv"),stringsAsFactors=FALSE)
strains$Full_name <- factor(strains$Full_name, levels=strains$Full_name,ordered=TRUE)

## Just to double check the data is right, look at BY/JL's other analysis and make sure the IDs match
other_titre_dat <- read.csv(paste0(fluscape_wd,"manuscripts/paired_titer/data/HI_titers_paired_R56Pilot.csv"),stringsAsFactors = FALSE)
ids_should_have <- unique(other_titre_dat$Participant_ID)


###################################################################
## 1. LOAD IN VISIT DATA
###################################################################
message("Load in visit data")
na_strings <- c("\\N","999")
#na_strings< -"\\N"
v1 <- read.table(paste0(fluscape_wd, "data/Participants_V1.csv"),header=TRUE,
                 stringsAsFactors = FALSE,sep=",",na.strings=na_strings)
v2 <- read.table(paste0(fluscape_wd, "data/Participants_V2.csv"),header=TRUE,
                 stringsAsFactors = FALSE,sep=",",na.strings=na_strings)
v3 <- read.table(paste0(fluscape_wd, "data/Participants_V3.csv"),header=TRUE,
                 stringsAsFactors = FALSE,sep=",",na.strings=na_strings)
v4 <- read.table(paste0(fluscape_wd, "data/Participants_V4.csv"),header=TRUE,
                 stringsAsFactors = FALSE,sep=",",na.strings=na_strings)

## If no blood sample, use question time as visit time
v1[is.na(v1$PART_SAMPLE_TIME),"PART_SAMPLE_TIME"] <- v1[is.na(v1$PART_SAMPLE_TIME),"PART_QuestTimeEnd"]
v2[is.na(v2$PART_SAMPLE_TIME),"PART_SAMPLE_TIME"] <- v2[is.na(v2$PART_SAMPLE_TIME),"PART_QuestTimeEnd"]
v3[is.na(v3$PART_SAMPLE_TIME),"PART_SAMPLE_TIME"] <- v3[is.na(v3$PART_SAMPLE_TIME),"PART_QuestTimeEnd"]
v4[is.na(v4$PART_SAMPLE_TIME),"PART_SAMPLE_TIME"] <- v4[is.na(v4$PART_SAMPLE_TIME),"PART_QuestTimeEnd"]
v1[is.na(v1$PART_SAMPLE_TIME),"PART_SAMPLE_TIME"] <- v1[is.na(v1$PART_SAMPLE_TIME),"PART_QuestTimeStart"]
v2[is.na(v2$PART_SAMPLE_TIME),"PART_SAMPLE_TIME"] <- v2[is.na(v2$PART_SAMPLE_TIME),"PART_QuestTimeStart"]
v3[is.na(v3$PART_SAMPLE_TIME),"PART_SAMPLE_TIME"] <- v3[is.na(v3$PART_SAMPLE_TIME),"PART_QuestTimeStart"]
v4[is.na(v4$PART_SAMPLE_TIME),"PART_SAMPLE_TIME"] <- v4[is.na(v4$PART_SAMPLE_TIME),"PART_QuestTimeStart"]

v1$Part_ID <- sprintf("L%02dH%02dP%02d",v1$LOC_ID,v1$HH_ID,v1$PARTICIPANT_ID)
v2$Part_ID <- sprintf("L%02dH%02dP%02d",v2$LOC_ID,v2$HH_ID,v2$PARTICIPANT_ID)
v3$Part_ID <- sprintf("L%02dH%02dP%02d",v3$LOC_ID,v3$HH_ID,v3$PARTICIPANT_ID)
v4$Part_ID <- sprintf("L%02dH%02dP%02d",v4$LOC_ID,v4$HH_ID,v4$PARTICIPANT_ID)

## Check which participant IDs are missing from the visit data, but present in the ID data
## Have all IDs
## Pass sanity check
titre_ids <- unique(fluscape_dat1$Participant_ID)
any(!(titre_ids %in% ids_should_have))
all_ids <- unique(c(v1$Part_ID, v2$Part_ID, v3$Part_ID, v4$Part_ID))
missing_ids <- titre_ids[which(!(titre_ids %in% all_ids))]
if(length(missing_ids) > 0){
  message("These IDs are missing from the visit data, but have titre data:")
  message(cat(missing_ids,sep="\t"))
}
###################################################################

###################################################################
## 2. TIDY AND COMBINE VISIT DATA
###################################################################
message("Tidy and combine visit data")
needed_names <- c("PARTICIPANT_ID","HH_ID","LOC_ID","PART_SAMPLE_TIME","PART_GENDER","PART_OCC_STATUS")
###################################################################
## V1 data

## V1 vars
## Gender; ILI in past month; Vaccine never/this year/last year/ 2-5 years/ 5+ years/unsure
extra_vars_v1 <- c("PART_FEV_ILI", "PART_LAST_VAC")
needed_names_v1 <- needed_names[!(needed_names %in% c("PART_BIRTH_YEAR","PART_BIRTH_MONTH"))]
v1 <- v1[,c(needed_names_v1,extra_vars_v1)]
v1$PART_SAMPLE_TIME <- as.Date(v1$PART_SAMPLE_TIME)
v1 <- melt(v1, id.vars=needed_names_v1)
v1$visit <- "V1"
###################################################################
###################################################################
## V2 data
## V2 vars
## Gender; ILI in past year; ILI since last visit;  ILI in last month; 
## Vaccine never/this year/last year/ 2-5 years/ 5+ years/unsure;
## Vcacine since last visit; Smoke last month; Vaccine last year
extra_vars_v2 <- c("PART_FEV_ILI_LAST_YR", "PART_FEV_ILI_LAST_VISIT","PART_LAST_VAC",
                   "PART_VAC_LAST_VISIT","PART_VAC_LAST_YR", "PART_SMOKE_LAST_MO")

v2 <- v2[,c(needed_names,extra_vars_v2)]
v2$PART_SAMPLE_TIME <- as.Date(v2$PART_SAMPLE_TIME)
v2 <- melt(v2, needed_names)
v2$visit <- "V2"
###################################################################
###################################################################
## V3/V4 data
## V3/V4 vars
## Gender; Height; Weight; 
## Heart disease; Diabetes; Hypertension; Cancer; Chronic Lung disease; Allergies;
## ILI fev in last month; multiple flu symptoms in last visit;
## Vaccine never/this year/last year/ 2-5 years/ 5+ years/unsure;
## Vcacine since last visit; Smoke last month; Vaccine last year
## Smoked in last month; Touch chickens; Touch ducks/geese; Touch pigs
extra_vars_v4 <- extra_vars_v3 <- c("PART_HEIGHT","PART_WEIGHT",
                                    "PART_DX_HEART_DISEASE" , "PART_DX_DIABETES", "PART_DX_HYPERTENSION", 
                                    "PART_DX_CANCER", "PART_DX_CHRONIC_LUNG","PART_DX_ALLERGIES",
                                    "PART_SYMPT_1YR_MULTIPLE",
                                    "PART_SYMPT_LSV_MULTIPLE",
                                    "PART_SMOKE_LAST_MO",
                                    "PART_CONT_CHICKEN_FREQ","PART_CONT_H2OFOWL_FREQ","PART_CONT_PIG_FREQ",
                                    "PART_LAST_VAC",
                                    "PART_VAC_LAST_VISIT","PART_VAC_LAST_YR")

v3 <- v3[,c(needed_names,extra_vars_v3)]
v3$PART_SAMPLE_TIME <- as.Date(v3$PART_SAMPLE_TIME)
v3 <- melt(v3, id.vars=needed_names)
v3$visit <- "V3"

## V4 is basically the same
needed_names_v4 <- c(needed_names,"PART_BIRTH_YEAR","PART_BIRTH_MONTH")
v4 <- v4[,c(needed_names_v4,extra_vars_v4)]
v4$PART_SAMPLE_TIME <- as.Date(v4$PART_SAMPLE_TIME)
#colnames(v4)[colnames(v4) %in% extra_vars_v4] <- paste0(colnames(v4)[colnames(v4) %in% extra_vars_v4],".V4")
v4 <- melt(v4, id.vars=needed_names_v4)
#colnames(v4)[colnames(v4) == "PART_SAMPLE_TIME"] <- "time"
v4$visit <- "V4"

#########################
## Now merge everything
DOBs <- unique(v4[,c("PARTICIPANT_ID", "HH_ID", "LOC_ID","PART_BIRTH_YEAR","PART_BIRTH_MONTH")])
v4 <- v4[,!(colnames(v4) %in% c("PART_BIRTH_YEAR","PART_BIRTH_MONTH"))]

part_info <- merge(v1, v2, id.vars=c("PARTICIPANT_ID","HH_ID","LOC_ID","PART_GENDER"),all=TRUE)
part_info <- merge(part_info,v3, id.vars=c("PARTICIPANT_ID","HH_ID","LOC_ID","PART_GENDER"),all=TRUE)
part_info <- merge(part_info,v4, id.vars=c("PARTICIPANT_ID","HH_ID","LOC_ID","PART_GENDER"),all=TRUE)
part_info <- merge(part_info, DOBs,all=TRUE)
PARTICIPANT_ID <- sprintf("L%02dH%02dP%02d",part_info$LOC_ID,part_info$HH_ID,part_info$PARTICIPANT_ID)
part_info$Participant_ID <- PARTICIPANT_ID

## Occupation status is whatever their status was at latest sample time
part_info1 <- ddply(part_info, ~Participant_ID, function(x){
  occs <- unique(x$PART_OCC_STATUS)
  x$PART_OCC_STATUS <- occs[length(occs)]
})
colnames(part_info1)[2] <- "TRUE_OCC_STATUS"
part_info <- merge(part_info, part_info1)

######
message("Fix gender entry")
## Get unique gender, as changes a bit due to entry error. Take the first gender recorded or NA
genders <- ddply(part_info, ~Participant_ID, function(x) {
  gender <- unique(x$PART_GENDER)
  gender <- gender[!is.na(gender)]
  res <- NA
  if(length(gender) != 0){
    res <- gender[1]
  }
  res
})
colnames(genders)[2] <- "PART_GENDER"

for(indiv in unique(part_info$Participant_ID)){
  actual_gender <- unique(genders[genders$Participant_ID == indiv,"PART_GENDER"])
  part_info[part_info$Participant_ID == indiv, "PART_GENDER"] <- actual_gender
}

## Dcast by sample time and convert back to date so that each row has all visit dates
samp_times <- unique(part_info[,c("Participant_ID","PART_SAMPLE_TIME","visit")])
samp_times <- dcast(samp_times, Participant_ID ~ visit, value.var="PART_SAMPLE_TIME")
samp_times <- data.frame(samp_times)
samp_times$V1 <- as.Date(samp_times$V1, origin="1970-01-01")
samp_times$V2 <- as.Date(samp_times$V2, origin="1970-01-01")
samp_times$V3 <- as.Date(samp_times$V3, origin="1970-01-01")
samp_times$V4 <- as.Date(samp_times$V4, origin="1970-01-01")

## Merge back together
part_info <- part_info[,colnames(part_info) != "PART_SAMPLE_TIME"]
part_info <- merge(part_info, samp_times)

## For each visit, find time of the previous visit if there, NA otherwise
visits <- c(NA, "V1","V2","V3","V4")
last_visit_times <- ddply(part_info, .(Participant_ID, visit), function(x){
  which_visit <- which(visits == unique(x$visit))
  last_visit <- visits[which_visit-1]
  last_visit_time <- NA
  if(!is.na(last_visit)){
    last_visit_time <- unique(x[,colnames(x) == last_visit])
  }
  last_visit_time
} )
last_visit_times$V1 <- as.Date(last_visit_times$V1, origin="1970-01-01")
colnames(last_visit_times)[3] <- "last_visit_time"
part_info <- merge(part_info, last_visit_times)


## Get number of months
if(get_datetime){
  if(buckets == 1){
    part_info$V1 <-  as.numeric(format(part_info$V1, "%Y"))
    part_info$V2 <-  as.numeric(format(part_info$V2, "%Y"))
    part_info$V3 <-  as.numeric(format(part_info$V3, "%Y"))
    part_info$V4 <-  as.numeric(format(part_info$V4, "%Y"))
    part_info$last_visit_time <-  as.numeric(format(part_info$last_visit_time, "%Y"))
  } else if(buckets == 4){
    part_info$V1 <-  as.numeric(format(part_info$V1, "%Y"))*12 + as.numeric(format(part_info$V1,"%m"))-1
    part_info$V2 <-  as.numeric(format(part_info$V2, "%Y"))*12 + as.numeric(format(part_info$V2,"%m"))-1
    part_info$V3 <-  as.numeric(format(part_info$V3, "%Y"))*12 + as.numeric(format(part_info$V3,"%m"))-1
    part_info$V4 <-  as.numeric(format(part_info$V4, "%Y"))*12 + as.numeric(format(part_info$V4,"%m"))-1
    part_info$last_visit_time <-  as.numeric(format(part_info$last_visit_time, "%Y"))*12 + as.numeric(format(part_info$last_visit_time,"%m"))-1
    
    part_info$V1 <- floor(part_info$V1*(buckets/12))
    part_info$V2 <- floor(part_info$V2*(buckets/12))
    part_info$V3 <- floor(part_info$V3*(buckets/12))
    part_info$V4 <- floor(part_info$V4*(buckets/12))
    part_info$last_visit_time <- floor(part_info$last_visit_time*(buckets/12))
  }
}

###################################################################
## 3. FIND TIME OF EPI EVENTS AND MAKE SEPARATE DATA FRAMES
###################################################################
message("Find time of epi events and make separate data frames")
## Gives the following:
##     - ili_last_year
##     - smoking_all
##     - vaccine_last_year
##     - ili_last_visit
##     - vaccine_last_visit
##     - diseases
source("scripts/aux/get_epi_event_times.R")


###################################################################
## 4. REMOVE NOW REDUNDANT VARIABLES AND CONSOLIDATE
###################################################################
## Remove variables we now now longer need
part_info <- part_info[!(part_info$variable %in% c("PART_FEV_ILI","PART_LAST_VAC", "PART_FEV_ILI_LAST_YR",
                                                   "PART_SYMPT_1YR_MULTIPLE","PART_FEV_ILI_LAST_VISIT", 
                                                   "PART_VAC_LAST_VISIT","PART_SMOKE_LAST_MO",
                                                   "PART_VAC_LAST_YR")),]
message("Find BMI")
## Extract height and weight for BMI
tmp_heights <- unique(part_info[part_info$variable == "PART_HEIGHT",c("Participant_ID","value","visit")])
colnames(tmp_heights)[2] <- "PART_HEIGHT"
tmp_weights <- unique(part_info[part_info$variable == "PART_WEIGHT",c("Participant_ID","value","visit")])
colnames(tmp_weights)[2] <- "PART_WEIGHT"
tmp_size <- merge(tmp_heights, tmp_weights)
tmp_size$BMI <- tmp_size$PART_WEIGHT/(tmp_size$PART_HEIGHT/100)^2

## Remove height and weight now
part_info <- part_info[!(part_info$variable %in% c("PART_HEIGHT","PART_WEIGHT")),]
part_info <- merge(part_info, tmp_size)

## Remove disease status
part_info <- part_info[!(part_info$variable %in% c("PART_DX_HEART_DISEASE","PART_DX_DIABETES","PART_DX_HYPERTENSION",
                                            "PART_DX_CANCER","PART_DX_CHRONIC_LUNG","PART_DX_ALLERGIES")),]
part_info_bmi <- unique(part_info[,!(colnames(part_info) %in% c("value","variable"))])
part_info_unique <- unique(part_info_bmi[,!(colnames(part_info_bmi) %in% c("PART_HEIGHT","PART_WEIGHT",
                                                                                 "BMI","last_visit_time","visit"))])
part_info_unique <- melt(part_info_unique, id.vars=c("Participant_ID", "PARTICIPANT_ID",
                                                     "HH_ID", "LOC_ID", "PART_GENDER", "PART_BIRTH_YEAR","PART_OCC_STATUS",
                                                     "PART_BIRTH_MONTH","TRUE_OCC_STATUS"))
colnames(part_info_unique)[c(10,11)] <- c("raw_visit","samples")
#part_info_unique$samples <- as.numeric(format(part_info_unique$samples_date, "%Y"))*12 + 
#  as.numeric(format(part_info_unique$samples_date,"%m"))-1

#part_info_unique$samples <- floor(part_info_unique$samples*(buckets/12))

message("Merge with titre data")
## Merge participant info with the data frame used in fitting
fluscape_dat <- merge(fluscape_dat1, unique(part_info_unique[,c("Participant_ID","PARTICIPANT_ID","HH_ID","LOC_ID",
                                                                "raw_visit","samples","PART_GENDER",
                                                                "TRUE_OCC_STATUS")]),
                      by=c("Participant_ID","samples"),all.x=TRUE)

###################################################################
## 6. GET HOUSEHOLD DATA
###################################################################
message("Loading household data")
household_dat <- load.and.merge.household.data.V1.V2.V3.V4(topdir=fluscape_wd)

## Pull out latitude and longitude for each visit and the date at which they changed
use_colnames <- colnames(household_dat)[c(grep("HH_Lat", colnames(household_dat)), 
                                          grep("HH_Lon",colnames(household_dat)),
                                          grep("HH_VisitA1DateTime",colnames(household_dat)))]
use_names_dates <- colnames(household_dat)[grep("HH_VisitA1DateTime",colnames(household_dat))]
hh_loc_dat <- household_dat[,c("HH_ID","LOC_ID",use_colnames)]
hh_loc_dat <- hh_loc_dat %>% dplyr::mutate(HH_VisitA1DateTime.V1=as.Date(HH_VisitA1DateTime.V1),
                                    HH_VisitA1DateTime.V2=as.Date(HH_VisitA1DateTime.V2),
                                    HH_VisitA1DateTime.V3=as.Date(HH_VisitA1DateTime.V3),
                                    HH_VisitA1DateTime.V4=as.Date(HH_VisitA1DateTime.V4)) %>%
    dplyr::rename(V1 = HH_VisitA1DateTime.V1, 
           V2 = HH_VisitA1DateTime.V2, 
           V3 = HH_VisitA1DateTime.V3, 
           V4 = HH_VisitA1DateTime.V4)

## Convert visit dates to buckets to be enumerated
if(get_datetime){
    if(buckets == 1){
        hh_loc_dat$V1 <-  as.numeric(format(hh_loc_dat$V1, "%Y"))
        hh_loc_dat$V2 <-  as.numeric(format(hh_loc_dat$V2, "%Y"))
        hh_loc_dat$V3 <-  as.numeric(format(hh_loc_dat$V3, "%Y"))
        hh_loc_dat$V4 <-  as.numeric(format(hh_loc_dat$V4, "%Y"))
    } else if(buckets == 4){
        hh_loc_dat$V1 <-  as.numeric(format(hh_loc_dat$V1, "%Y"))*12 + as.numeric(format(hh_loc_dat$V1,"%m"))-1
        hh_loc_dat$V2 <-  as.numeric(format(hh_loc_dat$V2, "%Y"))*12 + as.numeric(format(hh_loc_dat$V2,"%m"))-1
        hh_loc_dat$V3 <-  as.numeric(format(hh_loc_dat$V3, "%Y"))*12 + as.numeric(format(hh_loc_dat$V3,"%m"))-1
        hh_loc_dat$V4 <-  as.numeric(format(hh_loc_dat$V4, "%Y"))*12 + as.numeric(format(hh_loc_dat$V4,"%m"))-1

        hh_loc_dat$V1 <- floor(hh_loc_dat$V1*(buckets/12))
        hh_loc_dat$V2 <- floor(hh_loc_dat$V2*(buckets/12))
        hh_loc_dat$V3 <- floor(hh_loc_dat$V3*(buckets/12))
        hh_loc_dat$V4 <- floor(hh_loc_dat$V4*(buckets/12))
    }
}

unique_hh <- hh_loc_dat %>% dplyr::select(HH_ID,LOC_ID) %>% distinct()

tmp <- hh_loc_dat %>% left_join(expand_grid(unique_hh, date=seq(min(hh_loc_dat$V1,na.rm=TRUE),max(hh_loc_dat$V4,na.rm=TRUE),by=1)))
tmp <- tmp %>% 
    mutate(HH_Lat = ifelse(date < V1, HH_Lat.V1, 
                           ifelse(date < V2, HH_Lat.V2, 
                                  ifelse(date < V3, HH_Lat.V3, HH_Lat.V4)))) %>% 
    mutate(HH_Long = ifelse(date < V1, HH_Long.V1, 
                           ifelse(date < V2, HH_Long.V2, 
                                  ifelse(date < V3, HH_Long.V3, HH_Long.V4)))) %>% 
    group_by(HH_ID, LOC_ID) %>% fill(HH_Lat, .direction="down") %>% fill(HH_Long, .direction="down") %>%
    dplyr::select(HH_ID, LOC_ID, date, HH_Lat, HH_Long)

## IMPORTANT -- find the household lat/long in each bucket used for infection history inference
hh_loc_dat <- tmp

hh_colnames <- c("LOC_ID","HH_ID","HH_Income.V1","HH_Income.V2",
                 "HH_Income.V3","HH_Income.V4",
                "HH_Num_Member","HH_Num_Member.V3","HH_Num_Member.V4",
                "HH_Type.V1","HH_Type.V2", "HH_Type.V3", "HH_Type.V4",
                "HH_RunningWater.V1", "HH_RunningWater.V2", "HH_RunningWater.V3", 
                "HH_RunningWater.V4")
household_dat <- household_dat[,hh_colnames]

## Sort HH V1
visits <- c("V1","V2","V3","V4")
hh_names <- c("HH_Income","HH_Num_Member","HH_Type","HH_RunningWater")
hh_all <- NULL
for(visit in visits){
  use_hh_names <- c("HH_ID","LOC_ID",paste0(hh_names,".",visit))
  which_use <- which(use_hh_names[3:length(use_hh_names)] %in% hh_colnames)
  use_hh_names <- use_hh_names[c(1,2,which_use+2)]
  hh_tmp <- household_dat[,use_hh_names]
  colnames(hh_tmp) <- c("HH_ID","LOC_ID",hh_names[which_use])
  hh_tmp <- melt(hh_tmp, id.vars=c("HH_ID","LOC_ID"))
  hh_tmp$visit <- visit
  hh_all <- rbind(hh_all, hh_tmp)
}

## Income: 999 = unsure; 888 = decline
## Type: 1 = family home; 2 = elderly care home; 3 = worker's dormitory; 4 = multi-family home; 5 = Other; 999 = Unknown
## Running water: 1 = yes; 0 = no; 888 = decline; 999 = unsure
## Num_member: integers
hh_all[hh_all$value %in% c(999,888),"value"] <- NA
hh_all <- hh_all[complete.cases(hh_all),]
household_dat_all <- dcast(hh_all, HH_ID + LOC_ID + visit ~ variable, value.var = "value")


###################################################################
## 6. GET LOCATION DATA
###################################################################
message("Loading location data")
loc_dat <- load.and.merge.locs.V1.V2.V3(topdir=fluscape_wd,trim=TRUE)
fluscape_dat <- merge(fluscape_dat, loc_dat, by="LOC_ID")

###################################################################
## 7. CALCULATE TITRE DATA QUANTITIES
###################################################################
message("Finding ages and age groups...")
## Calculate age as of last sampling time
fluscape_dat$age <- (max(fluscape_dat$samples) - fluscape_dat$DOB)/buckets
fluscape_dat$age_group <- cut(fluscape_dat$age,breaks=c(0,10,20,30,40,50, 60,100),include.lowest=TRUE)
fluscape_dat$birth_cohort <- cut(fluscape_dat$DOB/buckets,breaks=c(1900,1968,1978,1988,2000,2010, 2100),include.lowest=TRUE)
birth_cohort_levels <- c("[1.9e+03,1.97e+03]"="[-1968]","(1.97e+03,1.98e+03]"="(1968,1978]",
                         "(1.98e+03,1.99e+03]"="(1978,1988]","(1.99e+03,2e+03]"="(1988,2000]",
                         "(2e+03,2.01e+03]"="(2000,2010]","(2.01e+03,2.1e+03]" ="(2010-]")
fluscape_dat$birth_cohort <- birth_cohort_levels[as.character(fluscape_dat$birth_cohort)]
fluscape_dat$birth_cohort <- factor(fluscape_dat$birth_cohort, levels=rev(birth_cohort_levels),ordered=TRUE)


tmp <- unique(fluscape_dat[,c("individual","samples")])
tmp <- ddply(tmp, .(individual), function(x) cbind(x, "visit"=1:nrow(x)))
fluscape_dat <- base::merge(fluscape_dat[,colnames(fluscape_dat) != "visit"], tmp)

## Order by age overall
message("Order factor level by age")
tmp <- unique(fluscape_dat[,c("individual","DOB")])
tmp_ordered <- tmp[sort(tmp$DOB,index.return=TRUE,decreasing = TRUE)$ix,]
tmp_ordered$order <- 1:nrow(tmp_ordered)
fluscape_dat <- merge(fluscape_dat, tmp_ordered)
strains$virus <- strains$virus/4 * buckets
fluscape_dat <- merge(fluscape_dat, strains)

message("Get titre change between visits and titre medians")
## Get change in titre based on medians across visits
first_visit_titres <- fluscape_dat[fluscape_dat$visit == 1, c("individual","virus","titre","run")]
colnames(first_visit_titres)[colnames(first_visit_titres) == "titre"] <- "v1_titre"
first_visit_titres <- ddply(first_visit_titres, .(individual,virus), function(x) max(x$v1_titre))
colnames(first_visit_titres)[3] <- "v1_titre"

second_visit_titres <- fluscape_dat[fluscape_dat$visit == 2, c("individual","virus","titre","run")]
colnames(second_visit_titres)[colnames(second_visit_titres) == "titre"] <- "v2_titre"
second_visit_titres <- ddply(second_visit_titres, .(individual,virus), function(x) max(x$v2_titre))
colnames(second_visit_titres)[3] <- "v2_titre"

tmp <- merge(first_visit_titres, second_visit_titres,all=TRUE)
tmp$change <- tmp$v2_titre - tmp$v1_titre
fluscape_dat <- merge(fluscape_dat, tmp,all.x=TRUE)

## Adding a slight offset for plotting
strains1 <- strains[complete.cases(strains),]
strains1[strains1$Year == 2009, "virus"] <- strains1[strains1$Year == 2009, "virus"] - 2
strains1[strains1$Year == 2010, "virus"] <- strains1[strains1$Year == 2010, "virus"] + 2

## Generate age group breaks
tmp <- unique(fluscape_dat[,c("individual","age_group")])
age_group_counts <- c(0,cumsum(ddply(tmp, ~age_group, nrow)$V1))
age_group_counts <- age_group_counts[1:(length(age_group_counts)-1)]

fluscape_dat$titre <- as.factor(fluscape_dat$titre)
colnames(fluscape_dat)[colnames(fluscape_dat) == "titre"] <- "log titre"

fluscape_dat[fluscape_dat$visit == 1,"visit"] <- "First visit"
fluscape_dat[fluscape_dat$visit == "2","visit"] <- "Second visit"

#################################
## Get seroconversions
seroconv_by_year <- ddply(fluscape_dat, ~Year, function(x){
  second_samp_titres <- x[x$visit == "Second visit",]
  second_samp_titres_change <- ddply(second_samp_titres, ~individual, function(y){
    mean(y$change)
  })
  all_tests <- nrow(second_samp_titres_change)
  all_pos <- nrow(second_samp_titres_change[second_samp_titres_change$V1 >= 2,])
  c(all_pos, all_tests)
})
cbind(seroconv_by_year[,1],signif(Hmisc::binconf(seroconv_by_year[,2],seroconv_by_year[,3])*100,3))


seroconv_by_year_and_age <- ddply(fluscape_dat, .(Year,age_group), function(x){
  second_samp_titres <- x[x$visit == "Second visit",]
  second_samp_titres_change <- ddply(second_samp_titres, ~individual, function(y){
    mean(y$change)
  })
  all_tests <- nrow(second_samp_titres_change)
  all_pos <- nrow(second_samp_titres_change[second_samp_titres_change$V1 >= 2,])
  c(all_pos, all_tests)
})

#################################
## Generate table of participant info
part_info <- unique(fluscape_dat[,c("individual","DOB","LOC_ID","Participant_ID","HH_ID",
                                    "PARTICIPANT_ID",  "raw_visit","age","age_group","visit","samples")])
part_info_v2 <- part_info[part_info$visit=="Second visit",]
part_info_v1 <- part_info[part_info$visit == "First visit",]
part_info <- rbind(part_info_v2, 
                   part_info_v1[part_info_v1$Participant_ID %in% 
                                  setdiff(part_info_v1$Participant_ID, part_info_v2$Participant_ID),])
tmp <- part_info_bmi[,c("Participant_ID","visit","PART_GENDER","BMI")]
colnames(tmp)[2] <- "raw_visit"
tmp <- merge(part_info, tmp)
tmp$bmi_group <- cut(tmp$BMI,breaks=c(0,18.5,24.9,29.9,39.9,100),include.lowest=TRUE)

colnames(vaccine_last_visit)[2] <- "raw_visit"
vaccine_last_visit1 <- merge(unique(part_info[,c("Participant_ID","raw_visit")]),vaccine_last_visit)

colnames(vaccine_last_year)[2] <- "raw_visit"
vaccine_last_year1 <- merge(unique(part_info[,c("Participant_ID","raw_visit")]),vaccine_last_year)


## ILI data
colnames(ili_last_year)[2] <- "raw_visit"
ili_last_year1 <- merge(unique(fluscape_dat[,c("Participant_ID","raw_visit")]),ili_last_year)


colnames(ili_last_visit)[2] <- "raw_visit"
ili_last_visit1 <- merge(unique(fluscape_dat[,c("Participant_ID","raw_visit")]),ili_last_visit)

## Smoking
smoking_all1 <- merge(unique(fluscape_dat[,c("Participant_ID","raw_visit","visit")]),smoking_all)

## Household data
colnames(hh_all)[5] <- "raw_visit"
hh_used <- unique(part_info[,c("HH_ID","LOC_ID","raw_visit")])
hh_all1 <- merge(hh_used, hh_all)
hh_all1$hh_id <- paste0(hh_all1$HH_ID,"_",hh_all1$LOC_ID)

hh_unique <- unique(hh_all1[,c("hh_id","raw_visit")])
hh_unique$visit_no <- as.numeric(hh_unique$raw_visit)
tmp <- ddply(hh_unique, ~hh_id, function(x) x[which.max(x$visit_no),])
hh_all1 <- merge(tmp, hh_all1)

table(hh_all1[hh_all1$variable == "HH_Income","value"])
table(hh_all1[hh_all1$variable == "HH_Type","value"])
table(hh_all1[hh_all1$variable == "HH_Num_Member","value"])
table(hh_all1[hh_all1$variable == "HH_RunningWater","value"])

loc_dat <- unique(fluscape_dat[,c("LOC_ID","DIST_FRM_GZ","URBAN","dens.1","dens.9","trav.km","trav.min")])
loc_dat$dist_grp <- cut(loc_dat$DIST_FRM_GZ,breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6),include.lowest = TRUE)
loc_dat$trav_dist_grp <- cut(loc_dat$trav.km,breaks=c(0,20,40,60,80,100),include.lowest = TRUE)
loc_dat$trav_time_grp <- cut(loc_dat$trav.min,breaks=c(0,20,40,60,80,100),include.lowest = TRUE)
loc_dat$log_dens1 <- log(loc_dat$dens.1)
loc_dat$dens1_grp <- cut(loc_dat$log_dens1, breaks=c(0,2,4,6,8,10,100))
loc_dat$log_dens9 <- log(loc_dat$dens.9)
loc_dat$dens9_grp <- cut(loc_dat$log_dens9, breaks=c(0,2,4,6,8,10,100))
## Plot distribution of sampling times
tmp <- unique(fluscape_dat[,c("individual","visit","samples")])
ylims <- range(fluscape_dat$samples)
ylab <- seq(ylims[1], ylims[2],by=1)
ylab <- ylab/4
to_append <- rep(c("Q1","Q2","Q3","Q4"),5)
to_append <- c("Q4",to_append,"Q1","Q2")
ylabels <- paste0(to_append, "-",floor(ylab))
sample_counts <- tmp %>% dplyr::group_by(visit, samples) %>% tally() %>% dplyr::rename(Visit=visit)
#sample_counts <- ddply(tmp, .(visit, samples), count)
#colnames(sample_counts)[1] <- "Visit"
p_samples <- ggplot(sample_counts) + 
  geom_bar(aes(x=samples,y=n,fill=Visit),stat="identity",col="black") +
  theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(labels=ylabels,breaks=ylab*4) + 
  theme(axis.text.x=element_text(angle=35,hjust=1,colour="black",size=8),
        axis.text.y=element_text(family="sans",colour="black",size=8),
        legend.position=c(0.2,0.9))+
  ylab("Count") +
  xlab("Sampling time") 
ggsave_jah(p_samples, figure_wd, "sample_times",6,4)

## sample_counts %>% mutate(samples=samples/4) %>% write.csv(file="~/Documents/GitHub/fluscape_infection_histories/data/figure_data/FigS1.csv",row.names=FALSE)

#fluscape_dat <- fluscape_dat[!is.na(fluscape_dat$change),]
fluscape_dat <- fluscape_dat %>% select(virus, individual, DOB, samples, LOC_ID, group, `log titre`,run,raw_visit, LOC_Lat,LOC_Long, dens.1, dens.9, age, age_group,birth_cohort,visit,order,Virus,Full_name, Strain, Year, v1_titre, v2_titre, change)

fluscape_dat$DOB <- floor(fluscape_dat$DOB/4)*4
fluscape_dat$DOB <- pmax(1930*4, fluscape_dat$DOB)
fluscape_dat$age <- floor(fluscape_dat$age)
fluscape_dat$Participant_ID <- fluscape_dat$individual

