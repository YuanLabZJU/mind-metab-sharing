# MIND diet score and NMR metabolomics in instance 1 (2012-2013)
## Step 1: 'Single-metabolite' analysis - replicate in NHS/HPFS?
## Step 2: Metabolomic signature of MIND diet in 50k sample
## Step 3: Metabolomic signature in 110k sample - association with dementia and brain MR
rm(list=ls())
setwd("~/hchen/mind-mb")
library(tidyverse)
library(corrplot)
library(tableone)
# library(ggnewscale)
library(patchwork)
library(readxl)
library(survival)
library(survminer)

nmr_ls <- read_excel("Nightingale_biomarker_groups-harmo.xlsx") |> 
  mutate(var_name = paste0("p", field_id, "_i0"))

## Load MIND diet score
load("~/data_gen/ukb/diet/mind_by_gram.rda") # 210948 participant w/ valid diet record
## Load covariates
load("~/data_gen/ukb/covars/covars.rda")

## Load NMR metabolomics dataset
load("~/data_gen/ukb/nmr-metab/nmr-metab.rda")
## Load APOE genotype
load("~/data_gen/ukb/genotype/apoe.rda") #APOE4 Status
## Load outcome data
load("~/data_gen/ukb/ho/ho_algo-231107.rda") # Dementia defined by algorithm
load("~/data_gen/ukb/ho/ho_death-240111.rda") # Death in registry
load("~/data_gen/ukb/ho/ho-cmd-fo.rda") # CMD first occurrences

## Merging datasets - for instance 0
alldata <- left_join(mind_diet_s, nmr_i0, by = c("eid"="eid"))
alldata <- left_join(alldata, data_cov, by = c("eid"="eid"))
alldata <- left_join(alldata, apoe, by = c("eid"="eid"))
alldata <- left_join(alldata, data_ho_algo, by = c("eid"="eid")) # 210948 participants
alldata <- left_join(alldata, cmd_processed, by = c("eid"="eid")) # 210948 participants
alldata <- left_join(alldata, ho_death, by = c("eid"="eid")) # 210948 participants
alldata <- data.frame(alldata)

eof_dem <- max(alldata$algo_dem_date, na.rm=TRUE)

alldata$death_date
alldata <- alldata |> 
  mutate(dem_status = ifelse(is.na(algo_dem_date), 0, 1),
         dem_eof = ifelse(dem_status==1, as.numeric(algo_dem_date),
                          ifelse(death==1, death_date, eof_dem)), # defining end of follow-up - death/dementia/or eof
         dem_tol = dem_eof - as.numeric(date_ins0), # defining time of follow-up
         self_cvd = ifelse(self_angina + self_hrtatt + self_strk>=1, 1, 0), # 10.1038/s41598-022-22469-6
         fo_cvd_base = ifelse(is.na(date_fo_cvd), 0, 
                              ifelse(date_fo_cvd-as.numeric(date_ins0)<=0, 1, 0)), # 10.1001/jamacardio.2022.3191
         Dyslipid  = ifelse(is.na(date_fo_dyslipid), 0, 
                            ifelse(date_fo_dyslipid-as.numeric(date_ins0)<=0, 1, 0)),
  )

alldata <- filter(alldata, MIND0>=0) # 70034 participants with MIND diet score at baseline
alldata <- filter(alldata, p23400_i0>=0) # 36883 also had metabolome assay
alldata <- filter(alldata, !(dem_status==1 & dem_tol<=5*365.25)) # 36841 w/o dementia in 5 years

alldata$TDIc <- cut_number(alldata$tdi, 3)
alldata$BMIc <- cut(alldata$bmi, breaks = c(0, 18.5, 25, 30, 999))
alldata$E4_num <- impute(alldata$apoe_e4_count, "Missing")
alldata$MIND_PER3 <- alldata$MIND0/3

expo_list <- c(names(mind_diet_s)[c(166:180)],"MIND_PER3")
mb_list <- setdiff(names(nmr_i0), "eid")

inv_norm <- function(x){
  # Zero-imputation
  x[is.na(x)] <- 0
  # Inverse normal transformation
  inv_norm_x <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(inv_norm_x)
}

# Calculate missing rate
missing_rates <- alldata[,mb_list] |> is.na() |> colSums() |> data.frame()
names(missing_rates) <- "missing_n"
missing_rates$var_name <- rownames(missing_rates)
missing_rates$missing_perc_dis <- missing_rates$missing_n/nrow(alldata)
missing_rates <- left_join(nmr_ls[,c("var_name", "title", "Group_S", "Subgroup")], missing_rates, by = c("var_name"="var_name"))
missing_rates |> write_csv("missing_rates-dis.csv")

# Inv_Norm transformation
alldata[,mb_list] <- sapply(alldata[,mb_list], inv_norm)

alldata <- alldata %>%
  mutate(
    but_ave = -(but_ave), # reverse
    rMeat_ave = -(rMeat_ave), # reverse
    fff_ave = -(fff_ave), # reverse
    pas_ave = -(pas_ave), # reverse
    whChe_ave = -(whChe_ave), # reverse
  )

## Step 0
alldata$sex <- 1-alldata$sex
vars <- c('MIND0', 
          'age_rec', 'sex', 'kcal_ave', 'eth_white', 'edu_college', 
          'TDIc', 'smoking', 'pa_cat', 
          'BMIc', 'drinking',
          'hbp_baseline', 'dm_baseline', 'fo_cvd_base')

catvars <- c('sex', 'edu_college', 'eth_white', 'smoking', 'pa_cat', 'drinking', 'TDIc', 'BMIc',
             'hbp_baseline', 'fo_cvd_base', 'dm_baseline')

alldata[,vars] <- alldata[,vars] |> mice(seed = 123) |> complete() # MICE for covariates
table1_overall <- CreateTableOne(vars=vars, data=alldata, factorVars=catvars)
table1_overall <- print(table1_overall,showAllLevels = FALSE,
                        contDigits=1, catDigits=1, nonnormal="MIND0")

write.csv((table1_overall), "Table1_UKB.csv")
