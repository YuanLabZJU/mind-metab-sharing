library(mediation)

# Dataset preparation and merging
alldata1 <- alldata1 |> 
  mutate(dem_tol = dem_tol - (as.numeric(date_ins1) - as.numeric(date_ins0)), # redefine survival time for the prospective dataset
         age_rec = age_rec + (as.numeric(date_ins1) - as.numeric(date_ins0))/365.25 # redefine baseline age for the prospective dataset
         )

alldata_cox <- alldata |> 
  dplyr::select(eid, date_ins0, date_ins2, dem_tol, dem_status, metSig, MIND0, age_rec, sex, edu_college, pa_cat, smoking, eth_white, TDIc, BMIc, kcal_ave, hbp_baseline, dm_baseline, fo_cvd_base, E4_num, drinking) |> 
  rename("MIND_ave" = "MIND0")
alldata1_cox <- alldata1 |> 
  dplyr::select(eid, date_ins0, date_ins2, dem_tol, dem_status, metSig, MIND_ave, age_rec, sex, edu_college, pa_cat, smoking, eth_white, TDIc, BMIc, kcal_ave, hbp_baseline, dm_baseline, fo_cvd_base, E4_num, drinking)

alldata_cox_both <- rbind(alldata_cox, alldata1_cox) # 45906
45906 - 29458
# Filter again: Age >=55 for cognition and dementia analysis
alldata_cox_both <- dplyr::filter(alldata_cox_both, age_rec>=55) # 29458
alldata_cox_both_r1 <- alldata_cox_both
alldata_cox_both <- dplyr::filter(alldata_cox_both, !(dem_status==1 & age_rec+dem_tol/365.25<65)) # 29449, 9 excluded
sum(alldata_cox_both$dem_status) # 482 dementia cases

sum(alldata_cox_both$dem_tol)/nrow(alldata_cox_both)/365.25 # 1.3e+08 P-Ys


# MIND diet score and dementia - UKB
alldata_cox_both <- alldata_cox_both |> filter(dem_tol>0)
cph_mets <- (coxph(Surv(dem_tol, dem_status)~scale(metSig)+age_rec+sex+edu_college+factor(pa_cat)+factor(smoking)+eth_white+TDIc+BMIc+kcal_ave+dm_baseline+hbp_baseline+fo_cvd_base, data=alldata_cox_both))
cph_mind <- (coxph(Surv(dem_tol, dem_status)~scale(MIND_ave)+age_rec+sex+edu_college+factor(pa_cat)+factor(smoking)+eth_white+TDIc+BMIc+kcal_ave+dm_baseline+hbp_baseline+fo_cvd_base, data=alldata_cox_both))
# ShowRegTable(cph_mets)
med.fit <- lm(scale(metSig) ~ scale(MIND_ave), data=alldata_cox_both) # 0.34, <2e-16
med.fit <- lm((metSig) ~ (MIND_ave), data=alldata_cox_both)
out.fit <- survreg(Surv(dem_tol, dem_status)~metSig + MIND_ave+age_rec+sex+
                     edu_college+factor(pa_cat)+factor(smoking)+eth_white+
                     TDIc+BMIc+kcal_ave+dm_baseline+hbp_baseline+fo_cvd_base,
                   data=alldata_cox_both)
set.seed(20240112)

med.out<-mediate(med.fit,out.fit,treat = "MIND_ave",mediator = "metSig",
                 sims = 500)
med_sum<-summary(med.out)
med_sum # Prop. Mediated (control)    0.4585  0.036 * 

# Subgroup ananlyses
pint_cox <- function(expo_var, strata_vars){
  pint <- NULL
  alldata_cox_both$expo <- alldata_cox_both[,expo_var]
  for (var in strata_vars){
    alldata_cox_both$var <- alldata_cox_both[,var]
    main_model <- coxph(Surv(dem_tol, dem_status)~scale(expo)+age_rec+sex+edu_college+factor(pa_cat)+factor(smoking)+eth_white+TDIc+BMIc+kcal_ave+dm_baseline+hbp_baseline+fo_cvd_base, 
                        data=alldata_cox_both)
    int_model <- coxph(Surv(dem_tol, dem_status)~scale(expo)*var+age_rec+sex+edu_college+factor(pa_cat)+factor(smoking)+eth_white+TDIc+BMIc+kcal_ave+dm_baseline+hbp_baseline+fo_cvd_base, 
                        data=alldata_cox_both)
    aov <- anova(main_model, int_model)
    pint <- c(pint, aov$`P(>|Chi|)`[2])
  }
  pint_table <- data.frame(expo = expo_var, var_name = strata_vars, pint = pint)
  return(pint_table)
}

pint_lmer <- function(expo_var, strata_vars){
  pint <- NULL
  mi_long$expo <- mi_long[,expo_var]
  for (var in strata_vars){
    mi_long$var <- mi_long[,var]
    main_model <- lmer(zscore~scale(expo)+age_rec+sex+edu_college+
                         factor(pa_cat)+factor(smoking)+eth_white+TDIc+
                         BMIc+kcal_ave+dm_baseline+hbp_baseline+fo_cvd_base +
                         zscore_base + (1|eid), 
                        data=filter(mi_long, domain == "global"))
    int_model <- lmer(zscore~scale(expo)*var+age_rec+sex+edu_college+
                        factor(pa_cat)+factor(smoking)+eth_white+TDIc+
                        BMIc+kcal_ave+dm_baseline+hbp_baseline+fo_cvd_base +
                        zscore_base + (1|eid), 
                      data=filter(mi_long, domain == "global"))
    aov <- anova(main_model, int_model)
    pint <- c(pint, aov$`Pr(>Chisq)`[2])
  }
  pint_table <- data.frame(expo = expo_var, var_name = strata_vars, pint = pint)
  return(pint_table)
}

sg_cox <- function(strata_vars){
  sgRegTable <- NULL
  for (var in strata_vars){
    print(var)
    alldata_cox_both$var <- factor(alldata_cox_both[,var])
    for (varlv in levels(alldata_cox_both$var)){
      cph_model_1 <- coxph(Surv(dem_tol, dem_status)~scale(metSig)+age_rec+sex+edu_college+factor(pa_cat)+factor(smoking)+eth_white+TDIc+BMIc+kcal_ave+dm_baseline+hbp_baseline+fo_cvd_base, 
                         data=alldata_cox_both |> filter(var==varlv))
      regtable <- summary(cph_model_1)
      regcoef <- c(var = var, exposure = "metSig", varlv = varlv, regtable$coefficients[1,])
      sgRegTable <- rbind(sgRegTable, regcoef)
      
      cph_model_2 <- coxph(Surv(dem_tol, dem_status)~scale(MIND_ave)+age_rec+sex+edu_college+factor(pa_cat)+factor(smoking)+eth_white+TDIc+BMIc+kcal_ave+dm_baseline+hbp_baseline+fo_cvd_base, 
                         data=alldata_cox_both |> filter(var==varlv))
      regtable <- summary(cph_model_2)
      regcoef <- c(var = var, exposure = "MIND_ave", varlv = varlv, regtable$coefficients[1,])
      sgRegTable <- rbind(sgRegTable, regcoef)
    }
  }
  return(sgRegTable)
}

sg_lmm <- function(strata_vars){
  sgRegTable <- NULL
  for (var in strata_vars){
    print(var)
    mi_long$var <- factor(mi_long[,var])
    for (varlv in levels(mi_long$var)){
      lmm_model_1 <- lmer(zscore~scale(metSig)+age_rec+sex+edu_college+
                          factor(pa_cat)+factor(smoking)+eth_white+TDIc+
                          BMIc+kcal_ave+dm_baseline+hbp_baseline+fo_cvd_base +
                          zscore_base + (1|eid), 
                        data=filter(mi_long, domain == "global" & var==varlv))
      
      regtable <- summary(lmm_model_1)
      regcoef <- c(var = var, expo = "metSig", varlv = varlv, regtable$coefficients[2,])
      sgRegTable <- rbind(sgRegTable, regcoef)
      
      lmm_model_2 <- lmer(zscore~scale(MIND_ave)+age_rec+sex+edu_college+
                            factor(pa_cat)+factor(smoking)+eth_white+TDIc+
                            BMIc+kcal_ave+dm_baseline+hbp_baseline+fo_cvd_base +
                            zscore_base + (1|eid), 
                          data=filter(mi_long, domain == "global" & var==varlv))
      regtable <- summary(lmm_model_2)
      regcoef <- c(var = var, expo = "MIND_ave", varlv = varlv, regtable$coefficients[2,])
      sgRegTable <- rbind(sgRegTable, regcoef)
    }
  }
  return(sgRegTable)
}


alldata_cox_both$age_b <- ifelse(alldata_cox_both$age_rec<=60, "Y <=60", "O >60")
alldata_cox_both$smo_b <- ifelse(alldata_cox_both$smoking==0, "Never", "Ever")
alldata_cox_both$ovw_b <- ifelse(alldata_cox_both$BMIc %in% c("(25,30]", "(30,999]"), "Ovw", "NonOvw")
alldata_cox_both$cvddm <- ifelse(alldata_cox_both$fo_cvd_base==1 | alldata_cox_both$dm_baseline==1, 1, 0)
mi_long$age_b <- ifelse(mi_long$age_rec<60, "Y <60", "O >60")
mi_long$smo_b <- ifelse(mi_long$smoking==0, "Never", "Ever")
mi_long$ovw_b <- ifelse(mi_long$BMIc %in% c("(25,30]", "(30,999]"), "Ovw", "NonOvw")
mi_long$cvddm <- ifelse(mi_long$fo_cvd_base==1 | mi_long$dm_baseline==1, 1, 0)
mi_long <- data.frame(mi_long)

sg_vars <- c("age_b", "sex", "edu_college", 
             "eth_white", "ovw_b", "hbp_baseline", "cvddm")

sg_lmm(sg_vars)
sg_cox_table <- sg_cox(sg_vars)
sg_lmm_table <- sg_lmm(sg_vars)
pint_cox_table <- rbind(pint_cox("MIND_ave", sg_vars), pint_cox("metSig", sg_vars))
pint_lmm_table <- rbind(pint_lmer("MIND_ave", sg_vars), pint_lmer("metSig", sg_vars))
sg_cox_table_output <- sg_cox_table |> 
  data.frame() |> 
  left_join(pint_cox_table, by = c("exposure"= "expo", "var" = "var_name")) |> 
  mutate(
    coef=as.numeric(coef),
    se.coef.=as.numeric(se.coef.),
    HR = exp(coef),
    HRL = exp(coef-1.96*se.coef.),
    HRU = exp(coef+1.96*se.coef.),
    HRCI = paste0(round(HR,2), " (", round(HRL,2), ", ", round(HRU,2), ")"))

sg_lmm_table_output <- sg_lmm_table |> 
  data.frame() |> 
  left_join(pint_lmm_table, by = c("expo"= "expo", "var" = "var_name")) |> 
  mutate(
    coef=as.numeric(Estimate),
    se.coef.=as.numeric(Std..Error),
    cil = coef-1.96*se.coef.,
    ciu = coef+1.96*se.coef.,
    coefCI = paste0(round(coef,2), " (", round(cil,2), ", ", round(ciu,2), ")"))

sg_cox_table_output |> write.csv(file = "table-cph-ukb-subgroup.csv")
sg_lmm_table_output |> write.csv(file = "table-lmm-ukb-subgroup.csv")

cph_mets_a <- coxph(Surv(dem_tol, dem_status)~scale(metSig)+factor(drinking)+age_rec+sex+edu_college+factor(pa_cat)+factor(smoking)+eth_white+TDIc+BMIc+kcal_ave+dm_baseline+hbp_baseline+fo_cvd_base, 
                    data=alldata_cox_both)
cph_mind_a <- coxph(Surv(dem_tol, dem_status)~scale(MIND_ave)+factor(drinking)+age_rec+sex+edu_college+factor(pa_cat)+factor(smoking)+eth_white+TDIc+BMIc+kcal_ave+dm_baseline+hbp_baseline+fo_cvd_base, 
                    data=alldata_cox_both)

table_cph <- rbind(cph_mets=ShowRegTable(cph_mets, printToggle = FALSE)[1,],
                   cph_mind=ShowRegTable(cph_mind, printToggle = FALSE)[1,],
                   cph_mets_a = ShowRegTable(cph_mets_a, printToggle = FALSE)[1,],
                   cph_mind_a = ShowRegTable(cph_mind_a, printToggle = FALSE)[1,]
                   )
table_cph |> write.csv(file = "table-cph-ukb.csv")

rm(apoe, cmd_processed)
save(list = ls(), file = "mind-mb-proj-data.RData")

