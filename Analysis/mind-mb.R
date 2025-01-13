## Step 1: identifying metabolites associated with MIND score
mb_sig <- c(NULL)
for (expo in expo_list){
  for (mb in mb_list){
    eval(parse(text = paste0("alldata$",mb, ".z=scale(alldata$", mb, ")")))
    formula_model <- paste0("mb.model <- lm(",
                            mb,
                            ".z~",
                            expo
                            ,"+age_rec+sex+edu_college+factor(pa_cat)+factor(smoking)+eth_white+TDIc+BMIc+kcal_ave+hbp_baseline+dm_baseline+fo_cvd_base,data=alldata)")
    eval(parse(text = formula_model))
    a <- summary(mb.model)
    mb_sig <- rbind(mb_sig,c(expo, mb, as.numeric(a$coefficients[, "Estimate"][2]), 
                             as.numeric(a$coefficients[, "Std. Error"][2]), 
                             as.numeric(a$coefficients[, "Pr(>|t|)"][2])))
  }
}

mb_sig_summary <- data.frame(mb_sig)
names(mb_sig_summary) <- c("Exposure", "Variable", "Estimate", "se", "P")
write_csv(mb_sig_summary, "mb_sig.csv")
mb_sig_summary <- read_csv("mb_sig.csv")

mb_sig_summary <- mb_sig_summary %>% 
  group_by(Exposure) %>% 
  mutate(Pfdr = p.adjust(P, method ="BH"),
         Pbon = p.adjust(P, method ="bonferroni"))

sum(mb_sig_summary$Pbon<=0.05)
sum(mb_sig_summary$Pfdr<=0.05)

mb_sig_summary <- left_join(mb_sig_summary, nmr_ls[,c("var_name", "title", "Group", "Subgroup")], by = c("Variable"="var_name"))
write_csv(mb_sig_summary, "mb_sig.csv")

