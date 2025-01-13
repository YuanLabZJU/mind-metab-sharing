
## Step 2: ElasticNet to construct metabolomic signature of MIND diet
library(glmnet)
require(doMC)
registerDoMC(cores = 8)
set.seed(2022)
alldata <- data.frame(alldata)
X_dat <- scale(as.matrix(alldata[,mb_list]))
y_dat <- scale(alldata$MIND0)
cv.fit <- cv.glmnet(X_dat, y_dat, alpha=0.5, parallel = TRUE,
                    type.measure = "mse", nfolds = 10, standardize=TRUE)
plot(cv.fit)
loads <- coef(cv.fit, s = "lambda.1se")
loads <- loads %>% 
  as.matrix() %>% 
  as.data.frame()
loads$var <- rownames(loads)
loads <- left_join(loads, nmr_ls, by = c("var" = "var_name")) |> filter(s1!=0)
loads <- loads %>% 
  mutate(var = rownames(loads),
         var = as.numeric(str_sub(var, 3, 7))) %>% 
  # left_join(nmr_ls, by = c("var" = "field_id")) %>% 
  filter(s1!=0)
write_csv(loads, "EN_Loads.csv")
metSig <- predict(cv.fit, X_dat, s="lambda.1se")
alldata$metSig <- metSig
text <- paste("Pearson's r =", round(cor(alldata$metSig, alldata$MIND0),2))

scatter.plot <- ggplot(data=alldata, aes(x = MIND0, y = metSig)) + 
  geom_point(alpha=0.5) + 
  geom_smooth(method = 'lm') + 
  theme_minimal() + 
  geom_text(aes(x = 2, y = 1.5, label = text)) + 
  xlab("aMIND diet score") + 
  ylab("MIND Metabolomic signature score")

fit_plot <- cbind(loglambda = log(cv.fit$lambda), 
                  cvm = cv.fit$cvm, 
                  cvup = cv.fit$cvup, 
                  cvlo = cv.fit$cvlo,
                  n = cv.fit$nzero) %>% 
  data.frame()

length(mb_list)
scale_n <- max(cv.fit$cvup) - min(cv.fit$cvlo)
from_n <- min(cv.fit$cvlo)
fit.plot <- ggplot(data = fit_plot, aes(x = loglambda, y = cvm)) + 
  geom_point(alpha = 0.5) + 
  geom_line(data = fit_plot, aes(x = loglambda, y = n/168*scale_n + from_n), linewidth=0.3, colour="#0077FF") +
  geom_errorbar(aes(ymin = cvlo, ymax = cvup), alpha = 0.5) + 
  scale_y_continuous(sec.axis = sec_axis(~(.-from_n) *168/scale_n, name = "Number of metabolites")) + 
  geom_vline(aes(xintercept = log(cv.fit$lambda.1se)), colour="#990000", linetype="dashed") + 
  ylab("Cross-Validated Mean Square Error") + 
  theme_minimal()

pdf("UKB_MetSig-development.pdf", width = 10, height = 5)
fit.plot + scatter.plot
dev.off()

