# Figure 2A

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

ggsave("Figure2A.pdf", fit.plot, width = 4, height = 4)