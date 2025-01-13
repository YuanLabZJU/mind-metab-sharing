# Figure 2C-E. Correlation between aMIND and MIND-MetS
library(patchwork)

scatter.plot1 <- ggplot(data=alldata_ukb_train, aes(x = MIND0, y = metSig)) + 
  geom_point(alpha=0.5) + 
  geom_smooth(method = 'lm') + 
  theme_minimal() + 
  geom_text(aes(x = 2, y = 1.5, label = text)) + 
  xlab("aMIND diet score") + 
  ylab("MIND Metabolomic signature score")

scatter.plot2 <- ggplot(data=alldata_ukb_val, aes(x = MIND, y = metSig)) + 
  geom_point(alpha=0.5) + 
  geom_smooth(method = 'lm') + 
  theme_minimal() + 
  geom_text(aes(x = 2, y = 1.5, label = text)) + 
  xlab("aMIND diet score") + 
  ylab("MIND Metabolomic signature score")

scatter.plot3 <- ggplot(data=alldata_whii, aes(x = MIND, y = metSig)) + 
  geom_point(alpha=0.5) + 
  geom_smooth(method = 'lm') + 
  theme_minimal() + 
  geom_text(aes(x = 2, y = 1.5, label = text)) + 
  xlab("aMIND diet score") + 
  ylab("MIND Metabolomic signature score")

plot_scatter <- scatter.plot1 + scatter.plot2 + scatter.plot3

ggsave("Fig2CDE.pdf", plot_scatter, width = 12, height = 5)