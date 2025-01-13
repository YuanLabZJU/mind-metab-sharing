# Figure 1B - Associations between aMIND and metabolites

library(readxl)
library(tidyverse)
library(ggsci)
library(patchwork)

groups <- c("Amino acids", 
            "Cholesterol", "Cholesterol esters", "Fatty acids", "Fluid balance",
            "Free cholesterol", "Glycolysis-related metabolites", "Inflammation",
            "Ketone bodies", "Lipids and lipoproteins")

## Figure 1. Metabolomic associations
# Load the data
mind_mb <- read_excel("05_Supplementary Tables.xlsx", sheet = "ST04",
                      skip=1) %>% 
  data.frame() %>%
  mutate(Group = ifelse(Group %in% groups, Group, "Lipids and lipoproteins") %>% 
           factor(levels = groups),
         Metabolite = str_remove(Metabolite, "Concentration of "),
         Metabolite = str_replace(Metabolite, "Extremely Large", "XL"),
         Metabolite = str_replace(Metabolite, "Chylomicrons", "CM"),
  ) %>% 
  filter((Estimate...5 > 0 & Estimate...10 > 0 & Estimate...16 > 0) | 
           (Estimate...5 < 0 & Estimate...10 < 0 & Estimate...16 < 0)) %>%
  filter(P.FDR...8<0.05 & P.FDR...13<0.05 & P.FDR...19<0.05) %>% 
  arrange(Group, Estimate...5)
# mind_mb <- read_excel("05_Supplementary Tables.xlsx", sheet = "ST04",
#                       skip=1) %>% 
#   filter((Estimate...5 > 0 & Estimate...10 > 0) | 
#            (Estimate...5 < 0 & Estimate...10 < 0)) %>%
#   filter(`P-FDR...8`<0.05 & `P-FDR...13`<0.05)
mind_mb_sig <- mind_mb %>% 
  select(Metabolite,
         Group,
         Whitehall.II.variable.name, 
         Estimate...5,
         Estimate...10,
         Estimate...16,) %>%
  pivot_longer(cols = c(Estimate...5,
                        Estimate...10,
                        Estimate...16,)) %>% 
  mutate(Cohort = ifelse(str_detect(name, "5"), "UKB Discovery", 
                         ifelse(str_detect(name, "10"), "UKB Validation", 
                                "WHII")),
         CohortID = ifelse(str_detect(name, "5"), 1, 
                           ifelse(str_detect(name, "10"), 2, 
                                  3)),
         Metabolite = factor(Metabolite, levels = mind_mb$Metabolite),
         Coef = sprintf("%.2f", value),
         panel = ifelse(row_number() < nrow(.)/2, "A", "B"))

plot_heat <- mind_mb_sig %>% 
  ggplot(aes(x = Cohort, y = Metabolite, fill = value)) +
  geom_tile() + 
  # geom_text(aes(label = Coef)) +
  coord_fixed(ratio = 0.75) +
  # ylim(0, 2) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        # axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) + 
  scale_fill_gradient2(low = "#2A6EBB", high = "#CD202C") + 
  scale_y_discrete(position = "right")

table(mind_mb_sig$Group)

plot_cat <- mind_mb_sig %>% 
  ggplot(aes(x = 1, y = Metabolite, fill = Group)) +
  geom_tile() + 
  coord_fixed(ratio = 1) +
  scale_fill_jco() + 
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

plot_mind_mb <- plot_cat + plot_heat + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom",
        legend.direction = "vertical")
ggsave("fig1_mind_mb.pdf", plot_mind_mb, width = 4, height = 12)
