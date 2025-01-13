# Figure 1C. Food - Metabolomic associations

fg_mb <- read_excel("05_Supplementary Tables.xlsx", sheet = "ST05") %>% 
  data.frame() %>% 
  filter(P.Bonferroni<0.05) %>%
  mutate(Group = ifelse(Group %in% groups, Group, "Lipids and lipoproteins") %>% 
           factor(levels = groups)) %>% 
  select(Metabolite, Exposure, Group, Dataset)

fg_mb_intersect <- fg_mb %>% 
  group_by(Metabolite, Exposure, Group) %>%
  summarise(n = n()) %>%
  filter(n==3) %>% 
  ungroup() %>% 
  mutate(Dataset = "Consistent significance across datasets") %>% 
  select(-n)

fg_mb_all <- bind_rows(fg_mb, fg_mb_intersect)

plot_fg <- ggplot(fg_mb_all, aes(y = Exposure, fill = Group)) +
  geom_bar() +
  facet_wrap(~Dataset) +
  # scale_fill_aaas() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) + 
  scale_fill_jco() + 
  theme(legend.position = "bottom")
plot_fg
fig_1 <- (plot_spacer() / plot_fg) | plot_mind_mb
# A4 format
ggsave("fig1_fg_mb.pdf", plot_fg, width = 7, height = 7)
