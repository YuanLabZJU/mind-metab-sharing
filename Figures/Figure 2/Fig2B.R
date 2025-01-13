# Figure 2B. Weights of metabolites in the elastic net model

ms_wgts <- read_excel("05_Supplementary Tables.xlsx", sheet = "ST06") %>%
  mutate(
    Metabolite = str_remove(Metabolite, "Concentration of "),
    Metabolite = str_replace(Metabolite, "Extremely Large", "XL"),
    Metabolite = str_replace(Metabolite, "Chylomicrons", "CM"),
    Metabolite = str_replace(Metabolite, "Free Cholesterol", "FC"),
    Metabolite = str_replace(Metabolite, "Cholesteryl Ester", "CE"),
    Weight = `Weight (for standardized variables)`) %>%
  select(Metabolite, Weight) %>%
  filter(Metabolite != "(Intercept)") %>% 
  arrange(Weight)

mat_wgts <- ms_wgts %>% 
  select(Weight) %>% 
  as.matrix()
rownames(mat_wgts) <- ms_wgts$Metabolite

library(circlize)
library(ComplexHeatmap)
mycol <- colorRamp2(c(-0.25, 0, 0.25), c("#2A6EBB", "white", "#CD202C"))

lg <- Legend(title = "Weight",
             col_fun = mycol,
             direction = c("vertical"),
             grid_height = unit(1, "cm"),
             grid_width = unit(0.5, "cm"))

pdf("fig2_ms.pdf", width = 6, height = 6)
circos.clear()
circos.heatmap(mat_wgts, col = mycol,
               rownames.side = "outside",
               rownames.col = "#333333",
               cell.border = "#D4D4D5",
               cell_width = rep(0.2, nrow(mat_wgts)),
               cluster = FALSE)
grid.draw(lg)
circos.clear()
dev.off()
