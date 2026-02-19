
library(ggplot2)
ggplot(hb_stage_2a, aes(x = cell_type, y = new_ratio, fill = cell_type)) + 
   geom_boxplot() +
  labs(x = "Cell Type", y = "New Ratio", title = "Cell Type vs Ratio Distributions") + 
  scale_fill_manual(values = hb_pal) +   
  theme_classic() +                   
  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "none")


ggsave("cell_type_vs_ratioo.png", width = 8, height = 6, dpi = 300)

library(readxl)
hb_stage_2b <- read_excel("C:/Users/HOME/OneDrive/Documents/HACKBIO_DATAVIZ/Stage Two/Visualization/hb_stage_2.xlsx", sheet = "b")
library(ggplot2)
library(ggrepel)

# 1. Define the exact thresholds from your image
half_cut <- 2.5 
alpha_cut <- -3.5 
hb_stage_2b_log2_df <- as.data.frame(hb_stage_2b)
View(hb_stage_2b_log2_df)

hb_stage_2b_log2_df$log2_half_life <- log2(hb_stage_2b_log2_df$half_life) 
hb_stage_2b_log2_df$log2_alpha <- log2(hb_stage_2b_log2_df$alpha)

 ggplot(hb_stage_2b_log2_df, aes(x = log2_half_life, y = log2_alpha)) + 
 
  geom_point(aes(color = interaction(log2_half_life > half_cut, log2_alpha > alpha_cut)), 
             alpha = 0.5, size = 0.8) + 
  
  geom_vline(xintercept = half_cut, linetype = "dashed") + 
  geom_hline(yintercept = alpha_cut, linetype = "dashed") + 
  
  scale_color_manual(values = c(
    "FALSE.FALSE" = "grey50",   
    "TRUE.FALSE"  = "royalblue", 
    "FALSE.TRUE"  = "green3",    
    "TRUE.TRUE"   = "red2"      
  )) + 
  
  geom_text_repel(
    data = subset(hb_stage_2b_log2_df, cell %in% c("Camp", "Ccr2")), 
    aes(label = cell),
    size = 7,               
    fontface = "bold", 
    nudge_x = 0.5,           
    point.padding = 0.5
  ) + 
  theme_classic() + 
  theme(legend.position = "none") + 
  labs(x = "log2(Half Life)", y = "log2(Alpha Life)")
ggsave("scatter_plot_half-life_vs_alpha_life.png", width = 8, height = 6, dpi = 300)

dev.off()

library(readxl)
hb_stage_2c <- read_excel("C:/Users/HOME/OneDrive/Documents/HACKBIO_DATAVIZ/Stage Two/Visualization/hb_stage_2.xlsx", sheet = "c")
View(hb_stage_2c)
rownames(hb_stage_2c) <- hb_stage_2c$genes
hb_stage_2c$genes <- NULL
hb_2c_matrix <- as.matrix(hb_stage_2c)


library(dplyr)
library(stringr)

col_info <- data.frame(Sample = colnames(hb_2c_matrix)) %>% mutate(CellType = str_extract(Sample, "^[^0-9]+"), CellType = str_replace(CellType, "n$", ""), Time = str_extract(Sample, "[0-9]+"))


rownames(col_info) <- col_info$Sample
head(col_info)

library(pheatmap)
ann_col <- col_info[, c("CellType", "Time")]


current_times <- unique(ann_col$Time)
current_cells <- unique(ann_col$CellType)

ann_colors_list <- list(
  Time = setNames(hb_pal[1:length(current_times)], current_times),
  CellType = setNames(hb_pal[11:(10 + length(current_cells))], current_cells)
)

library(pheatmap)




pheatmap(hb_2c_matrix,
         annotation_col = ann_col,
         annotation_colors = ann_colors_list, # Uses your existing list
         color = colorRampPalette(c("white", "aliceblue", "dodgerblue4"))(100),
         scale = "row",            
         cluster_rows = TRUE,      
         cluster_cols = FALSE,    
         show_rownames = FALSE, 
         show_colnames = FALSE,
         filename = "heatmap_across_cell_types_and_time.png", width = 10, height = 12, units = "in", dpi = 300)

library(readxl)
hb_stage_2d_1 <- read_excel("C:/Users/HOME/OneDrive/Documents/HACKBIO_DATAVIZ/Stage Two/Visualization/hb_stage_2.xlsx", sheet = "d_1")
View(hb_stage_2d_1)
hb_2d_1_matrix <- as.matrix(hb_stage_2d_1[,2:8])
View(hb_2d_1_matrix)

library(pheatmap)
p_e_h <- pheatmap(hb_2d_1_matrix, border_color = NA, color = colorRampPalette(c("firebrick", "#f0f9ff", "#084594"))(100), cluster_rows = F, cluster_cols = F, angle_col = 90)
ggsave("pathway_enrichment_heatmap.png", plot = p_e_h, width = 8, height = 6, dpi = 300)

hb_stage_2e <- read_excel("C:/Users/HOME/OneDrive/Documents/HACKBIO_DATAVIZ/Stage Two/Visualization/hb_stage_2.xlsx", sheet = "e")
View(hb_stage_2e)
 library(ggplot2)
BB <- ggplot(hb_stage_2e, aes(x = half_life, y = alpha, colour = stage, size = count)) + geom_point(alpha = 0.7) + scale_color_manual(values = c("72h" = "green", "6h" = "navy")) +
labs(x = "Half Life", y = "Alpha Life", color = NULL, size = NULL) + theme_classic() + scale_x_continuous(limits = c(0, 80)) +
 theme(
 legend.position = "inside",
 legend.position.inside = c(0.9, 0.5),
legend.spacing.y = unit(2, "cm"), legend.box = "vertical", legend.background = element_rect(color = "black", linetype = "solid"))

ggsave("Bubble_plot_of_kinetic_regimes.png", plot = BB, width = 8, height = 6, dpi = 300)

hb_stage_2f <- read_excel("C:/Users/HOME/OneDrive/Documents/HACKBIO_DATAVIZ/Stage Two/Visualization/hb_stage_2.xlsx", sheet = "f")
View(hb_stage_2f)


library(dplyr)

subset_stage <- hb_stage_2f %>% filter(stage %in% c("s00h", "s72h"))

SB <- ggplot(subset_stage, aes(x = stage, y = proportion, fill = cell_type)) + 
  geom_col(width = 0.5, position = "stack") + # Explicitly tell it to stack
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + # Ensures bars start at 0
  scale_fill_manual(values = c("B" = "pink", "Plasma" = "navy")) + 
  theme_classic() + theme(legend.position = "inside", legend.justification.inside = c(0.75, 0.7), legend.box.just = "right")
ggsave("stacked_bar_proportion.png", plot = SB, width = 8, height = 6, dpi = 300)

hb_stage_2g <- read_excel("C:/Users/HOME/OneDrive/Documents/HACKBIO_DATAVIZ/Stage Two/Visualization/hb_stage_2.xlsx", sheet = "g")
View(hb_stage_2g)

hb_stage_2g_df <- as.data.frame(hb_stage_2g)


rownames(hb_stage_2g_df) <- hb_stage_2g_df[, 1]

hb_2g_df_numeric <- hb_stage_2g_df[, -1]
View(hb_2g_df_numeric)

adj_matrix <- as.matrix(hb_2g_df_numeric)
library(igraph)
net <- graph_from_adjacency_matrix(adj_matrix, 
                                   mode = "directed", 
                                   weighted = TRUE, 
                                   diag = FALSE)
plot(net)
E(net)$width <- E(net)$weight * 7

 par(mar=c(1,1,1,1))
 V(net)$size <- 25
 V(net)$color <- "pink"
V(net)$label.color <- "black"
E(net)$color <- "darkgrey" 
E(net)$width <- 2 

 png("Directed_cell_interactions.png", width = 1000, height = 1000, res = 150)

  plot(net, layout = layout_nicely)

   dev.off()
 
 