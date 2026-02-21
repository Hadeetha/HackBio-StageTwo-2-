# Step 1: Install and load packages
install.packages(c("readxl", "ggplot2", "dplyr", "pheatmap", "ggrepel", "igraph", "tidyr"))
library(readxl)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(ggrepel)
library(igraph)
library(tidyr)

# Step 2: Define colors
hb_pal <- c("#4e79a7","#8cd17d","#e15759","#fabfd2","#a0cbe8",
            "#59a14f","#b07aa1","#ff9d9a","#f28e2b","#f1ce63",
            "#79706e","#d4a6c8","#e9e9e9","#ffbe7d","#bab0ac",
            "#9d7660","#d37295","#86bcb6","#362a39","#cd9942")

# Step 3: Set file path
file_path <- "C:/Users/Dell/Documents/HACKBIO_DATAVIZ/Stage Two/hb_stage_2.xlsx"

# --- Plot A: Boxplot ---
data_a <- read_excel(file_path, sheet = "a")
ggplot(data_a, aes(x = cell_type, y = new_ratio, fill = cell_type)) +
  geom_boxplot() +
  labs(title = "Cell Type vs Ratio", x = "Cell Type", y = "New Ratio") +
  scale_fill_manual(values = hb_pal) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none")
ggsave("plot_a_boxplot.png", width = 8, height = 6)

# --- Plot B: Scatter Plot ---
data_b <- read_excel(file_path, sheet = "b")
data_b$log2_half_life <- log2(data_b$half_life)
data_b$log2_alpha <- log2(data_b$alpha)

ggplot(data_b, aes(x = log2_half_life, y = log2_alpha)) +
  geom_point(aes(color = interaction(log2_half_life > 2.5, log2_alpha > -3.5)),
             alpha = 0.5, size = 1) +
  geom_vline(xintercept = 2.5, linetype = "dashed") +
  geom_hline(yintercept = -3.5, linetype = "dashed") +
  scale_color_manual(values = c("FALSE.FALSE" = "grey", "TRUE.FALSE" = "blue",
                                "FALSE.TRUE" = "green", "TRUE.TRUE" = "red")) +
  geom_text_repel(data = subset(data_b, cell %in% c("Camp", "Ccr2")),
                  aes(label = cell), size = 5, fontface = "bold") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "log2(Half Life)", y = "log2(Alpha)")
ggsave("plot_b_scatter.png", width = 8, height = 6)

# --- Plot C: Heatmap with hb_pal colors ---
data_c <- read_excel(file_path, sheet = "c")

# Use gene names as row names
rownames(data_c) <- data_c$genes
data_c$genes <- NULL
matrix_c <- as.matrix(data_c)

# Create annotation info for columns
col_info <- data.frame(Sample = colnames(matrix_c)) %>%
  mutate(CellType = stringr::str_extract(Sample, "^[^0-9]+"),
         Time = stringr::str_extract(Sample, "[0-9]+"))
rownames(col_info) <- col_info$Sample

# Get unique categories
cell_types <- unique(col_info$CellType)
times <- unique(col_info$Time)

# Assign colors from hb_pal automatically
ann_colors <- list(
  CellType = setNames(hb_pal[1:length(cell_types)], cell_types),
  Time = setNames(hb_pal[(length(cell_types)+1):(length(cell_types)+length(times))], times)
)

# Draw heatmap
pheatmap(matrix_c,
         annotation_col = col_info[, c("CellType", "Time")],
         annotation_colors = ann_colors,
         color = colorRampPalette(c("white", "lightblue", "darkblue"))(100),
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         filename = "plot_c_heatmap.png", width = 10, height = 12)

# --- Plot D: Pathway Heatmap ---
data_d <- read_excel(file_path, sheet = "d_1")

# Use pathway names as row names
pathway_names <- data_d$pathway
length(pathway_names)
matrix_d <- as.matrix(data_d[, -1])
rownames(matrix_d) <- pathway_names



# Draw heatmap
pheatmap(matrix_d,
         border_color = NA,
         color = colorRampPalette(c("red", "white", "blue"))(100), # red-white-blue gradient
         scale = "row",              # normalize each pathway row
         cluster_rows = FALSE,       # keep pathways in given order
         cluster_cols = FALSE,       # keep time points in given order
         angle_col = 90,             # rotate column labels
         show_rownames = TRUE,       # show pathway names
         show_colnames = TRUE,       # show time labels
         fontsize_row = 8,
         filename = "plot_d_heatmap.png",
         width = 10, height = 8)

# Create the plot object
plot_e <- ggplot(data_e, aes(x = half_life, y = alpha, 
                             color = stage, size = count)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("72h" = "green", "6h" = "blue")) +
  scale_size_continuous(range = c(2, 10)) +
  labs(x = "Half Life", y = "Alpha Life", color = "Stage", size = "Count") +
  theme_classic() +
  scale_x_continuous(limits = c(0, 80)) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.9, 0.5),
    legend.box = "vertical",
    legend.spacing.y = unit(2, "cm"),
    legend.background = element_rect(color = "black", linetype = "solid")
  )

# Show in plot pane
print(plot_e)

# Save to file
ggsave("plot_e_scatter.png", plot = plot_e, width = 8, height = 6, dpi = 300)

# --- Plot F: Stacked Bar Chart ---
data_f <- read_excel(file_path, sheet = "f")

# Keep only the two stages we want
subset_f <- filter(data_f, stage %in% c("s00h", "s72h"))

# Plot stacked bars
plot_f <- ggplot(subset_f, aes(x = stage, y = proportion, fill = cell_type)) +
  geom_col(width = 0.6, position = "stack") +   # stacked bars
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + # start at 0
  scale_fill_manual(values = c("B" = "pink", "Plasma" = "navy")) + # match example colors
  labs(x = "Stage", y = "Proportion", fill = "Cell Type") +
  theme_classic() +
  theme(
    legend.position = "inside",                  # legend inside plot
    legend.position.inside = c(0.65, 0.85),      # adjust location (x,y)
    legend.box = "vertical",                     # stack legend items vertically
    legend.background = element_rect(color = "black", linetype = "solid")
  )

# Show in plot pane
print(plot_f)

# Save to file
ggsave("plot_f_bar.png", plot = plot_f, width = 8, height = 6, dpi = 300)

# --- Plot G: Network Graph ---
data_g <- read_excel(file_path, sheet = "g")
rownames(data_g) <- data_g[[1]]
data_g <- data_g[, -1]

net <- graph_from_adjacency_matrix(as.matrix(data_g),
                                   mode = "directed", weighted = TRUE, diag = FALSE)

png("plot_g_network.png", width = 1000, height = 1000, res = 150)
plot(net, layout = layout_nicely,
     vertex.size = 25, vertex.color = "pink", vertex.label.color = "black",
     edge.color = "grey", edge.width = 2)
dev.off()