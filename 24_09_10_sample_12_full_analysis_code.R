#Step1: Create a csv file with cell . clone, clone count clone class, cell ttpe, pixel info asample name and NG names
library(dplyr)
require(stringr)
require(tidyverse)
require(xlsx)
require(Seurat)
library(STutility)
#added the count number along with the clone in this version


se2<- readRDS("~/Downloads/Everything Oxford/Lab Work/Spatial Transcriptomics/SpatialVDJ_forZenodo/data/tonsil/0_integrated/tonsil_se_PartIV.rds")

head(se2@meta.data)

st <- GetStaffli(se2)
head(st@meta.data)


library(Matrix)

# Assuming your matrix is stored in a variable named `clone_matrix`
clone_matrix <- se2@assays$CLON # replace with your actual matrix

clone_matrix <- as.data.frame(clone_matrix@counts)
# Calculate the number of clones per cell
clone_count_per_cell <- colSums(clone_matrix > 0)

# Identify cells with more than one clone
cells_with_multiple_clones <- colnames(clone_matrix)[clone_count_per_cell > 1]

# Extract both clone names and their corresponding counts for each cell
clones_and_counts_in_each_cell <- apply(clone_matrix[, cells_with_multiple_clones, drop = FALSE], 2, function(x) {
  # Get the clone names where the count is greater than 0
  clones <- rownames(clone_matrix)[x > 0]
  
  # Get the counts for those clones
  counts <- x[x > 0]
  
  # Return a data frame with clone names and their counts for the current cell
  data.frame(Clone = clones, Clone_Count = counts)
})

# Combine the list of data frames into a single data frame
result <- do.call(rbind, lapply(names(clones_and_counts_in_each_cell), function(cell) {
  # For each cell, create a data frame with Cell name, Clone, and Clone_Count
  data.frame(Cell = cell, clones_and_counts_in_each_cell[[cell]])
}))

# Print the table
print(result)

# Separate IG and TR clones
result$clone_class <- ifelse(grepl("^IG", result$Clone), "IG", "TR")
result[1:10,]


# Identify cells with both IG and TR clones
cells_with_ig_tr <- result %>%
  group_by(Cell) %>%
  summarize(ig_count = sum(clone_class == "IG"),
            tr_count = sum(clone_class == "TR")) 



#Match the names of the spatial barcode with the name of cell


result$cell_type <- se2@meta.data[result$Cell, "clustergroup"]


# Add cell type, pixel location, and sample information to the result
result$cell_type <- se2@meta.data[result$Cell, "clustergroup"]
result$pixel_x <- st@meta.data[result$Cell, "pixel_x"]
result$pixel_y <- st@meta.data[result$Cell, "pixel_y"]
result$sample <- se2@meta.data[result$Cell, "section"]
head(result)
View(result)
View(result)

x<- result[order(result$Clone_Count, decreasing = TRUE),]
x
x[x$clone_class=="TR",]
x[x$Clone== "IGLCclone58",]
table(duplicated(x[x$clone_class=="TR",]$Clone))

#Lets focus on sample 12
library(stringr)
result$isotype <- substr(result$Clone,1,4)
result12 <- result[result$sample == "12",]
head(result12)
NH_clust <- read.csv("~/Downloads/Everything Oxford/Lab Work/second project immune centrality/try codes and data/sample 12/Neighbour_cell_type.csv")
head(NH_clust)

rownames(NH_clust) <- NH_clust$Cell
result12$NH_names <- NH_clust[result12$Cell,"NH_names"]
head(result12)
table(result12$NH_names)
write.csv(result12,"section12_analysis.csv")


#separate the BCR and TCR and write it in different file

result12_bcr <- result12[result12$clone_class =="IG",]
write.csv(result12_bcr, "section12_bcr.csv")
result12_tcr <- result12[result12$clone_class =="TR",]
write.csv(result12_tcr, "section12_tcr.csv")


#Divide between heavy and light chain
result12_bcr$heavy_light <- substr(result12_bcr$Clone,1,3)
result12_bcr$heavy_light <- ifelse(result12_bcr$heavy_light == "IGH", "heavy","light")
table(result12_bcr$heavy_light)
#step 2: calculating percentage per isotype for each of the NH_sections and then plotting boxplot based on cell_type and the doing manova

library(MASS)

#Start with heave chain
result12_bcrh<- result12_bcr[result12_bcr$heavy_light =="heavy",]
result12_bcrl<- result12_bcr[!result12_bcr$heavy_light =="heavy",]
#Calculate the percentage of clones for each isotype in each NH_section

isotype_percentages <- result12_bcrh %>%
  group_by(NH_names, cell_type) %>%
  mutate(total_clones_per_section = n()) %>%  # Calculate total clones in each section per cell type
  group_by(NH_names, isotype, cell_type, total_clones_per_section) %>%  # Now group by isotype
  summarize(clone_count = n()) %>%  # Count clones for each isotype in each section
  mutate(percentage = (clone_count / total_clones_per_section) * 100)  # Calculate percentage
head(isotype_percentages)

pdf("heatmap of percentage heavy chain of BCR.pdf", height = 11, width =7)
ggplot(isotype_percentages, aes(x = isotype, y = NH_names, fill = percentage)) +
  geom_tile(color = "white") +  # Use geom_tile to create the heatmap tiles
  
  # Use viridis color scale for better publication quality
  scale_fill_viridis(option = "magma", direction = -1, name = "Percentage") +
  
  labs(title = "Heatmap of Heavychain Isotype Percentages",
       x = "Isotype",
       y = "Neighborhood") +
  
  # Minimal theme with adjustments for publication-ready aesthetics
  theme_minimal(base_size = 14) +
  
  # Further customization for clean visuals
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Rotate x-axis labels
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 14, face = "bold"),  # Bold axis titles
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Center the plot title
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid = element_blank(),  # Remove grid lines
    panel.border = element_blank()  # No border around the plot
  )
dev.off()


library(multcompView)
data <- isotype_percentages

# Convert isotype and cell_type to factors
data$isotype <- as.factor(data$isotype)
data$cell_type <- as.factor(data$cell_type)

# Perform ANOVA for each isotype and extract p-values and significance
anova_results <- list()
p_values <- data.frame(isotype = character(), p_value = numeric(), stringsAsFactors = FALSE)
significant_labels <- data.frame()

for (iso in levels(data$isotype)) {
  subset_data <- subset(data, isotype == iso)
  
  # Perform ANOVA
  anova_fit <- aov(percentage ~ cell_type, data = subset_data)
  summary_fit <- summary(anova_fit)
  
  # Extract the p-value from the ANOVA summary
  p_value <- summary_fit[[1]][["Pr(>F)"]][1]
  p_values <- rbind(p_values, data.frame(isotype = iso, p_value = p_value))
  anova_results[[iso]] <- summary_fit
  
  # Perform Tukey's HSD post-hoc test for pairwise comparisons
  tukey_fit <- TukeyHSD(anova_fit)
  
  # Get significance groups (compact letter display)
  tukey_summary <- tukey_fit$cell_type
  significant_pairs <- multcompLetters4(anova_fit, tukey_fit)$cell_type
  
  # Add the significant labels to the data
  labels <- data.frame(
    cell_type = names(significant_pairs$Letters),
    label = significant_pairs$Letters,
    isotype = iso
  )
  significant_labels <- rbind(significant_labels, labels)
}

# Print ANOVA p-values and Tukey's HSD results for each isotype
print(p_values)
print(significant_labels)

# Merge significance labels with the original data
data_with_labels <- merge(data, significant_labels, by = c("cell_type", "isotype"))

library(ggplot2)
library(scales)  # For controlling scales and color schemes

# Define a custom color palette with subtle, nature-friendly colors
nature_colors <- c("#4F81BD", "#C0504D", "#9BBB59", "#8064A2", "#F79646")

pdf("heavy chain percentage distribution acorss cell type.pdf", height = 7, width = 8)
# Plot percentage by cell_type for each isotype using ggplot2
ggplot(data_with_labels, aes(x = cell_type, y = percentage, fill = cell_type)) +
  geom_boxplot(outlier.shape = NA, color = "black", size = 0.75, alpha = 0.6) +  # Black outline for boxplot
  geom_jitter(width = 0.2, size = 1.5, aes(color = cell_type), alpha = 0.8) +  # Subtle jittered points
  scale_fill_manual(values = nature_colors) +  # Use muted color palette
  scale_color_manual(values = nature_colors) + # Apply same palette for jitter points
  facet_wrap(~ isotype, nrow = 2) +  # Split by isotype
  labs(title = "Heavy Chain Isotype Percentage Distribution Across Cell Types",
       x = "Cell Type",
       y = "Percentage (%)") +
  theme_minimal(base_size = 14, base_family = "Helvetica") +  # Clean font and size
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10, color = "black"),  # Rotate x-axis labels
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),  # Centered title
    strip.text = element_text(size = 12, face = "bold"),  # Facet labels
    panel.grid.major = element_blank(),  # Remove major gridlines for cleaner look
    panel.grid.minor = element_line(color = "grey90", size = 0.25),  # Subtle minor gridlines for reference
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),  # Thin black border around panels
    legend.position = "top",  # Place legend on top
    legend.title = element_blank(),  # Remove legend title for simplicity
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8, 'lines')  # Make legend key smaller
  ) +
  # Annotate p-values on the plot, positioned in each facet
  geom_text(data = p_values, aes(label = paste("p =", round(p_value, 5))),
            x = Inf, y = Inf, hjust = 1.1, vjust = 2, size = 4, color = "black", inherit.aes = FALSE)
dev.off()
#Light chain percentage distribution
# Assuming the dataset has columns: 'Clone', 'isotype', 'NH_names', 'cell_type'
isotype_percentages <- result12_bcrl %>%
  group_by(NH_names, cell_type) %>%
  mutate(total_clones_per_section = n()) %>%  # Calculate total clones in each section per cell type
  group_by(NH_names, isotype, cell_type, total_clones_per_section) %>%  # Now group by isotype
  summarize(clone_count = n()) %>%  # Count clones for each isotype in each section
  mutate(percentage = (clone_count / total_clones_per_section) * 100)  # Calculate percentage
head(isotype_percentages)

pdf("heatmap of percentage light chain of BCR.pdf", height = 11, width =7)
ggplot(isotype_percentages, aes(x = isotype, y = NH_names, fill = percentage)) +
  geom_tile(color = "white") +  # Use geom_tile to create the heatmap tiles
  
  # Use viridis color scale for better publication quality
  scale_fill_viridis(option = "magma", direction = -1, name = "Percentage") +
  
  labs(title = "Heatmap of Lightchain Isotype Percentages",
       x = "Isotype",
       y = "Neighborhood") +
  
  # Minimal theme with adjustments for publication-ready aesthetics
  theme_minimal(base_size = 14) +
  
  # Further customization for clean visuals
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Rotate x-axis labels
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 14, face = "bold"),  # Bold axis titles
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Center the plot title
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid = element_blank(),  # Remove grid lines
    panel.border = element_blank()  # No border around the plot
  )
dev.off()



library(multcompView)
data <- isotype_percentages

# Convert isotype and cell_type to factors
data$isotype <- as.factor(data$isotype)
data$cell_type <- as.factor(data$cell_type)

# Perform ANOVA for each isotype and extract p-values and significance
anova_results <- list()
p_values <- data.frame(isotype = character(), p_value = numeric(), stringsAsFactors = FALSE)
significant_labels <- data.frame()

for (iso in levels(data$isotype)) {
  subset_data <- subset(data, isotype == iso)
  
  # Perform ANOVA
  anova_fit <- aov(percentage ~ cell_type, data = subset_data)
  summary_fit <- summary(anova_fit)
  
  # Extract the p-value from the ANOVA summary
  p_value <- summary_fit[[1]][["Pr(>F)"]][1]
  p_values <- rbind(p_values, data.frame(isotype = iso, p_value = p_value))
  anova_results[[iso]] <- summary_fit
  
  # Perform Tukey's HSD post-hoc test for pairwise comparisons
  tukey_fit <- TukeyHSD(anova_fit)
  
  # Get significance groups (compact letter display)
  tukey_summary <- tukey_fit$cell_type
  significant_pairs <- multcompLetters4(anova_fit, tukey_fit)$cell_type
  
  # Add the significant labels to the data
  labels <- data.frame(
    cell_type = names(significant_pairs$Letters),
    label = significant_pairs$Letters,
    isotype = iso
  )
  significant_labels <- rbind(significant_labels, labels)
}

# Print ANOVA p-values and Tukey's HSD results for each isotype
print(p_values)
print(significant_labels)

# Merge significance labels with the original data
data_with_labels <- merge(data, significant_labels, by = c("cell_type", "isotype"))

library(ggplot2)
library(scales)  # For controlling scales and color schemes

# Define a custom color palette with subtle, nature-friendly colors
nature_colors <- c("#4F81BD", "#C0504D", "#9BBB59", "#8064A2", "#F79646")

pdf("light chain percentage distribution acorss cell type.pdf", height = 7, width = 8)
# Plot percentage by cell_type for each isotype using ggplot2
ggplot(data_with_labels, aes(x = cell_type, y = percentage, fill = cell_type)) +
  geom_boxplot(outlier.shape = NA, color = "black", size = 0.75, alpha = 0.6) +  # Black outline for boxplot
  geom_jitter(width = 0.2, size = 1.5, aes(color = cell_type), alpha = 0.8) +  # Subtle jittered points
  scale_fill_manual(values = nature_colors) +  # Use muted color palette
  scale_color_manual(values = nature_colors) + # Apply same palette for jitter points
  facet_wrap(~ isotype, nrow = 2) +  # Split by isotype
  labs(title = "Light Chain Isotype Percentage Distribution Across Cell Types",
       x = "Cell Type",
       y = "Percentage (%)") +
  theme_minimal(base_size = 14, base_family = "Helvetica") +  # Clean font and size
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10, color = "black"),  # Rotate x-axis labels
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),  # Centered title
    strip.text = element_text(size = 12, face = "bold"),  # Facet labels
    panel.grid.major = element_blank(),  # Remove major gridlines for cleaner look
    panel.grid.minor = element_line(color = "grey90", size = 0.25),  # Subtle minor gridlines for reference
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),  # Thin black border around panels
    legend.position = "top",  # Place legend on top
    legend.title = element_blank(),  # Remove legend title for simplicity
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8, 'lines')  # Make legend key smaller
  ) +
  # Annotate p-values on the plot, positioned in each facet
  geom_text(data = p_values, aes(label = paste("p =", round(p_value, 5))),
            x = Inf, y = Inf, hjust = 1.1, vjust = 2, size = 4, color = "black", inherit.aes = FALSE)
dev.off()




#Tcell receptor percentage distribution
#as some of the data has clon instaed. need to remove them
table(result12_tcr$isotype)
result12_tcr <- result12_tcr[result12_tcr$isotype != "clon",]
# Assuming the dataset has columns: 'Clone', 'isotype', 'NH_names', 'cell_type'
isotype_percentages <- result12_tcr %>%
  group_by(NH_names, cell_type) %>%
  mutate(total_clones_per_section = n()) %>%  # Calculate total clones in each section per cell type
  group_by(NH_names, isotype, cell_type, total_clones_per_section) %>%  # Now group by isotype
  summarize(clone_count = n()) %>%  # Count clones for each isotype in each section
  mutate(percentage = (clone_count / total_clones_per_section) * 100)  # Calculate percentage
head(isotype_percentages)

# Create the heatmap using viridis color scale
pdf("heatmap of percentage TCR.pdf", height = 11, width =7)
ggplot(isotype_percentages, aes(x = isotype, y = NH_names, fill = percentage)) +
  geom_tile(color = "white") +  # Use geom_tile to create the heatmap tiles
  
  # Use viridis color scale for better publication quality
  scale_fill_viridis(option = "magma", direction = -1, name = "Percentage") +
  
  labs(title = "Heatmap of TCR Isotype Percentages",
       x = "Isotype",
       y = "Neighborhood") +
  
  # Minimal theme with adjustments for publication-ready aesthetics
  theme_minimal(base_size = 14) +
  
  # Further customization for clean visuals
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Rotate x-axis labels
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 14, face = "bold"),  # Bold axis titles
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Center the plot title
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid = element_blank(),  # Remove grid lines
    panel.border = element_blank()  # No border around the plot
  )
dev.off()

library(multcompView)
data <- isotype_percentages

# Convert isotype and cell_type to factors
data$isotype <- as.factor(data$isotype)
data$cell_type <- as.factor(data$cell_type)

# Perform ANOVA for each isotype and extract p-values and significance
anova_results <- list()
p_values <- data.frame(isotype = character(), p_value = numeric(), stringsAsFactors = FALSE)
significant_labels <- data.frame()

for (iso in levels(data$isotype)) {
  subset_data <- subset(data, isotype == iso)
  
  # Perform ANOVA
  anova_fit <- aov(percentage ~ cell_type, data = subset_data)
  summary_fit <- summary(anova_fit)
  
  # Extract the p-value from the ANOVA summary
  p_value <- summary_fit[[1]][["Pr(>F)"]][1]
  p_values <- rbind(p_values, data.frame(isotype = iso, p_value = p_value))
  anova_results[[iso]] <- summary_fit
  
  # Perform Tukey's HSD post-hoc test for pairwise comparisons
  tukey_fit <- TukeyHSD(anova_fit)
  
  # Get significance groups (compact letter display)
  tukey_summary <- tukey_fit$cell_type
  significant_pairs <- multcompLetters4(anova_fit, tukey_fit)$cell_type
  
  # Add the significant labels to the data
  labels <- data.frame(
    cell_type = names(significant_pairs$Letters),
    label = significant_pairs$Letters,
    isotype = iso
  )
  significant_labels <- rbind(significant_labels, labels)
}

# Print ANOVA p-values and Tukey's HSD results for each isotype
print(p_values)
print(significant_labels)

# Merge significance labels with the original data
data_with_labels <- merge(data, significant_labels, by = c("cell_type", "isotype"))

library(ggplot2)
library(scales)  # For controlling scales and color schemes

# Define a custom color palette with subtle, nature-friendly colors
nature_colors <- c("#4F81BD", "#C0504D", "#9BBB59", "#8064A2", "#F79646")

pdf("TCR percentage distribution acorss cell type.pdf", height = 7, width = 8)
# Plot percentage by cell_type for each isotype using ggplot2
ggplot(data_with_labels, aes(x = cell_type, y = percentage, fill = cell_type)) +
  geom_boxplot(outlier.shape = NA, color = "black", size = 0.75, alpha = 0.6) +  # Black outline for boxplot
  geom_jitter(width = 0.2, size = 1.5, aes(color = cell_type), alpha = 0.8) +  # Subtle jittered points
  scale_fill_manual(values = nature_colors) +  # Use muted color palette
  scale_color_manual(values = nature_colors) + # Apply same palette for jitter points
  facet_wrap(~ isotype, nrow = 2) +  # Split by isotype
  labs(title = "Light Chain Isotype Percentage Distribution Across Cell Types",
       x = "Cell Type",
       y = "Percentage (%)") +
  theme_minimal(base_size = 14, base_family = "Helvetica") +  # Clean font and size
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10, color = "black"),  # Rotate x-axis labels
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),  # Centered title
    strip.text = element_text(size = 12, face = "bold"),  # Facet labels
    panel.grid.major = element_blank(),  # Remove major gridlines for cleaner look
    panel.grid.minor = element_line(color = "grey90", size = 0.25),  # Subtle minor gridlines for reference
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),  # Thin black border around panels
    legend.position = "top",  # Place legend on top
    legend.title = element_blank(),  # Remove legend title for simplicity
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8, 'lines')  # Make legend key smaller
  ) +
  # Annotate p-values on the plot, positioned in each facet
  geom_text(data = p_values, aes(label = paste("p =", round(p_value, 5))),
            x = Inf, y = Inf, hjust = 1.1, vjust = 2, size = 4, color = "black", inherit.aes = FALSE)
dev.off()


# Load necessary libraries
library(ggplot2)
library(viridis)  # For the viridis color palette

# Assuming your dataset is named 'isotype_percentages'

#Jaccard Overlap index for bcr

data<- result12_bcr
# Step 1: Initialize a symmetric matrix with NA values for Jaccard Index
neighborhood_clones <- data %>%
  group_by(NH_names) %>%
  summarize(clones = list(unique(Clone)))

# Create an empty matrix to hold Jaccard indices (symmetric)
nh_names <- neighborhood_clones$NH_names
jaccard_matrix <- matrix(NA, nrow = length(nh_names), ncol = length(nh_names), 
                         dimnames = list(nh_names, nh_names))

# Step 2: Calculate Jaccard Index and fill symmetric matrix
for (i in 1:nrow(jaccard_matrix)) {
  for (j in i:nrow(jaccard_matrix)) {
    clones1 <- neighborhood_clones$clones[[i]]
    clones2 <- neighborhood_clones$clones[[j]]
    jaccard <- jaccard_index(clones1, clones2)
    
    # Fill both sides of the matrix to ensure symmetry
    jaccard_matrix[i, j] <- jaccard
    jaccard_matrix[j, i] <- jaccard
  }
}

# Step 3: Ensure diagonal is filled with 1s (self-similarity)
diag(jaccard_matrix) <- NA

# Step 4: Only keep the upper triangle (or lower triangle) to show half the plot
jaccard_matrix[lower.tri(jaccard_matrix)] <- NA  # Optional: Keep upper triangle only

# Step 5: Melt the matrix for ggplot2
jaccard_long <- melt(jaccard_matrix, na.rm = TRUE)  # Remove NA values for plotting

midpoint <- (range(jaccard_long$value)[1]+range(jaccard_long$value)[2])/2
# Step 6: Plot heatmap using ggplot2 with a blue-white-red gradient
pdf("BCR jaccard overlap_dff_color.pdf", height =16, width = 16)
ggplot(jaccard_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +  # White grid lines
  scale_fill_gradient2(low = "white", mid = "blue", high = "red", 
                       midpoint = midpoint, name = "Jaccard") +  # Blue-White-Red gradient
  labs(
    title = "BCR Jaccard Overlap Heatmap",
    x = "Neighborhood 1",
    y = "Neighborhood 2"
  ) +
  scale_x_discrete(position = "top") +  # Move x-axis to top
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 1, margin = margin(t = 10)),  # Adjust margin
    axis.title.x = element_text(margin = margin(b = 20)),  # Add space between x-axis labels and plot
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.title.y = element_blank(),  # Optional: remove y-axis label if needed
    plot.margin = margin(t = 30, r = 10, b = 10, l = 10)  # Adjust plot margins
  ) +
  coord_fixed()  # Ensure square tiles
dev.off()


#jaccard overlap matrix across different cell types
# Assuming your dataset is stored in a variable called `jaccard_matrix`
library(reshape2)

# Step 2: Convert the matrix to a long format using `melt`
long_format <- melt(jaccard_matrix, varnames = c("NH1", "NH2"), value.name = "Jaccard_Index")

# Step 3: Remove rows where the Jaccard Index is NA or where NH1 == NH2 (if you don't want to include self-comparisons)
long_format <- na.omit(long_format)  # Remove NA values
long_format <- long_format[long_format$NH1 != long_format$NH2, ]  # Remove self-comparisons

# View the long format data
head(long_format)

# Step 1: Extract the part before the underscore in NH1 and NH2
long_format <- long_format %>%
  mutate(
    NH1_clean = sub("_.*", "", NH1),  # Remove everything after the underscore in NH1
    NH2_clean = sub("_.*", "", NH2)   # Remove everything after the underscore in NH2
  )

# Step 2: Create a new column `celltype` by concatenating the cleaned NH1 and NH2 with a dash
long_format <- long_format %>%
  mutate(celltype = paste(NH1_clean, NH2_clean, sep = "-"))

# Step 3: Remove intermediate columns if needed
long_format <- long_format[,c("NH1", "NH2", "Jaccard_Index", "celltype")]

# View the final dataset
table(long_format$celltype)

# Load necessary libraries
library(ggplot2)
library(ggpubr)  # For p-value annotations
library(dplyr)

# Assuming `data` is the dataset with 'celltype', 'Jaccard_Index'
# Step 1: Perform ANOVA for each celltype
anova_result <- aov(Jaccard_Index ~ celltype, data = long_format)
summary(anova_result)

# Extract p-value from ANOVA result
p_value <- summary(anova_result)[[1]]$`Pr(>F)`[1]

# Step 2: Define a custom color palette with 15 distinct colors
custom_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
                   "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", 
                   "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")

# Step 3: Create the boxplot with exact p-value annotation and 15 different colors
pdf("jaccard overlap index of BCR.pdf", height = 10, width =10)
ggplot(long_format, aes(x = celltype, y = Jaccard_Index, fill = celltype)) +
  geom_boxplot(outlier.shape = NA, color = "black", lwd = 0.7) +  # Black border, slightly thicker line
  geom_jitter(width = 0.1, alpha = 0.7, color = "black", size = 2) +  # Jitter points with black outlines for visibility
  scale_fill_manual(values = custom_colors) +  # Apply custom color palette
  labs(
    title = "Jaccard Overlap Index of BCR",
    x = "Celltype",
    y = "Jaccard Index"
  ) +
  # Display the p-value in the top right corner
  annotate("text", x = Inf, y = Inf, label = paste("p-value:", signif(p_value, digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 6, color = "black", fontface = "bold") +
  theme_minimal() +  # Minimal theme for a clean layout
  theme(
    text = element_text(family = "Helvetica", size = 14),  # Standard clean font
    axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 10)),  # Bold x-axis title
    axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 10)),  # Bold y-axis title
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 12),  # Rotated x-axis labels for readability
    axis.text.y = element_text(size = 12),  # Y-axis label size
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Center the plot title
    plot.margin = margin(15, 15, 15, 15),  # Add margin around the plot
    panel.grid.major = element_blank(),  # Remove major grid lines for clarity
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(size = 0.5),  # Add axis lines for clarity
    legend.position = "none"  # Remove legend since we are using x-axis labels for celltype
  )
dev.off()


#Jaccard overlap index for TCR

data<- result12_tcr
# Step 1: Initialize a symmetric matrix with NA values for Jaccard Index
neighborhood_clones <- data %>%
  group_by(NH_names) %>%
  summarize(clones = list(unique(Clone)))

# Create an empty matrix to hold Jaccard indices (symmetric)
nh_names <- neighborhood_clones$NH_names
jaccard_matrix <- matrix(NA, nrow = length(nh_names), ncol = length(nh_names), 
                         dimnames = list(nh_names, nh_names))

# Step 2: Calculate Jaccard Index and fill symmetric matrix
for (i in 1:nrow(jaccard_matrix)) {
  for (j in i:nrow(jaccard_matrix)) {
    clones1 <- neighborhood_clones$clones[[i]]
    clones2 <- neighborhood_clones$clones[[j]]
    jaccard <- jaccard_index(clones1, clones2)
    
    # Fill both sides of the matrix to ensure symmetry
    jaccard_matrix[i, j] <- jaccard
    jaccard_matrix[j, i] <- jaccard
  }
}

# Step 3: Ensure diagonal is filled with 1s (self-similarity)
diag(jaccard_matrix) <- NA

# Step 4: Only keep the upper triangle (or lower triangle) to show half the plot
jaccard_matrix[lower.tri(jaccard_matrix)] <- NA  # Optional: Keep upper triangle only

# Step 5: Melt the matrix for ggplot2
jaccard_long <- melt(jaccard_matrix, na.rm = TRUE)  # Remove NA values for plotting

midpoint <- (range(jaccard_long$value)[1]+range(jaccard_long$value)[2])/2
# Step 6: Plot heatmap using ggplot2 with a blue-white-red gradient
pdf("TCR jaccard overlap.pdf", height =16, width = 16)
ggplot(jaccard_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +  # White grid lines
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = midpoint, name = "Jaccard") +  # Blue-White-Red gradient
  labs(
    title = "TCR Jaccard Overlap Heatmap",
    x = "Neighborhood 1",
    y = "Neighborhood 2"
  ) +
  scale_x_discrete(position = "top") +  # Move x-axis to top
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 1, margin = margin(t = 10)),  # Adjust margin
    axis.title.x = element_text(margin = margin(b = 20)),  # Add space between x-axis labels and plot
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.title.y = element_blank(),  # Optional: remove y-axis label if needed
    plot.margin = margin(t = 30, r = 10, b = 10, l = 10)  # Adjust plot margins
  ) +
  coord_fixed()  # Ensure square tiles
dev.off()

#jaccard overlap matrix across different cell types
# Assuming your dataset is stored in a variable called `jaccard_matrix`
library(reshape2)

# Step 2: Convert the matrix to a long format using `melt`
long_format <- melt(jaccard_matrix, varnames = c("NH1", "NH2"), value.name = "Jaccard_Index")

# Step 3: Remove rows where the Jaccard Index is NA or where NH1 == NH2 (if you don't want to include self-comparisons)
long_format <- na.omit(long_format)  # Remove NA values
long_format <- long_format[long_format$NH1 != long_format$NH2, ]  # Remove self-comparisons

# View the long format data
head(long_format)

# Step 1: Extract the part before the underscore in NH1 and NH2
long_format <- long_format %>%
  mutate(
    NH1_clean = sub("_.*", "", NH1),  # Remove everything after the underscore in NH1
    NH2_clean = sub("_.*", "", NH2)   # Remove everything after the underscore in NH2
  )

# Step 2: Create a new column `celltype` by concatenating the cleaned NH1 and NH2 with a dash
long_format <- long_format %>%
  mutate(celltype = paste(NH1_clean, NH2_clean, sep = "-"))

# Step 3: Remove intermediate columns if needed
long_format <- long_format[,c("NH1", "NH2", "Jaccard_Index", "celltype")]

# View the final dataset
table(long_format$celltype)

# Load necessary libraries
library(ggplot2)
library(ggpubr)  # For p-value annotations
library(dplyr)

# Assuming `data` is the dataset with 'celltype', 'Jaccard_Index'
# Step 1: Perform ANOVA for each celltype
anova_result <- aov(Jaccard_Index ~ celltype, data = long_format)
summary(anova_result)

# Extract p-value from ANOVA result
p_value <- summary(anova_result)[[1]]$`Pr(>F)`[1]

# Step 2: Define a custom color palette with 15 distinct colors
custom_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
                   "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", 
                   "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")

# Step 3: Create the boxplot with exact p-value annotation and 15 different colors
pdf("jaccard overlap index of TCR.pdf", height = 10, width =10)
ggplot(long_format, aes(x = celltype, y = Jaccard_Index, fill = celltype)) +
  geom_boxplot(outlier.shape = NA, color = "black", lwd = 0.7) +  # Black border, slightly thicker line
  geom_jitter(width = 0.1, alpha = 0.7, color = "black", size = 2) +  # Jitter points with black outlines for visibility
  scale_fill_manual(values = custom_colors) +  # Apply custom color palette
  labs(
    title = "Jaccard Overlap Index of TCR",
    x = "Celltype",
    y = "Jaccard Index"
  ) +
  # Display the p-value in the top right corner
  annotate("text", x = Inf, y = Inf, label = paste("p-value:", signif(p_value, digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 3, color = "black", fontface = "bold") +
  theme_minimal() +  # Minimal theme for a clean layout
  theme(
    text = element_text(family = "Helvetica", size = 14),  # Standard clean font
    axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 10)),  # Bold x-axis title
    axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 10)),  # Bold y-axis title
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 12),  # Rotated x-axis labels for readability
    axis.text.y = element_text(size = 12),  # Y-axis label size
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Center the plot title
    plot.margin = margin(15, 15, 15, 15),  # Add margin around the plot
    panel.grid.major = element_blank(),  # Remove major grid lines for clarity
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(size = 0.5),  # Add axis lines for clarity
    legend.position = "none"  # Remove legend since we are using x-axis labels for celltype
  )
dev.off()

####################################################
#merging mixcr output with the result12 data

main_vdj <- read.csv("~/Downloads/Everything Oxford/Lab Work/second project immune centrality/try codes and data/sample 32/VDJ analysis across zone/MIXCR_pt086_Pacbio_run86_Tonsil-Mixcr_all.tsv", sep = "\t")

#remove the the ones that are empty
main_vdj <- main_vdj[main_vdj$allCHitsWithScore != "",]
#as i want to mathc the main_cdj file with the subsamples file with neighbours i need to make the names comparable
main_vdj$Clone <- paste(str_sub(main_vdj$allCHitsWithScore, 1,4),"clone",main_vdj$cloneId, sep = "")
library(dplyr)

nh_section <- unique(result12$NH_names)
nh_section

# Initialize an empty dataframe before the loop
all_final_res <- data.frame()

for (i in 1:length(nh_section)) {
  current_section <- nh_section[i]
  
  # Subset the result32 dataframe
  res_sub <- result12[result12$NH_names == current_section,]
  res_sub <- res_sub[, c("Clone", "Clone_Count", "clone_class", "NH_names")]
  
  # Summarize the clone count by group
  res_sub_summed <- res_sub %>%
    group_by(Clone, NH_names, clone_class) %>%
    summarise(Clone_Count = sum(Clone_Count)) %>%
    as.data.frame()
  
  # Merge with main_vdj dataframe
  final_res <- merge(res_sub_summed, main_vdj, by = "Clone")
  
  # Remove unwanted columns
  final_res$cloneCount <- NULL
  final_res$cloneFraction <- NULL
  
  # Recalculate cloneFraction
  final_res$cloneFraction <- final_res$Clone_Count / sum(final_res$Clone_Count)
  
  # Rename column
  colnames(final_res)[colnames(final_res) == "Clone_Count"] <- "cloneCount"
  
  # Add NH_section column
  final_res$NH_section <- current_section
  
  # Append the current final_res to the all_final_res dataframe
  all_final_res <- rbind(all_final_res, final_res)
  
}

View(all_final_res)


# Step 1: Extract the part before the underscore in NH1 and NH2
all_final_res <- all_final_res %>%
  mutate(
    cell_type = sub("_.*", "", NH_names))


# Extract clone type
all_final_res$isotype <- substr(all_final_res$allCHitsWithScore, 1, 4)


# a. Overall Isotype Usages
data<- all_final_res
# Calculate total counts for each clone_type within each NH_section
isotype_usage <- data %>%
  group_by(NH_section, clone_type) %>%
  summarise(total_count = sum(cloneCount, na.rm = TRUE)) %>%
  ungroup()

# Calculate total counts for each NH_section
total_counts_per_section <- isotype_usage %>%
  group_by(NH_section) %>%
  summarise(total_section_count = sum(total_count))

# Join with isotype_usage to calculate proportions
isotype_usage <- isotype_usage %>%
  left_join(total_counts_per_section, by = "NH_section") %>%
  mutate(proportion = total_count / total_section_count)

# Custom dark color palette for 8 clone types
dark_colors <- c("#1B4F72", "#196F3D", "#943126", "#6C3483", "#117864", "#D35400", "#7D6608", "#2C3E50")

# Plot Isotype Usage
pdf("stacked proportion across neighborhood.pdf", height = 6, width =12)
ggplot(isotype_usage, aes(x = NH_section, y = proportion, fill = clone_type)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.3) +  # Thin black outline for bars
  scale_fill_manual(values = dark_colors) +  # Apply custom dark color palette
  labs(title = "Isotype Usage by Neighborhood", x = "Neighborhood", y = "Proportions") +
  theme_minimal(base_size = 14) +  # Larger base font size for readability
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),  # Rotate x-axis labels by 90 degrees
    axis.text.y = element_text(size = 8),  # Adjust font size for y-axis
    axis.title = element_text(size = 14, face = "bold"),  # Bold axis titles
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Centered and bold title
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_blank(),  # Remove major gridlines for a cleaner look
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)  # Thin black border around the plot
  )
dev.off()

#CDR3 length distribution across different isotype

# Load required libraries
library(dplyr)
library(ggplot2)
library(multcompView)  # For compact letter display (Tukey's HSD)
library(scales)        # For color scales

# Define a custom color palette with subtle, nature-friendly colors
nature_colors <- c("#4F81BD", "#C0504D", "#9BBB59", "#8064A2", "#F79646")

# Assuming the dataset is named 'all_final_res' and contains the columns: NH_section, cell_type, aaSeqCDR3, isotype
data <- all_final_res

# Step 1: Calculate the mean CDR3 length for each isotype in each neighborhood
data <- data %>%
  mutate(CDR3_length = nchar(as.character(aaSeqCDR3))) %>%  # Calculate CDR3 length
  group_by(NH_section, cell_type, isotype) %>%  # Group by NH_section, cell_type, and isotype
  summarise(mean_CDR3_length = mean(CDR3_length, na.rm = TRUE)) %>%
  ungroup()

# Step 2: Perform ANOVA for each isotype and extract p-values and significance labels
anova_results <- list()
p_values <- data.frame(isotype = character(), p_value = numeric(), stringsAsFactors = FALSE)
significant_labels <- data.frame()

for (iso in levels(factor(data$isotype))) {
  subset_data <- subset(data, isotype == iso)
  
  # Perform ANOVA
  anova_fit <- aov(mean_CDR3_length ~ cell_type, data = subset_data)
  summary_fit <- summary(anova_fit)
  
  # Extract the p-value from the ANOVA summary
  p_value <- summary_fit[[1]][["Pr(>F)"]][1]
  p_values <- rbind(p_values, data.frame(isotype = iso, p_value = p_value))
  anova_results[[iso]] <- summary_fit
  
  # Perform Tukey's HSD post-hoc test for pairwise comparisons
  tukey_fit <- TukeyHSD(anova_fit)
  
  # Get significance groups (compact letter display)
  significant_pairs <- multcompLetters4(anova_fit, tukey_fit)$cell_type
  
  # Add the significant labels to the data
  labels <- data.frame(
    cell_type = names(significant_pairs$Letters),
    label = significant_pairs$Letters,
    isotype = iso
  )
  significant_labels <- rbind(significant_labels, labels)
}

# Print ANOVA p-values and Tukey's HSD results for each isotype
print(p_values)
print(significant_labels)

# Merge significance labels with the original data
data_with_labels <- merge(data, significant_labels, by = c("cell_type", "isotype"))

# Step 3: Plot mean CDR3 length by cell type for each isotype using ggplot2
pdf("Mean_CDR3_Length_Across_Cell_Types.pdf", height = 8, width = 10)

ggplot(data_with_labels, aes(x = cell_type, y = mean_CDR3_length, fill = cell_type)) +
  geom_boxplot(outlier.shape = NA, color = "black", size = 0.75, alpha = 0.6) +  # Boxplot with black outline
  geom_jitter(width = 0.2, size = 1.5, aes(color = cell_type), alpha = 0.8) +  # Subtle jittered points
  scale_fill_manual(values = nature_colors) +  # Use muted color palette
  scale_color_manual(values = nature_colors) + # Apply same palette for jitter points
  facet_wrap(~ isotype, nrow = 2) +  # Split by isotype
  labs(title = "Mean CDR3 Length Across Cell Types",
       x = "Cell Type",
       y = "Mean CDR3 Length") +
  theme_minimal(base_size = 14, base_family = "Helvetica") +  # Clean font and size
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10, color = "black"),  # Rotate x-axis labels
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),  # Centered title
    strip.text = element_text(size = 12, face = "bold"),  # Facet labels
    panel.grid.major = element_blank(),  # Remove major gridlines for cleaner look
    panel.grid.minor = element_line(color = "grey90", size = 0.25),  # Subtle minor gridlines for reference
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),  # Thin black border around panels
    legend.position = "top",  # Place legend on top
    legend.title = element_blank(),  # Remove legend title for simplicity
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8, 'lines')  # Make legend key smaller
  ) +
  # Annotate p-values on the plot, positioned in each facet
  geom_text(data = p_values, aes(label = paste("p =", round(p_value, 5))),
            x = Inf, y = Inf, hjust = 1.1, vjust = 2, size = 4, color = "black", inherit.aes = FALSE) 
# Close the PDF device
dev.off()

#heatmap of mean CDR3 length across different neighbourhood 

# Load required libraries
library(ggplot2)
library(scales)  # For scale adjustments
library(RColorBrewer)  # For a good color palette

# Define a custom color palette that transitions from white to blue to red
color_palette <- colorRampPalette(c("white", "blue", "red"))(100)

# Define a custom color palette for the heatmap (Blue-White-Red style)
color_palette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)  # Reverse for blue-red scheme

# Generate the heatmap using ggplot2
pdf("CDR3_Length_Heatmap.pdf", height = 7, width = 14)

ggplot(data_with_labels, aes(x = NH_section, y = isotype, fill = mean_CDR3_length)) +
  geom_tile(color = "white", size = 0.2) +  # Add white borders for each cell
  scale_fill_gradientn(colours = color_palette, name = "Mean CDR3 Length", 
                       limits = c(min(data$mean_CDR3_length), max(data$mean_CDR3_length)),
                       oob = squish) +  # Handle out-of-bound values gracefully
  labs(title = "Mean CDR3 Length by NH Section and Isotype",
       x = "NH Section",
       y = "Isotype") +
  theme_minimal(base_size = 14) +  # Minimal theme for a clean look
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),  # Title styling
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 8),  # Rotate x-axis labels
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),  # Remove gridlines for a cleaner heatmap
    panel.border = element_blank(),  # Remove outer border
    legend.position = "right",  # Place the legend on the right
    legend.key.size = unit(1, "cm"),  # Increase the size of the legend key
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )

# Close the PDF device to save the plot
dev.off()

#separate data between TCR and BCR
all_final_res_bcr <- all_final_res[all_final_res$clone_class == "IG",]

#only_BCR
data<- all_final_res_bcr
# Step 1: Calculate the frequency of each CDR3 sequence within each neighborhood
cdr3_counts <- data %>%
  group_by(NH_section, aaSeqCDR3) %>%
  summarize(clone_count = n()) %>%
  ungroup()

# Step 2: Classify sequences into unique, small, medium, large, and hyperexpanded categories
cdr3_counts <- cdr3_counts %>%
  mutate(clone_type = case_when(
    clone_count == 1 ~ "Unique",
    clone_count >= 2 & clone_count < 3 ~ "Expanded",
    clone_count >= 3 ~ "Hyperexpanded"
  ))

# Step 3: Calculate the proportion of each clone type within each neighborhood
isotype_usage <- cdr3_counts %>%
  group_by(NH_section, clone_type) %>%
  summarize(clone_count = sum(clone_count)) %>%  # Total clone count per type per neighborhood
  mutate(proportion = clone_count / sum(clone_count)) %>%  # Calculate proportions
  ungroup()

# Step 4: Create a stacked bar plot
# Load required libraries
library(ggplot2)

# Create a more elegant color palette for clone_type
nature_palette <- c("#4F81BD", "#70AD47", "#ED7D31", "#A5A5A5", "#FFC000")  # Blue, Green, Orange, Gray, Yellow

# Generate the plot
pdf("BCR_clonal_distribution_across_neighbours_nature.pdf", width = 14, height = 8)

ggplot(isotype_usage, aes(x = NH_section, y = proportion, fill = clone_type)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.3) +  # Stacked bar with black outlines
  scale_fill_manual(values = nature_palette, name = "Clone Type") +  # Custom color palette
  labs(title = "BCR Clonal Distribution Across Neighborhoods", 
       x = "Neighborhood", 
       y = "Proportion (%)") +  # Title and axis labels
  theme_minimal(base_size = 14) +  # Minimal theme for clarity
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 10)),  # Centered bold title with margin
    axis.title.x = element_text(size = 12, face = "bold"),  # Bold x-axis title
    axis.title.y = element_text(size = 12, face = "bold"),  # Bold y-axis title
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10, face = "plain", color = "black"),  # X-axis text
    axis.text.y = element_text(size = 10, color = "black"),  # Y-axis text size
    legend.position = "top",  # Legend on top
    legend.title = element_text(size = 14, face = "bold"),  # Bold legend title
    legend.text = element_text(size = 12),  # Adjust legend text size
    legend.key.size = unit(1, 'cm'),  # Increase legend key size for visibility
    panel.grid.major = element_blank(),  # Remove major gridlines for a cleaner look
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8),  # Black border around the plot
    plot.margin = unit(c(1, 1, 1, 1), "cm")  # Add margins around the plot for better readability
  )

# Close the PDF device
dev.off()


isotype_usage <- isotype_usage %>%
  mutate(
    cell_type = sub("_.*", "", NH_section))
head(isotype_usage)

data<- isotype_usage



# Assuming the dataset is called 'data'
# Step 1: Perform ANOVA for each clone_type across different cell_types
anova_results <- list()
p_values <- data.frame(clone_type = character(), p_value = numeric(), stringsAsFactors = FALSE)

# Loop through each unique clone_type and perform ANOVA
for (clone in unique(data$clone_type)) {
  
  # Subset data for the specific clone_type
  subset_data <- subset(data, clone_type == clone)
  
  # Perform ANOVA: Test if proportion differs across cell_type for each clone_type
  anova_fit <- aov(proportion ~ cell_type, data = subset_data)
  
  # Extract the p-value from ANOVA summary
  anova_pvalue <- summary(anova_fit)[[1]]["Pr(>F)"][1,1]
  
  # Store the p-value in p_values
  p_values <- rbind(p_values, data.frame(clone_type = clone, p_value = anova_pvalue))
  print(p_values)
  
  # Store the ANOVA result for each clone_type in the list
  anova_results[[clone]] <- summary(anova_fit)
}

# Print p-values
print(p_values)

# Step 2: Plot the boxplot for each clone_type, showing proportions across cell_type

# Define a custom color palette
color_palette <- c("#4F81BD", "#C0504D", "#9BBB59", "#8064A2", "#F79646")

# Create the boxplot and annotate with p-values
pdf("Proportion of BCR clone expansion across cell type.pdf", height = 7, width =10)
ggplot(data, aes(x = cell_type, y = proportion, fill = cell_type)) +
  geom_boxplot(outlier.shape = NA, color = "black", size = 0.75, alpha = 0.6) +  # Boxplot with light alpha
  geom_jitter(width = 0.2, size = 1.5, aes(color = cell_type), alpha = 0.8) +  # Subtle jittered points
  scale_fill_manual(values = color_palette) +  # Custom color palette for fill
  scale_color_manual(values = color_palette) +  # Same palette for jitter points
  facet_wrap(~ clone_type, nrow = 2, scales = "free_y") +  # Facet by clone_type
  labs(title = "Proportion of BCR Clone expansion Types Across Cell Types",
       x = "Cell Type",
       y = "Proportion (%)") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10, color = "black"),  # Rotate x-axis labels
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),  # Centered title
    strip.text = element_text(size = 12, face = "bold"),  # Facet labels
    panel.grid.major = element_blank(),  # Remove major gridlines for cleaner look
    panel.grid.minor = element_line(color = "grey90", size = 0.25),  # Subtle minor gridlines for reference
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),  # Thin black border around panels
    legend.position = "top",  # Place legend on top
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8, 'lines')  # Make legend key smaller
  ) +
  # Annotate p-values on the plot, positioned in each facet (top-right corner)
  geom_text(data = p_values, aes(label = paste("p =", round(p_value, 5))),
            x = Inf, y = Inf, hjust = 1.1, vjust = 2, size = 2, color = "black", inherit.aes = FALSE)
dev.off()

#separate data between TCR and BCR
all_final_res_tcr <- all_final_res[!all_final_res$clone_class == "IG",]

#only_TCR
data<- all_final_res_tcr
# Step 1: Calculate the frequency of each CDR3 sequence within each neighborhood
cdr3_counts <- data %>%
  group_by(NH_section, aaSeqCDR3) %>%
  summarize(clone_count = n()) %>%
  ungroup()

# Step 2: Classify sequences into unique, small, medium, large, and hyperexpanded categories
cdr3_counts <- cdr3_counts %>%
  mutate(clone_type = case_when(
    clone_count == 1 ~ "Unique",
    clone_count >= 2 & clone_count < 3 ~ "Expanded",
    clone_count >= 3 ~ "Hyperexpanded"
  ))

# Step 3: Calculate the proportion of each clone type within each neighborhood
isotype_usage <- cdr3_counts %>%
  group_by(NH_section, clone_type) %>%
  summarize(clone_count = sum(clone_count)) %>%  # Total clone count per type per neighborhood
  mutate(proportion = clone_count / sum(clone_count)) %>%  # Calculate proportions
  ungroup()

# Step 4: Create a stacked bar plot
# Load required libraries
library(ggplot2)

# Create a more elegant color palette for clone_type
nature_palette <- c("#4F81BD","#ED7D31", "#70AD47", "#A5A5A5", "#FFC000")  # Blue, Green, Orange, Gray, Yellow

# Generate the plot
pdf("TCR_clonal_distribution_across_neighbours_nature.pdf", width = 14, height = 8)

ggplot(isotype_usage, aes(x = NH_section, y = proportion, fill = clone_type)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.3) +  # Stacked bar with black outlines
  scale_fill_manual(values = nature_palette, name = "Clone Type") +  # Custom color palette
  labs(title = "TCR Clonal Distribution Across Neighborhoods", 
       x = "Neighborhood", 
       y = "Proportion (%)") +  # Title and axis labels
  theme_minimal(base_size = 14) +  # Minimal theme for clarity
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 10)),  # Centered bold title with margin
    axis.title.x = element_text(size = 12, face = "bold"),  # Bold x-axis title
    axis.title.y = element_text(size = 12, face = "bold"),  # Bold y-axis title
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10, face = "plain", color = "black"),  # X-axis text
    axis.text.y = element_text(size = 10, color = "black"),  # Y-axis text size
    legend.position = "top",  # Legend on top
    legend.title = element_text(size = 14, face = "bold"),  # Bold legend title
    legend.text = element_text(size = 12),  # Adjust legend text size
    legend.key.size = unit(1, 'cm'),  # Increase legend key size for visibility
    panel.grid.major = element_blank(),  # Remove major gridlines for a cleaner look
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8),  # Black border around the plot
    plot.margin = unit(c(1, 1, 1, 1), "cm")  # Add margins around the plot for better readability
  )

# Close the PDF device
dev.off()

isotype_usage <- isotype_usage %>%
  mutate(
    cell_type = sub("_.*", "", NH_section))
head(isotype_usage)

data<- isotype_usage



# Assuming the dataset is called 'data'
# Step 1: Perform ANOVA for each clone_type across different cell_types
anova_results <- list()
p_values <- data.frame(clone_type = character(), p_value = numeric(), stringsAsFactors = FALSE)

# Loop through each unique clone_type and perform ANOVA
for (clone in unique(data$clone_type)) {
  
  # Subset data for the specific clone_type
  subset_data <- subset(data, clone_type == clone)
  
  # Perform ANOVA: Test if proportion differs across cell_type for each clone_type
  anova_fit <- aov(proportion ~ cell_type, data = subset_data)
  
  # Extract the p-value from ANOVA summary
  anova_pvalue <- summary(anova_fit)[[1]]["Pr(>F)"][1,1]
  
  # Store the p-value in p_values
  p_values <- rbind(p_values, data.frame(clone_type = clone, p_value = anova_pvalue))
  print(p_values)
  
  # Store the ANOVA result for each clone_type in the list
  anova_results[[clone]] <- summary(anova_fit)
}

# Print p-values
print(p_values)

# Step 2: Plot the boxplot for each clone_type, showing proportions across cell_type

# Define a custom color palette
color_palette <- c("#4F81BD", "#C0504D", "#9BBB59", "#8064A2", "#F79646")

# Create the boxplot and annotate with p-values
pdf("Proportion of TCR clone expansion across cell type.pdf", height = 7, width =10)
ggplot(data, aes(x = cell_type, y = proportion, fill = cell_type)) +
  geom_boxplot(outlier.shape = NA, color = "black", size = 0.75, alpha = 0.6) +  # Boxplot with light alpha
  geom_jitter(width = 0.2, size = 1.5, aes(color = cell_type), alpha = 0.8) +  # Subtle jittered points
  scale_fill_manual(values = color_palette) +  # Custom color palette for fill
  scale_color_manual(values = color_palette) +  # Same palette for jitter points
  facet_wrap(~ clone_type, nrow = 2, scales = "free_y") +  # Facet by clone_type
  labs(title = "Proportion of TCR Clone expansion Types Across Cell Types",
       x = "Cell Type",
       y = "Proportion (%)") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10, color = "black"),  # Rotate x-axis labels
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),  # Centered title
    strip.text = element_text(size = 12, face = "bold"),  # Facet labels
    panel.grid.major = element_blank(),  # Remove major gridlines for cleaner look
    panel.grid.minor = element_line(color = "grey90", size = 0.25),  # Subtle minor gridlines for reference
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),  # Thin black border around panels
    legend.position = "top",  # Place legend on top
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8, 'lines')  # Make legend key smaller
  ) +
  # Annotate p-values on the plot, positioned in each facet (top-right corner)
  geom_text(data = p_values, aes(label = paste("p =", round(p_value, 5))),
            x = Inf, y = Inf, hjust = 1.1, vjust = 2, size = 2, color = "black", inherit.aes = FALSE)
dev.off()
