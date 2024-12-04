# Load necessary libraries
library(gdata)
library(dplyr)
library(readr)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(writexl)
#library(plotly)  # For 3D PCA Plot
# Install plotly if not installed
if (!requireNamespace("plotly", quietly = TRUE)) {
  install.packages("plotly", repos = "http://cran.us.r-project.org")
}

# Load plotly
library(plotly)

# Continue with the rest of your code...


# Install missing packages
# Install reshape2 if not installed
if (!requireNamespace("reshape2", quietly = TRUE)) {
  install.packages("reshape2", repos = "http://cran.us.r-project.org")
}

# Load reshape2 package
library(reshape2)

# Perform alternative clustering method: k-means clustering
# clusterGenomics not available, so using k-means as a substitute for PART clustering
if (!requireNamespace("clusterGenomics", quietly = TRUE)) {
  message("clusterGenomics package is not available, using k-means clustering instead.")
}

# Read the counts data from CSV
project_folder = "C:/Users/dhruv/OneDrive/Desktop/work space d/MSC/garden city university material/2nd sem/minor project/DGE candida/"
counts_data <- read.csv(paste0(project_folder, "data/data_for_dge_csv.csv"), header = TRUE, row.names = 1)
cts_integer <- mutate_all(counts_data, function(x) as.integer(x))

# Save the cts_integer data to .csv and .xlsx files
write.csv(cts_integer, "C:/Users/dhruv/OneDrive/Desktop/cts_integer.csv", row.names = TRUE)
write_xlsx(cts_integer, "C:/Users/dhruv/OneDrive/Desktop/cts_integer.xlsx")

# Load the metadata
metadata <- read.table("C:/Users/dhruv/OneDrive/Desktop/work space d/MSC/garden city university material/2nd sem/minor project/DGE candida/data/DGE_FILE_DATA.txt", header = TRUE, row.names = 1, sep = "\t")

# Convert metadata columns to factors
metadata$type <- as.factor(metadata$type)
metadata$timepoint <- as.factor(metadata$timepoint)

# Read the counts file into counts object
counts <- read.csv("C:/Users/dhruv/OneDrive/Desktop/work space d/MSC/garden city university material/2nd sem/minor project/DGE candida/cts_integer.csv", row.names = 1)

# Create a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ type)

# Normalize the data
dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized = TRUE)

# Log transformation
log_norm_counts <- log2(normalized_counts + 1)

# Get top expressed genes
average_expression <- rowMeans(log_norm_counts)
top_genes <- head(order(average_expression, decreasing = TRUE), 20)

# Save differential expression results
res <- results(dds)
resOrdered <- res[order(res$padj),]
top_diff_expressed <- head(resOrdered, 20)
write.csv(as.data.frame(top_diff_expressed), file = "C:/Users/dhruv/OneDrive/Desktop/top_diff_expressed_genes.csv")

# ------------------ K-means Clustering (Alternative to PART Clustering) ------------------
# Perform K-means clustering on the log-transformed normalized counts
kmeans_result <- kmeans(log_norm_counts, centers = 5)  # Choose 5 clusters as an example
clusters <- kmeans_result$cluster  # Extract cluster assignments

# ---------------------- Visualization ----------------------

# Heatmap of Top 20 Highly Expressed Genes
pheatmap(log_norm_counts[top_genes, ],
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         scale = "row",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "Heatmap of Top 20 Highly Expressed Genes")

# Boxplot of Top 20 Genes
par(mar = c(8, 4, 4, 2) + 0.1, las = 2, cex.axis = 0.7)
boxplot(log_norm_counts[top_genes, ], 
        las = 2, 
        col = "lightblue", 
        ylab = "Expression Level",
        xlab = "Genes",  # Label for X-axis
        main = "Top 20 Highly Expressed Genes")

# Barplot of Mean Expression Levels using ggplot2
gene_expression_df <- data.frame(
  Gene = rownames(log_norm_counts)[top_genes],
  Mean_Expression = rowMeans(log_norm_counts[top_genes, ])
)
ggplot(gene_expression_df, aes(x = reorder(Gene, -Mean_Expression), y = Mean_Expression, fill = Gene)) +
  geom_bar(stat = "identity", color = "black", show.legend = FALSE) +
  geom_text(aes(label = round(Mean_Expression, 2)), vjust = -0.5, size = 3) +
  labs(x = "Gene", y = "Mean Expression Level", title = "Top 20 Highly Expressed Genes") +
  theme_minimal() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set3")

# Volcano Plot
data <- read.csv("C:/Users/dhruv/OneDrive/Desktop/top_diff_expressed_genes.csv")
data$log10padj <- -log10(data$padj)
data$expression <- ifelse(data$log2FoldChange > 1, "Highly Expressed", "Others")
volcano_plot <- ggplot(data, aes(x=log2FoldChange, y=log10padj, color=expression)) +
  geom_point(alpha=0.8) + 
  scale_color_manual(values = c("Highly Expressed" = "red", "Others" = "blue")) +
  theme_minimal() + 
  labs(title = "Volcano Plot of Differential Gene Expression", 
       x = "Log2 Fold Change", 
       y = "-Log10 Adjusted P-value") +
  theme(plot.title = element_text(hjust = 0.5))
print(volcano_plot)

#linear plot to visualize the gene expression
# install.packages("ggplot2")

# Load necessary library
library(ggplot2)

# Load the data
data <- read.csv("C:/Users/dhruv/OneDrive/Desktop/work space d/MSC/garden city university material/2nd sem/minor project/DGE candida/results/top_diff_expressed_genes.csv", 
                 header = TRUE, encoding = "ISO-8859-1")

# Create the plot
plot <- ggplot(data, aes(x = Gene, y = baseMean, group = 1)) +  # group = 1 ensures all points are treated as part of the same group
  geom_line(color = "blue", size = 1) +  # Adds the line graph
  geom_point(color = "red", size = 3) +  # Adds points on the graph
  theme_minimal() +                      # Clean theme
  labs(title = "Linear Plot for Top Expressed Genes",  # Add title
       x = "Gene",                        # Label for x-axis
       y = "Expression Level (baseMean)") + # Label for y-axis
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) # Rotate x-axis labels

# View the plot
print(plot)




# ---------------------- MDS Plot ----------------------
# Calculate the distance between samples based on normalized counts
dist_matrix <- dist(t(log_norm_counts))
mds <- cmdscale(dist_matrix)

# Plot the MDS plot
plot(mds[,1], mds[,2], 
     xlab = "MDS1", 
     ylab = "MDS2", 
     main = "MDS Plot", 
     pch = 19, col = as.factor(metadata$type))

# Add sample labels to the plot
text(mds[,1], mds[,2], labels = rownames(log_norm_counts), pos = 3, cex = 0.6)

# Add a legend to the plot
legend("bottomleft", legend = unique(metadata$type), col = 1:length(unique(metadata$type)), pch = 19)

# ---------------------- PCA Plot ----------------------
# ---------------------- PCA Plot ----------------------
# Load the necessary libraries
library(plotly)

# Load the CSV file containing top differentially expressed genes
top_genes_data <- read.csv("C:/Users/dhruv/OneDrive/Desktop/work space d/MSC/garden city university material/2nd sem/minor project/DGE candida/results/top_diff_expressed_genes.csv", 
                           encoding = 'ISO-8859-1')

# Assuming the CSV file has a column called 'Gene_name' that lists the gene names of interest
# Extract the top 20 highly expressed genes
top_20_genes <- top_genes_data$Gene_name[1:20]  # Adjust if necessary

# Assuming the PCA data is already calculated and stored in 'pca' and 'clusters'
# Adjust the clusters to match the number of PCA samples (if necessary)
clusters <- clusters[1:18]

# Create the PCA data frame with sample names
pca_data_3d <- data.frame(Sample = rownames(pca$x), 
                          PC1 = pca$x[, 1], 
                          PC2 = pca$x[, 2], 
                          PC3 = pca$x[, 3],  # Third component for 3D plot
                          Cluster = as.factor(clusters))

# Create a column in pca_data_3d for labels of samples that are among the top 20 genes
pca_data_3d$Label <- ifelse(pca_data_3d$Sample %in% top_genes_data$Gene_name[1:20], 
                            pca_data_3d$Sample, "")

# Create the 3D PCA plot with sample names for top genes
pca_3d_plot <- plot_ly(data = pca_data_3d, 
                       x = ~PC1, 
                       y = ~PC2, 
                       z = ~PC3, 
                       type = 'scatter3d', 
                       mode = 'text+markers',  # Ensure both markers and text are shown
                       text = ~Sample,  # Use Sample names for labels
                       textposition = 'top right',  # Adjust label position
                       marker = list(size = 5, 
                                     color = ~Cluster, 
                                     colorscale = 'Viridis', 
                                     showscale = TRUE),
                       hoverinfo = 'text',  # Ensure hovering shows text
                       textfont = list(size = 5, color = 'black')) %>%  # Set text properties
  layout(title = "3D PCA Plot of K-means Clustering",
         scene = list(
           xaxis = list(title = "PC1"),
           yaxis = list(title = "PC2"),
           zaxis = list(title = "PC3")),
         legend = list(title = list(text = 'Cluster')))

# Display the plot
pca_3d_plot

