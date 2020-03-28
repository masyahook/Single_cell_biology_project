# Loading the libraries for the analysis
library(Seurat)
library(dplyr)
library(ggplot2)

# Reading the file
raw_counts <- read.table(file='/Users/MaximKryukov/Documents/GitHub/Single_cell_biology_project/data/MOUSE_BRAIN_DATASET_2_COUNTS.tsv')
mouse_rawdata <- CreateSeuratObject(counts = raw_counts, project = 'Mouse_cell_atlas')
sprintf('The total number of genes: %d. The total number of cells/samples: %d', dim(mouse_rawdata)[1], dim(mouse_rawdata)[2])

### DATA MANAGEMENT

# Plotting number of UMI counts for each cell
min_UMI_count = 600
UMI_counts = sort(mouse_rawdata$nCount_RNA, decreasing = TRUE)
cell_nums = 1:length(UMI_counts)
plot(cell_nums, UMI_counts, main = "Cells",
     xlab = "Barcodes (cell num)", ylab = "UMI count", col = 'green', log='y')
points(cell_nums[UMI_counts<=min_UMI_count], UMI_counts[UMI_counts<=min_UMI_count], main = "Cells",
       xlab = "Barcodes (cell num)", ylab = "UMI count", col = 'black')
legend("right", legend=c(sprintf("UMI counts > %d", min_UMI_count), sprintf("UMI counts < %d", min_UMI_count)),
       col=c("green", "black"), pch=c(1, 1), bty='n')

mouse_rawdata_filtered <- subset(mouse_rawdata, subset = nCount_RNA > min_UMI_count)

# Getting some statistics about our dataset
sprintf('The median number of detected genes per cell: %f. The mean number of detected genes per cell: %f',
        median(mouse_rawdata_filtered$nFeature_RNA), mean(mouse_rawdata_filtered$nFeature_RNA))
sprintf('The median number of UMI counts per cell: %f. The mean number of UMI counts per cell: %f', 
        median(mouse_rawdata_filtered$nCount_RNA), mean(mouse_rawdata_filtered$nCount_RNA))

### FILTERING

### Now let's plot some statistic about our dataset

# Creating the metadata parameter for percentage of mitochnondrial genes
mouse_rawdata_filtered[["percent.mt"]] <- PercentageFeatureSet(mouse_rawdata_filtered, pattern = "^mt-")

# Plotting the distribution of detected genes, UMI counts, and percentage of mt genes
VlnPlot(mouse_rawdata_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE)

# Now let's plot the dependance between metadata features
plot1 <- FeatureScatter(mouse_rawdata_filtered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mouse_rawdata_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# Now we can take a subset of the initial dataset that we will further work with
mt_percent_max = 5
nFeature_count_min = 300
nFeature_count_max = 4000

mouse_data <- subset(mouse_rawdata_filtered, subset = (nFeature_RNA > nFeature_count_min) & (nFeature_RNA < nFeature_count_max) & (percent.mt < mt_percent_max))
sprintf('The total number of genes: %d. The total number of cells/samples: %d', dim(mouse_data)[1], dim(mouse_data)[2])

### NORMALIZATION

# We normalize data with by correcting the sequencing depth and log transformatting
mouse_data <- NormalizeData(mouse_data, normalization.method = "LogNormalize", scale.factor = 10000)

# We search for variable features
mouse_data <- FindVariableFeatures(mouse_data, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10_variable_genes <- head(VariableFeatures(mouse_data), 10)

# Plot variable features with and without labels
variable_feature_plot <- VariableFeaturePlot(mouse_data)
LabelPoints(plot = variable_feature_plot, points = top10_variable_genes, repel = TRUE)

# Do the standardization (z-scoring) of the data
all_genes <- rownames(mouse_data)
mouse_data <- ScaleData(mouse_data, features = all_genes)

# DIMENSIONALITY REDUCTION

# Now we can find the principal components to find subtypes of cells
mouse_data <- RunPCA(mouse_data, features = VariableFeatures(object = mouse_data))

# Plot the explained variance of the PCs
ElbowPlot(mouse_data)

# Get the genes that contribute the most to the first 5 PCs
print(mouse_data[["pca"]], dims = 1:5, nfeatures = 5)

# Visualize the top genes associated with first 2 PCs
VizDimLoadings(mouse_data, dims = 1:2, reduction = "pca")

# Visualize the PCA with the first 2 PCs
DimPlot(mouse_data, reduction = "pca")

### Additional evaluation of the PCs
mouse_data <- JackStraw(mouse_data, num.replicate = 100)
mouse_data <- ScoreJackStraw(mouse_data, dims = 1:20)
JackStrawPlot(mouse_data, dims = 1:20)

### CLUSTERING

# Now we can cluster cells using Phenograph 
mouse_data <- FindNeighbors(mouse_data, dims = 1:12)
mouse_data <- FindClusters(mouse_data, resolution = 0.5)

# Plotting the clustered cells with UMAP
mouse_data <- RunUMAP(mouse_data, dims = 1:12)
DimPlot(mouse_data, reduction = "umap")

### TECHNICAL AND BIOLOGICAL BIAS DETECTION

# Detection for 3 different technical biases:
# 1) the sequencing depth (number of counts for each cell)
# 2) the cell cycle bias
# 3) the drouplet effect
# 4) the mitochondrial genes (dying cells)

# Plotting UMAP with 3 feautre coloring
FeaturePlot(mouse_data, features = c('nFeature_RNA', 'nCount_RNA', "percent.mt"), ncol=3)

# Getting cell-cycle canonical gene markers
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.gene

# Getting a subset that is present in mouse_data from initial list of gene markers
s_genes_present <- intersect(s_genes, all_genes)
g2m_genes_present <- intersect(g2m_genes, all_genes)

# The arrays above are empty -> no cell-cycle effect!

### DIFFERENTIAL EXPRESSION ANALYSIS

# Find markers for every cluster compared to all remaining cells, report only the positive ones
mouse_data.markers <- FindAllMarkers(mouse_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use='t')
de_genes <- data.frame(mouse_data.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))

# Heatmap expression plot for variable genes for each cluster
top10 <- mouse_data.markers  %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(mouse_data, features = top10$gene) + NoLegend()

# Let us try to do the same but with Wilcoxon test
mouse_data.markers_wilcoxon <- FindAllMarkers(mouse_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
de_genes_wilcoxon <- data.frame(mouse_data.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))

