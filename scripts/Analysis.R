# Packages 
library(Seurat)
#BiocManager::install("SingleR")
library(SingleR)
#BiocManager::install("celldex")
library(celldex)
#BiocManager::install("scrapper")
library(scrapper)
#BiocManager::install("dittoSeq")
library(dittoSeq)
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
#BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
# To make analyses faster 
#devtools::install_github('immunogenomics/presto')
library(presto)
library(DESeq2)
#BiocManager::install("apeglm")
library(apeglm)
library(dplyr)


# Load in Seurat object
data <- readRDS("../data/seurat_ass4.rds")

# inspect
data@meta.data

# Calculate mtDNA percentage
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")

Idents(data) <- "time"

# Plot QC metrics, looks a little crazy
before_filter <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1, group.by = "time", raster = F)

# Filter the data as per Kazer et al. (2024) for consistency of results
data <- subset(data, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & percent.mt < 15 & nCount_RNA > 750 & nCount_RNA < 100000)

#QC Metric plots after filtering

after_filter <- VlnPlot(data, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0.1, ncol = 3, group.by = "time")

before_filter / after_filter
#ggsave("../figures/elbow.png", plot = after_filter)

# Plot the Cell Counts per sample
cell.counts.trimmed = data.frame(table(data$orig.ident))
ggplot(cell.counts.trimmed, aes(x=Var1, y=Freq)) + geom_bar(aes(fill = Var1), stat = "identity") + labs(x = "", y = "Cell Count") + guides(fill = 'none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))


plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Normalisation with log normalisation using VST method, choose 5000 features
data_norm <- NormalizeData(data, normalization.method = "LogNormalize")
data_norm <- FindVariableFeatures(data_norm, selection.method = "vst", nfeatures = 5000) 

# Scaling the data

all_genes <- rownames(data_norm)
data_scaled <- ScaleData(data_norm, features = all_genes)

# Linear dimensional reduction
  
data_scaled_PC <- RunPCA(data_scaled, features = VariableFeatures(object = data_scaled))

elbow <- ElbowPlot(data_scaled_PC) # looks like 17 PCs
ggsave("../figures/elbow.png", plot = elbow, bg = "white")

# Visualise first 4 PCs
VizDimLoadings(data_scaled_PC, dims = 1:4, reduction = "pca")

# Plot PCA and heatmap
DimPlot(data_scaled_PC, reduction = "pca") + NoLegend()
DimHeatmap(data_scaled_PC, dims = 1:10, cells = 500, balanced = TRUE)

# Clustering using the first 17 PCs and 0.5 resolution
data_scaled_PC <- FindNeighbors(data_scaled_PC, dims = 1:17)

data_scaled_PC <- FindClusters(data_scaled_PC, resolution = 0.5)

# Check cluster IDs
head(Idents(data_scaled_PC), 5)

# UMAP on the data
data_UMAP <- RunUMAP(data_scaled_PC, dims = 1:17)

# Plot clusters
umap <- DimPlot(data_UMAP, reduction = "umap", label = T)
ggsave("../figures/umap_unlabelled.png", plot = umap, width = 10, height = 8, dpi = 300, bg = "white")

# Check number of clusters
length(unique(Idents(data_UMAP)))

# Annotating the clusters with singleR
# Loading in reference data
ref <- MouseRNAseqData()

# Get RNA counts from data
data_counts <- GetAssayData(data_UMAP, layer = "data")

# use singler

# Annotating the data with the fine scale data labels for more clarity
annotated_data <- SingleR(test = data_counts, ref = ref, clusters = Idents(data_UMAP), labels = ref$label.fine)

# Checking it worked
annotated_data
rownames(annotated_data)
levels(Idents(data_UMAP))


# Correct the labels- use pruned labels instead
cluster_labs <- annotated_data$pruned.labels
names(cluster_labs) <- rownames(annotated_data)
cell_clusters <- as.character(Idents(data_UMAP))
cell_labels <- cluster_labs[cell_clusters]
names(cell_labels) <- colnames(data_UMAP)

# Now add to metadata
data_UMAP$singleR_clusterlabels <- cell_labels

# Plot the clusters
UMAP_clusters <- DimPlot(data_UMAP, reduction = 'umap', group.by = 'singleR_clusterlabels', label = T) + ggtitle("UMAP Cell Type Clusters")

ggsave("../figures/clusters_UMAP.png", plot = UMAP_clusters, dpi = 300)

# have to specify idents
Idents(data_UMAP) <- "singleR_clusterlabels"

markers <- FindAllMarkers(data_UMAP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
length(unique(markers$cluster))
names(markers)

# identify what cluster 0 is annotated as 
cluster_zero <- as.data.frame(data_UMAP@meta.data)
cluster_zero %>%
  filter(seurat_clusters == "0") # Neurons

# Export top 20 markers for neurons
markers_for_neurons <- markers %>%
  filter(cluster == "Neurons") %>%
  arrange(desc(abs(avg_log2FC))) %>%
  slice_head(n = 20)
write.csv(markers_for_neurons, "../outputs/top20_markers_neurons.csv")
  

# Check prediction quality by looking at the scores
annotated_data$scores
plotScoreHeatmap(annotated_data)

# Clearly each cluster contains many different types of cells, so plot stacked bar plot
# Plotting cluster distribution by sample
tissues <- dittoBarPlot(data_UMAP, var = "singleR_clusterlabels", group.by = "organ_custom", main = "Cell Type Clusters in each Tissue Type", xlab = "Tissue")
ggsave("../figures/clusters_tissues.png", plot = tissues, width = 10, height = 8, dpi = 300)

# And by time
cluster_barplot <- dittoBarPlot(data_UMAP, var = "singleR_clusterlabels", group.by = "time", main = "Cell Type Clusters at Each Time after Infection", xlab = "Time after infection")

ggsave("../figures/clusters_timepoints.png", plot = cluster_barplot, width = 10, height = 8, dpi = 300)

fig3 <- cluster_barplot + tissues

ggsave("../figures/fig3.png", plot = fig3, width = 10, height = 8, dpi = 300)

# Neurons = top cluster- I will focus on this for further analysis

# Differential expression analysis

# Subset neurons
neurons <- subset(data_UMAP, subset = singleR_clusterlabels == "Neurons")

# need to put time in 
table(data_UMAP@meta.data$time)

#D02, D05, D08, D14, Naive

#add in 
Idents(neurons) <- "time"

# Pseudobulk samples 
pseudo_neurons <- AggregateExpression(neurons, assays = "RNA", return.seurat = T, group.by = c("time", "mouse_id", "organ_custom"))

# Differential expression with DESeq2
# Get counts data
neurons_counts <- GetAssayData(pseudo_neurons, assay = "RNA", layer = "counts")

# Extract metadata 
meta_data <- pseudo_neurons@meta.data

# Set the reference levels so that everything is compared to naive or RM, make sure everything is factors
meta_data$time <- factor(meta_data$time, levels = c("Naive", "D02", "D05", "D08", "D14"))
meta_data$organ_custom <- factor(meta_data$organ_custom, levels = c("RM", "OM", "LNG"))

#?DESeqDataSetFromMatrix

# include organs?
#random batch effects?

neurons@meta.data

# DE analysis 
de_obj_time <- DESeqDataSetFromMatrix(neurons_counts, colData = meta_data, design = ~ time)
de_obj_time <- DESeq(de_obj_time)
#de_obj_time <- results(de_obj_time)
 
# Check dispersion
plotDispEsts(de_obj_time)

# Shrink the coefs of each contrast
D02 <- lfcShrink(de_obj_time, coef = "time_D02_vs_Naive", type = "apeglm")
D02 <- as.data.frame(D02[!is.na(D02$padj), ]) # get rid of NAs

D05 <- lfcShrink(de_obj_time, coef = "time_D05_vs_Naive", type = "apeglm")
D05 <- as.data.frame(D05[!is.na(D05$padj), ])

D08 <- lfcShrink(de_obj_time, coef = "time_D08_vs_Naive", type = "apeglm")
D08 <- as.data.frame(D08[!is.na(D08$padj), ])

D14 <- lfcShrink(de_obj_time, coef = "time_D14_vs_Naive", type = "apeglm")
D14 <- as.data.frame(D14[!is.na(D14$padj), ])

# Write out, filtering for significant only and top 20
D14_signif <- D14 %>%
  filter(padj < 0.05) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = 20)
write.csv(D14_signif, "../outputs/Naive_vs_D14_significant.csv")
  
D02_signif <- D02 %>%
  filter(padj < 0.05) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = 20)
write.csv(D02_signif, "../outputs/Naive_vs_D02_significant.csv")

D05_signif <- D05 %>%
  filter(padj < 0.05) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = 20)
write.csv(D05_signif, "../outputs/Naive_vs_D05_significant.csv")

D08_signif <- D08 %>%
  filter(padj < 0.05) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = 20)
write.csv(D08_signif, "../outputs/Naive_vs_D08_significant.csv")

# Plot
# Volcano plot #-log10(padj)
D14_plot <- D14 %>%
  mutate(significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1,
                              "Significant", "Not significant"))

volcano_plot <- ggplot(D14_plot) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = significant)) +
  ggtitle("Volcano Plot of Neurons at Day 14 of Infection vs. Day 0") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_colour_manual(values = c("Significant" = "red", 
                                 "Not significant" = "blue"), name = "Significance (padj < 0.05)") + # red is signif
  scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "right", legend.title = element_text(size = 8),
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))
# what's up or down regulated at D14 compared to naive

ggsave("../figures/DE_d14_vs_naive_volcanoplot.png", plot = volcano_plot, width = 10, height = 8, dpi = 300)

# Gene Set Enrichment Analysis
# Plot
VlnPlot(neurons, features="nFeature_RNA")

# Get gene names and rank by LFC
# Rank the list
genes_D14_vs_naive <- D14$log2FoldChange
names(genes_D14_vs_naive) <- row.names(D14)
genes_D14_vs_naive <- sort(genes_D14_vs_naive, decreasing = T)

gsea_d14_vs_naive <- gseGO(genes_D14_vs_naive, ont="BP", keyType = "SYMBOL", OrgDb=org.Mm.eg.db)
dotplot <- dotplot(gsea_d14_vs_naive, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave("../figures/gsea_d14_vs_naive_dotplot.png", plot = dotplot, width = 10, height = 8, dpi = 300)

signif_results <- as.data.frame(gsea_d14_vs_naive[gsea_d14_vs_naive$p.adjust < 0.05, ])
dim(signif_results) # 527 categories are enriched

# Export
write.csv(signif_results, "../outputs/signif_results_GSEA_d14_vs_naive.csv")









