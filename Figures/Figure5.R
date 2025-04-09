#-------------------------------------------------------------------------------
# MULTIOME PAPER 1 - REVISION 1 - FIGURE 5
#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
# Loading Required Libraries
#-------------------------------------------------------------------------------

# Setting path for libraries
# .libPaths("/opt/R/4.1.1/lib/R/library")
.libPaths("/home/priya/R/x86_64-pc-linux-gnu-library/4.1/")

library(Signac) # for ATAC-seq analysis
library(Seurat) # for multimodal analysis
# library(EnsDb.Hsapiens.v86) # for human genome annotation
library(BSgenome.Hsapiens.UCSC.hg38) # human genome
library(SeuratDisk) # for storing metadata related to Seurat assay
library(stringr) # for string processing
library(BiocParallel) # for parallel processing
library(scDblFinder) # for doublet detection
library(pheatmap) # for mapping clusters to cell annotations
library(dplyr) # for stacked barplot
library(tidyverse) # for frequency of sample/cluster
library(MAST) # for DE analysis
library(DESeq2) # for DE analysis
library(ggplot2) # for bar chart
library(ggforce) # for function geom_sina()
library(plotly) # for plotting
library(TFBSTools) # for motif analysis
library(JASPAR2020) # for motif analysis
library(chromVAR) # for motif analysis
library(rhdf5) # to write .h5 file
library(ggpubr) # for box plot with p-value
library(data.table) # for function tstrsplit()
library(RColorBrewer) # for color palette
library(ComplexUpset) # for UpSet plot
library(bedr) # for peak distribution 
# library(ChIPseeker) # for peak distribution 
# library(TxDb.Hsapiens.UCSC.hg19.knownGene) # for peak distribution 
library(clusterProfiler) # for peak distribution 
library(annotables) # for peak distribution 
library(org.Hs.eg.db) # for peak distribution 
library(Matrix) # for function writemm()
library(AnnotationDbi) # for converting gene symbols to ENSEMBL IDs
library(rhdf5) # to write .h5 file
library(scater)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(magrittr)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(apeglm)
library(png)
library(openxlsx)
library(bedr)





#-------------------------------------------------------------------------------
# SCT1 Based UMAP - Grouped By Donor Status
#-------------------------------------------------------------------------------

# Loading data
seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothCD4TCMTregClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")


# Defining variables
active_ident = "bothcd4tcmtreg_SCT1_snn_res.0.1"


# Removing cluster 5
Idents(seurat_obj) = active_ident
seurat_obj_subset = subset(seurat_obj, idents = 5, invert = TRUE)


# Defining color palette
col_pal = c("T1DM" = "#d7191c",
            "AAB+" = "#fdae61",
            "HD" = "#2c7bb6")

# UMAP
DimPlot(seurat_obj_subset, repel = TRUE, reduction = "bothcd4tcmtreg.rna.umap",
        group.by = "status", cols = col_pal) +
  ggtitle("UMAP RNA - Status - \n CD4 TCM + Treg - PLNo + PLNn + Spleen") +
  xlab("UMAP1") + ylab("UMAP2")





#-------------------------------------------------------------------------------
# SCT1 Based UMAP - Grouped By SCT1 Res 0.1
#-------------------------------------------------------------------------------

# Loading data
seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothCD4TCMTregClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")


# Defining variables
active_ident = "bothcd4tcmtreg_SCT1_snn_res.0.1"


# Removing cluster 5
Idents(seurat_obj) = active_ident
seurat_obj_subset = subset(seurat_obj, idents = 5, invert = TRUE)


# Defining colors for each cluster
col_pal = c("0" = "#78c679",
            "1" = "#31a354",
            "2" = "#f768a1",
            "3" = "#c51b8a",
            "4" = "#feb24c",
            "5" = "#ff7f50",
            "6" = "#2b8cbe")

# UMAP
DimPlot(seurat_obj_subset, repel = TRUE, reduction = "bothcd4tcmtreg.rna.umap",
        group.by = active_ident, cols = col_pal) +
  ggtitle("UMAP RNA - Cluster SCT1 Res 0.1 - \n CD4 TCM + Treg - PLNo + PLNn + Spleen") +
  xlab("UMAP1") + ylab("UMAP2")





#-------------------------------------------------------------------------------
# Stacked Bar Plot - Donor Composition Per Cluster
#-------------------------------------------------------------------------------

# Loading data
seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothCD4TCMTregClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")



# Defining variables
active_ident = "bothcd4tcmtreg_SCT1_snn_res.0.1"



# Removing cluster 5
Idents(seurat_obj) = active_ident
seurat_obj_subset = subset(seurat_obj, idents = 5, invert = TRUE)



# Stacked barplot depicting donor status composition for each cluster
tmp1 = as.data.frame(table(seurat_obj_subset@meta.data[[active_ident]], seurat_obj_subset$status))
colnames(tmp1) = c("Cluster", "Status", "Freq")
tmp1 = tmp1[-which(tmp1$Cluster == 5), ]

# Converting cell type distribution per sample into percentage format
tmp2 = as.data.frame(tmp1 %>% pivot_wider(names_from = Status, values_from = Freq))
rownames(tmp2) = tmp2[, 1]
tmp2 = tmp2[, -1]

tmp3 = (tmp2/rowSums(tmp2)) * 100

tmp4 = cbind.data.frame(rownames(tmp3), tmp3)
colnames(tmp4) = c("Cluster", colnames(tmp3))

tmp5 = reshape2::melt(tmp4, id = "Cluster")
colnames(tmp5) = c("Cluster", "Status", "Percent")

tmp5$Status = factor(tmp5$Status, levels = c("T1DM", "AAB+", "HD"))


# Defining color plaette
col_pal = c("T1DM" = "#d7191c", "AAB+" = "#fdae61", "HD" = "#2c7bb6")

# Stacked bar plot depicting donor composition per cluster
stacked_bp = ggplot(data = tmp5, aes(x = Cluster, y = Percent, fill = Status)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlab("Cluster") +
  ylab("Percent") +
  scale_fill_manual(values = col_pal)





#-------------------------------------------------------------------------------
# Donut Plot Depicting Donor Composition Per Cluster
#-------------------------------------------------------------------------------

# Loading data
seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothCD4TCMTregClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")



# Defining variables
active_ident = "bothcd4tcmtreg_SCT1_snn_res.0.1"



# Subsetting Seurat object
Idents(seurat_obj) = active_ident
seurat_obj_subset = subset(seurat_obj, idents = 1)


# Donut plot
donutplot_df = data.frame(table(seurat_obj_subset$donor_ids_umap))
colnames(donutplot_df) = c("donor.id", "freq")

# Computing percentages
donutplot_df$fraction = donutplot_df$freq / sum(donutplot_df$freq)

# Computing cumulative percentages (top of each rectangle)
donutplot_df$ymax = cumsum(donutplot_df$fraction)

# Computing bottom of each rectangle
donutplot_df$ymin = c(0, head(donutplot_df$ymax, n = -1))

# Compute label position
donutplot_df$label.position = (donutplot_df$ymax + donutplot_df$ymin) / 2

# Computing good label
donutplot_df$label = paste0(donutplot_df$donor.id, "\n value: ", donutplot_df$freq)

# Defining color palette
col_pal = c("T1D1" = "#d94801",
            "T1D2" = "#ef3b2c",
            "T1D4" = "#cb181d",
            "T1D6" = "#e31a1c",
            "T1D10" = "#bd0026",
            "T1D12" = "#a50f15",
            "T1D15" = "#800026",
            "T1D17" = "#67000d",
            "AAB+1" = "#ffffcc",
            "AAB+2" = "#ffeda0",
            "AAB+3" = "#fdd0a2",
            "AAB+4" = "#fed976",
            "AAB+5" = "#fdae6b",
            "AAB+6" = "#feb24c",
            "AAB+7" = "#fd8d3c",
            "AAB+8" = "#fd8d3c",
            "HD1" = "#9ecae1",
            "HD3" = "#6baed6",
            "HD6" = "#4eb3d3",
            "HD10" = "#4292c6",
            "HD11" = "#2b8cbe",
            "HD13" = "#2171b5",
            "HD14" = "#08519c",
            "HD16" = "#08306b")

# Plotting
donut_plot = ggplot(donutplot_df, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = donor.id), col) +
  geom_rect() +
  geom_label(x = 3.5, aes(y = label.position, label = label), size = 2) +
  # geom_text(x = 2, aes(y = label.position, label = label, color = donor.id), size = 1) + # x here controls label position (inner/outer)
  scale_fill_manual(values = col_pal, name = "Donor IDs") +
  coord_polar(theta = "y") +
  xlim(c(2, 4)) +
  # xlim(c(-1, 4)) +
  theme_void()
# Ordering donor ids
donut_plot$data$donor.id = factor(x = donut_plot$data$donor.id,
                                  levels = c("T1D1", "T1D2", "T1D4", "T1D6", "T1D10", "T1D12", "T1D15", "T1D17",
                                             "AAB+1", "AAB+2", "AAB+3", "AAB+4", "AAB+5", "AAB+6", "AAB+7", "AAB+8",
                                             "HD1", "HD3", "HD6", "HD10", "HD11", "HD13", "HD14", "HD16"))
donut_plot





#-------------------------------------------------------------------------------
# Violin Plot - Gene Expression Per Cluster
#-------------------------------------------------------------------------------

# Loading data
seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothCD4TCMTregClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")


# Defining variables
active_ident = "bothcd4tcmtreg_SCT1_snn_res.0.1"


# Removing cluster 5
Idents(seurat_obj) = active_ident
seurat_obj_subset = subset(seurat_obj, idents = 5, invert = TRUE)


# Defining color palette
col_pal = c("0" = "#78c679",
            "1" = "#31a354",
            "2" = "#f768a1",
            "3" = "#c51b8a",
            "4" = "#feb24c",
            "6" = "#2b8cbe")

# Violin plot depicting gene expression per cluster
VlnPlot(seurat_obj_subset, features = "NFKB1", assay = "SCT1",
        slot = "data", group.by = active_ident, cols = col_pal)

VlnPlot(seurat_obj_subset, features = "BACH2", assay = "SCT1",
        slot = "data", group.by = active_ident, cols = col_pal)

VlnPlot(seurat_obj_subset, features = "TCF7", assay = "SCT1",
        slot = "data", group.by = active_ident, cols = col_pal)

VlnPlot(seurat_obj_subset, features = "IRF1", assay = "SCT1",
        slot = "data", group.by = active_ident, cols = col_pal)

VlnPlot(seurat_obj_subset, features = "FOXP3", assay = "SCT1",
        slot = "data", group.by = active_ident, cols = col_pal)





#-------------------------------------------------------------------------------
# Dot Plot - Gene Expression Per Cluster
#-------------------------------------------------------------------------------

# Loading data
seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothCD4TCMTregClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")


# Defining variables
active_ident = "bothcd4tcmtreg_SCT1_snn_res.0.1"


# Removing cluster 5
Idents(seurat_obj) = active_ident
seurat_obj_subset = subset(seurat_obj, idents = 5, invert = TRUE)


# Dot plot depicting gene expression per cluster
# features = c("NFKB1", "BACH2", "TCF7", "IRF1", "FOXP3")
features = c("FOXP3")

DotPlot(seurat_obj_subset, 
        features = features, 
        assay = "SCT1", 
        group.by = active_ident)





#-------------------------------------------------------------------------------
# SCT1 UMAP - Grouped by Tissue
#-------------------------------------------------------------------------------

# Loading data
seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothCD4TCMTregClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")


# Defining variables
active_ident = "bothcd4tcmtreg_SCT1_snn_res.0.1"


# Removing cluster 5
Idents(seurat_obj) = active_ident
seurat_obj_subset = subset(seurat_obj, idents = 5, invert = TRUE)


# Defining color palette
# col_pal = c("PLNo" = "#7b3294",
#             "PLNn" = "#f20adb",
#             "S" = "#008837")
col_pal = c("PLNo" = "#7b3294",
            "PLNn" = "#7b3294",
            "S" = "#008837")

# UMAP
DimPlot(seurat_obj_subset, repel = TRUE, reduction = "bothcd4tcmtreg.rna.umap",
        group.by = "tissue", cols = col_pal) +
  ggtitle("UMAP RNA - Tissue - \n CD4 TCM + Treg - PLNo + PLNn + Spleen") +
  xlab("UMAP1") + ylab("UMAP2")





#-------------------------------------------------------------------------------
# Feature Plot of Cytokines
#-------------------------------------------------------------------------------

# Loading data
seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothCD4TCMTregClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")



# Defining variables
fig_dir = "/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Figures/BothCD4TCMTreg_PLNoPLNnSpleen_MedianAge_SubsetPLN8AAB.8T1D.8HD/"
active_ident = "bothcd4tcmtreg_SCT1_snn_res.0.1"
assay = "SCT1"
max_pc = 27
res = 0.1
cell_type = c("CD4 TCM", "Treg")
tissue = "PLNo_PLNn_Spleen_MedianAge_Subset24Donors"



# Subsetting Seurat object to remove cluster 5
Idents(seurat_obj) = active_ident
seurat_obj_subset = subset(seurat_obj, idents = 5, invert = TRUE)



# Feature plot and violin plot of cytokines
DefaultAssay(seurat_obj_subset) = assay

features = c("IFNG", "TBX21", "GATA3", "IL4", "IL5", "IL15", "BCL6", "IL21", 
             "ICOS", "IL2", "CXCR5", "RORC")

# Defining colors for each cluster
col_pal = c("0" = "#78c679",
            "1" = "#31a354",
            "2" = "#f768a1",
            "3" = "#c51b8a",
            "4" = "#feb24c",
            "6" = "#2b8cbe")


for(i in 1:length(features))
{

  pdf(paste0(fig_dir, "FeaturePlot_UMAPNormSCT_PC", max_pc, "_Algo3Res", res, "_", gsub(" ", "", features[i]), "_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".pdf"),
      width = 6,
      height = 4)
  fp = FeaturePlot(seurat_obj_subset, reduction = "bothcd4tcmtreg.rna.umap", features = features[i])
  print(fp)
  dev.off()

  pdf(paste0(fig_dir, "ViolinPlot_NormSCT_PC", max_pc, "_Algo3Res", res, "_", gsub(" ", "", features[i]), "_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".pdf"),
      width = 8,
      height = 4)
  vp = VlnPlot(seurat_obj_subset, group.by = active_ident, cols = col_pal, features = features[i])
  print(vp)
  dev.off()

}





# #-----------------------------------------------------------------------------------------------
# # Cluster 1 (Treg) - Identifying Differential Features For T1D v/s HD, T1D v/s AAB+, AAB+ v/s HD
# #-----------------------------------------------------------------------------------------------
# 
# # Loading data
# seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothCD4TCMTregClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")
# 
# 
# 
# # Defining variables
# active_ident = "bothcd4tcmtreg_SCT1_snn_res.0.1"
# assay = "SCT1"
# 
# 
# 
# # Subsetting cluster 1
# Idents(seurat_obj) = active_ident
# seurat_obj_subset = subset(seurat_obj, idents = 1)
# 
# 
# DefaultAssay(seurat_obj_subset) = assay
# Idents(seurat_obj_subset) = active_ident
# 
# 
# # Declaring list to store  differential features
# markers_per_clust = list()
# 
# # Identifying differential features between T1D and HD
# markers = FindMarkers(seurat_obj_subset,
#                       ident.1 = "T1DM",
#                       ident.2 = "HD",
#                       slot = "data",
#                       logfc.threshold = 0,
#                       min.pct = 0)
# markers$features = rownames(markers)
# markers_per_clust[["clust1_T1DvsHD"]] = markers
# 
# # Identifying differential features between T1D and AAB+
# markers = FindMarkers(seurat_obj_subset,
#                       ident.1 = "T1DM",
#                       ident.2 = "AAB+",
#                       slot = "data",
#                       logfc.threshold = 0,
#                       min.pct = 0)
# markers$features = rownames(markers)
# markers_per_clust[["clust1_T1DvsAAB"]] = markers
# 
# # Identifying differential features between AAB+ and HD
# markers = FindMarkers(seurat_obj_subset,
#                       ident.1 = "AAB+",
#                       ident.2 = "HD",
#                       slot = "data",
#                       logfc.threshold = 0,
#                       min.pct = 0)
# markers$features = rownames(markers)
# markers_per_clust[["clust1_AABvsHD"]] = markers
# 
# 
# # Saving differential features
# ########
# # CHANGE
# ########
# openxlsx::write.xlsx(markers_per_clust,
#                      "/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Data/DEGenes_NormSCT_PC27_LSI15_ClustAlgo3Res0.1_DiffWilcox_Clust1BothCD4TCMTreg_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.xlsx")





#-------------------------------------------------------------------------------
# Subcluster of Cluster 1 (Treg) - Identifying Differential Features For
# T1D v/s HD, T1D v/s AAB+, AAB+ v/s HD
#-------------------------------------------------------------------------------

# Loading data
# seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothCD4TCMTregClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")
seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothCD4TCMTreg_SCT1Res0.1Clust1SubClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")




# # Defining variables
# active_ident = "bothcd4tcmtreg_SCT1_snn_res.0.1"




# # Sub-clustering cluster 1
# Idents(seurat_obj) = active_ident
# seurat_obj = FindSubCluster(seurat_obj,
#                             cluster = 1,
#                             graph.name = "SCT1_snn",
#                             resolution = 0.1,
#                             algorithm = 3,
#                             subcluster.name = "sub.cluster_sct1.res0.1.clust1")
# saveRDS(seurat_obj, file = "/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothCD4TCMTreg_SCT1Res0.1Clust1SubClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")



# # Subsetting Seurat object to remove cluster 5
# Idents(seurat_obj) = "sub.cluster_sct1.res0.1.clust1"
# seurat_obj_subset = subset(seurat_obj, idents = 5, invert = TRUE)
# 
# 
# 
# # SCT1 UMAP - grouped by cell type
# col_pal = c("CD4 TCM" = "#fc0fc0",
#             "Treg" = "#77b5fe")
# 
# DimPlot(seurat_obj_subset, repel = TRUE, reduction = "bothcd4tcmtreg.rna.umap",
#         group.by = "predicted.id.v1", cols = col_pal) +
#   ggtitle("UMAP RNA - Cell Type - CD4 TCM + Treg - PLNo + PLNn + Spleen") +
#   xlab("UMAP1") + ylab("UMAP2")
# 
# 
# 
# # Stacked barplot depicting tissue composition for each cluster
# tissue_v1 = seurat_obj_subset$tissue
# tissue_v1 = gsub("PLNo|PLNn", "PLN", tissue_v1)
# seurat_obj_subset$tissue_v1 = tissue_v1
# 
# 
# tmp1 = as.data.frame(table(seurat_obj_subset@meta.data[[active_ident]], seurat_obj_subset$tissue_v1))
# colnames(tmp1) = c("Cluster", "Tissue", "Freq")
# 
# # Converting cell type distribution per sample into percentage format
# tmp2 = as.data.frame(tmp1 %>% pivot_wider(names_from = Tissue, values_from = Freq))
# rownames(tmp2) = tmp2[, 1]
# tmp2 = tmp2[, -1]
# 
# tmp3 = (tmp2/rowSums(tmp2)) * 100
# 
# tmp4 = cbind.data.frame(rownames(tmp3), tmp3)
# colnames(tmp4) = c("Cluster", colnames(tmp3))
# 
# tmp5 = reshape2::melt(tmp4, id = "Cluster")
# colnames(tmp5) = c("Cluster", "Tissue", "Percent")
# 
# cols = c("PLN" = "#7b3294", "S" = "#008837")
# stacked_bp = ggplot(data = tmp5[complete.cases(tmp5), ], aes(x = Cluster, y = Percent, fill = Tissue)) +
#   geom_bar(stat = "identity", width = 0.7) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#   xlab("Cluster") +
#   ylab("Percent") +
#   scale_fill_manual(values = alpha(cols, .8)) +
#   theme_classic()
# stacked_bp
# 
# 
# 
# # SCT1 UMAP - grouped by cluster 1 sub-clusters
# col_pal = c("0" = "#bdbdbd",
#             "1_0" = "#78c679",
#             "1_1" = "#addd8e",
#             "1_2" = "#d9f0a3",
#             "2" = "#bdbdbd",
#             "3" = "#bdbdbd",
#             "4" = "#bdbdbd",
#             "6" = "#bdbdbd")
# 
# DimPlot(seurat_obj_subset, repel = TRUE, reduction = "bothcd4tcmtreg.rna.umap",
#         group.by = "sub.cluster_sct1.res0.1.clust1", cols = col_pal) +
#   ggtitle("UMAP RNA - SCT1 Res 0.1 Cluster 1 Subclusters - \n CD4 TCM+Treg - PLNo+PLNn+Spleen") +
#   xlab("UMAP1") + ylab("UMAP2")
# 
# 
# 
# # # Feature plot depicting markers
# # DefaultAssay(seurat_obj_subset) = "SCT1"
# # FeaturePlot(seurat_obj_subset,
# #             features = c("FOXP3", "CTLA4", "IL2RA"),
# #             reduction = "bothcd4tcmtreg.rna.umap")
# 
# 
# # Subsetting Seurat object to get cluster 1
# Idents(seurat_obj) = "sub.cluster_sct1.res0.1.clust1"
# seurat_obj_clust1 = subset(seurat_obj, idents = c(0, 2, 3, 4, 5, 6), invert = TRUE)
# 
# # Violin plot to check marker expression in sub cluster
# col_pal = c("1_0" = "#78c679",
#             "1_1" = "#addd8e",
#             "1_2" = "#d9f0a3")
# VlnPlot(seurat_obj_clust1,
#         features = c("FOXP3", "CTLA4", "IL2RA", "TOX", "TCF7", "PDCD1", "NFKB1"),
#         assay = "SCT1",
#         group.by = "sub.cluster_sct1.res0.1.clust1",
#         cols = col_pal)




# # Performing differential expression analysis
# Idents(seurat_obj) = "sub.cluster_sct1.res0.1.clust1"
# seurat_obj_treg = subset(seurat_obj, idents = c("0", "1_1", "2", "3", "4", "5", "6"), invert = TRUE)
# 
# 
# DefaultAssay(seurat_obj_treg) = "SCT1"
# Idents(seurat_obj_treg) = "status"
# 
# # Declaring list to store  differential features
# markers_per_clust = list()
# 
# # Identifying differential features for T1D v/s HD
# print("T1D v/s HD")
# markers = FindMarkers(seurat_obj_treg, ident.1 = "T1DM", ident.2 = "HD", slot = "data", logfc.threshold = 0, min.pct = 0)
# markers$features = rownames(markers)
# markers_per_clust[["T1D_vs_HD"]] = markers
# 
# # Identifying differential features for T1D v/s AAB+
# print("T1D v/s AAB+")
# markers = FindMarkers(seurat_obj_treg, ident.1 = "T1DM", ident.2 = "AAB+", slot = "data", logfc.threshold = 0, min.pct = 0)
# markers$features = rownames(markers)
# markers_per_clust[["T1D_vs_AAB"]] = markers
# 
# # Identifying differential features for T1D v/s HD
# print("AAB+ v/s HD")
# markers = FindMarkers(seurat_obj_treg, ident.1 = "AAB+", ident.2 = "HD", slot = "data", logfc.threshold = 0, min.pct = 0)
# markers$features = rownames(markers)
# markers_per_clust[["AAB_vs_HD"]] = markers
# 
# 
# # Saving differential features
# openxlsx::write.xlsx(markers_per_clust, "/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Data/DEGenes_NormSCT_PC27_LSI15_Algo3Res0.1Clust1SubClust_DiffWilcox_BothCD4TCMTreg_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.xlsx")
# 
# 
# 
# 
# # Performing differential accessibility analysis
# Idents(seurat_obj) = "sub.cluster_sct1.res0.1.clust1"
# seurat_obj_treg = subset(seurat_obj, idents = c("0", "1_1", "2", "3", "4", "5", "6"), invert = TRUE)
# 
# 
# DefaultAssay(seurat_obj_treg) = "peaks"
# Idents(seurat_obj_treg) = "status"
# 
# # Declaring list to store  differential features
# markers_per_clust = list()
# 
# # Identifying differential features for T1D v/s HD
# print("T1D v/s HD")
# markers = FindMarkers(seurat_obj_treg, ident.1 = "T1DM", ident.2 = "HD", slot = "data")
# markers$features = rownames(markers)
# markers_per_clust[["T1D_vs_HD"]] = markers
# 
# # Identifying differential features for T1D v/s AAB+
# print("T1D v/s AAB+")
# markers = FindMarkers(seurat_obj_treg, ident.1 = "T1DM", ident.2 = "AAB+", slot = "data")
# markers$features = rownames(markers)
# markers_per_clust[["T1D_vs_AAB"]] = markers
# 
# # Identifying differential features for T1D v/s HD
# print("AAB+ v/s HD")
# markers = FindMarkers(seurat_obj_treg, ident.1 = "AAB+", ident.2 = "HD", slot = "data")
# markers$features = rownames(markers)
# markers_per_clust[["AAB_vs_HD"]] = markers
# 
# 
# # Saving differential features
# openxlsx::write.xlsx(markers_per_clust, "/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Data/DAPeaks_NormSCT_PC27_LSI15_Algo3Res0.1Clust1SubClust_DiffWilcox_BothCD4TCMTreg_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.xlsx")




# Performing differential expression analysis after removing cells from spleen
Idents(seurat_obj) = "sub.cluster_sct1.res0.1.clust1"
seurat_obj_treg = subset(seurat_obj, idents = c("0", "1_1", "2", "3", "4", "5", "6"), invert = TRUE)

Idents(seurat_obj_treg) = "tissue"
seurat_obj_treg = subset(seurat_obj_treg, idents = c("S"), invert = TRUE)


DefaultAssay(seurat_obj_treg) = "SCT1"
Idents(seurat_obj_treg) = "status"

# Considering genes which are not ribosomal for differential analysis
genes.use = grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA", 
                 rownames(seurat_obj_treg), 
                 value = TRUE, 
                 invert = TRUE)  

# Declaring list to store  differential features
markers_per_clust = list()

# Identifying differential features for T1D v/s HD
print("T1D v/s HD")
# markers = FindMarkers(seurat_obj_treg, ident.1 = "T1DM", ident.2 = "HD", slot = "data", logfc.threshold = 0, min.pct = 0)
markers = FindMarkers(seurat_obj_treg, ident.1 = "T1DM", ident.2 = "HD", slot = "data", logfc.threshold = 0, min.pct = 0, features = genes.use)
markers$features = rownames(markers)
markers_per_clust[["T1D_vs_HD"]] = markers

# Identifying differential features for T1D v/s AAB+
print("T1D v/s AAB+")
# markers = FindMarkers(seurat_obj_treg, ident.1 = "T1DM", ident.2 = "AAB+", slot = "data", logfc.threshold = 0, min.pct = 0)
markers = FindMarkers(seurat_obj_treg, ident.1 = "T1DM", ident.2 = "AAB+", slot = "data", logfc.threshold = 0, min.pct = 0, features = genes.use)
markers$features = rownames(markers)
markers_per_clust[["T1D_vs_AAB"]] = markers

# Identifying differential features for T1D v/s HD
print("AAB+ v/s HD")
# markers = FindMarkers(seurat_obj_treg, ident.1 = "AAB+", ident.2 = "HD", slot = "data", logfc.threshold = 0, min.pct = 0)
markers = FindMarkers(seurat_obj_treg, ident.1 = "AAB+", ident.2 = "HD", slot = "data", logfc.threshold = 0, min.pct = 0, features = genes.use)
markers$features = rownames(markers)
markers_per_clust[["AAB_vs_HD"]] = markers


# Saving differential features
# openxlsx::write.xlsx(markers_per_clust, "/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Data/DEGenes_NormSCT_PC27_LSI15_Algo3Res0.1Clust1SubClust_DiffWilcox_BothCD4TCMTreg_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.xlsx")
# openxlsx::write.xlsx(markers_per_clust, "/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Data/DEGenes_NormSCT_PC27_LSI15_Algo3Res0.1Clust1SubClust_DiffWilcox_BothCD4TCMTreg_PLNo_PLNn_MedianAge_Subset24Donors.xlsx")
openxlsx::write.xlsx(markers_per_clust, "/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Data/NonRibosomalDEGenes_NormSCT_PC27_LSI15_Algo3Res0.1Clust1SubClust_DiffWilcox_BothCD4TCMTreg_PLNo_PLNn_MedianAge_Subset24Donors.xlsx")




# Performing differential accessibility analysis after removing cells from spleen
Idents(seurat_obj) = "sub.cluster_sct1.res0.1.clust1"
seurat_obj_treg = subset(seurat_obj, idents = c("0", "1_1", "2", "3", "4", "5", "6"), invert = TRUE)

Idents(seurat_obj_treg) = "tissue"
seurat_obj_treg = subset(seurat_obj_treg, idents = c("S"), invert = TRUE)


DefaultAssay(seurat_obj_treg) = "peaks"
Idents(seurat_obj_treg) = "status"

# Declaring list to store  differential features
markers_per_clust = list()

# Identifying differential features for T1D v/s HD
print("T1D v/s HD")
markers = FindMarkers(seurat_obj_treg, ident.1 = "T1DM", ident.2 = "HD", slot = "data")
markers$features = rownames(markers)
markers_per_clust[["T1D_vs_HD"]] = markers

# Identifying differential features for T1D v/s AAB+
print("T1D v/s AAB+")
markers = FindMarkers(seurat_obj_treg, ident.1 = "T1DM", ident.2 = "AAB+", slot = "data")
markers$features = rownames(markers)
markers_per_clust[["T1D_vs_AAB"]] = markers

# Identifying differential features for T1D v/s HD
print("AAB+ v/s HD")
markers = FindMarkers(seurat_obj_treg, ident.1 = "AAB+", ident.2 = "HD", slot = "data")
markers$features = rownames(markers)
markers_per_clust[["AAB_vs_HD"]] = markers


# Saving differential features
# openxlsx::write.xlsx(markers_per_clust, "/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Data/DAPeaks_NormSCT_PC27_LSI15_Algo3Res0.1Clust1SubClust_DiffWilcox_BothCD4TCMTreg_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.xlsx")
openxlsx::write.xlsx(markers_per_clust, "/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Data/DAPeaks_NormSCT_PC27_LSI15_Algo3Res0.1Clust1SubClust_DiffWilcox_BothCD4TCMTreg_PLNo_PLNn_MedianAge_Subset24Donors.xlsx")





#--------------------------------------------------------------------------------------
# Converting Cluster 2 and Cluster 3 Cells To .H5 For Cellarium Cell Annotation Service
#--------------------------------------------------------------------------------------

# Loading Seurat object
seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothCD4TCMTregClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")



# Subsetting Seurat object
Idents(seurat_obj) = "bothcd4tcmtreg_SCT1_snn_res.0.1"
seurat_obj_subset = subset(seurat_obj, idents = 3)



# Saving RNA-seq data as .h5 file
sct_counts = as.matrix(seurat_obj_subset@assays$SCT1@counts)


# Extracting feature/gene names
genes_sct_counts = as.data.frame(rownames(sct_counts))
colnames(genes_sct_counts) = c("Genes")

# Converting gene symbols to ENSEMBL IDs
genes_sct_counts$EnsemblID = mapIds(org.Hs.eg.db, 
                                    keys = genes_sct_counts$Genes, 
                                    column = "ENSEMBL", 
                                    keytype = "SYMBOL", 
                                    multiVals = "first")

na_indices = which(is.na(genes_sct_counts$EnsemblID) == TRUE)
sct_counts1 = sct_counts[-na_indices, ]
genes_sct_counts1 = genes_sct_counts[-na_indices, ]


# Saving data after removing genes whose ENSEMBL ID is not available
h5createFile("BothCD4TCMClustering_NormSCTAlgo3Res0.1Clust3_SCTCounts_V1.h5")
h5write(sct_counts1, "/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/CellAnnotation_Of_DiseaseCluster/BothCD4TCMClustering_NormSCTAlgo3Res0.1Clust3_SCTCounts_V1.h5", "RNA")

# Saving gene names after removing genes whose ENSEMBL ID is not available
write.csv(genes_sct_counts1, file = "/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/CellAnnotation_Of_DiseaseCluster/BothCD4TCMClustering_NormSCTAlgo3Res0.1Clust3_SCTCounts_Genes.csv", row.names = FALSE)

# Saving Ensembl ID 
feature_metadata = data.frame(EnsemblID = genes_sct_counts1[, -1])
rownames(feature_metadata) = genes_sct_counts1$Genes
write.csv(feature_metadata, file = "/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/CellAnnotation_Of_DiseaseCluster/BothCD4TCMClustering_NormSCTAlgo3Res0.1Clust3_SCTCounts_GenesMetadata.csv")

# Saving cell barcodes
cells_sct_counts = as.data.frame(colnames(sct_counts1))
colnames(cells_sct_counts) = c("Cells")
write.csv(cells_sct_counts, file = "/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/CellAnnotation_Of_DiseaseCluster/BothCD4TCMClustering_NormSCTAlgo3Res0.1Clust3_SCTCounts_Cells.csv", row.names = FALSE)





#-------------------------------------------------------------------------------
# Heatmap - Cytokines and Receptors, Transcription Factors, Cell Markers, and 
# Top DE Genes
#-------------------------------------------------------------------------------

# Loading data
seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothCD4TCMTregClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")


# Defining variables
active_ident = "bothcd4tcmtreg_SCT1_snn_res.0.1"

col_pal = c("0" = "#78c679",
            "1" = "#31a354",
            "2" = "#f768a1",
            "3" = "#c51b8a",
            "4" = "#feb24c",
            "6" = "#2b8cbe")


# Subsetting Seurat object to remove cluster 5
Idents(seurat_obj) = active_ident
seurat_obj_subset = subset(seurat_obj, idents = 5, invert = TRUE)


# Heatmap - cytokines and receptors
DoHeatmap(object = seurat_obj_subset, 
          features = c("IL2", "IL4", "IL5", "IL10", "IL17A", "IL17F", "IL6", 
                       "IL15", "IL21", "IL21R", "IL13", "IL18", "IL33", "IFNG", 
                       "IFNB1", "IFNAR2", "TNF", "IL23A", "IL2RB", "IL2RA", 
                       "IL7R", "TNFRSF1B", "TNFRSF1A", "TNFRSF25"), 
          group.by = active_ident, 
          group.colors = col_pal, 
          size = 3, 
          assay = "SCT1") + 
  # scale_fill_gradientn(colors = c("white", "gray", "yellow", "orange", "red"))
  # scale_fill_gradient2(low = "#ffffff", mid = "#fdf5e6", high = "#c51b8a")
  scale_fill_stepsn(colors = c("#ffffff", "#dee2e6", "#ffba08", "#ff477e", "#ff0a54"))

# Heatmap - transcription factors
DoHeatmap(object = seurat_obj_subset, 
          features = c("TBX21", "GATA3", "RORC", "RORA", "BACH2", "TCF7", "TOX", 
                       "MAF", "IKZF1", "IKZF2", "ZEB1", "ZEB2","FOXP3", "NFKB1", 
                       "NFKB2", "REL", "RELB", "NFAT5", "PRDM1", "BCL6", "NR4A1", 
                       "NR4A2", "IRF1"), 
          group.by = active_ident, 
          group.colors = col_pal, 
          size = 3, 
          assay = "SCT1") + 
  # scale_fill_gradientn(colors = c("white", "gray", "yellow", "orange", "red"))
  # scale_fill_gradient2(low = "#ffffff", mid = "#fdf5e6", high = "#c51b8a")
  scale_fill_stepsn(colors = c("#ffffff", "#dee2e6", "#ffba08", "#ff477e", "#ff0a54"))

# Heatmap - cell markers
DoHeatmap(object = seurat_obj_subset, 
          features = c("CD40LG", "CD19", "IL7R", "CD69", "CD99", "SELL", "CCR7", 
                       "CD27", "CR2", "CXCR4", "CXCR5", "TIGIT", "PDCD1", "LAG3", 
                       "CTLA4", "CD28", "ICOS", "CD44", "CD55"), 
          group.by = active_ident, 
          group.colors = col_pal, 
          size = 3, 
          assay = "SCT1") + 
  # scale_fill_gradientn(colors = c("white", "gray", "yellow", "orange", "red"))
  # scale_fill_gradient2(low = "#ffffff", mid = "#fdf5e6", high = "#c51b8a")
  scale_fill_stepsn(colors = c("#ffffff", "#dee2e6", "#ffba08", "#ff477e", "#ff0a54"))

# Heatmap - cluster 3 DE genes
DoHeatmap(object = seurat_obj_subset, 
          features = c("EZR", "PTPN22", "LMNA", "EZH2", "TXNIP", "SKIL", "SIK3", 
                       "HNRNPU", "ATP2B1", "GPR183"), 
          group.by = active_ident, 
          group.colors = col_pal, 
          size = 3, 
          assay = "SCT1") + 
  # scale_fill_gradientn(colors = c("white", "gray", "yellow", "orange", "red"))
  # scale_fill_gradient2(low = "#ffffff", mid = "#fdf5e6", high = "#c51b8a")
  scale_fill_stepsn(colors = c("#ffffff", "#dee2e6", "#ffba08", "#ff477e", "#ff0a54"))





#-------------------------------------------------------------------------------
# Lollipop Plot - (Number of CD4 TCM + Tregs) / (Donor)
#-------------------------------------------------------------------------------

# Loading data
seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothCD4TCMTregClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")


df1 = as.data.frame(table(seurat_obj$donor_ids_umap))
colnames(df1) = c("DonorIDs", "Freq")

df1$DonorIDs = factor(df1$DonorIDs, levels = rev(c("T1D1", "T1D2", "T1D4", "T1D6", "T1D10", "T1D12", "T1D15", "T1D17", 
                                                   "AAB+1", "AAB+2", "AAB+3", "AAB+4", "AAB+5", "AAB+6", "AAB+7", "AAB+8", 
                                                   "HD1", "HD3", "HD6", "HD10", "HD11", "HD13", "HD14", "HD16")))

ggplot(df1, aes(x = Freq, y = DonorIDs, label = Freq)) +
  geom_segment(aes(x = 0, y = DonorIDs, xend = Freq, yend = DonorIDs)) +
  geom_point(size = 8, color = "#feb24c") + 
  geom_text(size = 2) +
  labs(x = "Frequency", y = "Donor IDs")





#-------------------------------------------------------------------------------
# Coverage Plot
#-------------------------------------------------------------------------------

# Loading data
seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothCD4TCMTreg_PLNo_MedianAge_Subset24Donors_WithCD4TCMTregSCT1Algo3Res01Clusters.rds")




# Defining variables
features = c("CD4", "CD8A", "CD3DE", "TCF7", "LEF1", "BACH2")




# Removing cluster 5 from Seurat object
Idents(seurat_obj) = "cd4tcm.treg_sct1.algo3.res01_clusters"
seurat_obj = subset(seurat_obj, idents = 5, invert = TRUE)

levels(seurat_obj) = c(0, 1, 2, 3, 4, 6)



# Linking peaks to genes
DefaultAssay(seurat_obj) = "peaks"
seurat_obj = RegionStats(seurat_obj, genome = BSgenome.Hsapiens.UCSC.hg38, assay = "peaks")
seurat_obj = LinkPeaks(seurat_obj, 
                       peak.assay = "peaks", 
                       peak.slot = "data", 
                       expression.assay = "SCT", 
                       expression.slot = "data", 
                       genes.use = features)


# LOC finally used to generate coverage plot for paper figures

# CoveragePlot(seurat_obj, region = "EBF1", extend.upstream = 1000, 
#              extend.downstream = 1000, annotation = TRUE, peaks = FALSE, 
#              assay = "peaks")

p1 = CoveragePlot(seurat_obj, region = "CD4", annotation = TRUE, peaks = TRUE, assay = "peaks", links = TRUE, extend.upstream = 5000, extend.downstream = 5000)
p2 = CoveragePlot(seurat_obj, region = "CD8A", annotation = TRUE, peaks = TRUE, assay = "peaks", links = TRUE, extend.upstream = 5000, extend.downstream = 5000)
# p3 = CoveragePlot(seurat_obj, region = "CD3DE", annotation = TRUE, peaks = TRUE, assay = "peaks", links = TRUE, extend.upstream = 5000, extend.downstream = 5000)
p4 = CoveragePlot(seurat_obj, region = "TCF7", annotation = TRUE, peaks = TRUE, assay = "peaks", links = TRUE, extend.upstream = 5000, extend.downstream = 5000)
p5 = CoveragePlot(seurat_obj, region = "LEF1", annotation = TRUE, peaks = TRUE, assay = "peaks", links = TRUE, extend.upstream = 5000, extend.downstream = 5000)
p6 = CoveragePlot(seurat_obj, region = "BACH2", annotation = TRUE, peaks = TRUE, assay = "peaks", links = TRUE, extend.upstream = 5000, extend.downstream = 5000)
cowplot::plot_grid(p1, p2, p4, p5, p6, nrow = 2)





#-------------------------------------------------------------------------------
# Difference In Cytokine Expression Among T1D, AAB+, and Healthy Donors in 
# Cluster 0
#-------------------------------------------------------------------------------

#----- Loading data -----
seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothCD4TCMTregClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")




#----- Defining variables -----
active_ident = "bothcd4tcmtreg_SCT1_snn_res.0.1"

features = c("IFNA2", "IFNA1", "IFNAR1", "IFNAR2", "IFI27", "IFI30", "IFI35", 
             "IFI44", "IFI44L", "IFIH1", "IFIT2", "IFIT3", "IFITM1", "IFITM2", 
             "IFITM3", "IRF1", "IRF2", "IRF7", "IRF9", "IFNG", "CCL2", "CCL5", 
             "CCL7", "CXCL10", "CXCL11", "CXCL9", "IFIT1", "IRF4", "IRF5", "IRF8")




#----- Subsetting Seurat object -----
Idents(seurat_obj) = active_ident
seurat_obj_subset = subset(seurat_obj, idents = c(0))
Idents(seurat_obj_subset) = "tissue"
seurat_obj_subset = subset(seurat_obj_subset, idents = c("S"), invert = TRUE)


#----- Interferon expression across T1D, AAB+, and HD -----
col_pal = c("T1DM" = "#d7191c",
            "AAB+" = "#fdae61",
            "HD" = "#2c7bb6")

VlnPlot(seurat_obj_subset, features = features[21:30], ncol = 3, group.by = "status", 
        assay = "SCT1", col = col_pal)

# Idents(seurat_obj_subset) = "status"
# DoHeatmap(subset(seurat_obj_subset, downsample = 3000), features = features, 
#           size = 3, group.by = "status", group.colors = col_pal, assay = "SCT1")




#----- Performing differential expression analysis after removing cells from spleen -----
Idents(seurat_obj) = active_ident
seurat_obj_subset = subset(seurat_obj, idents = c("0"))

Idents(seurat_obj_subset) = "tissue"
seurat_obj_subset = subset(seurat_obj_subset, idents = c("S"), invert = TRUE)


DefaultAssay(seurat_obj_subset) = "SCT1"
Idents(seurat_obj_subset) = "status"

# Declaring list to store  differential features
markers_per_clust = list()

# Identifying differential features for T1D v/s HD
print("T1D v/s HD")
markers = FindMarkers(seurat_obj_subset, ident.1 = "T1DM", ident.2 = "HD", slot = "data", logfc.threshold = 0, min.pct = 0)
markers$features = rownames(markers)
markers_per_clust[["T1D_vs_HD"]] = markers

# Identifying differential features for T1D v/s AAB+
print("T1D v/s AAB+")
markers = FindMarkers(seurat_obj_subset, ident.1 = "T1DM", ident.2 = "AAB+", slot = "data", logfc.threshold = 0, min.pct = 0)
markers$features = rownames(markers)
markers_per_clust[["T1D_vs_AAB"]] = markers

# Identifying differential features for T1D v/s HD
print("AAB+ v/s HD")
markers = FindMarkers(seurat_obj_subset, ident.1 = "AAB+", ident.2 = "HD", slot = "data", logfc.threshold = 0, min.pct = 0)
markers$features = rownames(markers)
markers_per_clust[["AAB_vs_HD"]] = markers


# Saving differential features
openxlsx::write.xlsx(markers_per_clust, "/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Data/DEGenes_NormSCT_PC27_LSI15_Algo3Res0.1Clust0_DiffWilcox_BothCD4TCMTreg_PLNo_PLNn_MedianAge_Subset24Donors.xlsx")




#----- Performing differential accessibility analysis after removing cells from spleen -----
Idents(seurat_obj) = active_ident
seurat_obj_subset = subset(seurat_obj, idents = c("0"))

Idents(seurat_obj_subset) = "tissue"
seurat_obj_subset = subset(seurat_obj_subset, idents = c("S"), invert = TRUE)


DefaultAssay(seurat_obj_subset) = "peaks"
Idents(seurat_obj_subset) = "status"

# Declaring list to store  differential features
markers_per_clust = list()

# Identifying differential features for T1D v/s HD
print("T1D v/s HD")
markers = FindMarkers(seurat_obj_subset, ident.1 = "T1DM", ident.2 = "HD", slot = "data")
markers$features = rownames(markers)
markers_per_clust[["T1D_vs_HD"]] = markers

# Identifying differential features for T1D v/s AAB+
print("T1D v/s AAB+")
markers = FindMarkers(seurat_obj_subset, ident.1 = "T1DM", ident.2 = "AAB+", slot = "data")
markers$features = rownames(markers)
markers_per_clust[["T1D_vs_AAB"]] = markers

# Identifying differential features for T1D v/s HD
print("AAB+ v/s HD")
markers = FindMarkers(seurat_obj_subset, ident.1 = "AAB+", ident.2 = "HD", slot = "data")
markers$features = rownames(markers)
markers_per_clust[["AAB_vs_HD"]] = markers


# Saving differential features
openxlsx::write.xlsx(markers_per_clust, "/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Data/DAPeaks_NormSCT_PC27_LSI15_Algo3Res0.1Clust0_DiffWilcox_BothCD4TCMTreg_PLNo_PLNn_MedianAge_Subset24Donors.xlsx")





#-------------------------------------------------------------------------------
# Psuedobulk Differential Analysis - Subcluster of Cluster 1 (Treg) - 
# T1D v/s HD, T1D v/s AAB+, AAB+ v/s HD
#-------------------------------------------------------------------------------

#-------------
# Loading data
#-------------
seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothCD4TCMTreg_SCT1Res0.1Clust1SubClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")




#-------------------
# Defining variables
#-------------------

fig_dir = "/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Figures/BothCD4TCMTreg_PLNoPLNnSpleen_MedianAge_SubsetPLN8AAB.8T1D.8HD/"
data_dir = "/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Data/"

padj_cutoff = 0.05
log2fc_cutoff = 0.5

# Declaring list(s) to store DE genes
de_genes = list() # all
de_sig_genes = list() # significant

group = c("T1DM", "HD")




#-------------------------
# Subsetting Seurat object
#-------------------------
Idents(seurat_obj) = "sub.cluster_sct1.res0.1.clust1"
seurat = subset(seurat_obj, idents = c("0", "1_1", "2", "3", "4", "5", "6"), invert = TRUE)

Idents(seurat) = "tissue"
seurat = subset(seurat, idents = c("S"), invert = TRUE)

Idents(seurat) = "status"
seurat = subset(seurat, idents = group)




#---------------------------------------
# Creating Single Cell Experiment object
#---------------------------------------
# Extracting raw counts and metadata to create SingleCellExperiment (SCE) object
########
# CHANGE
########
counts = seurat@assays$SCT@counts
# counts = seurat@assays$peaks@counts
metadata = seurat@meta.data


# Setting up metadata as desired for aggregation and DE analysis
metadata$sample_id = seurat$donor_ids_umap
metadata$sample_id = factor(metadata$sample_id)

metadata$group_id = seurat$status


# Creating SCE object
sce = SingleCellExperiment(assays = list(counts = counts), colData = metadata)




#---------------------------
# Generating psuedobulk data
#---------------------------
# Extracting unique donor_ids which is equal to levels of sample_id factor variable
sample_names = levels(colData(sce)$sample_id)

# Identifying groups for aggregation of counts
groups = colData(sce)[, c("sample_id")]
head(groups)

# Aggregating across cluster-sample groups
# Transposing so that rownames of sce object matches rownames of groups
aggr_counts = aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum") 
dim(aggr_counts)

# Removing ribosomal genes
rp_genes = grep("^RP[LS]", colnames(aggr_counts), value = TRUE)
aggr_counts = aggr_counts[, !colnames(aggr_counts) %in% rp_genes]
dim(aggr_counts)




#--------------------------------------------------------
# Performing psuedobulk differential analysis using edgeR
#--------------------------------------------------------
# Transposing aggregated matrix to have genes as rows and samples as columns
counts = t(aggr_counts)
########
# CHANGE
########
group = gsub("(?<=HD|T1D)\\d+", "", colnames(counts), perl = TRUE)


# Constructing DGEList object
dge = DGEList(counts = counts)

# Normalizing counts matrix using TMM (trimmed means of m values) method
dge = calcNormFactors(dge)

design = model.matrix(~group)

dge = estimateDisp(dge, design)

fit = glmFit(dge, design)
lrt = glmLRT(fit)



# Table of results for all genes
tmp1 = topTags(lrt, n = nrow(counts), sort.by = "logFC", p.value = 1)$table
tmp1$gene = rownames(tmp1)
de_genes[[paste0(group[1], "_vs_", group[2])]] = tmp1

########
# CHANGE
########
comparison_check = grepl("T1D", topTags(lrt, n = nrow(counts), sort.by = "logFC", p.value = 1)$comparison, fixed = TRUE)



# Table of results for significant genes
# Subsetting the significant results
sig_res = topTags(lrt, n = nrow(counts), sort.by = "logFC", p.value = padj_cutoff)$table
sig_res$gene = rownames(sig_res)


if(!is.null(sig_res))
{
  sig_res = dplyr::filter(sig_res, abs(logFC) >= log2fc_cutoff)
  
  if(!is.null(sig_res))
  {
    de_sig_genes[[paste0(group[1], "_vs_", group[2])]] = sig_res$gene
  }
  
  # Extracting significantly up/down regulated genes
  n_sig_up = dplyr::filter(sig_res, logFC >= log2fc_cutoff) 
  n_sig_dn = dplyr::filter(sig_res, logFC <= -log2fc_cutoff)
}



# Heatmap of significant genes
if(!is.null(sig_res))
{
  if (nrow(sig_res) > 1)
  {
    
    # Extracting normalized counts from edgeR object
    normalized_counts = cpm(dge, log = FALSE)
    
    
    # Extracting normalized counts for significant genes only
    tmp1 = as.data.frame(normalized_counts[rownames(normalized_counts) %in% sig_res$gene, ])
    
    tmp2 = transpose(tmp1)
    rownames(tmp2) = colnames(tmp1)
    colnames(tmp2) = rownames(tmp1)
    
    tmp3 = cbind.data.frame(tmp2, cluster_metadata)
    tmp3 = tmp3[order(tmp3$group_id, decreasing = TRUE), ]
    
    tmp4 = transpose(tmp3)
    rownames(tmp4) = colnames(tmp3)
    colnames(tmp4) = rownames(tmp3)
    
    tmp5 = tmp4[1:nrow(tmp1), ]
    
    sig_counts = as.data.frame(sapply(tmp5, as.numeric))
    rownames(sig_counts) = rownames(tmp5)
    
    cluster_metadata1 = tmp3[, (ncol(tmp2) + 1):ncol(tmp3)]
    if("AAB+" %in% cluster_metadata1$group_id)
    {
      cluster_metadata1$group_id = replace(cluster_metadata1$group_id, cluster_metadata1$group_id == "AAB+", "AAB")
    }
    
    
    # Setting a color-blind friendly palette
    heat_colors = rev(brewer.pal(11, "PuOr"))
    
    # Setting color palette for donor status
    ########
    # CHANGE
    ########
    annot_colors = list(group_id = c(T1DM = "#ff0000", HD = "#9acd32"))
    # annot_colors = list(group_id = c(T1DM = "#ff0000", AAB = "#ff8c00"))
    # annot_colors = list(group_id = c(AAB = "#ff8c00", HD = "#9acd32"))
    
    
    ########
    # CHANGE
    ########
    pdf(paste0(fig_dir, "Heatmap_SubClust10And12edgeRPsuedobulkDEGenes.AdjpVal005.FC05_T1D.HD_CD4TCM.Treg_AgeMatchedPLNoPLNn.pdf"),
        width = 16,
        height = 10)
    
    # pheatmap using the metadata data frame for the annotation
    pheatmap(sig_counts,
             color = heat_colors,
             cluster_rows = TRUE,
             show_rownames = FALSE,
             annotation = cluster_metadata1[, c("group_id")], 
             annotation_colors = annot_colors, 
             border_color = NA,
             fontsize = 10,
             scale = "row",
             fontsize_row = 10,
             height = 20, 
             cluster_cols = FALSE)
    
    dev.off()
    
  }
}



# Saving results
openxlsx::write.xlsx(de_genes,
                     file = paste0(data_dir, "AllDEGs.SubClust10And12_edgeRPsuedobulk_T1D.HD_CD4TCM.Treg_AgeMatchedPLNoPLNn.xlsx"))
openxlsx::write.xlsx(de_sig_genes,
                     file = paste0(data_dir, "SignificantDEGsAdjpVal005FC05.SubClust10And12_edgeRPsuedobulk_T1D.HD_CD4TCM.Treg_AgeMatchedPLNoPLNn.xlsx"))





#-------------------------------------------------------------------------------
# Feature Plot, Violin Plot And Dot Plot - FOXP3 And TCF7
#-------------------------------------------------------------------------------

# Loading data
seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothCD4TCMTregClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")



# Defining variables
fig_dir = "/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Figures/BothCD4TCMTreg_PLNoPLNnSpleen_MedianAge_SubsetPLN8AAB.8T1D.8HD/"
active_ident = "bothcd4tcmtreg_SCT1_snn_res.0.1"
assay = "SCT1"
max_pc = 27
res = 0.1
cell_type = c("CD4 TCM", "Treg")
tissue = "PLNo_PLNn_Spleen_MedianAge_Subset24Donors"



# Subsetting Seurat object to remove cluster 5
Idents(seurat_obj) = active_ident
seurat_obj_subset = subset(seurat_obj, idents = 5, invert = TRUE)



# Feature plot and violin plot
DefaultAssay(seurat_obj_subset) = assay

features = c("RTKN2", "FOXP3", "AC133644.2", "CD4", "IL2RA", "TIGIT", "CTLA4", 
             "FCRL3", "LAIR2", "IKZF2")

# Defining colors for each cluster
col_pal = c("0" = "#78c679",
            "1" = "#31a354",
            "2" = "#f768a1",
            "3" = "#c51b8a",
            "4" = "#feb24c",
            "6" = "#2b8cbe")


# Feature plot and violin plot
for(i in 1:length(features))
{
  
  pdf(paste0(fig_dir, "FeaturePlot_UMAPNormSCT_PC", max_pc, "_Algo3Res", res, "_", gsub(" ", "", features[i]), "_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".pdf"),
      width = 6,
      height = 4)
  fp = FeaturePlot(seurat_obj_subset, reduction = "bothcd4tcmtreg.rna.umap", features = features[i])
  # fp = FeaturePlot(seurat_obj_subset, reduction = "bothcd4tcmtreg.rna.umap", features = features[i], max.cutoff = "q1")
  print(fp)
  dev.off()
  
  pdf(paste0(fig_dir, "ViolinPlot_NormSCT_PC", max_pc, "_Algo3Res", res, "_", gsub(" ", "", features[i]), "_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".pdf"),
      width = 8,
      height = 4)
  vp = VlnPlot(seurat_obj_subset, group.by = active_ident, cols = col_pal, features = features[i])
  print(vp)
  dev.off()
  
}



# Dot plot
DotPlot(seurat_obj_subset, 
        features = c("RTKN2", "FOXP3", "AC133644.2", "CD4", "IL2RA", "TIGIT", 
                     "CTLA4", "FCRL3", "LAIR2", "IKZF2"), 
        assay = "SCT1", 
        group.by = active_ident, 
        cols = "RdBu")










