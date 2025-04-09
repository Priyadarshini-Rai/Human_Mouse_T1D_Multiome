#-------------------------------------------------------------------------------
# MULTIOME PAPER 1 - REVISION 1 - FIGURE 6
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





# #-------------------------------------------------------------------------------
# # SCT1 Based UMAP - Grouped By Donor Status
# #-------------------------------------------------------------------------------
# 
# # Loading data
# seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothBintermediateBmemoryClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")
# 
# 
# # Defining color palette
# col_pal = c("T1DM" = "#d7191c",
#             "AAB+" = "#fdae61",
#             "HD" = "#2c7bb6")
# 
# # UMAP
# DimPlot(seurat_obj, repel = TRUE, reduction = "bothbintermediatebmemory.rna.umap", 
#         group.by = "status", cols = col_pal, raster = FALSE) +
#   ggtitle("UMAP RNA - Status - \n Bint + Bmem - PLNo + PLNn + Spleen") +
#   xlab("UMAP1") + ylab("UMAP2")
# 
# 
# 
# 
# 
# #-------------------------------------------------------------------------------
# # SCT1 Based UMAP - Grouped By SCT1 Res 0.1
# #-------------------------------------------------------------------------------
# 
# # Loading data
# seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothBintermediateBmemoryClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")
# 
# 
# # Defining variables
# active_ident = "bothbintermediatebmemory_SCT1_snn_res.0.1"
# 
# 
# # Defining colors for each cluster
# col_pal = c("0" = "#78c679",
#             "1" = "#31a354",
#             "2" = "#f768a1",
#             "3" = "#c51b8a",
#             "4" = "#ff7f50", 
#             "5" = "#1e90ff",
#             "6" = "#2b8cbe", 
#             "7" = "#ff0090", 
#             "8" = "#808000")
# 
# # UMAP
# DimPlot(seurat_obj, repel = TRUE, reduction = "bothbintermediatebmemory.rna.umap",
#         group.by = active_ident, cols = col_pal, raster = FALSE) +
#   ggtitle("UMAP RNA - Cluster SCT1 Res 0.1 - \n Bint + Bmem - PLNo + PLNn + Spleen") +
#   xlab("UMAP1") + ylab("UMAP2")
# 
# 
# 
# 
# 
# #-------------------------------------------------------------------------------
# # Stacked Bar Plot - Donor Composition Per Cluster
# #-------------------------------------------------------------------------------
# 
# # Loading data
# seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothBintermediateBmemoryClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")
# 
# 
# 
# # Defining variables
# active_ident = "bothbintermediatebmemory_SCT1_snn_res.0.1"
# 
# 
# 
# # Stacked barplot depicting donor status composition for each cluster
# tmp1 = as.data.frame(table(seurat_obj@meta.data[[active_ident]], seurat_obj$status))
# colnames(tmp1) = c("Cluster", "Status", "Freq")
# 
# # Converting cell type distribution per sample into percentage format
# tmp2 = as.data.frame(tmp1 %>% pivot_wider(names_from = Status, values_from = Freq))
# rownames(tmp2) = tmp2[, 1]
# tmp2 = tmp2[, -1]
# 
# tmp3 = (tmp2/rowSums(tmp2)) * 100
# 
# tmp4 = cbind.data.frame(rownames(tmp3), tmp3)
# colnames(tmp4) = c("Cluster", colnames(tmp3))
# 
# tmp5 = reshape2::melt(tmp4, id = "Cluster")
# colnames(tmp5) = c("Cluster", "Status", "Percent")
# 
# tmp5$Status = factor(tmp5$Status, levels = c("T1DM", "AAB+", "HD"))
# 
# 
# # Defining color plaette
# col_pal = c("T1DM" = "#d7191c", "AAB+" = "#fdae61", "HD" = "#2c7bb6")
# 
# # Stacked bar plot depicting donor composition per cluster
# stacked_bp = ggplot(data = tmp5, aes(x = Cluster, y = Percent, fill = Status)) +
#   geom_bar(stat = "identity") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#   xlab("Cluster") +
#   ylab("Percent") +
#   scale_fill_manual(values = col_pal)
# 
# 
# 
# 
# 
# #-------------------------------------------------------------------------------
# # Donut Plot Depicting Donor Composition Per Cluster
# #-------------------------------------------------------------------------------
# 
# # Loading data
# seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothBintermediateBmemoryClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")
# 
# 
# 
# # Defining variables
# active_ident = "bothbintermediatebmemory_SCT1_snn_res.0.1"
# 
# 
# 
# # Subsetting Seurat object
# Idents(seurat_obj) = active_ident
# #-------
# # CHANGE
# #-------
# seurat_obj_subset = subset(seurat_obj, idents = c(3, 7))
# 
# 
# # Donut plot
# donutplot_df = data.frame(table(seurat_obj_subset$donor_ids_umap))
# colnames(donutplot_df) = c("donor.id", "freq")
# 
# # Computing percentages
# donutplot_df$fraction = donutplot_df$freq / sum(donutplot_df$freq)
# 
# # Computing cumulative percentages (top of each rectangle)
# donutplot_df$ymax = cumsum(donutplot_df$fraction)
# 
# # Computing bottom of each rectangle
# donutplot_df$ymin = c(0, head(donutplot_df$ymax, n = -1))
# 
# # Compute label position
# donutplot_df$label.position = (donutplot_df$ymax + donutplot_df$ymin) / 2
# 
# # Computing good label
# donutplot_df$label = paste0(donutplot_df$donor.id, "\n value: ", donutplot_df$freq)
# 
# # Defining color palette
# col_pal = c("T1D1" = "#d94801",
#             "T1D2" = "#ef3b2c",
#             "T1D4" = "#cb181d",
#             "T1D6" = "#e31a1c",
#             "T1D10" = "#bd0026",
#             "T1D12" = "#a50f15",
#             "T1D15" = "#800026",
#             "T1D17" = "#67000d",
#             "AAB+1" = "#ffffcc",
#             "AAB+2" = "#ffeda0",
#             "AAB+3" = "#fdd0a2",
#             "AAB+4" = "#fed976",
#             "AAB+5" = "#fdae6b",
#             "AAB+6" = "#feb24c",
#             "AAB+7" = "#fd8d3c",
#             "AAB+8" = "#fd8d3c",
#             "HD1" = "#9ecae1",
#             "HD3" = "#6baed6",
#             "HD6" = "#4eb3d3",
#             "HD10" = "#4292c6",
#             "HD11" = "#2b8cbe",
#             "HD13" = "#2171b5",
#             "HD14" = "#08519c",
#             "HD16" = "#08306b")
# 
# # Plotting
# donut_plot = ggplot(donutplot_df, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = donor.id), col) +
#   geom_rect() +
#   geom_label(x = 3.5, aes(y = label.position, label = label), size = 2) +
#   # geom_text(x = 2, aes(y = label.position, label = label, color = donor.id), size = 1) + # x here controls label position (inner/outer)
#   scale_fill_manual(values = col_pal, name = "Donor IDs") +
#   coord_polar(theta = "y") +
#   xlim(c(2, 4)) +
#   # xlim(c(-1, 4)) +
#   theme_void()
# # Ordering donor ids
# donut_plot$data$donor.id = factor(x = donut_plot$data$donor.id,
#                                   levels = c("T1D1", "T1D2", "T1D4", "T1D6", "T1D10", "T1D12", "T1D15", "T1D17",
#                                              "AAB+1", "AAB+2", "AAB+3", "AAB+4", "AAB+5", "AAB+6", "AAB+7", "AAB+8",
#                                              "HD1", "HD3", "HD6", "HD10", "HD11", "HD13", "HD14", "HD16"))
# donut_plot
# 
# 
# 
# 
# 
# #-------------------------------------------------------------------------------
# # Violin Plot - Gene Expression Per Cluster
# #-------------------------------------------------------------------------------
# 
# # Loading data
# seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothBintermediateBmemoryClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")
# 
# 
# # Defining variables
# active_ident = "bothbintermediatebmemory_SCT1_snn_res.0.1"
# 
# 
# # Defining color palette
# col_pal = c("0" = "#78c679",
#             "1" = "#31a354",
#             "2" = "#f768a1",
#             "3" = "#c51b8a",
#             "4" = "#ff7f50", 
#             "5" = "#1e90ff",
#             "6" = "#2b8cbe", 
#             "7" = "#ff0090", 
#             "8" = "#808000")
# 
# # Violin plot depicting gene expression per cluster
# VlnPlot(seurat_obj, features = "NFKB1", assay = "SCT1", slot = "data", 
#         group.by = active_ident, cols = col_pal)
# 
# 
# 
# 
# 
# #-------------------------------------------------------------------------------
# # SCT1 UMAP - Grouped by Tissue
# #-------------------------------------------------------------------------------
# 
# # Loading data
# seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothBintermediateBmemoryClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")
# 
# 
# # Defining variables
# active_ident = "bothbintermediatebmemory_SCT1_snn_res.0.1"
# 
# 
# # Defining color palette
# # col_pal = c("PLNo" = "#f20adb",
# #             "PLNn" = "#7b3294",
# #             "S" = "#008837")
# col_pal = c("PLNo" = "#7b3294",
#             "PLNn" = "#7b3294",
#             "S" = "#008837")
# 
# # UMAP
# DimPlot(seurat_obj, repel = TRUE, reduction = "bothbintermediatebmemory.rna.umap", 
#         group.by = "tissue", cols = col_pal, raster = FALSE) +
#   ggtitle("UMAP RNA - Tissue - \n Bint + Bmem - PLNo + PLNn + Spleen") +
#   xlab("UMAP1") + ylab("UMAP2")
# 
# 
# 
# 
# 
# #-------------------------------------------------------------------------------
# # Feature Plot of Cytokines
# #-------------------------------------------------------------------------------
# 
# # Loading data
# seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothBintermediateBmemoryClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")
# 
# 
# 
# # Defining variables
# fig_dir = "/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Figures/BothBintBmem_PLNoPLNnSpleen_MedianAge_SubsetPLN8AAB.8T1D.8HD/"
# active_ident = "bothbintermediatebmemory_SCT1_snn_res.0.1"
# assay = "SCT1"
# max_pc = 27
# res = 0.1
# cell_type = c("B intermediate", "B memory")
# tissue = "PLNo_PLNn_Spleen_MedianAge_Subset24Donors"
# 
# 
# 
# # Feature plot and violin plot of cytokines
# DefaultAssay(seurat_obj) = assay
# 
# features = c("IFNG", "TBX21", "GATA3", "IL4", "IL5", "IL15", "BCL6", "IL21", 
#              "ICOS", "IL2", "CXCR5", "RORC")
# 
# # Defining colors for each cluster
# col_pal = c("0" = "#78c679",
#             "1" = "#31a354",
#             "2" = "#f768a1",
#             "3" = "#c51b8a",
#             "4" = "#ff7f50", 
#             "5" = "#1e90ff",
#             "6" = "#2b8cbe", 
#             "7" = "#ff0090", 
#             "8" = "#808000")
# 
# 
# for(i in 1:length(features))
# {
#   
#   pdf(paste0(fig_dir, "FeaturePlot_UMAPNormSCT_PC", max_pc, "_Algo3Res", res, "_", gsub(" ", "", features[i]), "_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".pdf"),
#       width = 6,
#       height = 4)
#   fp = FeaturePlot(seurat_obj, reduction = "bothbintermediatebmemory.rna.umap", features = features[i])
#   print(fp)
#   dev.off()
#   
#   pdf(paste0(fig_dir, "ViolinPlot_NormSCT_PC", max_pc, "_Algo3Res", res, "_", gsub(" ", "", features[i]), "_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".pdf"),
#       width = 8,
#       height = 4)
#   vp = VlnPlot(seurat_obj, group.by = active_ident, cols = col_pal, features = features[i])
#   print(vp)
#   dev.off()
#   
# }
# 
# 
# 
# 
# 
# #-------------------------------------------------------------------------------
# # Heatmap - Cytokines and Receptors, Transcription Factors, Cell Markers, and 
# # Top DE Genes
# #-------------------------------------------------------------------------------
# 
# # Loading data
# seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothBintermediateBmemoryClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")
# 
# 
# # Defining variables
# active_ident = "bothbintermediatebmemory_SCT1_snn_res.0.1"
# 
# col_pal = c("0" = "#78c679",
#             "1" = "#31a354",
#             "2" = "#f768a1",
#             "3" = "#c51b8a",
#             "4" = "#ff7f50", 
#             "5" = "#1e90ff",
#             "6" = "#2b8cbe", 
#             "7" = "#ff0090", 
#             "8" = "#808000")
# 
# 
# # Heatmap - cytokines and receptors
# DoHeatmap(object = seurat_obj, 
#           features = c("BACH2", "IRF4", "PRDM1", "BCL6", "XBP1", "EBF1", "PAX5", 
#                        "MYC", "BATF", "NFKB1", "REL", "BCL2"), 
#           group.by = active_ident, 
#           group.colors = col_pal, 
#           size = 3, 
#           assay = "SCT1") + 
#   # scale_fill_gradientn(colors = c("white", "gray", "yellow", "orange", "red"))
#   # scale_fill_gradient2(low = "#ffffff", mid = "#fdf5e6", high = "#c51b8a")
#   scale_fill_stepsn(colors = c("#ffffff", "#dee2e6", "#ffba08", "#ff477e", "#ff0a54"))
# 
# # Heatmap - transcription factors
# DoHeatmap(object = seurat_obj, 
#           features = c("ELK4", "ELK1", "NFKB1", "NFKB2", "MEF2D", "FLI1", "LRF", 
#                        "EBF", "BMYB", "GABPA", "ZAC1", "EBF2", "NFAT", "ERRG", 
#                        "ETV4", "BACH2", "MEF2B", "ELF1", "ETS", "EWS", "MYB", 
#                        "FOXD3", "GATA", "MYNN", "E2F6", "Olig2"), 
#           group.by = active_ident, 
#           group.colors = col_pal, 
#           size = 3, 
#           assay = "SCT1") + 
#   # scale_fill_gradientn(colors = c("white", "gray", "yellow", "orange", "red"))
#   # scale_fill_gradient2(low = "#ffffff", mid = "#fdf5e6", high = "#c51b8a")
#   scale_fill_stepsn(colors = c("#ffffff", "#dee2e6", "#ffba08", "#ff477e", "#ff0a54"))
# 
# # Heatmap - cell markers
# DoHeatmap(object = seurat_obj, 
#           features = c("CD80", "CD86", "CD10", "CD19", "MS4A1", "CD21", "CD38", 
#                        "CD40", "CD93", "CD95", "CD148", "HLA-DR", "TACI"), 
#           group.by = active_ident, 
#           group.colors = col_pal, 
#           size = 3, 
#           assay = "SCT1") + 
#   # scale_fill_gradientn(colors = c("white", "gray", "yellow", "orange", "red"))
#   # scale_fill_gradient2(low = "#ffffff", mid = "#fdf5e6", high = "#c51b8a")
#   scale_fill_stepsn(colors = c("#ffffff", "#dee2e6", "#ffba08", "#ff477e", "#ff0a54"))
# 
# # Heatmap - cluster 3 DE genes
# DoHeatmap(object = seurat_obj, 
#           features = c("HSPD1", "CD83", "EIF2AK3", "ZNF331", "KLF6", "IER5", 
#                        "JUN", "REL", "EZR", "NFKB1"), 
#           group.by = active_ident, 
#           group.colors = col_pal, 
#           size = 3, 
#           assay = "SCT1") + 
#   # scale_fill_gradientn(colors = c("white", "gray", "yellow", "orange", "red"))
#   # scale_fill_gradient2(low = "#ffffff", mid = "#fdf5e6", high = "#c51b8a")
#   scale_fill_stepsn(colors = c("#ffffff", "#dee2e6", "#ffba08", "#ff477e", "#ff0a54"))





#-------------------------------------------------------------------------------
# Differences in Gene Expression of Non-Naive B Cells Between PLN and Spleen of 
# T1D, HD, and AAB+
#-------------------------------------------------------------------------------

#--- Loading data ---
seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothBintermediateBmemoryClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")



#----- Performing differential expression analysis -----
# Declaring list to store  differential features
markers_per_clust = list()


Idents(seurat_obj) = "status"
seurat_obj_subset = subset(seurat_obj, idents = c("T1DM"))

DefaultAssay(seurat_obj_subset) = "SCT1"
Idents(seurat_obj_subset) = "tissue"
markers = FindMarkers(seurat_obj_subset, ident.1 = c("PLNo", "PLNn"), ident.2 = c("S"), slot = "data", logfc.threshold = 0, min.pct = 0)
markers$features = rownames(markers)
markers_per_clust[["PLN_vs_Spleen_T1D"]] = markers


Idents(seurat_obj) = "status"
seurat_obj_subset = subset(seurat_obj, idents = c("AAB+"))

DefaultAssay(seurat_obj_subset) = "SCT1"
Idents(seurat_obj_subset) = "tissue"
markers = FindMarkers(seurat_obj_subset, ident.1 = c("PLNo"), ident.2 = c("S"), slot = "data", logfc.threshold = 0, min.pct = 0)
markers$features = rownames(markers)
markers_per_clust[["PLN_vs_Spleen_AAB"]] = markers


Idents(seurat_obj) = "status"
seurat_obj_subset = subset(seurat_obj, idents = c("HD"))

DefaultAssay(seurat_obj_subset) = "SCT1"
Idents(seurat_obj_subset) = "tissue"
markers = FindMarkers(seurat_obj_subset, ident.1 = c("PLNo", "PLNn"), ident.2 = c("S"), slot = "data", logfc.threshold = 0, min.pct = 0)
markers$features = rownames(markers)
markers_per_clust[["PLN_vs_Spleen_HD"]] = markers


# Saving differential features
openxlsx::write.xlsx(markers_per_clust, "/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Data/DEGenes_PLNvsSpleenAcrossDonorStatus_DiffWilcox_BothBintermediateBmemory_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.xlsx")





#-------------------------------------------------------------------------------
# Dotplot
#-------------------------------------------------------------------------------

# Setting directory to save figures
setwd("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Figures/BothBintBmem_PLNoPLNnSpleen_MedianAge_SubsetPLN8AAB.8T1D.8HD/")

# Loading data
seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothBintermediateBmemoryClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")

# Defining variables
group = "bothbintermediatebmemory_SCT1_snn_res.0.1"

# Subsetting Seurat object
Idents(seurat_obj) = "bothbintermediatebmemory_SCT1_snn_res.0.1"
seurat_obj_subset = subset(seurat_obj, idents = c(0, 3, 6, 7))
DefaultAssay(seurat_obj_subset) = "SCT1"



#----- Cytokines -----
pdf("Dotplot_Cytokines_Bcells.pdf", width = 30, height = 10)
DotPlot(seurat_obj_subset, 
        features = c("IL2", "IL4", "IL5", "IL10", "IL17A", "IL17F", "IL6", 
                     "IL15", "IL21", "IL13", "IL18", "IL33", "IFNG", "IFNB1", 
                     "TNFRSF1B", "TNFRSF1A", "TNFRSF25", "TNF", "IL2RB", "IL2RA", 
                     "IL21R"),
        group.by = group, 
        cols = "RdBu", 
        dot.scale = 15)
dev.off()


#----- TFs -----
pdf("Dotplot_TFs_Bcells.pdf", width = 30, height = 10)
DotPlot(seurat_obj_subset, 
        features = c("EBF1", "PAX5", "MYC", "BATF", "RORA", "BACH2", "MAF", 
                     "IKZF1", "IKZF2", "ZEB1", "NFKB1", "NFKB2", "REL", "RELB", 
                     "NFAT5",  "NR4A1", "NR4A2", "IRF1", "IRF4", "FOS", "JUN", 
                     "JUNB", "TBX21", "PRDM1", "BCL6"), 
        group.by = group, 
        cols = "RdBu", 
        dot.scale = 15) 
dev.off()


#----- Cell surface markers -----
pdf("Dotmap_Markers_Bcells_V1_F.pdf", width = 30, height = 10)
DotPlot(seurat_obj_subset, 
        features = c("CD80", "CD86", "CD19", "MS4A1", "CD38", "CD40", "CD27", 
                     "CXCR4", "CD44", "CD55", "CXCR5", "CR2", "CD69", "TNFRSF13B"), 
        group.by = group, 
        cols = "RdBu", 
        dot.scale = 15) 
dev.off()


#----- Differential genes -----
pdf("Dotmap_DEG_Bcells.pdf", width = 30, height = 10)
DotPlot(seurat_obj_subset, 
        features = c("EZR", "HSPD1", "CD83", "EIF2AK3", "ZNF331", "KLF6", "IER5", 
                     "HES4", "CXCR4"), 
        group.by = group, 
        cols = "RdBu", 
        dot.scale = 15)
dev.off()





#-------------------------------------------------------------------------------
# Mapping Metadata "bothbintermediatebmemory_SCT1_snn_res.0.1" To PLNo And 
# Spleen Seurat Object
#-------------------------------------------------------------------------------

# Loading data 

# PLNo/Spleen Seurat object
########
# CHANGE
########
s1 = readRDS("/mnt/alvand/priya/PLN_Aggregate/R_WorkSpace/PLN_Seurat_Obj_AddedGeneticRiskScore_04.15.2024_V15.rds")
# s1 = readRDS("/mnt/alvand/priya/Spleen_Aggregate/R_WorkSpace/Spleen_Seurat_Obj_AddedGeneticRiskScore_06.07.2024_V10.rds")

s2 = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothBintermediateBmemoryClustering_NormSCTPC27_NormPeaksLSI15_BothSCTPeaks_PLNo_PLNn_Spleen_MedianAge_Subset24Donors.rds")




# Subsetting to get cells corresponding to age matched donors
Idents(s1) = "donor_ids_umap"

########
# CHANGE
########
s1_subset = subset(s1, idents = c("T1D1", "T1D2", "T1D4", "T1D6", "T1D10", "T1D12",
                                  "AAB+1", "AAB+2", "AAB+3", "AAB+4", "AAB+5", "AAB+6", "AAB+7", "AAB+8",
                                  "HD1", "HD3", "HD6", "HD10", "HD11", "HD13", "HD14")) # PLNo

# s1_subset = subset(s1, idents = c("T1D1", "T1D2", "T1D6", "T1D10", "T1D12", 
#                                   "AAB+1", "AAB+2", "AAB+3", "AAB+4", "AAB+5", "AAB+6", "AAB+7", "AAB+8", 
#                                   "HD3", "HD6", "HD10", "HD13", "HD14")) # Spleen


# Subsetting to get B intermediate and B memory cells 
Idents(s1_subset) = "predicted.id.v1"
s1_subset = subset(s1_subset, idents = c("B intermediate", "B memory"))

barcodes = s1_subset@meta.data[["barcodes"]]
########
# CHANGE
########
barcodes = paste0("PLNo_", barcodes)
# barcodes = paste0("S_", barcodes)
bint.bmem_sct1.algo3.res01_clusters = barcodes



Idents(s2) = "bothbintermediatebmemory_SCT1_snn_res.0.1"

s2_subset1 = subset(s2, idents = 0)
clust_barcodes1 = s2_subset1@assays[["SCT1"]]@data@Dimnames[[2]]
ind = match(clust_barcodes1, barcodes)
ind = ind[!is.na(ind)]
bint.bmem_sct1.algo3.res01_clusters[ind] = 0

s2_subset2 = subset(s2, idents = 1)
clust_barcodes2 = s2_subset2@assays[["SCT1"]]@data@Dimnames[[2]]
ind = match(clust_barcodes2, barcodes)
ind = ind[!is.na(ind)]
bint.bmem_sct1.algo3.res01_clusters[ind] = 1

s2_subset3 = subset(s2, idents = 2)
clust_barcodes3 = s2_subset3@assays[["SCT1"]]@data@Dimnames[[2]]
ind = match(clust_barcodes3, barcodes)
ind = ind[!is.na(ind)]
bint.bmem_sct1.algo3.res01_clusters[ind] = 2

s2_subset4 = subset(s2, idents = 3)
clust_barcodes4 = s2_subset4@assays[["SCT1"]]@data@Dimnames[[2]]
ind = match(clust_barcodes4, barcodes)
ind = ind[!is.na(ind)]
bint.bmem_sct1.algo3.res01_clusters[ind] = 3

s2_subset5 = subset(s2, idents = 4)
clust_barcodes5 = s2_subset5@assays[["SCT1"]]@data@Dimnames[[2]]
ind = match(clust_barcodes5, barcodes)
ind = ind[!is.na(ind)]
bint.bmem_sct1.algo3.res01_clusters[ind] = 4

s2_subset6 = subset(s2, idents = 5)
clust_barcodes6 = s2_subset6@assays[["SCT1"]]@data@Dimnames[[2]]
ind = match(clust_barcodes6, barcodes)
ind = ind[!is.na(ind)]
bint.bmem_sct1.algo3.res01_clusters[ind] = 5

s2_subset7 = subset(s2, idents = 6)
clust_barcodes7 = s2_subset7@assays[["SCT1"]]@data@Dimnames[[2]]
ind = match(clust_barcodes7, barcodes)
ind = ind[!is.na(ind)]
bint.bmem_sct1.algo3.res01_clusters[ind] = 6

s2_subset8 = subset(s2, idents = 7)
clust_barcodes8 = s2_subset8@assays[["SCT1"]]@data@Dimnames[[2]]
ind = match(clust_barcodes8, barcodes)
ind = ind[!is.na(ind)]
bint.bmem_sct1.algo3.res01_clusters[ind] = 7

s2_subset9 = subset(s2, idents = 8)
clust_barcodes9 = s2_subset9@assays[["SCT1"]]@data@Dimnames[[2]]
ind = match(clust_barcodes9, barcodes)
ind = ind[!is.na(ind)]
bint.bmem_sct1.algo3.res01_clusters[ind] = 8


s1_subset$bint.bmem_sct1.algo3.res01_clusters = bint.bmem_sct1.algo3.res01_clusters
########
# CHANGE
########
saveRDS(s1_subset, "/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothBintBmem_PLNo_MedianAge_Subset24Donors_WithBintBmemSCT1Algo3Res01Clusters.rds")
# saveRDS(s1_subset, "/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothBintBmem_Spleen_MedianAge_Subset24Donors_WithBintBmemSCT1Algo3Res01Clusters.rds")





#-------------------------------------------------------------------------------
# Coverage Plot - PLNo
#-------------------------------------------------------------------------------

#----- Setting working directory -----
setwd("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Figures/BothBintBmem_PLNoPLNnSpleen_MedianAge_SubsetPLN8AAB.8T1D.8HD/")



#----- Loading data -----
seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothBintBmem_PLNo_MedianAge_Subset24Donors_WithBintBmemSCT1Algo3Res01Clusters.rds")

# Subsetting Seurat object
Idents(seurat_obj) = "bint.bmem_sct1.algo3.res01_clusters"
seurat_obj_subset = subset(seurat_obj, idents = c(0, 3, 6, 7))

# Defining variables
col_pal = c("#da70d6", "#7b3294", "#9400d3", "#9966cc")



#----- Preparing data for coverage plot -----
DefaultAssay(seurat_obj_subset) = "peaks"
Idents(seurat_obj_subset) = "bint.bmem_sct1.algo3.res01_clusters"
seurat_obj_subset = RegionStats(seurat_obj_subset, 
                                genome = BSgenome.Hsapiens.UCSC.hg38, 
                                assay = "peaks")


#----- 1. BACH2 -----
seurat_obj_subset1 = LinkPeaks(seurat_obj_subset, 
                               peak.assay = "peaks", 
                               peak.slot = "counts", 
                               expression.assay = "SCT", 
                               expression.slot = "data", 
                               genes.use = "BACH2", 
                               distance = 1000000, 
                               score_cutoff = 0.03)

pdf("CoveragePlot_BintBmem_NormSCTAlgo3Res0.1.C0.C3.C6.C7_BACH2_PLNo_V2.pdf", width = 30, height = 10)
cp = CoveragePlot(seurat_obj_subset1, 
                  region = "BACH2", 
                  annotation = TRUE, 
                  peaks = TRUE, 
                  assay = "peaks", 
                  links = TRUE, 
                  extend.upstream = 400000, 
                  extend.downstream = 400000, 
                  expression.assay = "SCT", 
                  group.by = "bint.bmem_sct1.algo3.res01_clusters", 
                  ymax = 2)
cp & scale_fill_manual(values = col_pal)
print(cp)
dev.off()


#----- 2. NFKB1 -----
seurat_obj_subset1 = LinkPeaks(seurat_obj_subset, 
                               peak.assay = "peaks", 
                               peak.slot = "counts", 
                               expression.assay = "SCT", 
                               expression.slot = "data", 
                               genes.use = "NFKB1", 
                               distance = 1000000, 
                               score_cutoff = 0.01)

pdf("CoveragePlot_BintBmem_NormSCTAlgo3Res0.1.C0.C3.C6.C7_NFKB1_PLNo_V1.pdf", width = 30, height = 10)
cp = CoveragePlot(seurat_obj_subset1, 
                  region = "NFKB1", 
                  annotation = TRUE, 
                  peaks = TRUE, 
                  assay = "peaks", 
                  links = TRUE, 
                  extend.upstream = 200000, 
                  extend.downstream = 200000, 
                  expression.assay = "SCT", 
                  group.by = "bint.bmem_sct1.algo3.res01_clusters")
cp & scale_fill_manual(values = col_pal)
print(cp)
dev.off()


#----- 3. BCL2 -----
seurat_obj_subset1 = LinkPeaks(seurat_obj_subset, 
                               peak.assay = "peaks", 
                               peak.slot = "counts", 
                               expression.assay = "SCT", 
                               expression.slot = "data", 
                               genes.use = "BCL2", 
                               distance = 1000000, 
                               score_cutoff = 0.01)

pdf("CoveragePlot_BintBmem_NormSCTAlgo3Res0.1.C0.C3.C6.C7_BCL2_PLNo_V1.pdf", width = 30, height = 10)
cp = CoveragePlot(seurat_obj_subset1, 
                  region = "BCL2", 
                  annotation = TRUE, 
                  peaks = TRUE, 
                  assay = "peaks", 
                  links = TRUE, 
                  extend.upstream = 20000, 
                  extend.downstream = 20000, 
                  expression.assay = "SCT", 
                  expression.slot = "data", 
                  group.by = "bint.bmem_sct1.algo3.res01_clusters", 
                  ymax = 2)
cp & scale_fill_manual(values = col_pal)
print(cp)
dev.off()


#----- 4. EBF1 -----
seurat_obj_subset1 = LinkPeaks(seurat_obj_subset, 
                               peak.assay = "peaks", 
                               peak.slot = "counts", 
                               expression.assay = "SCT", 
                               expression.slot = "data", 
                               genes.use = "EBF1", 
                               distance = 1000000, 
                               score_cutoff = 0.01)

pdf("CoveragePlot_BintBmem_NormSCTAlgo3Res0.1.C0.C3.C6.C7_EBF1_PLNo_V1.pdf", width = 30, height = 10)
cp = CoveragePlot(seurat_obj_subset1, 
                  region = "EBF1", 
                  annotation = TRUE, 
                  peaks = TRUE, 
                  assay = "peaks", 
                  links = TRUE, 
                  extend.upstream = 20000, 
                  extend.downstream = 20000, 
                  expression.assay = "SCT", 
                  expression.slot = "data", 
                  group.by = "bint.bmem_sct1.algo3.res01_clusters", 
                  ymax = 2)
cp & scale_fill_manual(values = col_pal)
print(cp)
dev.off()





#-------------------------------------------------------------------------------
# Coverage Plot - Spleen
#-------------------------------------------------------------------------------

#----- Setting working directory -----
setwd("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Figures/BothBintBmem_PLNoPLNnSpleen_MedianAge_SubsetPLN8AAB.8T1D.8HD/")



#----- Loading data -----
seurat_obj = readRDS("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/BothBintBmem_Spleen_MedianAge_Subset24Donors_WithBintBmemSCT1Algo3Res01Clusters.rds")

# Subsetting Seurat object
Idents(seurat_obj) = "bint.bmem_sct1.algo3.res01_clusters"
seurat_obj_subset = subset(seurat_obj, idents = c("0", "3", "6", "7"))

# Defining variables
col_pal = c("#6b8e23", "#008837", "#8bbe1b", "#9dc209")



#----- Preparing data for coverage plot -----
DefaultAssay(seurat_obj_subset) = "peaks"
Idents(seurat_obj_subset) = "bint.bmem_sct1.algo3.res01_clusters"
seurat_obj_subset = RegionStats(seurat_obj_subset, 
                                genome = BSgenome.Hsapiens.UCSC.hg38, 
                                assay = "peaks")


#----- 1. BACH2 -----
seurat_obj_subset1 = LinkPeaks(seurat_obj_subset, 
                               peak.assay = "peaks", 
                               peak.slot = "counts", 
                               expression.assay = "SCT", 
                               expression.slot = "data", 
                               genes.use = "BACH2", 
                               distance = 1000000, 
                               score_cutoff = 0.03)

pdf("CoveragePlot_BintBmem_NormSCTAlgo3Res0.1.C0.C3.C6.C7_BACH2_S.pdf", width = 30, height = 10)
cp = CoveragePlot(seurat_obj_subset1, 
                  region = "BACH2", 
                  annotation = TRUE, 
                  peaks = TRUE, 
                  assay = "peaks", 
                  links = TRUE, 
                  extend.upstream = 400000, 
                  extend.downstream = 400000, 
                  expression.assay = "SCT", 
                  group.by = "bint.bmem_sct1.algo3.res01_clusters", 
                  ymin = 0, 
                  ymax = 120)
cp & scale_fill_manual(values = col_pal)
print(cp)
dev.off()


#----- 2. NFKB1 -----
seurat_obj_subset1 = LinkPeaks(seurat_obj_subset, 
                               peak.assay = "peaks", 
                               peak.slot = "counts", 
                               expression.assay = "SCT", 
                               expression.slot = "data", 
                               genes.use = "NFKB1", 
                               distance = 1000000, 
                               score_cutoff = 0.03)

pdf("CoveragePlot_BintBmem_NormSCTAlgo3Res0.1.C0.C3.C6.C7_NFKB1_S.pdf", width = 30, height = 10)
cp = CoveragePlot(seurat_obj_subset1, 
                  region = "NFKB1", 
                  annotation = TRUE, 
                  peaks = TRUE, 
                  assay = "peaks", 
                  links = TRUE, 
                  extend.upstream = 200000, 
                  extend.downstream = 200000, 
                  expression.assay = "SCT", 
                  group.by = "bint.bmem_sct1.algo3.res01_clusters", 
                  ymin = 0, 
                  ymax = 120)
cp & scale_fill_manual(values = col_pal)
print(cp)
dev.off()


#----- 3. BCL2 -----
seurat_obj_subset1 = LinkPeaks(seurat_obj_subset, 
                               peak.assay = "peaks", 
                               peak.slot = "counts", 
                               expression.assay = "SCT", 
                               expression.slot = "data", 
                               genes.use = "BCL2", 
                               distance = 1000000, 
                               score_cutoff = 0.02)

pdf("CoveragePlot_BintBmem_NormSCTAlgo3Res0.1.C0.C3.C6.C7_BCL2_S.pdf", width = 30, height = 10)
cp = CoveragePlot(seurat_obj_subset1, 
                  region = "BCL2", 
                  annotation = TRUE, 
                  peaks = TRUE, 
                  assay = "peaks", 
                  links = TRUE, 
                  extend.upstream = 20000, 
                  extend.downstream = 20000, 
                  expression.assay = "SCT", 
                  expression.slot = "data", 
                  group.by = "bint.bmem_sct1.algo3.res01_clusters", 
                  ymin = 0, 
                  ymax = 120)
cp & scale_fill_manual(values = col_pal)
print(cp)
dev.off()


#----- 4. EBF1 -----
seurat_obj_subset1 = LinkPeaks(seurat_obj_subset, 
                               peak.assay = "peaks", 
                               peak.slot = "counts", 
                               expression.assay = "SCT", 
                               expression.slot = "data", 
                               genes.use = "EBF1", 
                               distance = 1000000, 
                               score_cutoff = 0.02)

pdf("CoveragePlot_BintBmem_NormSCTAlgo3Res0.1.C0.C3.C6.C7_EBF1_S.pdf", width = 30, height = 10)
cp = CoveragePlot(seurat_obj_subset1, 
                  region = "EBF1", 
                  annotation = TRUE, 
                  peaks = TRUE, 
                  assay = "peaks", 
                  links = TRUE, 
                  extend.upstream = 20000, 
                  extend.downstream = 20000, 
                  expression.assay = "SCT", 
                  expression.slot = "data", 
                  group.by = "bint.bmem_sct1.algo3.res01_clusters", 
                  ymin = 0, 
                  ymax = 100)
cp & scale_fill_manual(values = col_pal)
print(cp)
dev.off()









