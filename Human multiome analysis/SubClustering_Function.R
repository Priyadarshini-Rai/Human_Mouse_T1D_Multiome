#--------------------------------------------------------------------------------------------------
# FUNCTION TO APPLY CLUSTERING, TOGETHER ON TWO CELL TYPES OF COMBINED PLN OLD, PLN NEW, AND SPLEEN
#--------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
# Loading Required Libraries
#-------------------------------------------------------------------------------

# Setting path for libraries
# .libPaths("/opt/R/4.1.1/lib/R/library")
.libPaths("/home/priya/R/x86_64-pc-linux-gnu-library/4.1/")

library(Seurat) # for multimodal analysis
library(Signac) # for ATAC-seq analysis
# library(EnsDb.Hsapiens.v86) # for human genome annotation
library(BSgenome.Hsapiens.UCSC.hg38) # human genome
library(SeuratDisk) # for storing metadata related to Seurat assay
library(stringr) # for string processing
library(BiocParallel) # for parallel processing
library(scDblFinder) # for doublet detection
# library(pheatmap) # for mapping clusters to cell annotations
library(dplyr) # for stacked barplot
library(tidyverse) # for frequency of sample/cluster
library(MAST) # for DE analysis
library(DESeq2) # for DE analysis
library(ggplot2) # for bar chart
# library(ggforce) # for function geom_sina()
# library(plotly) # for plotting
library(TFBSTools) # for motif analysis
library(JASPAR2020) # for motif analysis
library(chromVAR) # for motif analysis
library(rhdf5) # to write .h5 file
# library(ggpubr) # for box plot with p-value
library(data.table) # for function tstrsplit()
# library(RColorBrewer) # for color palette
# library(ComplexUpset) # for UpSet plot
library(bedr) # for peak distribution 
# library(ChIPseeker) # for peak distribution 
# library(TxDb.Hsapiens.UCSC.hg19.knownGene) # for peak distribution 
library(clusterProfiler) # for peak distribution 
library(annotables) # for peak distribution 
library(org.Hs.eg.db) # for peak distribution 
library(Matrix) # for function writemm()





#-------------------------------------------------------------------------------
# Setting Parameter(s) Globally
#-------------------------------------------------------------------------------

# Following parameter allows to plot all the labels on the two-dimensional plot
options(ggrepel.max.overlaps = Inf)

options(future.globals.maxSize = 8000 * 1024^2)





#-------------------------------------------------------------------------------
# Function To Combine PLN Old And PLN New For A Particular Cell Type
#-------------------------------------------------------------------------------

combine_tissues <- function(cell_type)
{
  
  # Loading data
  s1_o = readRDS("/mnt/alvand/priya/PLN_Aggregate/R_WorkSpace/PLN_Seurat_Obj_AddedGeneticRiskScore_04.15.2024_V15.rds") # PLNo
  tissue1 = "PLNo"
  
  s2_o = readRDS("/mnt/alvand/priya/PLNn/R_Files/PLNn_AddedAgeGenderGenericDonorIDHbA1c_09.13.2024_V1_F.rds")
  tissue2 = "PLNn"
  
  s3_o = readRDS("/mnt/alvand/priya/Spleen_Aggregate/R_WorkSpace/Spleen_Seurat_Obj_AddedGeneticRiskScore_06.07.2024_V10.rds")
  tissue3 = "S"
  
  
  
  # Subsetting Seurat objects for following donors so that median age for T1D, 
  # AAb+, and HD is comparable
  
  PLNo_donors = c("T1D1", "T1D2", "T1D4", "T1D6", "T1D10", "T1D12", 
                  "AAB+1", "AAB+2", "AAB+3", "AAB+4", "AAB+5", "AAB+6", "AAB+7", 
                  "AAB+8", 
                  "HD1", "HD3", "HD6", "HD10", "HD11", "HD13", "HD14")
  Idents(s1_o) = "donor_ids_umap"
  s1_o = subset(s1_o, idents = PLNo_donors)
  
  PLNn_donors = c("T1D15", "T1D17", 
                  "HD16")
  Idents(s2_o) = "donor_ids_umap"
  s2_o = subset(s2_o, idents = PLNn_donors)
  
  S_donors = c("T1D1", "T1D2", "T1D6", "T1D10", "T1D12", 
               "AAB+1", "AAB+2", "AAB+3", "AAB+4", "AAB+5", "AAB+6", "AAB+7", 
               "AAB+8", 
               "HD3", "HD6", "HD10", "HD13", "HD14")
  Idents(s3_o) = "donor_ids_umap"
  s3_o = subset(s3_o, idents = S_donors)
  
  
  
  # Subsetting Seurat objects for required cell types
  Idents(s1_o) = "predicted.id.v1"
  s1 = subset(s1_o, idents = cell_type)
  
  Idents(s2_o) = "predicted.id.v1"
  s2 = subset(s2_o, idents = cell_type)
  
  Idents(s3_o) = "predicted.id.v1"
  s3 = subset(s3_o, idents = cell_type)
  
  
  
  # Extracting "data" from SCT assay
  m_s1 = as.matrix(s1@assays[["SCT"]]@counts)
  m_s2 = as.matrix(s2@assays[["SCT"]]@counts)
  m_s3 = as.matrix(s3@assays[["SCT"]]@counts)
  
  # Identifying common features between tissue types
  m_s1_features = toupper(rownames(m_s1))
  m_s2_features = toupper(rownames(m_s2))
  m_s3_features = toupper(rownames(m_s3))
  common_features = intersect(intersect(m_s1_features, m_s2_features), m_s3_features)
  ind_s1 = match(common_features, m_s1_features)
  ind_s2 = match(common_features, m_s2_features)
  ind_s3 = match(common_features, m_s3_features)
  
  # Keeping only common features
  m_s1 = m_s1[ind_s1, ]
  colnames(m_s1) = paste0(tissue1, "_", colnames(m_s1))
  m_s2 = m_s2[ind_s2, ]
  colnames(m_s2) = paste0(tissue2, "_", colnames(m_s2))
  m_s3 = m_s3[ind_s3, ]
  colnames(m_s3) = paste0(tissue3, "_", colnames(m_s3))
  combined_sct = cbind(m_s1, m_s2, m_s3)
  
  
  
  # Extracting "data" from peaks assay
  m_s1 = as.matrix(s1@assays[["peaks"]]@counts)
  chrom.assay_s1 = CreateChromatinAssay(counts = m_s1, sep = c(":", "-"))
  chrom.assay_s1_seurat = CreateSeuratObject(chrom.assay_s1, assay = "peaks")
  chrom.assay_s1_seurat$dataset = tissue1
  
  m_s2 = as.matrix(s2@assays[["peaks"]]@counts)
  chrom.assay_s2 = CreateChromatinAssay(counts = m_s2, sep = c(":", "-"))
  chrom.assay_s2_seurat = CreateSeuratObject(chrom.assay_s2, assay = "peaks")
  chrom.assay_s2_seurat$dataset = tissue2
  
  m_s3 = as.matrix(s3@assays[["peaks"]]@counts)
  chrom.assay_s3 = CreateChromatinAssay(counts = m_s3, sep = c(":", "-"))
  chrom.assay_s3_seurat = CreateSeuratObject(chrom.assay_s3, assay = "peaks")
  chrom.assay_s3_seurat$dataset = tissue3
  
  
  combined = merge(x = chrom.assay_s1_seurat, 
                   y = list(chrom.assay_s2_seurat, chrom.assay_s3_seurat) , 
                   add.cell.ids = c(tissue1, tissue2, tissue3))
  combined_peaks = as.matrix(combined@assays[["peaks"]]@counts)
  
  
  
  # Generating metadata
  metadata_s1 = data.frame(predicted.id.v1 = s1$predicted.id.v1,
                           status = s1$status,
                           donor_ids = s1$donor_ids,
                           sample = s1$sample,
                           age = s1$age,
                           gender = s1$gender,
                           barcodes = s1$barcodes,
                           libcodes = s1$libcodes,
                           donor_ids_umap = s1$donor_ids_umap,
                           donor_ids_celldistribution.acrosssample = s1$donor_ids_celldistribution.acrosssample, 
                           tissue = rep(tissue1, ncol(m_s1)))
  rownames(metadata_s1) = paste0(tissue1, "_", rownames(metadata_s1))
  
  metadata_s2 = data.frame(predicted.id.v1 = s2$predicted.id.v1,
                           status = s2$status,
                           donor_ids = s2$donor_ids,
                           sample = s2$orig.ident,
                           age = s2$age,
                           gender = s2$gender,
                           barcodes = s2$barcodes,
                           libcodes = s2$libcodes,
                           donor_ids_umap = s2$donor_ids_umap,
                           donor_ids_celldistribution.acrosssample = s2$donor_ids_celldistribution.acrosssample, 
                           tissue = rep(tissue2, ncol(m_s2)))
  rownames(metadata_s2) = paste0(tissue2, "_", rownames(metadata_s2))
  
  metadata_s3 = data.frame(predicted.id.v1 = s3$predicted.id.v1,
                           status = s3$status,
                           donor_ids = s3$donor_ids,
                           sample = s3$sample,
                           age = s3$age,
                           gender = s3$gender,
                           barcodes = s3$barcodes,
                           libcodes = s3$libcodes,
                           donor_ids_umap = s3$donor_ids_umap,
                           donor_ids_celldistribution.acrosssample = s3$donor_ids_celldistribution.acrosssample, 
                           tissue = rep(tissue3, ncol(m_s3)))
  rownames(metadata_s3) = paste0(tissue3, "_", rownames(metadata_s3))
  
  metadata_combined = rbind(metadata_s1, metadata_s2, metadata_s3)
  
  
  
  # Combining tissues to create a single Seurat object
  seurat_combined = CreateSeuratObject(counts = combined_sct,
                                       meta.data = metadata_combined,
                                       names.delim = "-",
                                       assay = "SCT")
  
  seurat_combined[["peaks"]] = CreateChromatinAssay(counts = combined_peaks, sep = c(":", "-"))
  
  
  
  saveRDS(seurat_combined, file = paste0("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".rds"))
  
  
  
  # Returning Seurat object
  return(seurat_combined)
  
}





#-------------------------------------------------------------------------------
# Function To Perform Clustering On A Particular Cell Type
#-------------------------------------------------------------------------------

celltype.clustering <- function(cell_type, tissue, perform_rna_norm, fig_dir, resolution_rna, perform_atac_norm, resolution_atac, resolution_rna_atac) {
  
  cell.type_clust.info_gene = list()
  cell.type_clust.info_peaks = list()
  cell.type_clust.info_multimodal = list()
  
  # Combining PLN old and PLN new for a particular cell type
  seurat_obj_subset = combine_tissues(cell_type)
  
  
  
  
  #-----------------------
  # SCT - No Normalization
  #-----------------------
  
  if(toupper(perform_rna_norm) == "NO")
  {
    DefaultAssay(seurat_obj_subset) = "SCT"
    
    
    seurat_obj_subset = FindVariableFeatures(seurat_obj_subset)
    
    seurat_obj_subset = RunPCA(seurat_obj_subset, 
                               features = VariableFeatures(object = seurat_obj_subset), 
                               reduction.name = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".pca"), 
                               reduction.key = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), "pca_"))
    
    png(paste0(fig_dir, "ElbowPlotPCA_SCT_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".png"), 
        width = 600, 
        height = 400)
    elbow_plot = ElbowPlot(seurat_obj_subset,
                           ndims = 50,
                           reduction = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".pca"))
    print(elbow_plot)
    dev.off()
    
    variance = seurat_obj_subset[[paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".pca")]]@stdev^2
    cum_var = ((cumsum(variance)) / (sum(variance)))
    max_pc = as.numeric(which(cum_var >= 0.79)[1]) # retaining number of components based on a variance threshold (e.g., 90%)
    
    
    seurat_obj_subset = FindNeighbors(seurat_obj_subset, 
                                      dims = 1:max_pc, 
                                      reduction = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".pca"))
    
    seurat_obj_subset = FindClusters(seurat_obj_subset, 
                                     algorithm = 3, 
                                     resolution = c(0.1, 0.2, 0.3, 0.4, 0.5))
    for(r in seq(from = 0.1, to = 0.5, by = 0.1))
    {
      seurat_obj_subset@meta.data[[paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), "_SCT_snn_res.", r)]] = seurat_obj_subset@meta.data[[paste0("SCT_snn_res.", r)]]
    }
    
    seurat_obj_subset = RunUMAP(seurat_obj_subset, 
                                reduction = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".pca"), 
                                dims = 1:max_pc, 
                                reduction.name = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".rna.umap"), 
                                reduction.key = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), "rnaUMAP_"))
    
    for(res in seq(from = 0.1, to = 0.5, by = 0.1))
    {
      png(paste0(fig_dir, "UMAP_SCT_PC", max_pc, "_Algo3Res", res, "_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".png"), 
          width = 600, 
          height = 400)
      umap_sct = DimPlot(seurat_obj_subset, 
                         repel = TRUE, 
                         reduction = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".rna.umap"), 
                         label = TRUE, 
                         group.by = paste0("SCT_snn_res.", res))
      print(umap_sct)
      dev.off()
    }
    
    active_ident = paste0("SCT_snn_res.", resolution_rna)
    
    
    
    # Stacked barplot depicting donor status composition for each cluster
    tmp1 = as.data.frame(table(seurat_obj_subset@meta.data[[active_ident]], seurat_obj_subset$status))
    colnames(tmp1) = c("Cluster", "Status", "Freq")
    
    # Converting cell type distribution per sample into percentage format
    tmp2 = as.data.frame(tmp1 %>% pivot_wider(names_from = Status, values_from = Freq))
    rownames(tmp2) = tmp2[, 1]
    tmp2 = tmp2[, -1]
    
    tmp3 = (tmp2/rowSums(tmp2)) * 100
    
    tmp4 = cbind.data.frame(rownames(tmp3), tmp3)
    colnames(tmp4) = c("Cluster", colnames(tmp3))
    
    tmp5 = reshape2::melt(tmp4, id = "Cluster")
    colnames(tmp5) = c("Cluster", "Status", "Percent")
    
    
    cols = c("T1DM" = "#d7191c", "AAB+" = "#fdae61", "HD" = "#2c7bb6")
    png(paste0(fig_dir, "StackedBarPlot_DonorStatusCompositionPerClust_SCT_PC", max_pc, "_Algo3Res", resolution_rna, "_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".png"), 
        width = 600, 
        height = 400)
    stacked_bp = ggplot(data = tmp5, aes(x = Cluster, y = Percent, fill = Status)) + 
      geom_bar(stat = "identity") + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      xlab("Cluster") + 
      ylab("Percent") + 
      scale_fill_manual(values = cols)
    print(stacked_bp)
    dev.off()
    
    
    # Lollipop plot depicting number of cells per cluster
    tmp6 = as.data.frame(table(seurat_obj_subset@meta.data[[active_ident]]))
    colnames(tmp6) = c("Cluster", "Freq")
    png(paste0(fig_dir, "LollipopPlot_ClustFreq_SCT_PC", max_pc, "_Algo3Res", resolution_rna, "_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".png"), 
        width = 600, 
        height = 400)
    lollipop_plot = ggplot(tmp6, aes(x = Cluster, y = Freq)) +
      geom_segment(aes(x = Cluster, xend = Cluster, y = 0, yend = Freq), color = "gray", lwd = 1.5) +
      geom_point(size = 11.5, pch = 21, bg = 4, col = 1) +
      geom_text(aes(label = Freq), color = "white", size = 3) +
      coord_flip() +
      theme_minimal()
    print(lollipop_plot)
    dev.off()
    
    
    
    tmp4$T1D_AAb = tmp4$T1DM + tmp4$`AAB+`
    colnames(tmp4)[as.numeric(which(colnames(tmp4) == "AAB+"))] = "AAb"
    clust1 = as.numeric(tmp4$Cluster[which(tmp4$T1D_AAb == max(tmp4$T1D_AAb))])
    
    
    df1 = tmp2[, c("T1DM", "AAB+", "HD")]
    colnames(df1)[as.numeric(which(colnames(df1) == "AAB+"))] = "AAb"
    
    
    df2 = data.frame(table(seurat_obj_subset$status))
    colnames(df2) = c("Status", "Freq")
    rownames(df2) = df2$Status
    rownames(df2)[as.numeric(which(rownames(df2) == "AAB+"))] = "AAb"
    df2 = df2[c("T1DM", "AAb", "HD"), ]
    
    df3 = df1
    df3$T1DM = df3$T1DM/df2$Freq[which(df2$Status == "T1DM")]
    df3$AAb = df3$AAb/df2$Freq[which(df2$Status == "AAB+")]
    df3$HD = df3$HD/df2$Freq[which(df2$Status == "HD")]
    
    df4 = (df3/rowSums(df3)) * 100
    df4$Cluster = rownames(df4)
    
    clust2 = c()
    for(j in 1:nrow(df4))
    {
      if(df4$T1DM[j] > 70 | df4$AAb[j] > 70 | (df4$T1DM[j] + df4$AAb[j]) > 70)
      {
        clust2 = c(clust2, df4$Cluster[j])
      }
    }
    
    clust = union(clust1, as.numeric(clust2))
    
    
    cluster_info = data.frame(matrix(nrow = 0, ncol = 3))
    colnames(cluster_info) = c("DonorID", "Freq", "Cluster")
    
    for(k in 1:length(clust))
    {
      Idents(seurat_obj_subset) = active_ident
      seurat_obj_subset_clust = subset(seurat_obj_subset, idents = clust[k])
      tmp7 = data.frame(table(seurat_obj_subset_clust$donor_ids_umap))
      colnames(tmp7) = c("DonorID", "Freq")
      tmp7$Cluster = clust[k]
      tmp7 = tmp7[order(tmp7$Freq, decreasing = TRUE), ]
      cluster_info = rbind(cluster_info, tmp7)
    }
    
    cluster_info$BestClust = rep(clust1, nrow(cluster_info))
    cluster_info$MaxPC = rep(max_pc, nrow(cluster_info))
    cell.type_clust.info_gene[[paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))]] = cluster_info
    
  }
  
  
  
  
  #-----------------------------------------------------------------
  # SCT - Performing SCT normalization using "counts" of data subset
  #-----------------------------------------------------------------
  
  if(toupper(perform_rna_norm) == "YES")
  {
    DefaultAssay(seurat_obj_subset) = "SCT"
    
    
    seurat_obj_subset = SCTransform(seurat_obj_subset, assay = "SCT", new.assay.name = "SCT1")
    
    all_features = rownames(seurat_obj_subset)
    seurat_obj_subset = ScaleData(seurat_obj_subset, features = all_features)
    
    seurat_obj_subset = FindVariableFeatures(seurat_obj_subset)
    
    seurat_obj_subset = RunPCA(seurat_obj_subset, 
                               features = VariableFeatures(object = seurat_obj_subset), 
                               reduction.name = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".pca"), 
                               reduction.key = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), "pca_"))
    
    png(paste0(fig_dir, "ElbowPlotPCA_NormSCT_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".png"), 
        width = 600, 
        height = 400)
    elbow_plot = ElbowPlot(seurat_obj_subset,
                           ndims = 50,
                           reduction = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".pca"))
    print(elbow_plot)
    dev.off()
    
    variance = seurat_obj_subset[[paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".pca")]]@stdev^2
    cum_var = ((cumsum(variance)) / (sum(variance)))
    max_pc = as.numeric(which(cum_var >= 0.79)[1]) # retaining number of components based on a variance threshold (e.g., 90%)
    
    
    seurat_obj_subset = FindNeighbors(seurat_obj_subset, 
                                      dims = 1:max_pc, 
                                      reduction = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".pca"))
    
    seurat_obj_subset = FindClusters(seurat_obj_subset, 
                                     algorithm = 3, 
                                     resolution = c(0.1, 0.2, 0.3, 0.4, 0.5))
    for(r in seq(from = 0.1, to = 0.5, by = 0.1))
    {
      seurat_obj_subset@meta.data[[paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), "_SCT1_snn_res.", r)]] = seurat_obj_subset@meta.data[[paste0("SCT1_snn_res.", r)]]
    }
    
    seurat_obj_subset = RunUMAP(seurat_obj_subset, 
                                reduction = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".pca"), 
                                dims = 1:max_pc, 
                                reduction.name = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".rna.umap"), 
                                reduction.key = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), "rnaUMAP_"))
    
    for(res in seq(from = 0.1, to = 0.5, by = 0.1))
    {
      png(paste0(fig_dir, "UMAP_NormSCT_PC", max_pc, "_Algo3Res", res, "_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".png"),
          width = 600,
          height = 400)
      umap_sct = DimPlot(seurat_obj_subset,
                         repel = TRUE,
                         reduction = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".rna.umap"),
                         label = TRUE,
                         group.by = paste0("SCT1_snn_res.", res), 
                         raster = FALSE)
      print(umap_sct)
      dev.off()
    }
    
    active_ident = paste0("SCT1_snn_res.", resolution_rna)
    
    
    
    # Stacked barplot depicting donor status composition for each cluster
    tmp1 = as.data.frame(table(seurat_obj_subset@meta.data[[active_ident]], seurat_obj_subset$status))
    colnames(tmp1) = c("Cluster", "Status", "Freq")
    
    # Converting cell type distribution per sample into percentage format
    tmp2 = as.data.frame(tmp1 %>% pivot_wider(names_from = Status, values_from = Freq))
    rownames(tmp2) = tmp2[, 1]
    tmp2 = tmp2[, -1]
    
    tmp3 = (tmp2/rowSums(tmp2)) * 100
    
    tmp4 = cbind.data.frame(rownames(tmp3), tmp3)
    colnames(tmp4) = c("Cluster", colnames(tmp3))
    
    tmp5 = reshape2::melt(tmp4, id = "Cluster")
    colnames(tmp5) = c("Cluster", "Status", "Percent")
    
    
    cols = c("T1DM" = "#d7191c", "AAB+" = "#fdae61", "HD" = "#2c7bb6")
    png(paste0(fig_dir, "StackedBarPlot_DonorStatusCompositionPerClust_NormSCT_PC", max_pc, "_Algo3Res", resolution_rna, "_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".png"), 
        width = 600, 
        height = 400)
    stacked_bp = ggplot(data = tmp5, aes(x = Cluster, y = Percent, fill = Status)) + 
      geom_bar(stat = "identity") + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      xlab("Cluster") + 
      ylab("Percent") + 
      scale_fill_manual(values = cols)
    print(stacked_bp)
    dev.off()
    
    
    # Lollipop plot depicting number of cells per cluster
    tmp6 = as.data.frame(table(seurat_obj_subset@meta.data[[active_ident]]))
    colnames(tmp6) = c("Cluster", "Freq")
    png(paste0(fig_dir, "LollipopPlot_ClustFreq_NormSCT_PC", max_pc, "_Algo3Res", resolution_rna, "_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".png"), 
        width = 600, 
        height = 400)
    lollipop_plot = ggplot(tmp6, aes(x = Cluster, y = Freq)) +
      geom_segment(aes(x = Cluster, xend = Cluster, y = 0, yend = Freq), color = "gray", lwd = 1.5) +
      geom_point(size = 11.5, pch = 21, bg = 4, col = 1) +
      geom_text(aes(label = Freq), color = "white", size = 3) +
      coord_flip() +
      theme_minimal()
    print(lollipop_plot)
    dev.off()
    
    
    
    tmp4$T1D_AAb = tmp4$T1DM + tmp4$`AAB+`
    colnames(tmp4)[as.numeric(which(colnames(tmp4) == "AAB+"))] = "AAb"
    clust1 = as.numeric(tmp4$Cluster[which(tmp4$T1D_AAb == max(tmp4$T1D_AAb))])
    
    
    df1 = tmp2[, c("T1DM", "AAB+", "HD")]
    colnames(df1)[as.numeric(which(colnames(df1) == "AAB+"))] = "AAb"
    
    
    df2 = data.frame(table(seurat_obj_subset$status))
    colnames(df2) = c("Status", "Freq")
    rownames(df2) = df2$Status
    rownames(df2)[as.numeric(which(rownames(df2) == "AAB+"))] = "AAb"
    df2 = df2[c("T1DM", "AAb", "HD"), ]
    
    df3 = df1
    df3$T1DM = df3$T1DM/df2$Freq[which(df2$Status == "T1DM")]
    df3$AAb = df3$AAb/df2$Freq[which(df2$Status == "AAB+")]
    df3$HD = df3$HD/df2$Freq[which(df2$Status == "HD")]
    
    df4 = (df3/rowSums(df3)) * 100
    df4$Cluster = rownames(df4)
    
    clust2 = c()
    for(j in 1:nrow(df4))
    {
      if(df4$T1DM[j] > 70 | df4$AAb[j] > 70 | (df4$T1DM[j] + df4$AAb[j]) > 70)
      {
        clust2 = c(clust2, df4$Cluster[j])
      }
    }
    
    clust = union(clust1, as.numeric(clust2))
    
    
    cluster_info = data.frame(matrix(nrow = 0, ncol = 3))
    colnames(cluster_info) = c("DonorID", "Freq", "Cluster")
    
    for(k in 1:length(clust))
    {
      Idents(seurat_obj_subset) = active_ident
      seurat_obj_subset_clust = subset(seurat_obj_subset, idents = clust[k])
      tmp7 = data.frame(table(seurat_obj_subset_clust$donor_ids_umap))
      colnames(tmp7) = c("DonorID", "Freq")
      tmp7$Cluster = clust[k]
      tmp7 = tmp7[order(tmp7$Freq, decreasing = TRUE), ]
      cluster_info = rbind(cluster_info, tmp7)
    }
    
    cluster_info$BestClust = rep(clust1, nrow(cluster_info))
    cluster_info$MaxPC = rep(max_pc, nrow(cluster_info))
    cell.type_clust.info_gene[[paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))]] = cluster_info
    
  }
  
  
  
  
  print("SCT -- DONE")
  
  
  
  
  #-----------------
  # Peaks - No TFIDF
  #-----------------
  
  if(toupper(perform_atac_norm) == "NO")
  {
    
    DefaultAssay(seurat_obj_subset) = "peaks"
    
    
    seurat_obj_subset = FindTopFeatures(object = seurat_obj_subset, min.cutoff = 5)
    
    seurat_obj_subset = RunSVD(object = seurat_obj_subset, 
                               reduction.name = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".lsi"), 
                               reduction.key = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), "lsi_"))
    
    png(paste0(fig_dir, "ElbowPlotLSI_Peaks_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".png"), 
        width = 600, 
        height = 400)
    elbow_plot = ElbowPlot(seurat_obj_subset,
                           ndims = 50,
                           reduction = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".lsi"))
    print(elbow_plot)
    dev.off()
    
    # variance = seurat_obj_subset[[paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".lsi")]]@stdev^2
    # cum_var = ((cumsum(variance)) / (sum(variance)))
    # max_lsi = as.numeric(which(cum_var >= 0.79)[1]) # retaining number of components based on a variance threshold (e.g., 90%)
    max_lsi = 15
    
    seurat_obj_subset = FindNeighbors(seurat_obj_subset, 
                                      reduction = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".lsi"), 
                                      dims = 2:max_lsi)
    
    seurat_obj_subset = FindClusters(seurat_obj_subset, algorithm = 3, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5))
    for(r in seq(from = 0.1, to = 0.5, by = 0.1))
    {
      seurat_obj_subset@meta.data[[paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), "_peaks_snn_res.", r)]] = seurat_obj_subset@meta.data[[paste0("peaks_snn_res.", r)]]
    }
    
    seurat_obj_subset = RunUMAP(seurat_obj_subset, 
                                reduction = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".lsi"), 
                                dims = 2:max_lsi, 
                                reduction.name = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".atac.umap"), 
                                reduction.key = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), "atacUMAP_"))
    
    for(res in seq(from = 0.1, to = 0.5, by = 0.1))
    {
      png(paste0(fig_dir, "UMAP_Peaks_LSI", max_lsi, "_Algo3Res", res, "_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".png"), 
          width = 600, 
          height = 400)
      umap_peaks = DimPlot(seurat_obj_subset, 
                           repel = TRUE, 
                           reduction = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".atac.umap"), 
                           label = TRUE, 
                           group.by = paste0("peaks_snn_res.", res))
      print(umap_peaks)
      dev.off()
    }
    
    active_ident = paste0("peaks_snn_res.", resolution_atac)
    
    
    
    # Stacked barplot depicting donor status composition for each cluster
    tmp1 = as.data.frame(table(seurat_obj_subset@meta.data[[active_ident]], seurat_obj_subset$status))
    colnames(tmp1) = c("Cluster", "Status", "Freq")
    
    # Converting cell type distribution per sample into percentage format
    tmp2 = as.data.frame(tmp1 %>% pivot_wider(names_from = Status, values_from = Freq))
    rownames(tmp2) = tmp2[, 1]
    tmp2 = tmp2[, -1]
    
    tmp3 = (tmp2/rowSums(tmp2)) * 100
    
    tmp4 = cbind.data.frame(rownames(tmp3), tmp3)
    colnames(tmp4) = c("Cluster", colnames(tmp3))
    
    tmp5 = reshape2::melt(tmp4, id = "Cluster")
    colnames(tmp5) = c("Cluster", "Status", "Percent")
    
    
    cols = c("T1DM" = "#d7191c", "AAB+" = "#fdae61", "HD" = "#2c7bb6")
    png(paste0(fig_dir, "StackedBarPlot_DonorStatusCompositionPerClust_Peaks_LSI", max_lsi, "_Algo3Res", resolution_atac, "_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".png"), 
        width = 600, 
        height = 400)
    stacked_bp = ggplot(data = tmp5, aes(x = Cluster, y = Percent, fill = Status)) + 
      geom_bar(stat = "identity") + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      xlab("Cluster") + 
      ylab("Percent") + 
      scale_fill_manual(values = cols)
    print(stacked_bp)
    dev.off()
    
    
    # Lollipop plot depicting number of cells per cluster
    tmp6 = as.data.frame(table(seurat_obj_subset@meta.data[[active_ident]]))
    colnames(tmp6) = c("Cluster", "Freq")
    png(paste0(fig_dir, "LollipopPlot_ClustFreq_Peaks_LSI", max_lsi, "_Algo3Res", resolution_atac, "_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".png"), 
        width = 600, 
        height = 400)
    lollipop_plot = ggplot(tmp6, aes(x = Cluster, y = Freq)) +
      geom_segment(aes(x = Cluster, xend = Cluster, y = 0, yend = Freq), color = "gray", lwd = 1.5) +
      geom_point(size = 11.5, pch = 21, bg = 4, col = 1) +
      geom_text(aes(label = Freq), color = "white", size = 3) +
      coord_flip() +
      theme_minimal()
    print(lollipop_plot)
    dev.off()
    
    
    
    tmp4$T1D_AAb = tmp4$T1DM + tmp4$`AAB+`
    colnames(tmp4)[as.numeric(which(colnames(tmp4) == "AAB+"))] = "AAb"
    clust1 = as.numeric(tmp4$Cluster[which(tmp4$T1D_AAb == max(tmp4$T1D_AAb))])
    
    
    df1 = tmp2[, c("T1DM", "AAB+", "HD")]
    colnames(df1)[as.numeric(which(colnames(df1) == "AAB+"))] = "AAb"
    
    
    df2 = data.frame(table(seurat_obj_subset$status))
    colnames(df2) = c("Status", "Freq")
    rownames(df2) = df2$Status
    rownames(df2)[as.numeric(which(rownames(df2) == "AAB+"))] = "AAb"
    df2 = df2[c("T1DM", "AAb", "HD"), ]
    
    df3 = df1
    df3$T1DM = df3$T1DM/df2$Freq[which(df2$Status == "T1DM")]
    df3$AAb = df3$AAb/df2$Freq[which(df2$Status == "AAB+")]
    df3$HD = df3$HD/df2$Freq[which(df2$Status == "HD")]
    
    df4 = (df3/rowSums(df3)) * 100
    df4$Cluster = rownames(df4)
    
    clust2 = c()
    for(j in 1:nrow(df4))
    {
      if(df4$T1DM[j] > 70 | df4$AAb[j] > 70 | (df4$T1DM[j] + df4$AAb[j]) > 70)
      {
        clust2 = c(clust2, df4$Cluster[j])
      }
    }
    
    clust = union(clust1, as.numeric(clust2))
    
    
    cluster_info = data.frame(matrix(nrow = 0, ncol = 3))
    colnames(cluster_info) = c("DonorID", "Freq", "Cluster")
    
    for(k in 1:length(clust))
    {
      Idents(seurat_obj_subset) = active_ident
      seurat_obj_subset_clust = subset(seurat_obj_subset, idents = clust[k])
      tmp7 = data.frame(table(seurat_obj_subset_clust$donor_ids_umap))
      colnames(tmp7) = c("DonorID", "Freq")
      tmp7$Cluster = clust[k]
      tmp7 = tmp7[order(tmp7$Freq, decreasing = TRUE), ]
      cluster_info = rbind(cluster_info, tmp7)
    }
    
    cluster_info$BestClust = rep(clust1, nrow(cluster_info))
    cluster_info$MaxLSI = rep(max_lsi, nrow(cluster_info))
    cell.type_clust.info_peaks[[paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))]] = cluster_info
    
  }
  
  
  
  
  #----------------------------------------------
  # Peaks - Performing TFIDF again on data subset
  #----------------------------------------------
  
  if(toupper(perform_atac_norm) == "YES")
  {
    
    DefaultAssay(seurat_obj_subset) = "peaks"
    
    
    seurat_obj_subset = RunTFIDF(object = seurat_obj_subset)
    
    seurat_obj_subset = FindTopFeatures(object = seurat_obj_subset, min.cutoff = 5)
    
    seurat_obj_subset = RunSVD(object = seurat_obj_subset, 
                               reduction.name = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".lsi"), 
                               reduction.key = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), "lsi_"))
    
    png(paste0(fig_dir, "ElbowPlotLSI_NormPeaks_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".png"), 
        width = 600, 
        height = 400)
    elbow_plot = ElbowPlot(seurat_obj_subset,
                           ndims = 50,
                           reduction = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".lsi"))
    print(elbow_plot)
    dev.off()
    
    # variance = seurat_obj_subset[[paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".lsi")]]@stdev^2
    # cum_var = ((cumsum(variance)) / (sum(variance)))
    # max_lsi = as.numeric(which(cum_var >= 0.79)[1]) # retaining number of components based on a variance threshold (e.g., 90%)
    max_lsi = 15
    
    seurat_obj_subset = FindNeighbors(seurat_obj_subset, 
                                      reduction = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".lsi"), 
                                      dims = 2:max_lsi)
    
    seurat_obj_subset = FindClusters(seurat_obj_subset, algorithm = 3, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5))
    for(r in seq(from = 0.1, to = 0.5, by = 0.1))
    {
      seurat_obj_subset@meta.data[[paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), "_peaks_snn_res.", r)]] = seurat_obj_subset@meta.data[[paste0("peaks_snn_res.", r)]]
    }
    
    seurat_obj_subset = RunUMAP(seurat_obj_subset, 
                                reduction = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".lsi"), 
                                dims = 2:max_lsi, 
                                reduction.name = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".atac.umap"), 
                                reduction.key = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), "atacUMAP_"))
    
    for(res in seq(from = 0.1, to = 0.5, by = 0.1))
    {
      png(paste0(fig_dir, "UMAP_NormPeaks_LSI", max_lsi, "_Algo3Res", res, "_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".png"),
          width = 600,
          height = 400)
      umap_peaks = DimPlot(seurat_obj_subset,
                           repel = TRUE,
                           reduction = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".atac.umap"),
                           label = TRUE,
                           group.by = paste0("peaks_snn_res.", res), 
                           raster = FALSE)
      print(umap_peaks)
      dev.off()
    }
    
    active_ident = paste0("peaks_snn_res.", resolution_atac)
    
    
    
    # Stacked barplot depicting donor status composition for each cluster
    tmp1 = as.data.frame(table(seurat_obj_subset@meta.data[[active_ident]], seurat_obj_subset$status))
    colnames(tmp1) = c("Cluster", "Status", "Freq")
    
    # Converting cell type distribution per sample into percentage format
    tmp2 = as.data.frame(tmp1 %>% pivot_wider(names_from = Status, values_from = Freq))
    rownames(tmp2) = tmp2[, 1]
    tmp2 = tmp2[, -1]
    
    tmp3 = (tmp2/rowSums(tmp2)) * 100
    
    tmp4 = cbind.data.frame(rownames(tmp3), tmp3)
    colnames(tmp4) = c("Cluster", colnames(tmp3))
    
    tmp5 = reshape2::melt(tmp4, id = "Cluster")
    colnames(tmp5) = c("Cluster", "Status", "Percent")
    
    
    cols = c("T1DM" = "#d7191c", "AAB+" = "#fdae61", "HD" = "#2c7bb6")
    png(paste0(fig_dir, "StackedBarPlot_DonorStatusCompositionPerClust_NormPeaks_LSI", max_lsi, "_Algo3Res", resolution_atac, "_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".png"), 
        width = 600, 
        height = 400)
    stacked_bp = ggplot(data = tmp5, aes(x = Cluster, y = Percent, fill = Status)) + 
      geom_bar(stat = "identity") + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      xlab("Cluster") + 
      ylab("Percent") + 
      scale_fill_manual(values = cols)
    print(stacked_bp)
    dev.off()
    
    
    # Lollipop plot depicting number of cells per cluster
    tmp6 = as.data.frame(table(seurat_obj_subset@meta.data[[active_ident]]))
    colnames(tmp6) = c("Cluster", "Freq")
    png(paste0(fig_dir, "LollipopPlot_ClustFreq_NormPeaks_LSI", max_lsi, "_Algo3Res", resolution_atac, "_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".png"), 
        width = 600, 
        height = 400)
    lollipop_plot = ggplot(tmp6, aes(x = Cluster, y = Freq)) +
      geom_segment(aes(x = Cluster, xend = Cluster, y = 0, yend = Freq), color = "gray", lwd = 1.5) +
      geom_point(size = 11.5, pch = 21, bg = 4, col = 1) +
      geom_text(aes(label = Freq), color = "white", size = 3) +
      coord_flip() +
      theme_minimal()
    print(lollipop_plot)
    dev.off()
    
    
    
    tmp4$T1D_AAb = tmp4$T1DM + tmp4$`AAB+`
    colnames(tmp4)[as.numeric(which(colnames(tmp4) == "AAB+"))] = "AAb"
    clust1 = as.numeric(tmp4$Cluster[which(tmp4$T1D_AAb == max(tmp4$T1D_AAb))])
    
    
    df1 = tmp2[, c("T1DM", "AAB+", "HD")]
    colnames(df1)[as.numeric(which(colnames(df1) == "AAB+"))] = "AAb"
    
    
    df2 = data.frame(table(seurat_obj_subset$status))
    colnames(df2) = c("Status", "Freq")
    rownames(df2) = df2$Status
    rownames(df2)[as.numeric(which(rownames(df2) == "AAB+"))] = "AAb"
    df2 = df2[c("T1DM", "AAb", "HD"), ]
    
    df3 = df1
    df3$T1DM = df3$T1DM/df2$Freq[which(df2$Status == "T1DM")]
    df3$AAb = df3$AAb/df2$Freq[which(df2$Status == "AAB+")]
    df3$HD = df3$HD/df2$Freq[which(df2$Status == "HD")]
    
    df4 = (df3/rowSums(df3)) * 100
    df4$Cluster = rownames(df4)
    
    clust2 = c()
    for(j in 1:nrow(df4))
    {
      if(df4$T1DM[j] > 70 | df4$AAb[j] > 70 | (df4$T1DM[j] + df4$AAb[j]) > 70)
      {
        clust2 = c(clust2, df4$Cluster[j])
      }
    }
    
    clust = union(clust1, as.numeric(clust2))
    
    
    cluster_info = data.frame(matrix(nrow = 0, ncol = 3))
    colnames(cluster_info) = c("DonorID", "Freq", "Cluster")
    
    for(k in 1:length(clust))
    {
      Idents(seurat_obj_subset) = active_ident
      seurat_obj_subset_clust = subset(seurat_obj_subset, idents = clust[k])
      tmp7 = data.frame(table(seurat_obj_subset_clust$donor_ids_umap))
      colnames(tmp7) = c("DonorID", "Freq")
      tmp7$Cluster = clust[k]
      tmp7 = tmp7[order(tmp7$Freq, decreasing = TRUE), ]
      cluster_info = rbind(cluster_info, tmp7)
    }
    
    cluster_info$BestClust = rep(clust1, nrow(cluster_info))
    cluster_info$MaxLSI = rep(max_lsi, nrow(cluster_info))
    cell.type_clust.info_peaks[[paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))]] = cluster_info
    
  }
  
  
  
  
  print("Peaks -- DONE")
  
  
  
  
  #-----------------------------------
  # Multimodal - Without normalization
  #-----------------------------------
  
  if((toupper(perform_rna_norm) == "NO") & (toupper(perform_atac_norm) == "NO"))
  {
    
    seurat_obj_subset = FindMultiModalNeighbors(seurat_obj_subset, 
                                                reduction.list = list(paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".pca"), paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".lsi")), 
                                                dims.list = list(1:max_pc, 2:max_lsi), 
                                                modality.weight.name = c(paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".SCT.weight"), paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".peaks.weight")), 
                                                knn.graph.name = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".wknn"), 
                                                snn.graph.name = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".wsnn"), 
                                                weighted.nn.name = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".weighted.nn"))
    
    seurat_obj_subset = RunUMAP(seurat_obj_subset, 
                                nn.name = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".weighted.nn"), 
                                reduction.name = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".wnn.umap"), 
                                reduction.key = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), "wnnUMAP_"))
    
    seurat_obj_subset = FindClusters(seurat_obj_subset, 
                                     graph.name = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".wsnn"), 
                                     algorithm = 3, 
                                     resolution = c(0.1, 0.2, 0.3, 0.4, 0.5))
    
    for(res in seq(from = 0.1, to = 0.5, by = 0.1))
    {
      png(paste0(fig_dir, "UMAP_BothSCTPeaks_PC", max_pc, "_LSI", max_lsi, "_Algo3Res", res, "_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".png"), 
          width = 600, 
          height = 400)
      umap_sct_peaks = DimPlot(seurat_obj_subset, 
                               repel = TRUE, 
                               reduction = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".wnn.umap"), 
                               label = TRUE, 
                               group.by = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".wsnn_res.", res))
      print(umap_sct_peaks)
      dev.off()
    }
    
    active_ident = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".wsnn_res.", resolution_rna_atac)
    
    
    
    # Stacked barplot depicting donor status composition for each cluster
    tmp1 = as.data.frame(table(seurat_obj_subset@meta.data[[active_ident]], seurat_obj_subset$status))
    colnames(tmp1) = c("Cluster", "Status", "Freq")
    
    # Converting cell type distribution per sample into percentage format
    tmp2 = as.data.frame(tmp1 %>% pivot_wider(names_from = Status, values_from = Freq))
    rownames(tmp2) = tmp2[, 1]
    tmp2 = tmp2[, -1]
    
    tmp3 = (tmp2/rowSums(tmp2)) * 100
    
    tmp4 = cbind.data.frame(rownames(tmp3), tmp3)
    colnames(tmp4) = c("Cluster", colnames(tmp3))
    
    tmp5 = reshape2::melt(tmp4, id = "Cluster")
    colnames(tmp5) = c("Cluster", "Status", "Percent")
    
    
    cols = c("T1DM" = "#d7191c", "AAB+" = "#fdae61", "HD" = "#2c7bb6")
    png(paste0(fig_dir, "StackedBarPlot_DonorStatusCompositionPerClust_BothSCTPeaks_PC", max_pc, "_LSI", max_lsi, "_Algo3Res", resolution_rna_atac, "_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".png"), 
        width = 600, 
        height = 400)
    stacked_bp = ggplot(data = tmp5, aes(x = Cluster, y = Percent, fill = Status)) + 
      geom_bar(stat = "identity") + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      xlab("Cluster") + 
      ylab("Percent") + 
      scale_fill_manual(values = cols)
    print(stacked_bp)
    dev.off()
    
    
    # Lollipop plot depicting number of cells per cluster
    tmp6 = as.data.frame(table(seurat_obj_subset@meta.data[[active_ident]]))
    colnames(tmp6) = c("Cluster", "Freq")
    png(paste0(fig_dir, "LollipopPlot_ClustFreq_BothSCTPeaks_PC", max_pc, "_LSI", max_lsi, "_Algo3Res", resolution_rna_atac, "_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".png"), 
        width = 600, 
        height = 400)
    lollipop_plot = ggplot(tmp6, aes(x = Cluster, y = Freq)) +
      geom_segment(aes(x = Cluster, xend = Cluster, y = 0, yend = Freq), color = "gray", lwd = 1.5) +
      geom_point(size = 11.5, pch = 21, bg = 4, col = 1) +
      geom_text(aes(label = Freq), color = "white", size = 3) +
      coord_flip() +
      theme_minimal()
    print(lollipop_plot)
    dev.off()
    
    
    
    tmp4$T1D_AAb = tmp4$T1DM + tmp4$`AAB+`
    colnames(tmp4)[as.numeric(which(colnames(tmp4) == "AAB+"))] = "AAb"
    clust1 = as.numeric(tmp4$Cluster[which(tmp4$T1D_AAb == max(tmp4$T1D_AAb))])
    
    
    df1 = tmp2[, c("T1DM", "AAB+", "HD")]
    colnames(df1)[as.numeric(which(colnames(df1) == "AAB+"))] = "AAb"
    
    
    df2 = data.frame(table(seurat_obj_subset$status))
    colnames(df2) = c("Status", "Freq")
    rownames(df2) = df2$Status
    rownames(df2)[as.numeric(which(rownames(df2) == "AAB+"))] = "AAb"
    df2 = df2[c("T1DM", "AAb", "HD"), ]
    
    df3 = df1
    df3$T1DM = df3$T1DM/df2$Freq[which(df2$Status == "T1DM")]
    df3$AAb = df3$AAb/df2$Freq[which(df2$Status == "AAB+")]
    df3$HD = df3$HD/df2$Freq[which(df2$Status == "HD")]
    
    df4 = (df3/rowSums(df3)) * 100
    df4$Cluster = rownames(df4)
    
    clust2 = c()
    for(j in 1:nrow(df4))
    {
      if(df4$T1DM[j] > 70 | df4$AAb[j] > 70 | (df4$T1DM[j] + df4$AAb[j]) > 70)
      {
        clust2 = c(clust2, df4$Cluster[j])
      }
    }
    
    clust = union(clust1, as.numeric(clust2))
    
    
    cluster_info = data.frame(matrix(nrow = 0, ncol = 3))
    colnames(cluster_info) = c("DonorID", "Freq", "Cluster")
    
    for(k in 1:length(clust))
    {
      Idents(seurat_obj_subset) = active_ident
      seurat_obj_subset_clust = subset(seurat_obj_subset, idents = clust[k])
      tmp7 = data.frame(table(seurat_obj_subset_clust$donor_ids_umap))
      colnames(tmp7) = c("DonorID", "Freq")
      tmp7$Cluster = clust[k]
      tmp7 = tmp7[order(tmp7$Freq, decreasing = TRUE), ]
      cluster_info = rbind(cluster_info, tmp7)
    }
    
    cluster_info$BestClust = rep(clust1, nrow(cluster_info))
    cluster_info$MaxPC = rep(max_pc, nrow(cluster_info))
    cluster_info$MaxLSI = rep(max_lsi, nrow(cluster_info))
    cell.type_clust.info_multimodal[[paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))]] = cluster_info 
    
  }
  
  
  
  
  #--------------------------------
  # Multimodal - With normalization
  #--------------------------------
  
  if((toupper(perform_rna_norm) == "YES") & (toupper(perform_atac_norm) == "YES"))
  {
    
    seurat_obj_subset = FindMultiModalNeighbors(seurat_obj_subset, 
                                                reduction.list = list(paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".pca"), paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".lsi")), 
                                                dims.list = list(1:max_pc, 2:max_lsi), 
                                                modality.weight.name = c(paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".SCT.weight"), paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".peaks.weight")), 
                                                knn.graph.name = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".wknn"), 
                                                snn.graph.name = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".wsnn"), 
                                                weighted.nn.name = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".weighted.nn"))
    
    seurat_obj_subset = RunUMAP(seurat_obj_subset, 
                                nn.name = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".weighted.nn"), 
                                reduction.name = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".wnn.umap"), 
                                reduction.key = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), "wnnUMAP_"))
    
    seurat_obj_subset = FindClusters(seurat_obj_subset, 
                                     graph.name = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".wsnn"), 
                                     algorithm = 3, 
                                     resolution = c(0.1, 0.2, 0.3, 0.4, 0.5))
    
    for(res in seq(from = 0.1, to = 0.5, by = 0.1))
    {
      png(paste0(fig_dir, "UMAP_BothSCTPeaksNorm_PC", max_pc, "_LSI", max_lsi, "_Algo3Res", res, "_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".png"),
          width = 600,
          height = 400)
      umap_sct_peaks = DimPlot(seurat_obj_subset,
                               repel = TRUE,
                               reduction = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".wnn.umap"),
                               label = TRUE,
                               group.by = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".wsnn_res.", res), 
                               raster = FALSE)
      print(umap_sct_peaks)
      dev.off()
    }
    
    active_ident = paste0(tolower(paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))), ".wsnn_res.", resolution_rna_atac)
    
    
    
    # Stacked barplot depicting donor status composition for each cluster
    tmp1 = as.data.frame(table(seurat_obj_subset@meta.data[[active_ident]], seurat_obj_subset$status))
    colnames(tmp1) = c("Cluster", "Status", "Freq")
    
    # Converting cell type distribution per sample into percentage format
    tmp2 = as.data.frame(tmp1 %>% pivot_wider(names_from = Status, values_from = Freq))
    rownames(tmp2) = tmp2[, 1]
    tmp2 = tmp2[, -1]
    
    tmp3 = (tmp2/rowSums(tmp2)) * 100
    
    tmp4 = cbind.data.frame(rownames(tmp3), tmp3)
    colnames(tmp4) = c("Cluster", colnames(tmp3))
    
    tmp5 = reshape2::melt(tmp4, id = "Cluster")
    colnames(tmp5) = c("Cluster", "Status", "Percent")
    
    
    cols = c("T1DM" = "#d7191c", "AAB+" = "#fdae61", "HD" = "#2c7bb6")
    png(paste0(fig_dir, "StackedBarPlot_DonorStatusCompositionPerClust_BothSCTPeaksNorm_PC", max_pc, "_LSI", max_lsi, "_Algo3Res", resolution_rna_atac, "_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".png"), 
        width = 600, 
        height = 400)
    stacked_bp = ggplot(data = tmp5, aes(x = Cluster, y = Percent, fill = Status)) + 
      geom_bar(stat = "identity") + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      xlab("Cluster") + 
      ylab("Percent") + 
      scale_fill_manual(values = cols)
    print(stacked_bp)
    dev.off()
    
    
    # Lollipop plot depicting number of cells per cluster
    tmp6 = as.data.frame(table(seurat_obj_subset@meta.data[[active_ident]]))
    colnames(tmp6) = c("Cluster", "Freq")
    png(paste0(fig_dir, "LollipopPlot_ClustFreq_BothSCTPeaksNorm_PC", max_pc, "_LSI", max_lsi, "_Algo3Res", resolution_rna_atac, "_", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "_", tissue, ".png"), 
        width = 600, 
        height = 400)
    lollipop_plot = ggplot(tmp6, aes(x = Cluster, y = Freq)) +
      geom_segment(aes(x = Cluster, xend = Cluster, y = 0, yend = Freq), color = "gray", lwd = 1.5) +
      geom_point(size = 11.5, pch = 21, bg = 4, col = 1) +
      geom_text(aes(label = Freq), color = "white", size = 3) +
      coord_flip() +
      theme_minimal()
    print(lollipop_plot)
    dev.off()
    
    
    
    tmp4$T1D_AAb = tmp4$T1DM + tmp4$`AAB+`
    colnames(tmp4)[as.numeric(which(colnames(tmp4) == "AAB+"))] = "AAb"
    clust1 = as.numeric(tmp4$Cluster[which(tmp4$T1D_AAb == max(tmp4$T1D_AAb))])
    
    
    df1 = tmp2[, c("T1DM", "AAB+", "HD")]
    colnames(df1)[as.numeric(which(colnames(df1) == "AAB+"))] = "AAb"
    
    
    df2 = data.frame(table(seurat_obj_subset$status))
    colnames(df2) = c("Status", "Freq")
    rownames(df2) = df2$Status
    rownames(df2)[as.numeric(which(rownames(df2) == "AAB+"))] = "AAb"
    df2 = df2[c("T1DM", "AAb", "HD"), ]
    
    df3 = df1
    df3$T1DM = df3$T1DM/df2$Freq[which(df2$Status == "T1DM")]
    df3$AAb = df3$AAb/df2$Freq[which(df2$Status == "AAB+")]
    df3$HD = df3$HD/df2$Freq[which(df2$Status == "HD")]
    
    df4 = (df3/rowSums(df3)) * 100
    df4$Cluster = rownames(df4)
    
    clust2 = c()
    for(j in 1:nrow(df4))
    {
      if(df4$T1DM[j] > 70 | df4$AAb[j] > 70 | (df4$T1DM[j] + df4$AAb[j]) > 70)
      {
        clust2 = c(clust2, df4$Cluster[j])
      }
    }
    
    clust = union(clust1, as.numeric(clust2))
    
    
    cluster_info = data.frame(matrix(nrow = 0, ncol = 3))
    colnames(cluster_info) = c("DonorID", "Freq", "Cluster")
    
    for(k in 1:length(clust))
    {
      Idents(seurat_obj_subset) = active_ident
      seurat_obj_subset_clust = subset(seurat_obj_subset, idents = clust[k])
      tmp7 = data.frame(table(seurat_obj_subset_clust$donor_ids_umap))
      colnames(tmp7) = c("DonorID", "Freq")
      tmp7$Cluster = clust[k]
      tmp7 = tmp7[order(tmp7$Freq, decreasing = TRUE), ]
      cluster_info = rbind(cluster_info, tmp7)
    }
    
    cluster_info$BestClust = rep(paste(clust1, collapse = "_"), nrow(cluster_info))
    cluster_info$MaxPC = rep(max_pc, nrow(cluster_info))
    cluster_info$MaxLSI = rep(max_lsi, nrow(cluster_info))
    cell.type_clust.info_multimodal[[paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))]] = cluster_info
    
  }
  
  
  
  
  print("Multimodal -- DONE")
  
  
  
  
  # #---------------------
  # # Saving Seurat object
  # #---------------------
  # 
  # if((toupper(perform_rna_norm) == "NO") & (toupper(perform_atac_norm) == "NO"))
  # {
  #   saveRDS(seurat_obj_subset, file = paste0("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "Clustering_SCTPC", max_pc, "_PeaksLSI", max_lsi, "_BothSCTPeaks_", tissue, ".rds"))
  # }
  # 
  # if((toupper(perform_rna_norm) == "YES") & (toupper(perform_atac_norm) == "YES"))
  # {
  #   saveRDS(seurat_obj_subset, file = paste0("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "Clustering_NormSCTPC", max_pc, "_NormPeaksLSI", max_lsi, "_BothSCTPeaks_", tissue, ".rds"))
  # }
  
  
  
  
  #---------------------------
  # Saving cluster information
  #---------------------------
  if((toupper(perform_rna_norm) == "NO") & (toupper(perform_atac_norm) == "NO"))
  {
    
    cell_types = paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))
    
    openxlsx::write.xlsx(cell.type_clust.info_gene, file = paste0("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Data/ClusterInfo_SCT_", cell_types, "_", tissue, ".xlsx"))
    openxlsx::write.xlsx(cell.type_clust.info_peaks, file = paste0("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Data/ClusterInfo_Peaks_", cell_types, "_", tissue, ".xlsx"))
    openxlsx::write.xlsx(cell.type_clust.info_multimodal, file = paste0("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Data/ClusterInfo_BothSCTPeaks_", cell_types, "_", tissue, ".xlsx"))
    
  }
  
  if((toupper(perform_rna_norm) == "YES") & (toupper(perform_atac_norm) == "YES"))
  {
    
    cell_types = paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2]))
    
    openxlsx::write.xlsx(cell.type_clust.info_gene, file = paste0("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Data/ClusterInfo_NormSCT_", cell_types, "_", tissue, ".xlsx"))
    openxlsx::write.xlsx(cell.type_clust.info_peaks, file = paste0("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Data/ClusterInfo_NormPeaks_", cell_types, "_", tissue, ".xlsx"))
    openxlsx::write.xlsx(cell.type_clust.info_multimodal, file = paste0("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Data/ClusterInfo_BothSCTPeaksNorm_", cell_types, "_", tissue, ".xlsx"))
    
  }
  
  
  
  
  #---------------------
  # Saving Seurat object
  #---------------------
  
  if((toupper(perform_rna_norm) == "NO") & (toupper(perform_atac_norm) == "NO"))
  {
    saveRDS(seurat_obj_subset, file = paste0("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "Clustering_SCTPC", max_pc, "_PeaksLSI", max_lsi, "_BothSCTPeaks_", tissue, ".rds"))
  }
  
  if((toupper(perform_rna_norm) == "YES") & (toupper(perform_atac_norm) == "YES"))
  {
    saveRDS(seurat_obj_subset, file = paste0("/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Workspace/", paste0("Both", gsub(" ", "", cell_type[1]), gsub(" ", "", cell_type[2])), "Clustering_NormSCTPC", max_pc, "_NormPeaksLSI", max_lsi, "_BothSCTPeaks_", tissue, ".rds"))
  }
  
}





#-------------------------------------------------------------------------------
# Calling Function
#-------------------------------------------------------------------------------

# cell_type = c("CD4 TCM", "Treg")
cell_type = c("B intermediate", "B memory")
tissue = "PLNo_PLNn_Spleen_MedianAge_Subset24Donors"
perform_rna_norm = "Yes"
fig_dir = "/mnt/alvand/priya/MultiomePaper1_R1/CellTypeClustering_UsingFunction/Figures/PLNo_PLNn_Spleen_MedianAge_SubsetPLN8AAB.8T1D.8HD/"
resolution_rna = 0.1
perform_atac_norm = "Yes"
resolution_atac = 0.1
resolution_rna_atac = 0.1

t1 = celltype.clustering(cell_type, tissue, perform_rna_norm, fig_dir, resolution_rna, perform_atac_norm, resolution_atac, resolution_rna_atac)









