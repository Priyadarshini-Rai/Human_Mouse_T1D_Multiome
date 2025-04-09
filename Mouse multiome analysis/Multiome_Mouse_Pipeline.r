#-------------------------------------------------------------------------------
# MULTIOME DATA FILTERING AND PREPROCESSING PIPELINE
#-------------------------------------------------------------------------------




#-------------------------------------------------------------------------------
# Loading Required Libraries
#-------------------------------------------------------------------------------

library(Signac) # for ATAC-seq analysis
library(Seurat) # for multimodal analysis
library(EnsDb.Mmusculus.v79) # for mouse genome annotation
library(BSgenome.Mmusculus.UCSC.mm10) # mouse genome
library(SeuratDisk) # for storing metadata related to Seurat assay
library(stringr) # for string processing
library(BiocParallel) # for parallel processing
library(scDblFinder) # for doublet detection
library(rhdf5) # to write .h5 file




#-------------------------------------------------------------------------------
# Setting Data Directory And Working Directory
#-------------------------------------------------------------------------------

# Path to the data directory
Data_Dir = "/mnt/alvand/VahediLab/GenomicsData/Mouse/T_lymphocytes/Multiome/CD45/Merge/Multiome_aggregates_9_and_19_weeks_NOD_Ea16_09_13_24/outs/"

# # Setting working directory to save the .RData file
# setwd("/...../...../")




#-------------------------------------------------------------------------------
# Reading RNA And ATAC Data From Aggregate Analysis Using Cellranger aggr
#-------------------------------------------------------------------------------

print("Loading data")
ptm = proc.time()


# Loading RNA-seq and ATAC-seq data
Data = Read10X_h5(file.path(Data_Dir, "filtered_feature_bc_matrix.h5"))

# Path to fragments file required by Seurat
ATAC.Fragments = file.path(Data_Dir, "atac_fragments.tsv.gz")

# .CSV file to create metadata information
Info = read.csv(file.path(Data_Dir, "aggr.csv"), stringsAsFactors = F)

Annotation <- readRDS("/mnt/alvand/VahediLab/Databases/annotations/annotation_mm10.rds")
seqlevelsStyle(Annotation) <- "UCSC"
blacklist <- blacklist_mm10 

print("Loading data -- DONE")
print("Time taken:")
proc.time() - ptm




#-------------------------------------------------------------------------------
# Aggregating Multiome Data Into A Seurat Object
#-------------------------------------------------------------------------------

print("Creating Seurat object")
ptm = proc.time()


# Creating Seurat object
Seurat_Obj = CreateSeuratObject(
  counts = Data$`Gene Expression`,
  assay = "RNA",
  min.cells = 3
)


# Creating chromatin assay and adding it to Seurat object
Seurat_Obj[["ATAC"]] = CreateChromatinAssay(
  counts = Data$Peaks,
  sep = c(":", "-"),
  fragments = ATAC.Fragments,
  annotation = Annotation
)


print("Creating Seurat object -- DONE")
print("Time taken:")
proc.time() - ptm

#-------------------------------------------------------------------------------
# Assigning Sample IDs and Other Info. to the Cells
#-------------------------------------------------------------------------------
CellInfo = as.data.frame(Seurat_Obj@assays$RNA@counts@Dimnames[[2]])
colnames(CellInfo) = "barcodes"
rownames(CellInfo) = CellInfo$barcodes

CellInfo$libcodes = as.numeric(gsub(pattern = ".+-", replacement = "", CellInfo$barcodes))

CellInfo$samples = as.vector(Info$library_id[CellInfo$libcodes])

#-------------------------------------------------------------------------------
# Adding Metadata Information
#-------------------------------------------------------------------------------

Seurat_Obj$barcodes = CellInfo$barcodes
Seurat_Obj$libcodes = CellInfo$libcodes
Seurat_Obj$samples = CellInfo$samples


#-------------------------------------------------------------------------------
# Doublet Removal
#-------------------------------------------------------------------------------

print("Identifying doublets")
ptm = proc.time()


# Converting Seurat object into SingleCellExperiment object
SCE_Seurat_Obj = as.SingleCellExperiment(Seurat_Obj, assay = "RNA")

# Setting seed for reproducibility
set.seed(1234)

# Identifying doublets
#SCE_Seurat_Obj = scDblFinder(SCE_Seurat_Obj, clusters = FALSE, BPPARAM = MulticoreParam(10))
SCE_Seurat_Obj = scDblFinder(SCE_Seurat_Obj, clusters = FALSE, samples = "samples", BPPARAM = MulticoreParam(10))

# Saving R workspace
save.image("RWorkSpace_Till_DoubletRemoval.RData")


print("Identifying doublets -- DONE")
print("Time taken:")
proc.time() - ptm




#-------------------------------------------------------------------------------
# Quality Control
#-------------------------------------------------------------------------------

print("Filtering bad quality cells")
ptm = proc.time()



# By deafult, different functions will work on RNA assay
DefaultAssay(Seurat_Obj) = "RNA"

# # Identifying mitochondrial genes in the RNA assay
# Computing mitochondrial percentage for each cell type
#grep( "mt-", rownames(Seurat_Obj), value = T)
Seurat_Obj[["percent.mt"]] = PercentageFeatureSet(Seurat_Obj, pattern = "^mt-")


# By deafult, different functions will work on ATAC assay
DefaultAssay(Seurat_Obj) = "ATAC"

# Ratio of mononucleosomal (147-294 bp) to nucleosome-free (< 147 bp) fragments
# sequenced for the cell. Provides a value corresponding to each sample.
Seurat_Obj = NucleosomeSignal(Seurat_Obj)

# Ratio of mean number of Tn5 insertion events centered on TSS sites (+/- 500 bp
# of TSS) to mean number of Tn5 insertion events at TSS flanking regions (+/-
# 900 to 1000 bp of TSS). Provides a value corresponding to each sample.
Seurat_Obj = TSSEnrichment(Seurat_Obj)


# Violin plot depicting range of different parameters prior to quality control
VlnPlot(
  object = Seurat_Obj,
  features = c("nFeature_RNA", "nCount_RNA", "nFeature_ATAC", "nCount_ATAC", "percent.mt", "nucleosome_signal", "TSS.enrichment"),
  ncol = 4,
  pt.size = 0
)



# If a sample has a good coverage (>= minCov), then don't set a lower thresold for nCount, it's already pretty good.
Min_Cov = 1000

if(min(Seurat_Obj$nCount_RNA) >= Min_Cov)
{
  Count_Low = min(Seurat_Obj$nCount_RNA)
}else
{
  Count_Low = quantile(Seurat_Obj$nCount_RNA, prob = c(0.01))
}

Count_High = quantile(Seurat_Obj$nCount_RNA, prob = 0.99)

Feature_Low = quantile(Seurat_Obj$nFeature_RNA, prob = 0.01)

# Assigning doublet information to cell metadata
Seurat_Obj$scDblFinder.class = SCE_Seurat_Obj$scDblFinder.class

Seurat_Obj_Filtered = subset(Seurat_Obj,
                             subset = nFeature_RNA > Feature_Low &
                               nCount_RNA > Count_Low  & nCount_RNA < Count_High &
                               nucleosome_signal < 2 &
                               2 < TSS.enrichment &
                               percent.mt < 25 &
                               scDblFinder.class == "singlet")


# Violin plot depicting range of different parameters after quality control
VlnPlot(
  object = Seurat_Obj_Filtered,
  features = c("nFeature_RNA", "nCount_RNA", "nFeature_ATAC", "nCount_ATAC", "percent.mt", "nucleosome_signal", "TSS.enrichment"),
  ncol = 4,
  pt.size = 0
)



# Saving R workspace
save.image("RWorkSpace_Till_QualityControl.RData")



print("Filtering bad quality cells -- DONE")
print("Time taken:")
proc.time() - ptm




#-------------------------------------------------------------------------------
# Processing Gene Expression Data
#-------------------------------------------------------------------------------

print("Processing RNA-seq data")
ptm = proc.time()


# Normalizing gene expression data using SCTransform, and reducing its dimensionality using PCA

DefaultAssay(Seurat_Obj_Filtered) = "RNA"

Seurat_Obj_Filtered = SCTransform(Seurat_Obj_Filtered)
Seurat_Obj_Filtered = ScaleData(Seurat_Obj_Filtered)
Seurat_Obj_Filtered = FindVariableFeatures(Seurat_Obj_Filtered)
Seurat_Obj_Filtered = RunPCA(Seurat_Obj_Filtered, features = VariableFeatures(Seurat_Obj_Filtered))


# Saving R workspace
save.image("RWorkSpace_Till_RNA-Preprocessing.RData")


print("Processing RNA-seq data -- DONE")
print("Time taken:")
proc.time() - ptm



#-------------------------------------------------------------------------------
# Peak Calling
#-------------------------------------------------------------------------------

print("Applying MACS2")
ptm = proc.time()


# Calling peaks using MACS2

DefaultAssay(Seurat_Obj_Filtered) = "ATAC"

Peaks = CallPeaks(Seurat_Obj_Filtered,
                  macs2.path = "/mnt/alvand/apps/anaconda2/bin/macs2"
)

save.image("RWorkSpace_Till_Peaks.RData")


# Removing peaks on nonstandard chromosomes and in genomic blacklist regions
Peaks = keepStandardChromosomes(Peaks, pruning.mode = "coarse")
Peaks = subsetByOverlaps(x = Peaks, ranges = blacklist_mm10, invert = TRUE)

# Quantifying counts in each peak
MACS2_Counts = FeatureMatrix(
  fragments = Fragments(Seurat_Obj_Filtered),
  features = Peaks
)

save.image("RWorkSpace_Till_MACS2.RData")


# Creating a new assay using the MACS2 peak set and adding it to the Seurat object
Seurat_Obj_Filtered[["peaks"]] = CreateChromatinAssay(
  counts = MACS2_Counts,
  fragments = ATAC.Fragments,
  annotation = Annotation
)


# Saving R workspace
save.image("RWorkSpace_Till_PeakCalling.RData")


print("Applying MACS2 -- DONE")
print("Time taken:")
proc.time() - ptm




#-------------------------------------------------------------------------------
# Processing DNA Accessibility Data
#-------------------------------------------------------------------------------

print("Processing ATAC-seq data")
ptm = proc.time()


# Normalizing scATAC-seq dataset using latent semantic indexing (LSI), and reducing its dimensionality using SVD

DefaultAssay(Seurat_Obj_Filtered) = "peaks"

Seurat_Obj_Filtered = FindTopFeatures(Seurat_Obj_Filtered, min.cutoff = 5)
Seurat_Obj_Filtered = RunTFIDF(Seurat_Obj_Filtered)
Seurat_Obj_Filtered = RunSVD(Seurat_Obj_Filtered)


# Saving R workspace
save.image("RWorkSpace_Till_ATAC-Preprocessing.RData")


print("Processing ATAC-seq data -- DONE")
print("Time taken:")
proc.time() - ptm




#-------------------------------------------------------------------------------
# Multimodal Clustering using Seurat
#-------------------------------------------------------------------------------

print("Performing multimodal clustering")
ptm = proc.time()


DefaultAssay(Seurat_Obj_Filtered) = "SCT"

ElbowPlot(Seurat_Obj_Filtered, reduction = "pca")

ElbowPlot(Seurat_Obj_Filtered, reduction = "lsi")

# Seurat uses both assays to identify nearest neighbors
Seurat_Obj_Filtered = FindMultiModalNeighbors(Seurat_Obj_Filtered, 
                                              reduction.list = list("pca","lsi"), 
                                              dims.list = list(1:15, 2:10), 
                                              modality.weight.name = "RNA.weight")


# Applying UMAP for two-dimensional visualisation
Seurat_Obj_Filtered = RunUMAP(Seurat_Obj_Filtered, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

# Using nearest neighbor graph obtained using FindMultiModalNeighbors for identifying clusters
Seurat_Obj_Filtered = FindClusters(Seurat_Obj_Filtered, graph.name = "wsnn", algorithm = 3, resolution = c(0.4, 0.6, 0.8), verbose = FALSE)


# Saving workspace image
save.image("RWorkSpace_Till_MultimodalClustering_Algo3.RData")


print("Performing multimodal clustering -- DONE")
print("Time taken:")
proc.time() - ptm

#Data Visualization
DimPlot(Seurat_Obj_Filtered, repel = FALSE, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(Seurat_Obj_Filtered, repel = FALSE, reduction = "wnn.umap", group.by = "samples", label = FALSE)
DimPlot(Seurat_Obj_Filtered, repel = FALSE, reduction = "wnn.umap", group.by = "wsnn_res.0.4", label = TRUE)
DimPlot(Seurat_Obj_Filtered, repel = FALSE, reduction = "wnn.umap", group.by = "wsnn_res.0.6", label = TRUE)
DimPlot(Seurat_Obj_Filtered, repel = FALSE, reduction = "wnn.umap", group.by = "wsnn_res.0.8", label = TRUE)

Col = c("#ff8c00", "#ff0000", "#fa8072", "#00ffff", "#1e90ff", "#ff69b4", "#cc3366", "#fc0fc0", "#00bfff", "#9f00c5", "#9457eb", "#e0b0ff", "#0000ff", "#2f4f4f", "#4169e1", "#20b2aa", "#9acd32", "#ffd700", "#77b5fe", "#4a79b6", "#eb7d4a", "#e2643e")
DimPlot(Seurat_Obj_Filtered, repel = FALSE, reduction = "wnn.umap", label = TRUE, cols = Col, group.by = "seurat_clusters")

#Manual Cell Annotation
DefaultAssay(Seurat_Obj_Filtered)="SCT"
VlnPlot(Seurat_Obj_Filtered, features = c("Cd3d","Cd4","Cd8a", "Tcf7","Ccr7", "Sell", "Ebf1", "Cd19", "Ms4a1", "Cd14", "Fcgr3", "Ncam1", "Foxp3", "Ctla4", "Il2ra"), group.by = "seurat_clusters")


Idents(Seurat_Obj_Filtered)="seurat_clusters"
levels(Seurat_Obj_Filtered)
new.cluster.ids <- c("CD8 Naive", "B naive", "CD4 Naive", "Treg", "CD8 Naive", "B naive", "CD4 Naive", "Unknown", "CD4 Naive", "CD4 Naive", "Unknown", "B naive","B naive" ,"CD4 Naive", "CD4 TCM","CD4 TEM", "CD8 TCM", "Unknown", "CD8 TEM", "B memory", "cDC1", "CD4 Naive")
levels(Seurat_Obj_Filtered)
names(new.cluster.ids) <- levels(Seurat_Obj_Filtered)
Seurat_Obj_Filtered <- RenameIdents(Seurat_Obj_Filtered, new.cluster.ids)

Seurat_Obj_Filtered$manual_annot=Seurat_Obj_Filtered@active.ident
Final_CellTypes_Col_P = c("B intermediate" = "#ff8c00", "B memory" = "#ff0000",
                          "B naive" = "#fa8072", "CD14 Mono" = "#00ffff",
                          "CD16 Mono" = "#1e90ff", "CD4 Naive" = "#cc3366",
                          "CD4 TCM" = "#fc0fc0", "CD4 TEM" = "#00bfff",
                          "CD8 Naive" = "#9f00c5", "CD8 TCM" = "#9457eb",
                          "CD8 TEM" = "#e0b0ff", "gdT" = "#0000ff",
                          "ILC" = "#2f4f4f", "MAIT" = "#4169e1",
                          "NK" = "#20b2aa", "NK_CD56bright" = "#9acd32",
                          "Plasmablast" = "#ffd700", "Treg" = "#77b5fe",
                          "Unknown" = "#f2f3f4", "CD4 CTL" = "#ff69b4", 
                          "Platelet" = "#fff000", 
                          "CD4 Proliferating" = "#fc0bb0", 
                          "CD8 Proliferating" = "#d0a0ff", 
                          "NK Proliferating" = "#66cdaa", 
                          "ASDC" = "#00ff00", "cDC1" = "#adff2f", 
                          "cDC2" = "#7fff00", "pDC" = "#ccff00", 
                          "dnT" = "#BBBBBB", "Eryth" = "#CCCCCC", 
                          "HSPC" = "#DDDDDD")
						  
DimPlot(Seurat_Obj_Filtered, repel = FALSE, reduction = "wnn.umap", label = FALSE, group.by = "manual_annot", cols=Final_CellTypes_Col_P)