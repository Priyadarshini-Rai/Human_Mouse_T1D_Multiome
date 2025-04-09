#-------------------------------------------------------------------------------
# MULTIOME DATA FILTERING AND PREPROCESSING PIPELINE
#-------------------------------------------------------------------------------




#-------------------------------------------------------------------------------
# Loading Required Libraries
#-------------------------------------------------------------------------------

library(Signac) # for ATAC-seq analysis
library(Seurat) # for multimodal analysis
library(EnsDb.Hsapiens.v86) # for human genome annotation
library(BSgenome.Hsapiens.UCSC.hg38) # human genome
library(SeuratDisk) # for storing metadata related to Seurat assay
library(stringr) # for string processing
library(BiocParallel) # for parallel processing
library(scDblFinder) # for doublet detection
library(rhdf5) # to write .h5 file




#-------------------------------------------------------------------------------
# Setting Data Directory And Working Directory
#-------------------------------------------------------------------------------

# Path to the data directory
Data_Dir = "/mnt/alvand/VahediLab/GenomicsData/Human/CD45/Multiome/Merge/Multiome_aggregates_04_14_23_Spleen/outs/"

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

# Loading donors information
Our_Donors = read.csv("/mnt/alvand/VahediLab/Projects/Multiome_HPAP/Shared/Donor_information_Agg_04042023.CSV")

# Reference data for cell annotation
Reference_Data = LoadH5Seurat("/mnt/alvand/VahediLab/Databases/annotations/pbmc_multimodal.h5seurat")
#We’ll use an annotated PBMC reference dataset from Hao et al. (2020), available for download here: https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat


print("Loading data -- DONE")
print("Time taken:")
proc.time() - ptm




#-------------------------------------------------------------------------------
# Assigning Sample ID And Other Information To Each Cell
#-------------------------------------------------------------------------------

print("Creating cell information dataframe")
ptm = proc.time()


CellInfo = as.data.frame(Data$`Gene Expression`@Dimnames[[2]])
colnames(CellInfo) = "barcodes"
rownames(CellInfo) = CellInfo$barcodes

CellInfo$libcodes = as.numeric(gsub(pattern = ".+-", replacement = "", CellInfo$barcodes))

CellInfo$sample = as.vector(Info$library_id[CellInfo$libcodes])

Donor_ID = str_extract(CellInfo$sample, "[^_]+")
CellInfo$donor_id = as.vector(Donor_ID)
Unique_ID = unique(str_extract(CellInfo$sample, "[^_]+"))

CellInfo$gender = rep(0, ncol(Data$`Gene Expression`))

CellInfo$age = rep(0, ncol(Data$`Gene Expression`))

CellInfo$status = rep(0, ncol(Data$`Gene Expression`))

Our_Donors$donor_ID = gsub("-", "", Our_Donors$donor_ID)

for(i in 1:length(Unique_ID))
{
  Index1 = which(CellInfo$donor_id == Unique_ID[i])
  Index2 = as.numeric(match(Unique_ID[i], Our_Donors$donor_ID))
  CellInfo[Index1, 'gender'] = Our_Donors[Index2, 'gender']
  CellInfo[Index1, 'age'] = Our_Donors[Index2, 'age_years']
  CellInfo[Index1, 'status'] = Our_Donors[Index2, 'Status']
}


print("Creating cell information dataframe -- DONE")
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

# Getting gene annotations for hg38
Annotation = readRDS("/mnt/alvand/VahediLab/Databases/annotations/annotation_hg38.rds")

# Setting naming convention to UCSC.
seqlevelsStyle(Annotation) = "UCSC"

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
# Adding Metadata Information
#-------------------------------------------------------------------------------

print("Adding cell information to the Seurat object")
ptm = proc.time()


Seurat_Obj$barcodes = CellInfo$barcodes
Seurat_Obj$libcodes = CellInfo$libcodes
Seurat_Obj$sample = CellInfo$sample
Seurat_Obj$donor_ids = CellInfo$donor_id
Seurat_Obj$gender = CellInfo$gender
Seurat_Obj$age = CellInfo$age
Seurat_Obj$status = CellInfo$status


print("Adding cell information to the Seurat Object -- DONE")
print("Time taken:")
proc.time() - ptm




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
SCE_Seurat_Obj = scDblFinder(SCE_Seurat_Obj, clusters = FALSE, samples = "sample", BPPARAM = MulticoreParam(10))

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
Seurat_Obj[["percent.mt"]] = PercentageFeatureSet(Seurat_Obj, pattern = "^MT-")


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
save.image("RWorkSpace_Till_RNA-PreprocessingRData")


print("Processing RNA-seq data -- DONE")
print("Time taken:")
proc.time() - ptm




#-------------------------------------------------------------------------------
# Annotating Cell Types using Seurat
#-------------------------------------------------------------------------------

print("Performing cell annotation")
ptm = proc.time()


DefaultAssay(Seurat_Obj_Filtered) = "SCT"

# Transfering cell type labels from reference to query
Transfer_Anchors = FindTransferAnchors(
  reference = Reference_Data,
  query = Seurat_Obj_Filtered,
  normalization.method = "SCT",
  reference.reduction = "spca",
  recompute.residuals = FALSE,
  dims = 1:50
)

Predictions = TransferData(
  anchorset = Transfer_Anchors,
  refdata = Reference_Data$celltype.l2,
  weight.reduction = Seurat_Obj_Filtered[['pca']],
  dims = 1:50
)

Seurat_Obj_Filtered = AddMetaData(
  object = Seurat_Obj_Filtered,
  metadata = Predictions
)

# Setting the cell identities to the cell type predictions
Idents(Seurat_Obj_Filtered) = "predicted.id"


# Saving R workspace
save.image("RWorkSpace_Till_CellAnnotation.RData")


print("Performing cell annotation -- DONE")
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
                  macs2.path = "/mnt/alvand/apps/anaconda2/bin/macs2",
                  group.by = "predicted.id"
)

# Removing peaks on nonstandard chromosomes and in genomic blacklist regions
Peaks = keepStandardChromosomes(Peaks, pruning.mode = "coarse")
Peaks = subsetByOverlaps(x = Peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# Quantifying counts in each peak
MACS2_Counts = FeatureMatrix(
  fragments = Fragments(Seurat_Obj_Filtered),
  features = Peaks,
  cells = colnames(Seurat_Obj_Filtered)
)

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

# Seurat uses both assays to identify nearest neighbors
Seurat_Obj_Filtered = FindMultiModalNeighbors(Seurat_Obj_Filtered, 
                                              reduction.list = list("pca","lsi"), 
                                              dims.list = list(1:50, 2:40), 
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
