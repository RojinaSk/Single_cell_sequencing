######title:single cell sequencing of normal and fibroid lungs in mouse and human######
#install necessary packages and load them##using seurat V5
install.packages("Seurat")
#loading the library
library(Seurat)
library(tidyverse)
library(patchwork)


#read barcodes, features(genes) and matrix files

NML_1 <- Read10X(data.dir = "~/Desktop/Project_seurat_lungs/GSE132771_RAW/NML1/")
NML_2 <- Read10X(data.dir = "~/Desktop/Project_seurat_lungs/GSE132771_RAW/NML2/")
NML_3 <- Read10X(data.dir = "~/Desktop/Project_seurat_lungs/GSE132771_RAW/NML3/")


#create seurat objects

NML_1 <- CreateSeuratObject(counts = NML_1, project = "NML_1", min.cells = 3, min.features = 200)#cells expressing less than 200 genes is considered poor
NML_2 <- CreateSeuratObject(counts = NML_2, project = "NML_2", min.cells = 3, min.features = 200)
NML_3 <- CreateSeuratObject(counts = NML_3, project = "NML_3", min.cells = 3, min.features = 200)


##saving seurat object##

saveRDS(NML_1, file = "~/Desktop/Project_seurat_lungs/NML_1.RDS")
saveRDS(NML_2, file = "~/Desktop/Project_seurat_lungs/NML_2.RDS")
saveRDS(NML_3, file = "~/Desktop/Project_seurat_lungs/NML_3.RDS")


###merging all three seurat object#####

merged_NML <- merge(NML_1, y= c(NML_2, NML_3),
                    add.cell.ids =ls()[1:3],
                    project="Merged_NMLs")

#View(merged_NML@meta.data)
#saving the merged seurat object
saveRDS(merged_NML, file = "~/Desktop/Project_seurat_lungs/Merged_NML.RDS")

####standard preprocessing workflow########

#storing mitochindrial percentage in object meta data
merged_NML <- PercentageFeatureSet(merged_NML, pattern = "^MT-", 
                                   col.name = "percent.mt")
#saving the merged seurat object with percent
saveRDS(merged_NML, file = "~/Desktop/Project_seurat_lungs/Merged_NML_per.RDS")


####selecting cells for further analysis######

#visulaising qc metrics as a violin plot

VlnPlot(merged_NML, features = c( "nCount_RNA", "nFeature_RNA","percent.mt"), 
        ncol = 3)

#filtering out 
merged_NML_subset <- subset(merged_NML, subset = nFeature_RNA > 500 
                            & nFeature_RNA < 4000  & 
                              nCount_RNA < 20000 & percent.mt < 10)

#visualising qc metrics after filtering
VlnPlot(merged_NML_subset, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), 
        ncol=3)

#save seurat object
saveRDS(merged_NML_subset, file = "~/Desktop/Project_seurat_lungs/merged_NMLs_subset.RDS")

######DAtA normalization,feature selection and data scaling######
#data normalization

merged_NML_normalized <- NormalizeData(merged_NML_subset, assay = "RNA",
                                       normalization.method = "LogNormalize",
                                       scale.factor= 10000)

##identification of highly variable features(feature selection)
merged_NML_normalized <- FindVariableFeatures(merged_NML_normalized, 
                                              selection.method ="vst",
                                              nfeatures = 2000)

#creating a feature plot

VariableFeaturePlot(merged_NML_normalized)


###data scaling###
merged_NML_normalized <- ScaleData(merged_NML_normalized) #for the 2000 genes

all.genes <- row.names(merged_NML_normalized)

merged_NML_normalized_all <- ScaleData(merged_NML_normalized, features = all.genes) #for all the genes
head(merged_NML_normalized@assays[["RNA"]]@layers[["scale.data"]])

#saving normalized data
saveRDS(merged_NML_normalized, file = "~/Desktop/Project_seurat_lungs/merged_NML_norm.RDS")


###Performing PCA on the scaled data(linear dimensionality reduction)#######

merged_NML_normalized <- RunPCA(merged_NML_normalized)

###examine and visualise PCA results 
DimPlot(merged_NML_normalized, reduction = "pca", dims = c(1,2))
#DimPlot(merged_NML_normalized, reduction = "pca", dims = c(1,10))
#DimPlot(merged_NML_normalized, reduction = "pca", dims = c(1,20))

##determine the dimension of the dataset##
ElbowPlot(merged_NML_normalized)
ElbowPlot(merged_NML_normalized, ndims = 50, reduction = "pca")


##clustering the cells######
merged_NML_normalized <- FindNeighbors(merged_NML_normalized, dims = 1:20)
merged_NML_normalized <- FindClusters(merged_NML_normalized) 

##Running UMAP#######Running non-linear dimensionality reduction(UMAP/TSNE)####
merged_NML_normalized <- RunUMAP(merged_NML_normalized, dims = 1:20)
p1 <- DimPlot(merged_NML_normalized, reduction = "umap", label = TRUE, repel = TRUE)
p2 <- DimPlot(merged_NML_normalized, reduction = "umap", label = TRUE, repel = TRUE,
              group.by = "orig.ident", cols = c("green", "blue","red"))

grid.arrange(p1, p2, ncol=2, nrow=2)

merged_NML_normalized <- RunTSNE(merged_NML_normalized)
DimPlot(merged_NML_normalized, reduction = "tsne")




####integrating scrna data#####

##setting up Seurat object##
merged_NML <- readRDS('Merged_NML_norm.RDS')

merged_NML.list <- SplitObject(merged_NML_normalized, split.by = 'orig.ident')
View(merged_NML.list)


#normalize and identify variable features for each dataset independently
for (i in 1:length(merged_NML.list)) {
  merged_NML.list[[i]] <- NormalizeData(object = merged_NML.list[[i]])
  merged_NML.list[[i]] <- FindVariableFeatures(object = merged_NML.list[[i]])
}


#select integration features
features <- SelectIntegrationFeatures(object.list=merged_NML.list)
dim(features)
#find integration anchors
mergedNML_anchors <-  FindIntegrationAnchors(merged_NML.list)

#integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)

#scale data, run pca, umap and visualise integrated data(standard workflow)

seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(seurat.integrated)
#seurat.integrated <- FindNeighbors(seurat.integrated, reduction ='pca', dims=1:30)
#seurat.integrated <- FindClusters(seurat.integrated, resolution=0.3)
seurat.integrated <- RunUMAP(seurat.integrated, dims=1:30)

p3 <- DimPlot(seurat.integrated, reduction = 'umap', label= TRUE)
p4 <-  DimPlot(seurat.integrated, reduction = 'umap', group.by = 'orig,ident',
               cols = "red", "green","blue")


grid.arrange(p3, p4, ncol=2, nrow=2)