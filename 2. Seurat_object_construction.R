############Here we took the adrenal gland as an example to introduce seurat object construction.
setwd("")
getwd()
{library(Seurat)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(cowplot)
  library(EnsDb.Mmusculus.v79)
  library(Signac)
  library(patchwork)}
#filtered_feature_bc_matrix.h5 file loading and seurat object construction
FD_1_adrenal_counts<- Read10X_h5("FD_1/filtered_feature_bc_matrix.h5")
FD_1_adrenal <- CreateSeuratObject(counts = FD_1_adrenal_counts,project = "FD_1_adrenal")

FD_2_adrenal_counts<- Read10X_h5("FD_2/filtered_feature_bc_matrix.h5")
FD_2_adrenal <- CreateSeuratObject(counts = FD_2_adrenal_counts,project = "FD_2_adrenal")

FS_1_adrenal_counts<- Read10X_h5("FS_1/filtered_feature_bc_matrix.h5")
FS_1_adrenal <- CreateSeuratObject(counts = FS_1_adrenal_counts,project = "FS_1_adrenal")

FS_2_adrenal_counts<- Read10X_h5("FS_2/filtered_feature_bc_matrix.h5")
FS_2_adrenal <- CreateSeuratObject(counts = FS_2_adrenal_counts,project = "FS_2_adrenal")

MC_1_adrenal_counts<- Read10X_h5("MC_1/filtered_feature_bc_matrix.h5")
MC_1_adrenal <- CreateSeuratObject(counts = MC_1_adrenal_counts,project = "MC_1_adrenal")

MC_2_adrenal_counts<- Read10X_h5("MC_2/filtered_feature_bc_matrix.h5")
MC_2_adrenal <- CreateSeuratObject(counts = MC_2_adrenal_counts,project = "MC_2_adrenal")

MS_1_adrenal_counts<- Read10X_h5("MS_1/filtered_feature_bc_matrix.h5")
MS_1_adrenal <- CreateSeuratObject(counts = MS_1_adrenal_counts,project = "MS_1_adrenal")

MS_2_adrenal_counts<- Read10X_h5("MS_2/filtered_feature_bc_matrix.h5")
MS_2_adrenal <- CreateSeuratObject(counts = MS_2_adrenal_counts,project = "MS_2_adrenal")

#merge eight seurat objects into one seurat object
adrenal.combined <-merge(FD_1_adrenal, y = c(FD_2_adrenal,FS_1_adrenal,FS_2_adrenal,MC_1_adrenal,MC_2_adrenal,MS_1_adrenal,MS_2_adrenal), add.cell.ids = c("FD_1_adrenal", "FD_2_adrenal","FS_1_adrenal","FS_2_adrenal","MC_1_adrenal","MC_2_adrenal","MS_1_adrenal","MS_2_adrenal"), project = "adrenal.combined")

#calculate the percentage of mitochondrial genes of each cell
adrenal.combined[["percent.mt"]] <- PercentageFeatureSet(adrenal.combined, pattern = "^mt-")
VlnPlot(adrenal.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)

#cell filtering for downstream analysis
adrenal.combined.filter<-subset(adrenal.combined,subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 10)#50928

###variable features definition
adrenal.combined.filter.list <- SplitObject(adrenal.combined.filter, split.by = "orig.ident")
adrenal.combined.filter.list <- lapply(X = adrenal.combined.filter.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = adrenal.combined.filter.list)
adrenal.combined.filter.list <- lapply(X = adrenal.combined.filter.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

#######remove adrenal.combined.filter to release computation memory 
rm(adrenal.combined.filter)
gc()
anchors <- FindIntegrationAnchors(object.list = adrenal.combined.filter.list, reduction = "rpca", 
                                  dims = 1:50)
#This step would waste much time
adrenal.combined.filter.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
rm(adrenal.combined.filter.list)
gc()

DefaultAssay(adrenal.combined.filter.integrated) <- "integrated"
adrenal.combined.filter.integrated <- ScaleData(adrenal.combined.filter.integrated, verbose = FALSE)
adrenal.combined.filter.integrated <- RunPCA(adrenal.combined.filter.integrated, verbose = FALSE)
adrenal.combined.filter.integrated <- RunUMAP(adrenal.combined.filter.integrated, dims = 1:50)
adrenal.combined.filter.integrated <- FindNeighbors(adrenal.combined.filter.integrated, dims = 1:50)
adrenal.combined.filter.integrated <- FindClusters(adrenal.combined.filter.integrated, resolution = 0.20,verbose = TRUE)#23 clusters

table((adrenal.combined.filter.integrated$seurat_clusters))
Idents(adrenal.combined.filter.integrated)<-adrenal.combined.filter.integrated$seurat_clusters
table(adrenal.combined.filter.integrated@meta.data$orig.ident)
treatment<-c(rep("FD_adrenal",times=10892),rep("FS_adrenal",times=7857),rep("MC_adrenal",times=8891),rep("MS_adrenal",times=6223))    
adrenal.combined.filter.integrated$treatment <-treatment
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', "#b20000",'#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')

#visualization
DimPlot(adrenal.combined.filter.integrated, reduction = "umap",cols = my36colors,label = TRUE)
DimPlot(adrenal.combined.filter.integrated,reduction = "umap",cols = my36colors,split.by="orig.ident",label = TRUE)
DimPlot(adrenal.combined.filter.integrated,reduction = "umap",cols = my36colors,split.by="treatment",label = TRUE)

####Medulla cell
FeaturePlot(adrenal.combined.filter.integrated, features = c("Sox10"),reduction = "umap",cols= c("gray", "red"))
####Cortex cell
FeaturePlot(adrenal.combined.filter.integrated, features = c("Nr5a1"),reduction = "umap",cols= c("gray", "red"))

###############Find markers of each cell type
DefaultAssay(adrenal.combined.filter.integrated)<-"RNA"
adrenal.combined.filter.integrated.markers <- FindAllMarkers(adrenal.combined.filter.integrated,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(adrenal.combined.filter.integrated.markers,"adrenal.combined.filter.integrated.markers.csv")

###########Heatmap plot for top10 markers of each cell type
top10_fc_adrenal <- adrenal.combined.filter.integrated.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
library(ggplot2)
library(RColorBrewer)
DoHeatmap(subset(adrenal.combined.filter.integrated,downsample=100), group.colors=my36colors,features=top10_fc_adrenal$gene) + NoLegend()+scale_fill_gradientn(colors = rev(colorRampPalette(brewer.pal(7, "RdYlBu"))(100)))

###########Dotplot for representative markers of each cell type
markers.adrenal.cluster0<-c("Pecam1","Igfbp5")
markers.adrenal.cluster1<-c("Nr5a1","1500015O10Rik")
markers.adrenal.cluster2<-c("Hba-a1")
markers.adrenal.cluster3<-c("Cyp4b1")
markers.adrenal.cluster4<-c("Cyp11b1")
markers.adrenal.cluster5<-c("Flt1")
markers.adrenal.cluster6<-c("Aqp1")
markers.adrenal.cluster7<-c("C1qa")
markers.adrenal.cluster8<-c("Col1a1","Cfh")
markers.adrenal.cluster9<-c("Gzma")
markers.adrenal.cluster10<-c("Cd3e")
markers.adrenal.cluster11<-c("Mgp")
markers.adrenal.cluster12<-c("Igkc")
markers.adrenal.cluster13<-c("Mki67")
markers.adrenal.cluster14<-c("Gsn")
markers.adrenal.cluster15<-c("Lsp1")
markers.adrenal.cluster16<-c("Apod")
markers.adrenal.cluster17<-c("Acta2")
markers.adrenal.cluster18<-c("Akr1d1")
markers.adrenal.cluster19<-c("Cma1")
markers.adrenal.cluster20<-c("Krt19")
markers.adrenal.cluster21<-c("Chga")
markers.adrenal.cluster22<-c("Siglech")
markers.to.plot.adrenal <- c(markers.adrenal.cluster0,markers.adrenal.cluster1,markers.adrenal.cluster2,markers.adrenal.cluster3,markers.adrenal.cluster4,markers.adrenal.cluster5,markers.adrenal.cluster6,markers.adrenal.cluster7,
                             markers.adrenal.cluster8,markers.adrenal.cluster9,markers.adrenal.cluster10,markers.adrenal.cluster11,markers.adrenal.cluster12,
                             markers.adrenal.cluster13,markers.adrenal.cluster14,markers.adrenal.cluster15,markers.adrenal.cluster16,markers.adrenal.cluster17,
                             markers.adrenal.cluster18,markers.adrenal.cluster19,markers.adrenal.cluster20,markers.adrenal.cluster21,markers.adrenal.cluster22
)
DotPlot(adrenal.combined.filter.integrated, features = c(markers.to.plot.adrenal),dot.scale = 4, cols = c("lightgrey", "red")) +
  RotatedAxis()

#save RDS
saveRDS(adrenal.combined.filter.integrated,file = "adrenal.combined.filter.integrated.20211216.rds")

#read RDS
adrenal.combined.filter.integrated<-readRDS("adrenal.combined.filter.integrated.20211216.rds")



