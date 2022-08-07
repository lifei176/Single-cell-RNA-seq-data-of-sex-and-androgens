#############################################################
setwd("../Pathway_enrichment/")
deg<-read.csv('../DEG.csv')
head(deg)
#############################################################signalling pathway analysis
library(DOSE)
library(enrichplot)
library(clusterProfiler)
options(connectionObserver = NULL)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(biomaRt)

{library(RColorBrewer)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  color70 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[-c(6,17,18,19)]
  #37ç§?
  color37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")
  #pie(rep(1,37), col=c(color37, 37))
  #20ç§?
  color20<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
  #pie(rep(1,20), col=c(color20, 20))
  mycolors<-unique(c(color70,color37,color20))}
#human <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
#mouse <- useMart("ensembl",dataset="mmusculus_gene_ensembl")
#human_homo<- getLDS(attributes = c("hgnc_symbol"),filters="hgnc_symbol",values=up_shared$Detail,mart=human,attributesL = "mgi_symbol",martL=mouse,uniqueRows = TRUE)
#head(geneMm_deg_up_inter)
#geneMm_unique_deg_up_inter<-unique(geneMm_deg_up_inter$HGNC.ID)
#length(geneMm_unique_deg_up_inter)#215
#id<-gsub("HGNC:","",geneMm_unique_deg_up_inter) 
#################################Define significantly changed pathways for each organ
#################up_adipose
up_adipose <- bitr(deg[which(deg$Organ=="adipose"&deg$Change=="up"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                     toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                     OrgDb = org.Mm.eg.db)
head(up_adipose)
dim(up_adipose)
################GO
up_adipose_go <- enrichGO(gene      = up_adipose$ENTREZID,
                            OrgDb         = org.Mm.eg.db,
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.01,
                            qvalueCutoff  = 0.01,
                            readable      = TRUE)
dim(up_adipose_go)
barplot(up_adipose_go , showCategory=30)
cnetplot(up_adipose_go, categorySize="pvalue")
dotplot(up_adipose_go ,showCategory=10,font.size=10,title='KEGG')
heatplot(up_adipose_go)
up_adipose_go[up_adipose_go@result$Description=="positive regulation of cytokine production",]
#################up_adrenal
{up_adrenal <- bitr(deg[which(deg$Organ=="adrenal"&deg$Change=="up"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                     toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                     OrgDb = org.Mm.eg.db)
head(up_adrenal)
up_adrenal_go <- enrichGO(gene      = up_adrenal$ENTREZID,
                            OrgDb         = org.Mm.eg.db,
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.01,
                            qvalueCutoff  = 0.01,
                            readable      = TRUE)
dim(up_adrenal_go)}
#################up_bone_marrow
{up_bone_marrow <- bitr(deg[which(deg$Organ=="bone_marrow"&deg$Change=="up"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                     toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                     OrgDb = org.Mm.eg.db)
head(up_bone_marrow)
up_bone_marrow_go <- enrichGO(gene      = up_bone_marrow$ENTREZID,
                            OrgDb         = org.Mm.eg.db,
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.01,
                            qvalueCutoff  = 0.01,
                            readable      = TRUE)
dim(up_bone_marrow_go)}
#################up_brain
{up_brain <- bitr(deg[which(deg$Organ=="brain"&deg$Change=="up"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                          toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                          OrgDb = org.Mm.eg.db)
  head(up_brain)
  up_brain_go <- enrichGO(gene      = up_brain$ENTREZID,
                                  OrgDb         = org.Mm.eg.db,
                                  ont           = "BP",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff  = 0.01,
                                  qvalueCutoff  = 0.01,
                                  readable      = TRUE)
  dim(up_brain_go)}
#################up_colon
{up_colon <- bitr(deg[which(deg$Organ=="colon"&deg$Change=="up"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                    toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                    OrgDb = org.Mm.eg.db)
  head(up_colon)
  up_colon_go <- enrichGO(gene      = up_colon$ENTREZID,
                            OrgDb         = org.Mm.eg.db,
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.01,
                            qvalueCutoff  = 0.01,
                            readable      = TRUE)
  dim(up_colon_go)}
#################up_heart
{up_heart <- bitr(deg[which(deg$Organ=="heart"&deg$Change=="up"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                    toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                    OrgDb = org.Mm.eg.db)
  head(up_heart)
  up_heart_go <- enrichGO(gene      = up_heart$ENTREZID,
                            OrgDb         = org.Mm.eg.db,
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.01,
                            qvalueCutoff  = 0.01,
                            readable      = TRUE)
  dim(up_heart_go)}
#################up_intestine
{up_intestine <- bitr(deg[which(deg$Organ=="intestine"&deg$Change=="up"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                    toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                    OrgDb = org.Mm.eg.db)
  head(up_intestine)
  up_intestine_go <- enrichGO(gene      = up_intestine$ENTREZID,
                            OrgDb         = org.Mm.eg.db,
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.01,
                            qvalueCutoff  = 0.01,
                            readable      = TRUE)
  dim(up_intestine_go)}
#################up_kidney
{up_kidney <- bitr(deg[which(deg$Organ=="kidney"&deg$Change=="up"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                        toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                        OrgDb = org.Mm.eg.db)
  head(up_kidney)
  up_kidney_go <- enrichGO(gene      = up_kidney$ENTREZID,
                                OrgDb         = org.Mm.eg.db,
                                ont           = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.01,
                                readable      = TRUE)
  dim(up_kidney_go)}
#################up_lacrimal
{up_lacrimal <- bitr(deg[which(deg$Organ=="lacrimal"&deg$Change=="up"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                     toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                     OrgDb = org.Mm.eg.db)
  head(up_lacrimal)
  up_lacrimal_go <- enrichGO(gene      = up_lacrimal$ENTREZID,
                             OrgDb         = org.Mm.eg.db,
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.01,
                             readable      = TRUE)
  dim(up_lacrimal_go)}
#################up_liver
{up_liver <- bitr(deg[which(deg$Organ=="liver"&deg$Change=="up"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                     toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                     OrgDb = org.Mm.eg.db)
  head(up_liver)
  up_liver_go <- enrichGO(gene      = up_liver$ENTREZID,
                             OrgDb         = org.Mm.eg.db,
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.01,
                             readable      = TRUE)
  dim(up_liver_go)}
#################up_lung
{up_lung <- bitr(deg[which(deg$Organ=="lung"&deg$Change=="up"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                     toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                     OrgDb = org.Mm.eg.db)
  head(up_lung)
  up_lung_go <- enrichGO(gene      = up_lung$ENTREZID,
                             OrgDb         = org.Mm.eg.db,
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.01,
                             readable      = TRUE)
  dim(up_lung_go)}
#################up_pancreas
{up_pancreas <- bitr(deg[which(deg$Organ=="pancreas"&deg$Change=="up"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                     toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                     OrgDb = org.Mm.eg.db)
  head(up_pancreas)
  up_pancreas_go <- enrichGO(gene      = up_pancreas$ENTREZID,
                             OrgDb         = org.Mm.eg.db,
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.01,
                             readable      = TRUE)
  dim(up_pancreas_go)}
#################up_salivary
{up_salivary <- bitr(deg[which(deg$Organ=="salivary"&deg$Change=="up"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                     toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                     OrgDb = org.Mm.eg.db)
  head(up_salivary)
  up_salivary_go <- enrichGO(gene      = up_salivary$ENTREZID,
                             OrgDb         = org.Mm.eg.db,
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.01,
                             readable      = TRUE)
  dim(up_salivary_go)}
#################up_skeletalmuscle
{up_skeletalmuscle <- bitr(deg[which(deg$Organ=="skeletalmuscle"&deg$Change=="up"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                     toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                     OrgDb = org.Mm.eg.db)
  head(up_skeletalmuscle)
  up_skeletalmuscle_go <- enrichGO(gene      = up_skeletalmuscle$ENTREZID,
                             OrgDb         = org.Mm.eg.db,
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.01,
                             readable      = TRUE)
  dim(up_skeletalmuscle_go)}
#################up_spleen
{up_spleen <- bitr(deg[which(deg$Organ=="spleen"&deg$Change=="up"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                     toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                     OrgDb = org.Mm.eg.db)
  head(up_spleen)
  up_spleen_go <- enrichGO(gene      = up_spleen$ENTREZID,
                             OrgDb         = org.Mm.eg.db,
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.01,
                             readable      = TRUE)
  dim(up_spleen_go)
  }
#################up_stomach
{up_stomach <- bitr(deg[which(deg$Organ=="stomach"&deg$Change=="up"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                   toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                   OrgDb = org.Mm.eg.db)
  head(up_stomach)
  up_stomach_go <- enrichGO(gene      = up_stomach$ENTREZID,
                           OrgDb         = org.Mm.eg.db,
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.01,
                           readable      = TRUE)
  dim(up_stomach_go)}
#################up_thymus
{up_thymus <- bitr(deg[which(deg$Organ=="thymus"&deg$Change=="up"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                     toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                     OrgDb = org.Mm.eg.db)
  head(up_thymus)
  up_thymus_go <- enrichGO(gene      = up_thymus$ENTREZID,
                             OrgDb         = org.Mm.eg.db,
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.01,
                             readable      = TRUE)
  dim(up_thymus_go)}
organ <- c("adipose","adrenal","bone_marrow","brain","colon","heart","intestine","kidney","lacrimal",
           "liver","lung","pancreas","salivary","skeletalmuscle","spleen","stomach","thymus")

setwd("E:/lifei/androgen_project/20211216_Reanalyses/DEG/Pathway_enrichment/20220115/")

save(up_adipose_go,up_adrenal_go,up_brain_go,up_bone_marrow_go,up_colon_go,up_heart_go,up_intestine_go,
     up_kidney_go,up_lacrimal_go,up_liver_go,up_lung_go,up_pancreas_go,up_salivary_go,up_skeletalmuscle_go,up_spleen_go,
     up_stomach_go,up_thymus_go,file="up_go_filter.Rdata")

{up_adipose_go_data<-data.frame(up_adipose_go)
row.names(up_adipose_go_data)<-up_adipose_go_data$Description
up_adrenal_go_data<-data.frame(up_adrenal_go)
row.names(up_adrenal_go_data)<-up_adrenal_go_data$Description
up_bone_marrow_go_data<-data.frame(up_bone_marrow_go)
row.names(up_bone_marrow_go_data)<-up_bone_marrow_go_data$Description
up_brain_go_data<-data.frame(up_brain_go)
row.names(up_brain_go_data)<-up_brain_go_data$Description
up_colon_go_data<-data.frame(up_colon_go)
row.names(up_colon_go_data)<-up_colon_go_data$Description
up_heart_go_data<-data.frame(up_heart_go)
row.names(up_heart_go_data)<-up_heart_go_data$Description
up_intestine_go_data<-data.frame(up_intestine_go)
row.names(up_intestine_go_data)<-up_intestine_go_data$Description
up_kidney_go_data<-data.frame(up_kidney_go)
row.names(up_kidney_go_data)<-up_kidney_go_data$Description

up_lacrimal_go_data<-data.frame(up_lacrimal_go)
row.names(up_lacrimal_go_data)<-up_lacrimal_go_data$Description
up_liver_go_data<-data.frame(up_liver_go)
row.names(up_liver_go_data)<-up_liver_go_data$Description
up_lung_go_data<-data.frame(up_lung_go)
row.names(up_lung_go_data)<-up_lung_go_data$Description
up_pancreas_go_data<-data.frame(up_pancreas_go)
row.names(up_pancreas_go_data)<-up_pancreas_go_data$Description
up_salivary_go_data<-data.frame(up_salivary_go)
row.names(up_salivary_go_data)<-up_salivary_go_data$Description
up_skeletalmuscle_go_data<-data.frame(up_skeletalmuscle_go)
row.names(up_skeletalmuscle_go_data)<-up_skeletalmuscle_go_data$Description
up_spleen_go_data<-data.frame(up_spleen_go)
row.names(up_spleen_go_data)<-up_spleen_go_data$Description
up_stomach_go_data<-data.frame(up_stomach_go)
row.names(up_stomach_go_data)<-up_stomach_go_data$Description
up_thymus_go_data<-data.frame(up_thymus_go)
row.names(up_thymus_go_data)<-up_thymus_go_data$Description}

organ <- c("adipose","adrenal","bone_marrow","brain","colon","heart","intestine","kidney","lacrimal",
           "liver","lung","pancreas","salivary","skeletalmuscle","spleen","stomach","thymus")
organ_up_pathway<-data.frame(organ,c(nrow(up_adipose_go_data),nrow(up_adrenal_go_data),nrow(up_bone_marrow_go_data),nrow(up_brain_go_data),nrow(up_colon_go_data),nrow(up_heart_go_data),nrow(up_intestine_go_data),
nrow(up_kidney_go_data),nrow(up_lacrimal_go_data),nrow(up_liver_go_data),nrow(up_lung_go_data),nrow(up_pancreas_go_data),nrow(up_salivary_go_data),nrow(up_skeletalmuscle_go_data),nrow(up_spleen_go_data),nrow(up_stomach_go_data),nrow(up_thymus_go_data)))
colnames(organ_up_pathway)<-c("Organ","Enriched_pathway_number")
write.csv(organ_up_pathway,"organ_up_pathway.csv")

organ_up_pathway$Color<-(mycolors[c(1:17)])
library(ggpubr)
ggdotchart(organ_up_pathway, x = "Organ", y = "Enriched_pathway_number",
           color = "Organ",
           palette = organ_up_pathway[order(organ_up_pathway$Enriched_pathway_number,decreasing = F),]$Color,
           sorting = "descending",                        
           add = "segments",                             
           xlab="", 
           rotate = TRUE,
           dot.size = 6 ,
           label=organ_up_pathway$Enriched_pathway_number
)


#################################Construct a union for up_go data
up_union<-union(x=c(up_adipose_go_data$Description,up_adrenal_go_data$Description,up_bone_marrow_go_data$Description,up_brain_go_data$Description,
                 up_colon_go_data$Description,up_heart_go_data$Description,up_intestine_go_data$Description,up_kidney_go_data$Description), 
                 y=c(up_lacrimal_go_data$Description,up_liver_go_data$Description,up_lung_go_data$Description,up_pancreas_go_data$Description,
                     up_salivary_go_data$Description,up_skeletalmuscle_go_data$Description,up_spleen_go_data$Description,up_stomach_go_data$Description,up_thymus_go_data$Description))
length(up_union)

up_union1<-unique(c(up_adipose_go_data$Description,up_adrenal_go_data$Description,up_bone_marrow_go_data$Description,up_brain_go_data$Description,
                    up_colon_go_data$Description,up_heart_go_data$Description,up_intestine_go_data$Description,up_kidney_go_data$Description, 
                up_lacrimal_go_data$Description,up_liver_go_data$Description,up_lung_go_data$Description,up_pancreas_go_data$Description,
                    up_salivary_go_data$Description,up_skeletalmuscle_go_data$Description,up_spleen_go_data$Description,up_stomach_go_data$Description,up_thymus_go_data$Description))

length(up_union1)


{up_adipose_go_data_union<-up_adipose_go_data[up_union,]
row.names(up_adipose_go_data_union)<-up_union
up_adrenal_go_data_union<-up_adrenal_go_data[up_union,]
row.names(up_adrenal_go_data_union)<-up_union
up_bone_marrow_go_data_union<-up_bone_marrow_go_data[up_union,]
row.names(up_bone_marrow_go_data_union)<-up_union
up_brain_go_data_union<-up_brain_go_data[up_union,]
row.names(up_brain_go_data_union)<-up_union
up_colon_go_data_union<-up_colon_go_data[up_union,]
row.names(up_colon_go_data_union)<-up_union
up_heart_go_data_union<-up_heart_go_data[up_union,]
row.names(up_heart_go_data_union)<-up_union
up_intestine_go_data_union<-up_intestine_go_data[up_union,]
row.names(up_intestine_go_data_union)<-up_union
up_kidney_go_data_union<-up_kidney_go_data[up_union,]
row.names(up_kidney_go_data_union)<-up_union

up_lacrimal_go_data_union<-up_lacrimal_go_data[up_union,]
row.names(up_lacrimal_go_data_union)<-up_union
up_liver_go_data_union<-up_liver_go_data[up_union,]
row.names(up_liver_go_data_union)<-up_union
up_lung_go_data_union<-up_lung_go_data[up_union,]
row.names(up_lung_go_data_union)<-up_union
up_pancreas_go_data_union<-up_pancreas_go_data[up_union,]
row.names(up_pancreas_go_data_union)<-up_union
up_salivary_go_data_union<-up_salivary_go_data[up_union,]
row.names(up_salivary_go_data_union)<-up_union
up_skeletalmuscle_go_data_union<-up_skeletalmuscle_go_data[up_union,]
row.names(up_skeletalmuscle_go_data_union)<-up_union
up_spleen_go_data_union<-up_spleen_go_data[up_union,]
row.names(up_spleen_go_data_union)<-up_union
up_stomach_go_data_union<-up_stomach_go_data[up_union,]
row.names(up_stomach_go_data_union)<-up_union
up_thymus_go_data_union<-up_thymus_go_data[up_union,]
row.names(up_thymus_go_data_union)<-up_union}

nrow(up_stomach_go_data)
#1245

up_stomach_go_data_union<-up_stomach_go_data[up_union,]
write.csv(up_stomach_go_data_union,"up_stomach_go_data_union.csv")
row.names(up_stomach_go_data_union)<-up_union
nrow(up_stomach_go_data_union[which(up_stomach_go_data_union$pvalue<0.01),])
#1247
length(row.names(up_stomach_go_data_union))
up_stomach_go_data_union[]
length(intersect(row.names(up_stomach_go_data),up_union))
#1245

up_all_go_data_union<-cbind(up_adipose_go_data_union$pvalue,up_adrenal_go_data_union$pvalue,up_bone_marrow_go_data_union$pvalue,up_brain_go_data_union$pvalue,
                        up_colon_go_data_union$pvalue,up_heart_go_data_union$pvalue,up_intestine_go_data_union$pvalue,up_kidney_go_data_union$pvalue,
                        up_lacrimal_go_data_union$pvalue,up_liver_go_data_union$pvalue,up_lung_go_data_union$pvalue,up_pancreas_go_data_union$pvalue,
                        up_salivary_go_data_union$pvalue,up_skeletalmuscle_go_data_union$pvalue,up_spleen_go_data_union$pvalue,up_stomach_go_data_union$pvalue,up_thymus_go_data_union$pvalue)

row.names(up_all_go_data_union)<-up_union
colnames(up_all_go_data_union)<-organ

up_all_go_data_union_log<-(-log(up_all_go_data_union,10))
head(up_all_go_data_union_log)
pheatmap::pheatmap(up_all_go_data_union_log,cluster_rows =FALSE,cluster_cols = FALSE,color = colorRampPalette(c("firebrick3", "firebrick3"))(50),labels_row = FALSE)
write.csv(up_all_go_data_union_log,"up_all_go_data_union_log.csv")

up_all_go_data_union_log<-data.frame(up_all_go_data_union_log)
row.names(up_all_go_data_union_log[which(up_all_go_data_union_log$stomach>2),])
setdiff(row.names(up_all_go_data_union_log[which(up_all_go_data_union_log$stomach>2),]),up_stomach_go_data$Description)




Merge_count<-c()
for (i in 1:length(up_union)){
Merge_count1<-table((up_all_go_data_union_log[i,])>1)[[1]]
Merge_count<-c(Merge_count,Merge_count1)
}
length(Merge_count)

Median_pvalue<-c()
for (i in 1:length(up_union)){
  Median_pvalue1<-median(up_all_go_data_union_log[i,],na.rm = TRUE)
  Median_pvalue<-c(Median_pvalue,Median_pvalue1)
}
up_all_go_data_union_count<-as.data.frame(cbind(up_all_go_data_union_log,Merge_count,Median_pvalue))

head(up_all_go_data_union_count)
table(up_all_go_data_union_count$Merge_count)

library(ggpubr)
ggbarplot(data.frame(table(up_all_go_data_union_count$Merge_count)), x = "Var1", y = "Freq",
          color = "white",fill = "Var1",label=data.frame(table(up_all_go_data_union_count$Merge_count))$Freq,
          xlab="Shared sets", ylab="Number of pathways (GO BP)")+geom_vline(xintercept ="9",linetype="dashed")

up_all_go_data_union_count_order<-up_all_go_data_union_count[order(up_all_go_data_union_count$Median_pvalue,decreasing = T),]
write.csv(up_all_go_data_union_count_order,"up_all_go_data_union_count_order.csv")

############################visualization
############################visualization
############################visualization
############################visualization
############################visualization
############################visualization
library(DOSE)
library(enrichplot)
library(clusterProfiler)
options(connectionObserver = NULL)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(biomaRt)
library(viridis)
setwd("E:/20220314Êý¾Ý¿½±´/androgen_project/20211216_Reanalyses/DEG/Pathway_enrichment/20220115/")
load("up_go_filter.Rdata")
organ_up_pathway<-read.csv("organ_up_pathway.csv")
up_all_go_data_union_count_order<-read.csv("up_all_go_data_union_count_order.csv",row.names = 1)
up_all_go_data_union_log<-read.csv("up_all_go_data_union_log.csv",row.names = 1)

pheatmap::pheatmap(up_all_go_data_union_log,cluster_rows =FALSE,cluster_cols = FALSE,color = colorRampPalette(c("firebrick3", "firebrick3"))(50),show_rownames=FALSE,na_col = "lightgrey")
pheatmap::pheatmap(up_all_go_data_union_log,cluster_rows =FALSE,cluster_cols = FALSE,color = colorRampPalette(c("firebrick3", "firebrick3"))(50),show_rownames=FALSE,na_col = "white")

pheatmap::pheatmap(t(up_all_go_data_union_log),cluster_rows =FALSE,cluster_cols = FALSE,color = colorRampPalette(c("firebrick3", "firebrick3"))(50),show_colnames=FALSE,na_col = "lightgrey")
pheatmap::pheatmap(t(up_all_go_data_union_log),cluster_rows =FALSE,cluster_cols = FALSE,color = colorRampPalette(c("firebrick3", "firebrick3"))(50),show_colnames=FALSE,na_col = "white")

mat_up<-up_all_go_data_union_log
mat_up[is.na(mat_up)]<-0
#pheatmap::pheatmap(mat_up,labels_row = FALSE, color= inferno(10))
#pheatmap::pheatmap(mat_up,color=plasma(length(mat_breaks_up) - 1),breaks= mat_breaks_up,show_rownames=FALSE)
#pheatmap::pheatmap(up_all_go_data_union_log,cluster_rows =FALSE,cluster_cols = FALSE,color=plasma(length(mat_breaks_up) - 1),breaks= mat_breaks_up,show_rownames=FALSE)
mat_breaks_up<-c(0,2,4,6,8,10)
pheatmap::pheatmap(mat_up,cluster_rows =FALSE,cluster_cols = FALSE,color=plasma(length(mat_breaks_up) - 1),breaks= mat_breaks_up,labels_row = FALSE,show_rownames=FALSE)

{library(RColorBrewer)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  color70 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[-c(6,17,18,19)]
  #37ç§?
  color37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")
  #pie(rep(1,37), col=c(color37, 37))
  #20ç§?
  color20<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
  #pie(rep(1,20), col=c(color20, 20))
  mycolors<-unique(c(color70,color37,color20))}
organ_up_pathway$Color<-(mycolors[c(1:17)])
library(ggpubr)
ggdotchart(organ_up_pathway, x = "Organ", y = "Enriched_pathway_number",
           color = "Organ",
           palette = organ_up_pathway[order(organ_up_pathway$Enriched_pathway_number,decreasing = F),]$Color,
           sorting = "descending",                        
           add = "segments",                             
           xlab="", 
           rotate = TRUE,
           dot.size = 6 ,
           label=organ_up_pathway$Enriched_pathway_number
)


library(ggpubr)
ggbarplot(data.frame(table(up_all_go_data_union_count_order$Merge_count)), x = "Var1", y = "Freq",title="Enriched_pathways_overlapping (up_DEGs)",
          color = "white",fill = "#E64B35FF",label=data.frame(table(up_all_go_data_union_count_order$Merge_count))$Freq, 
          xlab="Shared sets", ylab="Number of pathways (GO BP)")+geom_vline(xintercept ="9",linetype="dashed")

#############Heatmap for showing overlapped pathways
#annotation_row_pvalue<-data.frame(row.names(up_all_go_data_union_count_order[which(up_all_go_data_union_count_order$Merge_count>8),]),
#                                  up_all_go_data_union_count_order[which(up_all_go_data_union_count_order$Merge_count>8),]$Median_pvalue,rep("up",nrow(up_all_go_data_union_count_order[which(up_all_go_data_union_count_order$Merge_count>8),])))
#row.names(annotation_row_pvalue)<-row.names(up_all_go_data_union_count_order[which(up_all_go_data_union_count_order$Merge_count>8),])
#annotation_row_pvalue<-annotation_row_pvalue[,-1]
#colnames(annotation_row_pvalue)<-c("Neg_log10_pvalue","DEGs")
#ann_colors = list(
#  Neg_log10_pvalue = c("white", "firebrick3"),DEGs =c(up="white"))
#pheatmap::pheatmap(up_all_go_data_union_count_order[which(up_all_go_data_union_count_order$Merge_count>8),][,-c(18,19)] ,annotation_row=annotation_row_pvalue,annotation_colors = ann_colors,
#                   color = colorRampPalette(c("pink", "firebrick3"))(100),cluster_cols = FALSE,na_col = "white",fontsize_row = 12)
#up_all_go_data_union_count_order_filter<-up_all_go_data_union_count_order[which(up_all_go_data_union_count_order$Merge_count>8),]
#pheatmap::pheatmap(up_all_go_data_union_count_order_filter[order(up_all_go_data_union_count_order_filter$Merge_count,decreasing = T),][,-c(18,19)]  ,annotation_row=annotation_row_pvalue,annotation_colors = ann_colors,color = colorRampPalette(c("pink", "firebrick3"))(100),
#                   cluster_cols = FALSE,cluster_rows = FALSE,na_col = "white",fontsize_row = 12)

up_all_go_data_union_count_order_filter<-up_all_go_data_union_count_order[which(up_all_go_data_union_count_order$Merge_count>8),]
pheatmap::pheatmap(up_all_go_data_union_count_order_filter[order(up_all_go_data_union_count_order_filter$Merge_count,decreasing = T),][,-c(18,19)],color = colorRampPalette(c("pink", "firebrick3"))(100),
                   cluster_cols = FALSE,cluster_rows = FALSE,na_col = "white",fontsize_row = 12)

pheatmap::pheatmap(up_all_go_data_union_count_order_filter[order(up_all_go_data_union_count_order_filter$Merge_count,decreasing = T),][,-c(18,19)],color = colorRampPalette(c("pink", "firebrick3"))(100),
                   cluster_cols = FALSE,cluster_rows = FALSE,na_col = "white",fontsize_row = 12)

##################Looking for overlapping genes in specific pathways
library(VennDetail)
library(ggplot2)
{up_adipose_go_data<-data.frame(up_adipose_go)
  row.names(up_adipose_go_data)<-up_adipose_go_data$Description
  up_adrenal_go_data<-data.frame(up_adrenal_go)
  row.names(up_adrenal_go_data)<-up_adrenal_go_data$Description
  up_bone_marrow_go_data<-data.frame(up_bone_marrow_go)
  row.names(up_bone_marrow_go_data)<-up_bone_marrow_go_data$Description
  up_brain_go_data<-data.frame(up_brain_go)
  row.names(up_brain_go_data)<-up_brain_go_data$Description
  up_colon_go_data<-data.frame(up_colon_go)
  row.names(up_colon_go_data)<-up_colon_go_data$Description
  up_heart_go_data<-data.frame(up_heart_go)
  row.names(up_heart_go_data)<-up_heart_go_data$Description
  up_intestine_go_data<-data.frame(up_intestine_go)
  row.names(up_intestine_go_data)<-up_intestine_go_data$Description
  up_kidney_go_data<-data.frame(up_kidney_go)
  row.names(up_kidney_go_data)<-up_kidney_go_data$Description
  
  up_lacrimal_go_data<-data.frame(up_lacrimal_go)
  row.names(up_lacrimal_go_data)<-up_lacrimal_go_data$Description
  up_liver_go_data<-data.frame(up_liver_go)
  row.names(up_liver_go_data)<-up_liver_go_data$Description
  up_lung_go_data<-data.frame(up_lung_go)
  row.names(up_lung_go_data)<-up_lung_go_data$Description
  up_pancreas_go_data<-data.frame(up_pancreas_go)
  row.names(up_pancreas_go_data)<-up_pancreas_go_data$Description
  up_salivary_go_data<-data.frame(up_salivary_go)
  row.names(up_salivary_go_data)<-up_salivary_go_data$Description
  up_skeletalmuscle_go_data<-data.frame(up_skeletalmuscle_go)
  row.names(up_skeletalmuscle_go_data)<-up_skeletalmuscle_go_data$Description
  up_spleen_go_data<-data.frame(up_spleen_go)
  row.names(up_spleen_go_data)<-up_spleen_go_data$Description
  up_stomach_go_data<-data.frame(up_stomach_go)
  row.names(up_stomach_go_data)<-up_stomach_go_data$Description
  up_thymus_go_data<-data.frame(up_thymus_go)
  row.names(up_thymus_go_data)<-up_thymus_go_data$Description}

top_pathway_up<-list(adipose=strsplit(up_adipose_go_data[which(up_adipose_go_data$Description=="reactive oxygen species metabolic process"),]$geneID,split="/")[[1]],
                          adrenal=strsplit(up_adrenal_go_data[which(up_adrenal_go_data$Description=="reactive oxygen species metabolic process"),]$geneID,split="/")[[1]],
                          bone_marrow=strsplit(up_bone_marrow_go_data[which(up_bone_marrow_go_data$Description=="reactive oxygen species metabolic process"),]$geneID,split="/")[[1]],
                          heart=strsplit(up_heart_go_data[which(up_heart_go_data$Description=="reactive oxygen species metabolic process"),]$geneID,split="/")[[1]],
                          intestine=strsplit(up_intestine_go_data[which(up_intestine_go_data$Description=="reactive oxygen species metabolic process"),]$geneID,split="/")[[1]],
                          liver=strsplit(up_liver_go_data[which(up_liver_go_data$Description=="reactive oxygen species metabolic process"),]$geneID,split="/")[[1]],
                          lung=strsplit(up_lung_go_data[which(up_lung_go_data$Description=="reactive oxygen species metabolic process"),]$geneID,split="/")[[1]],
                          pancreas=strsplit(up_pancreas_go_data[which(up_pancreas_go_data$Description=="reactive oxygen species metabolic process"),]$geneID,split="/")[[1]],
                          skeletalmuscle=strsplit(up_skeletalmuscle_go_data[which(up_skeletalmuscle_go_data$Description=="reactive oxygen species metabolic process"),]$geneID,split="/")[[1]],
                          stomach=strsplit(up_stomach_go_data[which(up_stomach_go_data$Description=="reactive oxygen species metabolic process"),]$geneID,split="/")[[1]]
)
venn_top_pathway_up<-(venndetail(top_pathway_up))
#plot(venn_top_pathway_up, type = "upset")
#plot(venn_top_pathway_up, type = "vennpie")
#result(venn_leukocyte_migration)
#head(venn_leukocyte_migration@wide)
#write.csv(venn_leukocyte_migration@wide,"pathway_overlap/up_venn_leukocyte_migration.csv")
table(venn_top_pathway_up@wide$SharedSets)

library(ggpubr)
ggbarplot(data.frame(table(venn_top_pathway_up@wide$SharedSets)), x = "Var1", y = "Freq",title = "reactive oxygen species metabolic process",
          color = "white",fill = "#E64B35FF",label=data.frame(table(venn_top_pathway_up@wide$SharedSets))$Freq,
          xlab="Shared sets", ylab="Number of Genes")

head(venn_top_pathway_up@wide)
p1<-cnetplot(up_adipose_go, showCategory="reactive oxygen species metabolic process",categorySize="pvalue")+ggtitle("up_adipose_go")
p2<-cnetplot(up_adrenal_go, showCategory="reactive oxygen species metabolic process",categorySize="pvalue")+ggtitle("up_adrenal_go")
p3<-cnetplot(up_bone_marrow_go, showCategory="reactive oxygen species metabolic process",categorySize="pvalue")+ggtitle("up_bone_marrow_go")
p4<-cnetplot(up_heart_go, showCategory="reactive oxygen species metabolic process",categorySize="pvalue")+ggtitle("up_heart_go")
p5<-cnetplot(up_intestine_go, showCategory="reactive oxygen species metabolic process",categorySize="pvalue")+ggtitle("up_heart_go")
p6<-cnetplot(up_liver_go, showCategory="reactive oxygen species metabolic process",categorySize="pvalue")+ggtitle("up_liver_go")
p7<-cnetplot(up_lung_go, showCategory="reactive oxygen species metabolic process",categorySize="pvalue")+ggtitle("up_lung_go")
p8<-cnetplot(up_pancreas_go, showCategory="reactive oxygen species metabolic process",categorySize="pvalue")+ggtitle("up_pancreas_go")
p9<-cnetplot(up_skeletalmuscle_go, showCategory="reactive oxygen species metabolic process",categorySize="pvalue")+ggtitle("up_skeletalmuscle_go")
p10<-cnetplot(up_stomach_go, showCategory="reactive oxygen species metabolic process",categorySize="pvalue")+ggtitle("up_stomach_go")

p1+p2+p3+p4+p5+p6+p7+p8+p9+p10

cnetplot(up_adipose_go, showCategory=10,categorySize="pvalue")+ggtitle("up_adipose_go")

########################Network(Top10)
setwd("E:/lifei/androgen_project/20211216_Reanalyses/DEG/Pathway_enrichment/20220115/")
load("up_go_filter.Rdata")
up_go_list_all<-list(up_adipose_go,up_adrenal_go,up_bone_marrow_go,up_brain_go,up_colon_go,up_heart_go,up_intestine_go,
                       up_kidney_go,up_lacrimal_go,up_liver_go,up_lung_go,up_pancreas_go,up_salivary_go,up_skeletalmuscle_go,up_spleen_go,
                       up_stomach_go,up_thymus_go)
organ <- c("adipose","adrenal","bone_marrow","brain","colon","heart","intestine","kidney","lacrimal",
           "liver","lung","pancreas","salivary","skeletalmuscle","spleen","stomach","thymus")

library(Cairo)
library(viridis)
i=1
emapplot(pairwise_termsim(up_go_list_all[[i]]),showCategory = 10,cex_line = 0.5)+ggtitle(paste0("top10_up_",organ[i],"_go"))+scale_color_viridis(option = "inferno")

###########################Analysis on specific pathway for cytoscape
up_go_list<-list(up_adipose_go,up_adrenal_go,up_bone_marrow_go,up_heart_go,up_intestine_go,
            up_liver_go,up_lung_go,up_pancreas_go,up_skeletalmuscle_go,up_stomach_go)

head(up_go_list[[1]])
up_go<-data.frame()
for (i in (1:length(up_go_list))) {
  up_go_1<-as.data.frame(up_go_list[[i]][up_go_list[[i]]$Description=="reactive oxygen species metabolic process",])
  up_go<-rbind(up_go,up_go_1)
}
up_go$Organ<-c("adipose","adrenal","bone_marrow","heart","intestine","liver","lung","pancreas","skeletalmuscle","stomach")
write.csv(up_go,"go_up_reactive oxygen species metabolic process.csv")

up_go_cytoscape<-data.frame()
for (i in (1:length(up_go_list))) {
  up_go_cytoscape1<-data.frame(strsplit(up_go$geneID[i],"/")[[1]],rep(up_go$Organ[i],length(strsplit(up_go$geneID[i],"/")[[1]])))
  up_go_cytoscape<-rbind(up_go_cytoscape,up_go_cytoscape1)
}
colnames(up_go_cytoscape)<-c("Gene","Organ")
write.csv(up_go_cytoscape,"go_up_cytoscape.csv")

###########################Analysis on overlapping pathway for cytoscape
organ <- c("adipose","adrenal","bone_marrow","brain","colon","heart","intestine","kidney","lacrimal",
           "liver","lung","pancreas","salivary","skeletalmuscle","spleen","stomach","thymus")

go_up_filter<-data.frame()
for (i in 1:length(organ)){
  go_up_filter1<-data.frame(rownames(up_all_go_data_union_count_order_filter),up_all_go_data_union_count_order_filter[,i],rep(organ[i],nrow(up_all_go_data_union_count_order_filter)))
  go_up_filter<-rbind(go_up_filter,go_up_filter1)
  
}
head(go_up_filter)
colnames(go_up_filter)<-c("Pathway","P_value","Organ")
go_up_filter<-na.omit(go_up_filter)
dim(go_up_filter)
table(go_up_filter$Organ)
write.csv(go_up_filter,"go_up_filter.csv")

###########################Analysis on all the pathways for cytoscape
go_up_all<-data.frame()
for (i in 1:length(organ)){
  go_up_all1<-data.frame(rownames(up_all_go_data_union_log),up_all_go_data_union_log[,i],rep(organ[i],nrow(up_all_go_data_union_log)))
  go_up_all<-rbind(go_up_all,go_up_all1)
  
}
head(go_up_all)
dim(go_up_all)
colnames(go_up_all)<-c("Pathway","P_value","Organ")
go_up_all<-na.omit(go_up_all)
dim(go_up_all)
table(go_up_all$Organ)
write.csv(go_up_all,"go_up_all.csv")

###########################Visualization on the enriched pathways of each organ
organ <- c("adipose","adrenal","bone_marrow","brain","colon","heart","intestine","kidney","lacrimal",
           "liver","lung","pancreas","salivary","skeletalmuscle","spleen","stomach","thymus")

up_go_list_all<-list(up_adipose_go,up_adrenal_go,up_bone_marrow_go,up_brain_go,up_colon_go,up_heart_go,up_intestine_go,
     up_kidney_go,up_lacrimal_go,up_liver_go,up_lung_go,up_pancreas_go,up_salivary_go,up_skeletalmuscle_go,up_spleen_go,
     up_stomach_go,up_thymus_go)
head(up_go_list_all[[4]])

up_go_all<-data.frame()
for (i in (1:length(up_go_list_all))) {
  if (nrow(up_go_list_all[[i]])>10) 
      up_go_all_1<-as.data.frame(up_go_list_all[[i]][c(1:10),]) else
        up_go_all_1<-as.data.frame(up_go_list_all[[i]])
  up_go_all_1$Organ<-rep(organ[i],nrow(up_go_all_1))
  up_go_all<-rbind(up_go_all,up_go_all_1)
}
table(up_go_all$Organ)
up_go_all$Pathway_organ<-paste0(up_go_all$Description,"_",up_go_all$Organ,sep="")
up_go_all$LogP<-round((-log(up_go_all$pvalue,10)),2)
head(up_go_all)

ggdotchart(up_go_all, x = "Pathway_organ", y = "LogP",
           color = "Organ",group = "Organ", title="Top10 enriched pathways of each organ (Up)",
           sorting = "descending", 
           palette = (mycolors[c(1:15,17)]),
           add = "segments",                             
           xlab="", 
           rotate = FALSE,
           dot.size = 3)+geom_vline(xintercept =c(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150),linetype="dashed")

ggdotchart(up_go_all, x = "Pathway_organ", y = "LogP",group = "Organ",color = "#E64B35FF",fill = "#E64B35FF",
           sorting = "descending",  title="Top10 enriched pathways of each organ (Up)",
           add = "segments",                             
           xlab="", 
           rotate = FALSE,
           dot.size = 3)+geom_vline(xintercept =c(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150),linetype="dashed")

#################################20220426
###############MSvsFS_up
load("E:/20220314Êý¾Ý¿½±´/androgen_project/20211216_Reanalyses/DEG/Pathway_enrichment/20220115/up_go_filter.Rdata")
organ <- c("adipose","adrenal","bone marrow","brain","colon","heart","intestine","kidney","lacrimal",
           "liver","lung","pancreas","salivary","skeletal muscle","spleen","stomach","thymus")

{up_adipose_go_data<-data.frame(up_adipose_go)
  row.names(up_adipose_go_data)<-up_adipose_go_data$Description
  up_adrenal_go_data<-data.frame(up_adrenal_go)
  row.names(up_adrenal_go_data)<-up_adrenal_go_data$Description
  up_bone_marrow_go_data<-data.frame(up_bone_marrow_go)
  row.names(up_bone_marrow_go_data)<-up_bone_marrow_go_data$Description
  up_brain_go_data<-data.frame(up_brain_go)
  row.names(up_brain_go_data)<-up_brain_go_data$Description
  up_colon_go_data<-data.frame(up_colon_go)
  row.names(up_colon_go_data)<-up_colon_go_data$Description
  up_heart_go_data<-data.frame(up_heart_go)
  row.names(up_heart_go_data)<-up_heart_go_data$Description
  up_intestine_go_data<-data.frame(up_intestine_go)
  row.names(up_intestine_go_data)<-up_intestine_go_data$Description
  up_kidney_go_data<-data.frame(up_kidney_go)
  row.names(up_kidney_go_data)<-up_kidney_go_data$Description
  
  up_lacrimal_go_data<-data.frame(up_lacrimal_go)
  row.names(up_lacrimal_go_data)<-up_lacrimal_go_data$Description
  up_liver_go_data<-data.frame(up_liver_go)
  row.names(up_liver_go_data)<-up_liver_go_data$Description
  up_lung_go_data<-data.frame(up_lung_go)
  row.names(up_lung_go_data)<-up_lung_go_data$Description
  up_pancreas_go_data<-data.frame(up_pancreas_go)
  row.names(up_pancreas_go_data)<-up_pancreas_go_data$Description
  up_salivary_go_data<-data.frame(up_salivary_go)
  row.names(up_salivary_go_data)<-up_salivary_go_data$Description
  up_skeletalmuscle_go_data<-data.frame(up_skeletalmuscle_go)
  row.names(up_skeletalmuscle_go_data)<-up_skeletalmuscle_go_data$Description
  up_spleen_go_data<-data.frame(up_spleen_go)
  row.names(up_spleen_go_data)<-up_spleen_go_data$Description
  up_stomach_go_data<-data.frame(up_stomach_go)
  row.names(up_stomach_go_data)<-up_stomach_go_data$Description
  up_thymus_go_data<-data.frame(up_thymus_go)
  row.names(up_thymus_go_data)<-up_thymus_go_data$Description}

up_go_data_list_all<-list(up_adipose_go_data,up_adrenal_go_data,up_bone_marrow_go_data,up_brain_go_data,up_colon_go_data,up_heart_go_data,up_intestine_go_data,
                          up_kidney_go_data,up_lacrimal_go_data,up_liver_go_data,up_lung_go_data,up_pancreas_go_data,up_salivary_go_data,up_skeletalmuscle_go_data,up_spleen_go_data,
                          up_stomach_go_data,up_thymus_go_data)
dim(up_skeletalmuscle_go_data)
up_go_data_list_final<-data.frame()
for (i in (1:17)){
  up_go_data_list_select<-up_go_data_list_all[[i]]
  up_go_data_list_select$Tissue<-rep(organ[i],nrow(up_go_data_list_select))
  up_go_data_list_select$Change<-rep("Up",nrow(up_go_data_list_select))
  up_go_data_list_final<-rbind(up_go_data_list_final,up_go_data_list_select)
}
dim(up_go_data_list_final)
table(up_go_data_list_final$Tissue)
length(table(up_go_data_list_final$Tissue))
#######################################MSvsFS_down
load("E:/20220314Êý¾Ý¿½±´/androgen_project/20211216_Reanalyses/DEG/Pathway_enrichment/20220115/down_go_filter.Rdata")
organ <- c("adipose","adrenal","bone marrow","brain","colon","heart","intestine","kidney","lacrimal",
           "liver","lung","pancreas","salivary","skeletal muscle","spleen","stomach","thymus")

{down_adipose_go_data<-data.frame(down_adipose_go)
  row.names(down_adipose_go_data)<-down_adipose_go_data$Description
  down_adrenal_go_data<-data.frame(down_adrenal_go)
  row.names(down_adrenal_go_data)<-down_adrenal_go_data$Description
  down_bone_marrow_go_data<-data.frame(down_bone_marrow_go)
  row.names(down_bone_marrow_go_data)<-down_bone_marrow_go_data$Description
  down_brain_go_data<-data.frame(down_brain_go)
  row.names(down_brain_go_data)<-down_brain_go_data$Description
  down_colon_go_data<-data.frame(down_colon_go)
  row.names(down_colon_go_data)<-down_colon_go_data$Description
  down_heart_go_data<-data.frame(down_heart_go)
  row.names(down_heart_go_data)<-down_heart_go_data$Description
  down_intestine_go_data<-data.frame(down_intestine_go)
  row.names(down_intestine_go_data)<-down_intestine_go_data$Description
  down_kidney_go_data<-data.frame(down_kidney_go)
  row.names(down_kidney_go_data)<-down_kidney_go_data$Description
  
  down_lacrimal_go_data<-data.frame(down_lacrimal_go)
  row.names(down_lacrimal_go_data)<-down_lacrimal_go_data$Description
  down_liver_go_data<-data.frame(down_liver_go)
  row.names(down_liver_go_data)<-down_liver_go_data$Description
  down_lung_go_data<-data.frame(down_lung_go)
  row.names(down_lung_go_data)<-down_lung_go_data$Description
  down_pancreas_go_data<-data.frame(down_pancreas_go)
  row.names(down_pancreas_go_data)<-down_pancreas_go_data$Description
  down_salivary_go_data<-data.frame(down_salivary_go)
  row.names(down_salivary_go_data)<-down_salivary_go_data$Description
  down_skeletalmuscle_go_data<-data.frame(down_skeletalmuscle_go)
  row.names(down_skeletalmuscle_go_data)<-down_skeletalmuscle_go_data$Description
  down_spleen_go_data<-data.frame(down_spleen_go)
  row.names(down_spleen_go_data)<-down_spleen_go_data$Description
  down_stomach_go_data<-data.frame(down_stomach_go)
  row.names(down_stomach_go_data)<-down_stomach_go_data$Description
  down_thymus_go_data<-data.frame(down_thymus_go)
  row.names(down_thymus_go_data)<-down_thymus_go_data$Description}

down_go_data_list_all<-list(down_adipose_go_data,down_adrenal_go_data,down_bone_marrow_go_data,down_brain_go_data,down_colon_go_data,down_heart_go_data,down_intestine_go_data,
                          down_kidney_go_data,down_lacrimal_go_data,down_liver_go_data,down_lung_go_data,down_pancreas_go_data,down_salivary_go_data,down_skeletalmuscle_go_data,down_spleen_go_data,
                          down_stomach_go_data,down_thymus_go_data)

down_go_data_list_final<-data.frame()
for (i in (1:17)){
  down_go_data_list_select<-down_go_data_list_all[[i]]
  down_go_data_list_select$Tissue<-rep(organ[i],nrow(down_go_data_list_select))
  down_go_data_list_select$Change<-rep("down",nrow(down_go_data_list_select))
  down_go_data_list_final<-rbind(down_go_data_list_final,down_go_data_list_select)
}
dim(down_go_data_list_final)
table(down_go_data_list_final$Tissue)

Biological_pathway_MSvsFS<-rbind(up_go_data_list_final,down_go_data_list_final)
write.csv(Biological_pathway_MSvsFS,"Biological_pathway_MSvsFS.csv")

rm(list=ls())
gc()

###############FDvsFS_up
load("../Pathway_enrichment/up_go_filter_FDVSFS.Rdata")
organ <- c("adipose","adrenal","bone marrow","brain","colon","heart","intestine","kidney","lacrimal",
           "liver","lung","pancreas","salivary","skeletal muscle","spleen","stomach","thymus")

{up_adipose_go_data<-data.frame(up_adipose_go)
  row.names(up_adipose_go_data)<-up_adipose_go_data$Description
  up_adrenal_go_data<-data.frame(up_adrenal_go)
  row.names(up_adrenal_go_data)<-up_adrenal_go_data$Description
  up_bone_marrow_go_data<-data.frame(up_bone_marrow_go)
  row.names(up_bone_marrow_go_data)<-up_bone_marrow_go_data$Description
  up_brain_go_data<-data.frame(up_brain_go)
  row.names(up_brain_go_data)<-up_brain_go_data$Description
  up_colon_go_data<-data.frame(up_colon_go)
  row.names(up_colon_go_data)<-up_colon_go_data$Description
  up_heart_go_data<-data.frame(up_heart_go)
  row.names(up_heart_go_data)<-up_heart_go_data$Description
  up_intestine_go_data<-data.frame(up_intestine_go)
  row.names(up_intestine_go_data)<-up_intestine_go_data$Description
  up_kidney_go_data<-data.frame(up_kidney_go)
  row.names(up_kidney_go_data)<-up_kidney_go_data$Description
  
  up_lacrimal_go_data<-data.frame(up_lacrimal_go)
  row.names(up_lacrimal_go_data)<-up_lacrimal_go_data$Description
  up_liver_go_data<-data.frame(up_liver_go)
  row.names(up_liver_go_data)<-up_liver_go_data$Description
  up_lung_go_data<-data.frame(up_lung_go)
  row.names(up_lung_go_data)<-up_lung_go_data$Description
  up_pancreas_go_data<-data.frame(up_pancreas_go)
  row.names(up_pancreas_go_data)<-up_pancreas_go_data$Description
  up_salivary_go_data<-data.frame(up_salivary_go)
  row.names(up_salivary_go_data)<-up_salivary_go_data$Description
  up_skeletalmuscle_go_data<-data.frame(up_skeletalmuscle_go)
  row.names(up_skeletalmuscle_go_data)<-up_skeletalmuscle_go_data$Description
  up_spleen_go_data<-data.frame(up_spleen_go)
  row.names(up_spleen_go_data)<-up_spleen_go_data$Description
  up_stomach_go_data<-data.frame(up_stomach_go)
  row.names(up_stomach_go_data)<-up_stomach_go_data$Description
  up_thymus_go_data<-data.frame(up_thymus_go)
  row.names(up_thymus_go_data)<-up_thymus_go_data$Description}

up_go_data_list_all<-list(up_adipose_go_data,up_adrenal_go_data,up_bone_marrow_go_data,up_brain_go_data,up_colon_go_data,up_heart_go_data,up_intestine_go_data,
                          up_kidney_go_data,up_lacrimal_go_data,up_liver_go_data,up_lung_go_data,up_pancreas_go_data,up_salivary_go_data,up_skeletalmuscle_go_data,up_spleen_go_data,
                          up_stomach_go_data,up_thymus_go_data)

up_go_data_list_final<-data.frame()
for (i in (1:17)){
  up_go_data_list_select<-up_go_data_list_all[[i]]
  up_go_data_list_select$Tissue<-rep(organ[i],nrow(up_go_data_list_select))
  up_go_data_list_select$Change<-rep("Up",nrow(up_go_data_list_select))
  up_go_data_list_final<-rbind(up_go_data_list_final,up_go_data_list_select)
}
dim(up_go_data_list_final)
table(up_go_data_list_final$Tissue)
length(table(up_go_data_list_final$Tissue))
#######################################FDvsFS_down
load("../Pathway_enrichment/down_go_filter_FDVSFS.Rdata")
organ <- c("adipose","adrenal","bone marrow","brain","colon","heart","intestine","kidney","lacrimal",
           "liver","lung","pancreas","salivary","skeletal muscle","spleen","stomach","thymus")

{down_adipose_go_data<-data.frame(down_adipose_go)
  row.names(down_adipose_go_data)<-down_adipose_go_data$Description
  down_adrenal_go_data<-data.frame(down_adrenal_go)
  row.names(down_adrenal_go_data)<-down_adrenal_go_data$Description
  down_bone_marrow_go_data<-data.frame(down_bone_marrow_go)
  row.names(down_bone_marrow_go_data)<-down_bone_marrow_go_data$Description
  down_brain_go_data<-data.frame(down_brain_go)
  row.names(down_brain_go_data)<-down_brain_go_data$Description
  down_colon_go_data<-data.frame(down_colon_go)
  row.names(down_colon_go_data)<-down_colon_go_data$Description
  down_heart_go_data<-data.frame(down_heart_go)
  row.names(down_heart_go_data)<-down_heart_go_data$Description
  down_intestine_go_data<-data.frame(down_intestine_go)
  row.names(down_intestine_go_data)<-down_intestine_go_data$Description
  down_kidney_go_data<-data.frame(down_kidney_go)
  row.names(down_kidney_go_data)<-down_kidney_go_data$Description
  
  down_lacrimal_go_data<-data.frame(down_lacrimal_go)
  row.names(down_lacrimal_go_data)<-down_lacrimal_go_data$Description
  down_liver_go_data<-data.frame(down_liver_go)
  row.names(down_liver_go_data)<-down_liver_go_data$Description
  down_lung_go_data<-data.frame(down_lung_go)
  row.names(down_lung_go_data)<-down_lung_go_data$Description
  down_pancreas_go_data<-data.frame(down_pancreas_go)
  row.names(down_pancreas_go_data)<-down_pancreas_go_data$Description
  down_salivary_go_data<-data.frame(down_salivary_go)
  row.names(down_salivary_go_data)<-down_salivary_go_data$Description
  down_skeletalmuscle_go_data<-data.frame(down_skeletalmuscle_go)
  row.names(down_skeletalmuscle_go_data)<-down_skeletalmuscle_go_data$Description
  down_spleen_go_data<-data.frame(down_spleen_go)
  row.names(down_spleen_go_data)<-down_spleen_go_data$Description
  down_stomach_go_data<-data.frame(down_stomach_go)
  row.names(down_stomach_go_data)<-down_stomach_go_data$Description
  down_thymus_go_data<-data.frame(down_thymus_go)
  row.names(down_thymus_go_data)<-down_thymus_go_data$Description}

down_go_data_list_all<-list(down_adipose_go_data,down_adrenal_go_data,down_bone_marrow_go_data,down_brain_go_data,down_colon_go_data,down_heart_go_data,down_intestine_go_data,
                            down_kidney_go_data,down_lacrimal_go_data,down_liver_go_data,down_lung_go_data,down_pancreas_go_data,down_salivary_go_data,down_skeletalmuscle_go_data,down_spleen_go_data,
                            down_stomach_go_data,down_thymus_go_data)

down_go_data_list_final<-data.frame()
for (i in (1:17)){
  down_go_data_list_select<-down_go_data_list_all[[i]]
  down_go_data_list_select$Tissue<-rep(organ[i],nrow(down_go_data_list_select))
  down_go_data_list_select$Change<-rep("down",nrow(down_go_data_list_select))
  down_go_data_list_final<-rbind(down_go_data_list_final,down_go_data_list_select)
}
dim(down_go_data_list_final)
table(down_go_data_list_final$Tissue)

Biological_pathway_FDvsFS<-rbind(up_go_data_list_final,down_go_data_list_final)
write.csv(Biological_pathway_FDvsFS,"Biological_pathway_FDvsFS.csv")

rm(list=ls())
gc()


###############MCvsMS_up
load("../Pathway_enrichment/up_go_filter_MCVSMS.Rdata")
organ <- c("adipose","adrenal","bone marrow","brain","colon","heart","intestine","kidney","lacrimal",
           "liver","lung","pancreas","salivary","skeletal muscle","spleen","stomach","thymus")

{up_adipose_go_data<-data.frame(up_adipose_go)
  row.names(up_adipose_go_data)<-up_adipose_go_data$Description
  up_adrenal_go_data<-data.frame(up_adrenal_go)
  row.names(up_adrenal_go_data)<-up_adrenal_go_data$Description
  up_bone_marrow_go_data<-data.frame(up_bone_marrow_go)
  row.names(up_bone_marrow_go_data)<-up_bone_marrow_go_data$Description
  up_brain_go_data<-data.frame(up_brain_go)
  row.names(up_brain_go_data)<-up_brain_go_data$Description
  up_colon_go_data<-data.frame(up_colon_go)
  row.names(up_colon_go_data)<-up_colon_go_data$Description
  up_heart_go_data<-data.frame(up_heart_go)
  row.names(up_heart_go_data)<-up_heart_go_data$Description
  up_intestine_go_data<-data.frame(up_intestine_go)
  row.names(up_intestine_go_data)<-up_intestine_go_data$Description
  up_kidney_go_data<-data.frame(up_kidney_go)
  row.names(up_kidney_go_data)<-up_kidney_go_data$Description
  
  up_lacrimal_go_data<-data.frame(up_lacrimal_go)
  row.names(up_lacrimal_go_data)<-up_lacrimal_go_data$Description
  up_liver_go_data<-data.frame(up_liver_go)
  row.names(up_liver_go_data)<-up_liver_go_data$Description
  up_lung_go_data<-data.frame(up_lung_go)
  row.names(up_lung_go_data)<-up_lung_go_data$Description
  up_pancreas_go_data<-data.frame(up_pancreas_go)
  row.names(up_pancreas_go_data)<-up_pancreas_go_data$Description
  up_salivary_go_data<-data.frame(up_salivary_go)
  row.names(up_salivary_go_data)<-up_salivary_go_data$Description
  up_skeletalmuscle_go_data<-data.frame(up_skeletalmuscle_go)
  row.names(up_skeletalmuscle_go_data)<-up_skeletalmuscle_go_data$Description
  up_spleen_go_data<-data.frame(up_spleen_go)
  row.names(up_spleen_go_data)<-up_spleen_go_data$Description
  up_stomach_go_data<-data.frame(up_stomach_go)
  row.names(up_stomach_go_data)<-up_stomach_go_data$Description
  up_thymus_go_data<-data.frame(up_thymus_go)
  row.names(up_thymus_go_data)<-up_thymus_go_data$Description}

up_go_data_list_all<-list(up_adipose_go_data,up_adrenal_go_data,up_bone_marrow_go_data,up_brain_go_data,up_colon_go_data,up_heart_go_data,up_intestine_go_data,
                          up_kidney_go_data,up_lacrimal_go_data,up_liver_go_data,up_lung_go_data,up_pancreas_go_data,up_salivary_go_data,up_skeletalmuscle_go_data,up_spleen_go_data,
                          up_stomach_go_data,up_thymus_go_data)

up_go_data_list_final<-data.frame()
for (i in (1:17)){
  up_go_data_list_select<-up_go_data_list_all[[i]]
  up_go_data_list_select$Tissue<-rep(organ[i],nrow(up_go_data_list_select))
  up_go_data_list_select$Change<-rep("Up",nrow(up_go_data_list_select))
  up_go_data_list_final<-rbind(up_go_data_list_final,up_go_data_list_select)
}
dim(up_go_data_list_final)
table(up_go_data_list_final$Tissue)
length(table(up_go_data_list_final$Tissue))
#######################################MCvsMS_down
load("../Pathway_enrichment/down_go_filter_MCVSMS.Rdata")
organ <- c("adipose","adrenal","bone marrow","brain","colon","heart","intestine","kidney","lacrimal",
           "liver","lung","pancreas","salivary","skeletal muscle","spleen","stomach","thymus")

{down_adipose_go_data<-data.frame(down_adipose_go)
  row.names(down_adipose_go_data)<-down_adipose_go_data$Description
  down_adrenal_go_data<-data.frame(down_adrenal_go)
  row.names(down_adrenal_go_data)<-down_adrenal_go_data$Description
  down_bone_marrow_go_data<-data.frame(down_bone_marrow_go)
  row.names(down_bone_marrow_go_data)<-down_bone_marrow_go_data$Description
  down_brain_go_data<-data.frame(down_brain_go)
  row.names(down_brain_go_data)<-down_brain_go_data$Description
  down_colon_go_data<-data.frame(down_colon_go)
  row.names(down_colon_go_data)<-down_colon_go_data$Description
  down_heart_go_data<-data.frame(down_heart_go)
  row.names(down_heart_go_data)<-down_heart_go_data$Description
  down_intestine_go_data<-data.frame(down_intestine_go)
  row.names(down_intestine_go_data)<-down_intestine_go_data$Description
  down_kidney_go_data<-data.frame(down_kidney_go)
  row.names(down_kidney_go_data)<-down_kidney_go_data$Description
  
  down_lacrimal_go_data<-data.frame(down_lacrimal_go)
  row.names(down_lacrimal_go_data)<-down_lacrimal_go_data$Description
  down_liver_go_data<-data.frame(down_liver_go)
  row.names(down_liver_go_data)<-down_liver_go_data$Description
  down_lung_go_data<-data.frame(down_lung_go)
  row.names(down_lung_go_data)<-down_lung_go_data$Description
  down_pancreas_go_data<-data.frame(down_pancreas_go)
  row.names(down_pancreas_go_data)<-down_pancreas_go_data$Description
  down_salivary_go_data<-data.frame(down_salivary_go)
  row.names(down_salivary_go_data)<-down_salivary_go_data$Description
  down_skeletalmuscle_go_data<-data.frame(down_skeletalmuscle_go)
  row.names(down_skeletalmuscle_go_data)<-down_skeletalmuscle_go_data$Description
  down_spleen_go_data<-data.frame(down_spleen_go)
  row.names(down_spleen_go_data)<-down_spleen_go_data$Description
  down_stomach_go_data<-data.frame(down_stomach_go)
  row.names(down_stomach_go_data)<-down_stomach_go_data$Description
  down_thymus_go_data<-data.frame(down_thymus_go)
  row.names(down_thymus_go_data)<-down_thymus_go_data$Description}

down_go_data_list_all<-list(down_adipose_go_data,down_adrenal_go_data,down_bone_marrow_go_data,down_brain_go_data,down_colon_go_data,down_heart_go_data,down_intestine_go_data,
                            down_kidney_go_data,down_lacrimal_go_data,down_liver_go_data,down_lung_go_data,down_pancreas_go_data,down_salivary_go_data,down_skeletalmuscle_go_data,down_spleen_go_data,
                            down_stomach_go_data,down_thymus_go_data)

down_go_data_list_final<-data.frame()
for (i in (1:17)){
  down_go_data_list_select<-down_go_data_list_all[[i]]
  down_go_data_list_select$Tissue<-rep(organ[i],nrow(down_go_data_list_select))
  down_go_data_list_select$Change<-rep("down",nrow(down_go_data_list_select))
  down_go_data_list_final<-rbind(down_go_data_list_final,down_go_data_list_select)
}
dim(down_go_data_list_final)
table(down_go_data_list_final$Tissue)

Biological_pathway_MCvsMS<-rbind(up_go_data_list_final,down_go_data_list_final)
write.csv(Biological_pathway_MCvsMS,"Biological_pathway_MCvsMS.csv")

rm(list=ls())
gc()

