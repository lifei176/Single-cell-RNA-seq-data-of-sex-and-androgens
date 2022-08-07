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
#human <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
#mouse <- useMart("ensembl",dataset="mmusculus_gene_ensembl")
#human_homo<- getLDS(attributes = c("hgnc_symbol"),filters="hgnc_symbol",values=down_shared$Detail,mart=human,attributesL = "mgi_symbol",martL=mouse,uniqueRows = TRUE)
#head(geneMm_deg_down_inter)
#geneMm_unique_deg_down_inter<-unique(geneMm_deg_down_inter$HGNC.ID)
#length(geneMm_unique_deg_down_inter)#215
#id<-gsub("HGNC:","",geneMm_unique_deg_down_inter) 
#################################Define significantly changed pathways for each organ
#################down_adipose
down_adipose <- bitr(deg[which(deg$Organ=="adipose"&deg$Change=="down"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                   toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                   OrgDb = org.Mm.eg.db)
head(down_adipose)
dim(down_adipose)
################GO
down_adipose_go <- enrichGO(gene      = down_adipose$ENTREZID,
                          OrgDb         = org.Mm.eg.db,
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.01,
                          readable      = TRUE)
dim(down_adipose_go)
barplot(down_adipose_go , showCategory=30)
cnetplot(down_adipose_go, showCategory="positive regulation of cytokine production",categorySize="pvalue")
dotplot(down_adipose_go ,showCategory=10,font.size=10,title='KEGG')
heatplot(down_adipose_go)
down_adipose_go[down_adipose_go@result$Description=="positive regulation of cytokine production",]
#################down_adrenal
{down_adrenal <- bitr(deg[which(deg$Organ=="adrenal"&deg$Change=="down"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                    toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                    OrgDb = org.Mm.eg.db)
  head(down_adrenal)
  down_adrenal_go <- enrichGO(gene      = down_adrenal$ENTREZID,
                            OrgDb         = org.Mm.eg.db,
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.01,
                            qvalueCutoff  = 0.01,
                            readable      = TRUE)
  dim(down_adrenal_go)}
#################down_bone_marrow
{down_bone_marrow <- bitr(deg[which(deg$Organ=="bone_marrow"&deg$Change=="down"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                        toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                        OrgDb = org.Mm.eg.db)
  head(down_bone_marrow)
  down_bone_marrow_go <- enrichGO(gene      = down_bone_marrow$ENTREZID,
                                OrgDb         = org.Mm.eg.db,
                                ont           = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.01,
                                readable      = TRUE)
  dim(down_bone_marrow_go)}
#################down_brain
{down_brain <- bitr(deg[which(deg$Organ=="brain"&deg$Change=="down"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                  toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                  OrgDb = org.Mm.eg.db)
  head(down_brain)
  down_brain_go <- enrichGO(gene      = down_brain$ENTREZID,
                          OrgDb         = org.Mm.eg.db,
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.01,
                          readable      = TRUE)
  dim(down_brain_go)}
#################down_colon
{down_colon <- bitr(deg[which(deg$Organ=="colon"&deg$Change=="down"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                  toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                  OrgDb = org.Mm.eg.db)
  head(down_colon)
  down_colon_go <- enrichGO(gene      = down_colon$ENTREZID,
                          OrgDb         = org.Mm.eg.db,
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.01,
                          readable      = TRUE)
  dim(down_colon_go)}
#################down_heart
{down_heart <- bitr(deg[which(deg$Organ=="heart"&deg$Change=="down"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                  toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                  OrgDb = org.Mm.eg.db)
  head(down_heart)
  down_heart_go <- enrichGO(gene      = down_heart$ENTREZID,
                          OrgDb         = org.Mm.eg.db,
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.01,
                          readable      = TRUE)
  dim(down_heart_go)}
#################down_intestine
{down_intestine <- bitr(deg[which(deg$Organ=="intestine"&deg$Change=="down"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                      toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                      OrgDb = org.Mm.eg.db)
  head(down_intestine)
  down_intestine_go <- enrichGO(gene      = down_intestine$ENTREZID,
                              OrgDb         = org.Mm.eg.db,
                              ont           = "BP",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.01,
                              qvalueCutoff  = 0.01,
                              readable      = TRUE)
  dim(down_intestine_go)}
#################down_kidney
{down_kidney <- bitr(deg[which(deg$Organ=="kidney"&deg$Change=="down"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                   toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                   OrgDb = org.Mm.eg.db)
  head(down_kidney)
  down_kidney_go <- enrichGO(gene      = down_kidney$ENTREZID,
                           OrgDb         = org.Mm.eg.db,
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.01,
                           readable      = TRUE)
  dim(down_kidney_go)}
#################down_lacrimal
{down_lacrimal <- bitr(deg[which(deg$Organ=="lacrimal"&deg$Change=="down"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                     toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                     OrgDb = org.Mm.eg.db)
  head(down_lacrimal)
  down_lacrimal_go <- enrichGO(gene      = down_lacrimal$ENTREZID,
                             OrgDb         = org.Mm.eg.db,
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.01,
                             readable      = TRUE)
  dim(down_lacrimal_go)}
#################down_liver
{down_liver <- bitr(deg[which(deg$Organ=="liver"&deg$Change=="down"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                  toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                  OrgDb = org.Mm.eg.db)
  head(down_liver)
  down_liver_go <- enrichGO(gene      = down_liver$ENTREZID,
                          OrgDb         = org.Mm.eg.db,
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.01,
                          readable      = TRUE)
  dim(down_liver_go)}
#################down_lung
{down_lung <- bitr(deg[which(deg$Organ=="lung"&deg$Change=="down"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                 toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                 OrgDb = org.Mm.eg.db)
  head(down_lung)
  down_lung_go <- enrichGO(gene      = down_lung$ENTREZID,
                         OrgDb         = org.Mm.eg.db,
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.01,
                         readable      = TRUE)
  dim(down_lung_go)}
#################down_pancreas
{down_pancreas <- bitr(deg[which(deg$Organ=="pancreas"&deg$Change=="down"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                     toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                     OrgDb = org.Mm.eg.db)
  head(down_pancreas)
  down_pancreas_go <- enrichGO(gene      = down_pancreas$ENTREZID,
                             OrgDb         = org.Mm.eg.db,
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.01,
                             readable      = TRUE)
  dim(down_pancreas_go)}
#################down_salivary
{down_salivary <- bitr(deg[which(deg$Organ=="salivary"&deg$Change=="down"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                     toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                     OrgDb = org.Mm.eg.db)
  head(down_salivary)
  down_salivary_go <- enrichGO(gene      = down_salivary$ENTREZID,
                             OrgDb         = org.Mm.eg.db,
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.01,
                             readable      = TRUE)
  dim(down_salivary_go)}
#################down_skeletalmuscle
{down_skeletalmuscle <- bitr(deg[which(deg$Organ=="skeletalmuscle"&deg$Change=="down"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                           toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                           OrgDb = org.Mm.eg.db)
  head(down_skeletalmuscle)
  down_skeletalmuscle_go <- enrichGO(gene      = down_skeletalmuscle$ENTREZID,
                                   OrgDb         = org.Mm.eg.db,
                                   ont           = "BP",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.01,
                                   readable      = TRUE)
  dim(down_skeletalmuscle_go)}
#################down_spleen
{down_spleen <- bitr(deg[which(deg$Organ=="spleen"&deg$Change=="down"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                   toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                   OrgDb = org.Mm.eg.db)
  head(down_spleen)
  down_spleen_go <- enrichGO(gene      = down_spleen$ENTREZID,
                           OrgDb         = org.Mm.eg.db,
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.01,
                           readable      = TRUE)
  dim(down_spleen_go)}
#################down_stomach
{down_stomach <- bitr(deg[which(deg$Organ=="stomach"&deg$Change=="down"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                    toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                    OrgDb = org.Mm.eg.db)
  head(down_stomach)
  down_stomach_go <- enrichGO(gene      = down_stomach$ENTREZID,
                            OrgDb         = org.Mm.eg.db,
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.01,
                            qvalueCutoff  = 0.01,
                            readable      = TRUE)
  dim(down_stomach_go)}
#################down_thymus
{down_thymus <- bitr(deg[which(deg$Organ=="thymus"&deg$Change=="down"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                   toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                   OrgDb = org.Mm.eg.db)
  head(down_thymus)
  down_thymus_go <- enrichGO(gene      = down_thymus$ENTREZID,
                           OrgDb         = org.Mm.eg.db,
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.01,
                           readable      = TRUE)
  dim(down_thymus_go)}

save(down_adipose_go,down_adrenal_go,down_brain_go,down_bone_marrow_go,down_colon_go,down_heart_go,down_intestine_go,
     down_kidney_go,down_lacrimal_go,down_liver_go,down_lung_go,down_pancreas_go,down_salivary_go,down_skeletalmuscle_go,down_spleen_go,
     down_stomach_go,down_thymus_go,file="down_go_filter.Rdata")

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


organ <- c("adipose","adrenal","bone_marrow","brain","colon","heart","intestine","kidney","lacrimal",
           "liver","lung","pancreas","salivary","skeletalmuscle","spleen","stomach","thymus")
organ_down_pathway<-data.frame(organ,c(nrow(down_adipose_go_data),nrow(down_adrenal_go_data),nrow(down_bone_marrow_go_data),nrow(down_brain_go_data),nrow(down_colon_go_data),nrow(down_heart_go_data),nrow(down_intestine_go_data),
                                     nrow(down_kidney_go_data),nrow(down_lacrimal_go_data),nrow(down_liver_go_data),nrow(down_lung_go_data),nrow(down_pancreas_go_data),nrow(down_salivary_go_data),nrow(down_skeletalmuscle_go_data),nrow(down_spleen_go_data),nrow(down_stomach_go_data),nrow(down_thymus_go_data)))
colnames(organ_down_pathway)<-c("Organ","Enriched_pathway_number")
write.csv(organ_down_pathway,"organ_down_pathway.csv")
library(ggpubr)
ggdotchart(organ_down_pathway, x = "Organ", y = "Enriched_pathway_number",
           color = "Organ",
           palette = rev(mycolors[c(1:17)]),
           sorting = "descending",                        
           add = "segments",                             
           xlab="", 
           rotate = TRUE,
           dot.size = 6 ,
           label=organ_down_pathway$Enriched_pathway_number)

#################################Construct a union for down_go data
down_union<-union(x=c(down_adipose_go_data$Description,down_adrenal_go_data$Description,down_bone_marrow_go_data$Description,down_brain_go_data$Description,
                    down_colon_go_data$Description,down_heart_go_data$Description,down_intestine_go_data$Description,down_kidney_go_data$Description), 
                y=c(down_lacrimal_go_data$Description,down_liver_go_data$Description,down_lung_go_data$Description,down_pancreas_go_data$Description,
                    down_salivary_go_data$Description,down_skeletalmuscle_go_data$Description,down_spleen_go_data$Description,down_stomach_go_data$Description,down_thymus_go_data$Description))
length(down_union)

{down_adipose_go_data_union<-down_adipose_go_data[down_union,]
  row.names(down_adipose_go_data_union)<-down_union
  down_adrenal_go_data_union<-down_adrenal_go_data[down_union,]
  row.names(down_adrenal_go_data_union)<-down_union
  down_bone_marrow_go_data_union<-down_bone_marrow_go_data[down_union,]
  row.names(down_bone_marrow_go_data_union)<-down_union
  down_brain_go_data_union<-down_brain_go_data[down_union,]
  row.names(down_brain_go_data_union)<-down_union
  down_colon_go_data_union<-down_colon_go_data[down_union,]
  row.names(down_colon_go_data_union)<-down_union
  down_heart_go_data_union<-down_heart_go_data[down_union,]
  row.names(down_heart_go_data_union)<-down_union
  down_intestine_go_data_union<-down_intestine_go_data[down_union,]
  row.names(down_intestine_go_data_union)<-down_union
  down_kidney_go_data_union<-down_kidney_go_data[down_union,]
  row.names(down_kidney_go_data_union)<-down_union
  
  down_lacrimal_go_data_union<-down_lacrimal_go_data[down_union,]
  row.names(down_lacrimal_go_data_union)<-down_union
  down_liver_go_data_union<-down_liver_go_data[down_union,]
  row.names(down_liver_go_data_union)<-down_union
  down_lung_go_data_union<-down_lung_go_data[down_union,]
  row.names(down_lung_go_data_union)<-down_union
  down_pancreas_go_data_union<-down_pancreas_go_data[down_union,]
  row.names(down_pancreas_go_data_union)<-down_union
  down_salivary_go_data_union<-down_salivary_go_data[down_union,]
  row.names(down_salivary_go_data_union)<-down_union
  down_skeletalmuscle_go_data_union<-down_skeletalmuscle_go_data[down_union,]
  row.names(down_skeletalmuscle_go_data_union)<-down_union
  down_spleen_go_data_union<-down_spleen_go_data[down_union,]
  row.names(down_spleen_go_data_union)<-down_union
  down_stomach_go_data_union<-down_stomach_go_data[down_union,]
  row.names(down_stomach_go_data_union)<-down_union
  down_thymus_go_data_union<-down_thymus_go_data[down_union,]
  row.names(down_thymus_go_data_union)<-down_union}

down_all_go_data_union<-cbind(down_adipose_go_data_union$pvalue,down_adrenal_go_data_union$pvalue,down_bone_marrow_go_data_union$pvalue,down_brain_go_data_union$pvalue,
                            down_colon_go_data_union$pvalue,down_heart_go_data_union$pvalue,down_intestine_go_data_union$pvalue,down_kidney_go_data_union$pvalue,
                            down_lacrimal_go_data_union$pvalue,down_liver_go_data_union$pvalue,down_lung_go_data_union$pvalue,down_pancreas_go_data_union$pvalue,
                            down_salivary_go_data_union$pvalue,down_skeletalmuscle_go_data_union$pvalue,down_spleen_go_data_union$pvalue,down_stomach_go_data_union$pvalue,down_thymus_go_data_union$pvalue)

row.names(down_all_go_data_union)<-down_union
colnames(down_all_go_data_union)<-organ

down_all_go_data_union_log<-(-log(down_all_go_data_union,10))
head(down_all_go_data_union_log)
pheatmap::pheatmap(down_all_go_data_union_log,cluster_rows =FALSE,cluster_cols = FALSE,color = colorRampPalette(c("firebrick3", "firebrick3"))(50),labels_row = FALSE)
write.csv(down_all_go_data_union_log,"down_all_go_data_union_log.csv")
Merge_count<-c()
for (i in 1:length(down_union)){
  Merge_count1<-table((down_all_go_data_union_log[i,])>1)[[1]]
  Merge_count<-c(Merge_count,Merge_count1)
}
length(Merge_count)

Median_pvalue<-c()
for (i in 1:length(down_union)){
  Median_pvalue1<-median(down_all_go_data_union_log[i,],na.rm = TRUE)
  Median_pvalue<-c(Median_pvalue,Median_pvalue1)
}
down_all_go_data_union_count<-as.data.frame(cbind(down_all_go_data_union_log,Merge_count,Median_pvalue))

head(down_all_go_data_union_count)
table(down_all_go_data_union_count$Merge_count)

library(ggpubr)
ggbarplot(data.frame(table(down_all_go_data_union_count$Merge_count)), x = "Var1", y = "Freq",
          color = "white",fill = "Var1",label=data.frame(table(down_all_go_data_union_count$Merge_count))$Freq,
          xlab="Shared sets", ylab="Number of pathways (GO BP)")+geom_vline(xintercept ="9",linetype="dashed")

down_all_go_data_union_count_order<-down_all_go_data_union_count[order(down_all_go_data_union_count$Median_pvalue,decreasing = T),]

annotation_row_pvalue<-data.frame(row.names(down_all_go_data_union_count_order[which(down_all_go_data_union_count_order$Merge_count>8),]),
                                  down_all_go_data_union_count_order[which(down_all_go_data_union_count_order$Merge_count>8),]$Median_pvalue,rep("down",nrow(down_all_go_data_union_count_order[which(down_all_go_data_union_count_order$Merge_count>8),])))
row.names(annotation_row_pvalue)<-row.names(down_all_go_data_union_count_order[which(down_all_go_data_union_count_order$Merge_count>8),])
annotation_row_pvalue<-annotation_row_pvalue[,-1]
colnames(annotation_row_pvalue)<-c("Neg_log10_pvalue","DEGs")
ann_colors = list(
  Neg_log10_pvalue = c("white", "firebrick3"),DEGs =c(down="white"))

colnames(down_all_go_data_union_count_order)<-c("adipose","adrenal","bone_marrow","brain","colon","heart","intestine","kidney","lacrimal",
                                              "liver","lung","pancreas","salivary","skeletalmuscle","spleen","stomach","thymus","Merge_count","Median_pvalue")

pheatmap::pheatmap(down_all_go_data_union_count_order[which(down_all_go_data_union_count_order$Merge_count>8),][,-c(18,19)] ,annotation_row=annotation_row_pvalue,annotation_colors = ann_colors,
                   color = colorRampPalette(c("pink", "firebrick3"))(100),cluster_cols = FALSE,na_col = "white",fontsize_row = 4)
pheatmap::pheatmap(down_all_go_data_union_count_order[which(down_all_go_data_union_count_order$Merge_count>8),][,-c(18,19)] ,annotation_row=annotation_row_pvalue,annotation_colors = ann_colors,
                   color = colorRampPalette(c("pink", "firebrick3"))(100),cluster_rows = FALSE,cluster_cols = FALSE,na_col = "white",fontsize_row = 4)

write.csv(down_all_go_data_union_count_order,"down_all_go_data_union_count_order.csv")
##################visualization for specific signalling pathway
cnetplot(down_adipose_go, showCategory="positive regulation of cytokine production",categorySize="pvalue")+ggtitle("down_adipose_go")
cnetplot(down_adrenal_go, showCategory="positive regulation of cytokine production",categorySize="pvalue")+ggtitle("down_adrenal_go")
cnetplot(down_bone_marrow_go, showCategory="positive regulation of cytokine production",categorySize="pvalue")+ggtitle("down_bone_marrow_go")
cnetplot(down_heart_go, showCategory="positive regulation of cytokine production",categorySize="pvalue")+ggtitle("down_heart_go")
cnetplot(down_lacrimal_go, showCategory="positive regulation of cytokine production",categorySize="pvalue")+ggtitle("down_lacrimal_go")
cnetplot(down_liver_go, showCategory="positive regulation of cytokine production",categorySize="pvalue")+ggtitle("down_liver_go")
cnetplot(down_lung_go, showCategory="positive regulation of cytokine production",categorySize="pvalue")+ggtitle("down_lung_go")
cnetplot(down_pancreas_go, showCategory="positive regulation of cytokine production",categorySize="pvalue")+ggtitle("down_pancreas_go")
cnetplot(down_salivary_go, showCategory="positive regulation of cytokine production",categorySize="pvalue")+ggtitle("down_salivary_go")
cnetplot(down_skeletalmuscle_go, showCategory="positive regulation of cytokine production",categorySize="pvalue")+ggtitle("down_skeletalmuscle_go")
cnetplot(down_stomach_go, showCategory="positive regulation of cytokine production",categorySize="pvalue")+ggtitle("down_stomach_go")

############################visualization
############################visualization
############################visualization
setwd("E:/lifei/androgen_project/20211216_Reanalyses/DEG/Pathway_enrichment/20220115/")
load("down_go_filter.Rdata")

library(DOSE)
library(enrichplot)
library(clusterProfiler)
options(connectionObserver = NULL)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(viridis)
organ_down_pathway<-read.csv("organ_down_pathway.csv")
down_all_go_data_union_count_order<-read.csv("down_all_go_data_union_count_order.csv",row.names = 1)
down_all_go_data_union_log<-read.csv("down_all_go_data_union_log.csv",row.names = 1)
pheatmap::pheatmap(down_all_go_data_union_log,cluster_rows =FALSE,cluster_cols = FALSE,color = colorRampPalette(c("#3C5488FF", "#3C5488FF"))(50),show_rownames=FALSE,na_col = "lightgrey")
pheatmap::pheatmap(down_all_go_data_union_log,cluster_rows =FALSE,cluster_cols = FALSE,color = colorRampPalette(c("#3C5488FF", "#3C5488FF"))(50),show_rownames=FALSE,na_col = "white")

pheatmap::pheatmap(t(down_all_go_data_union_log),cluster_rows =FALSE,cluster_cols = FALSE,color = colorRampPalette(c("#3C5488FF", "#3C5488FF"))(50),show_colnames=FALSE,na_col = "lightgrey")
pheatmap::pheatmap(t(down_all_go_data_union_log),cluster_rows =FALSE,cluster_cols = FALSE,color = colorRampPalette(c("#3C5488FF", "#3C5488FF"))(50),show_colnames=FALSE,na_col = "white")

c("#00A087FF","#DC0000FF","#3C5488FF")

#mat_down<-down_all_go_data_union_log
#mat_down[is.na(mat_down)]<-0
#pheatmap::pheatmap(mat_down,labels_row = FALSE, color= inferno(10))
#mat_breaks_down <- seq(min(mat), max(mat), length.out = 10)
mat_breaks_down<-c(0,2,4,6,8,10)
#pheatmap::pheatmap(mat_down,color=viridis(length(mat_breaks_down) - 1),breaks= mat_breaks_down,labels_row = FALSE,show_rownames=FALSE)
#pheatmap::pheatmap(down_all_go_data_union_log,cluster_rows =FALSE,cluster_cols = FALSE,color=viridis(length(mat_breaks_down) - 1),breaks= mat_breaks_down,show_rownames=FALSE)
pheatmap::pheatmap(mat_down,cluster_rows =FALSE,cluster_cols = FALSE,color=viridis(length(mat_breaks_down) - 1),breaks= mat_breaks_down,show_rownames=FALSE)

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
library(ggpubr)
organ_down_pathway$Color<-(mycolors[c(1:17)])
ggdotchart(organ_down_pathway, x = "Organ", y = "Enriched_pathway_number",
           color = "Organ",
           palette = organ_down_pathway[order(organ_down_pathway$Enriched_pathway_number,decreasing = F),]$Color,
           sorting = "descending",                        
           add = "segments",                             
           xlab="", 
           rotate = TRUE,
           dot.size = 6 ,
           label=organ_down_pathway$Enriched_pathway_number
)
library(ggpubr)
ggbarplot(data.frame(table(down_all_go_data_union_count_order$Merge_count)), x = "Var1", y = "Freq",title="Enriched_pathways_overlapping (down_DEGs)",
          color = "white",fill = "#3C5488FF",label=data.frame(table(down_all_go_data_union_count_order$Merge_count))$Freq,
          xlab="Shared sets", ylab="Number of pathways (GO BP)")+geom_vline(xintercept ="9",linetype="dashed")

#############Heatmap for showing overlapped pathways
#annotation_row_pvalue<-data.frame(row.names(down_all_go_data_union_count_order[which(down_all_go_data_union_count_order$Merge_count>8),]),
#                                  down_all_go_data_union_count_order[which(down_all_go_data_union_count_order$Merge_count>8),]$Median_pvalue,rep("down",nrow(down_all_go_data_union_count_order[which(down_all_go_data_union_count_order$Merge_count>8),])))
#row.names(annotation_row_pvalue)<-row.names(down_all_go_data_union_count_order[which(down_all_go_data_union_count_order$Merge_count>8),])
#annotation_row_pvalue<-annotation_row_pvalue[,-1]
#colnames(annotation_row_pvalue)<-c("Neg_log10_pvalue","DEGs")
#ann_colors = list(
#  Neg_log10_pvalue = c("white", "firebrick3"),DEGs =c(down="white"))
#pheatmap::pheatmap(down_all_go_data_union_count_order_filter[order(down_all_go_data_union_count_order_filter$Merge_count,decreasing = T),][,-c(18,19)] ,annotation_row=annotation_row_pvalue,annotation_colors = ann_colors,
#                   color = colorRampPalette(c("steelblue1", "steelblue4"))(100),cluster_cols = FALSE,cluster_rows = FALSE,na_col = "white",fontsize_row = 12)
down_all_go_data_union_count_order_filter<-down_all_go_data_union_count_order[which(down_all_go_data_union_count_order$Merge_count>8),]
pheatmap::pheatmap(down_all_go_data_union_count_order_filter[order(down_all_go_data_union_count_order_filter$Merge_count,decreasing = T),][,-c(18,19)] ,
                   color = colorRampPalette(c("steelblue1", "#3C5488FF"))(100),cluster_cols = FALSE,cluster_rows = FALSE,na_col = "white",fontsize_row = 10)

pheatmap::pheatmap(rbind(up_all_go_data_union_count_order_filter[order(up_all_go_data_union_count_order_filter$Merge_count,decreasing = T),][,-c(18,19)],-down_all_go_data_union_count_order_filter[order(down_all_go_data_union_count_order_filter$Merge_count,decreasing = T),][,-c(18,19)]) ,
                  cluster_cols = FALSE,cluster_rows = FALSE,na_col = "white",fontsize_row = 12)

library(RColorBrewer)
pheatmap::pheatmap(rbind(up_all_go_data_union_count_order_filter[order(up_all_go_data_union_count_order_filter$Merge_count,decreasing = T),][,-c(18,19)],-down_all_go_data_union_count_order_filter[order(down_all_go_data_union_count_order_filter$Merge_count,decreasing = T),][,-c(18,19)]) ,
                   cluster_cols = FALSE,cluster_rows = FALSE,na_col = "white",fontsize_row = 12,color=rev(brewer.pal(n = 11, name = "RdBu")))

##################Looking for overlapping genes in specifci pathways
library(VennDetail)
library(ggplot2)
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

#################################################
##################Looking for overlapping genes
library(VennDetail)
library(ggplot2)
top_pathway_down<-list(adipose=strsplit(down_adipose_go_data[which(down_adipose_go_data$Description=="positive regulation of cytokine production"),]$geneID,split="/")[[1]],
                          adrenal=strsplit(down_adrenal_go_data[which(down_adrenal_go_data$Description=="positive regulation of cytokine production"),]$geneID,split="/")[[1]],
                          bone_marrow=strsplit(down_bone_marrow_go_data[which(down_bone_marrow_go_data$Description=="positive regulation of cytokine production"),]$geneID,split="/")[[1]],
                          heart=strsplit(down_heart_go_data[which(down_heart_go_data$Description=="positive regulation of cytokine production"),]$geneID,split="/")[[1]],
                          lacrimal=strsplit(down_lacrimal_go_data[which(down_lacrimal_go_data$Description=="positive regulation of cytokine production"),]$geneID,split="/")[[1]],
                          liver=strsplit(down_liver_go_data[which(down_liver_go_data$Description=="positive regulation of cytokine production"),]$geneID,split="/")[[1]],
                          lung=strsplit(down_lung_go_data[which(down_lung_go_data$Description=="positive regulation of cytokine production"),]$geneID,split="/")[[1]],
                          salivary=strsplit(down_salivary_go_data[which(down_salivary_go_data$Description=="positive regulation of cytokine production"),]$geneID,split="/")[[1]],
                          skeletalmuscle=strsplit(down_skeletalmuscle_go_data[which(down_skeletalmuscle_go_data$Description=="positive regulation of cytokine production"),]$geneID,split="/")[[1]],
                          stomach=strsplit(down_stomach_go_data[which(down_stomach_go_data$Description=="positive regulation of cytokine production"),]$geneID,split="/")[[1]]
)
venn_top_pathway_down<-(venndetail(top_pathway_down))
#plot(venn_cytokine_production, type = "downset")
#plot(venn_cytokine_production, type = "vennpie")
#result(venn_cytokine_production)
#head(venn_cytokine_production@wide)
#write.csv(venn_cytokine_production@wide,"pathway_overlap/down_venn_cytokine_production.csv")
table(venn_top_pathway_down@wide$SharedSets)
head(venn_top_pathway_down@wide)

library(ggpubr)
ggbarplot(data.frame(table(venn_top_pathway_down@wide$SharedSets)), x = "Var1", y = "Freq",title = "positive regulation of cytokine production",
          color = "white",fill = "#3C5488FF",label=data.frame(table(venn_top_pathway_down@wide$SharedSets))$Freq,
          xlab="Shared sets", ylab="Number of Genes")

p1<-cnetplot(down_adipose_go, showCategory="positive regulation of cytokine production",categorySize="pvalue")+ggtitle("down_adipose_go")
p2<-cnetplot(down_adrenal_go, showCategory="positive regulation of cytokine production",categorySize="pvalue")+ggtitle("down_adrenal_go")
p3<-cnetplot(down_bone_marrow_go, showCategory="positive regulation of cytokine production",categorySize="pvalue")+ggtitle("down_bone_marrow_go")
p4<-cnetplot(down_heart_go, showCategory="positive regulation of cytokine production",categorySize="pvalue")+ggtitle("down_heart_go")
p5<-cnetplot(down_lacrimal_go, showCategory="positive regulation of cytokine production",categorySize="pvalue")+ggtitle("down_lacrimal_go")
p6<-cnetplot(down_liver_go, showCategory="positive regulation of cytokine production",categorySize="pvalue")+ggtitle("down_liver_go")
p7<-cnetplot(down_lung_go, showCategory="positive regulation of cytokine production",categorySize="pvalue")+ggtitle("down_lung_go")
p8<-cnetplot(down_salivary_go, showCategory="positive regulation of cytokine production",categorySize="pvalue")+ggtitle("down_salivary_go")
p9<-cnetplot(down_skeletalmuscle_go, showCategory="positive regulation of cytokine production",categorySize="pvalue")+ggtitle("down_skeletalmuscle_go")
p10<-cnetplot(down_stomach_go, showCategory="positive regulation of cytokine production",categorySize="pvalue")+ggtitle("down_stomach_go")
p1+p2+p3+p4+p5+p6+p7+p8+p9+p10


########################Network(Top10)
library(Cairo)
library(viridis)
load("down_go_filter.Rdata")
down_go_list_all<-list(down_adipose_go,down_adrenal_go,down_bone_marrow_go,down_brain_go,down_colon_go,down_heart_go,down_intestine_go,
                       down_kidney_go,down_lacrimal_go,down_liver_go,down_lung_go,down_pancreas_go,down_salivary_go,down_skeletalmuscle_go,down_spleen_go,
                       down_stomach_go,down_thymus_go)
organ <- c("adipose","adrenal","bone_marrow","brain","colon","heart","intestine","kidney","lacrimal",
           "liver","lung","pancreas","salivary","skeletalmuscle","spleen","stomach","thymus")

i=17
CairoPNG(file=paste0("E:/lifei/androgen_project/20211216_Reanalyses/DEG/Pathway_enrichment/Figure/Network_top10/top10_down_",organ[i],"_go.PNG",sep=""),width=400,height=400)
emapplot(pairwise_termsim(down_go_list_all[[i]]),showCategory = 10,cex_line = 0.5)+ggtitle(paste0("top10_up_",organ[i],"_go"))+scale_color_viridis(option = "viridis")
dev.off()

i=17
emapplot(pairwise_termsim(down_go_list_all[[i]]),showCategory = 10,cex_line = 0.5)+ggtitle(paste0("top10_down_",organ[i],"_go"))+scale_color_viridis(option = "viridis")

###########################Analysis on specific pathway for cytoscape
down_go_list<-list(down_adipose_go,down_adrenal_go,down_bone_marrow_go,down_heart_go,down_lacrimal_go,
                 down_liver_go,down_lung_go,down_salivary_go,down_skeletalmuscle_go,down_stomach_go)

head(down_go_list[[1]])
down_go<-data.frame()
for (i in (1:length(down_go_list))) {
  down_go_1<-as.data.frame(down_go_list[[i]][down_go_list[[i]]$Description=="positive regulation of cytokine production",])
  down_go<-rbind(down_go,down_go_1)
}

down_go$Organ<-c("adipose","adrenal","bone_marrow","heart","lacrimal","liver","lung","salivary","skeletalmuscle","stomach")
write.csv(down_go,"go_down.positive regulation of cytokine production.csv")

down_go_cytoscape<-data.frame()
for (i in (1:length(down_go_list))) {
  down_go_cytoscape1<-data.frame(strsplit(down_go$geneID[i],"/")[[1]],rep(down_go$Organ[i],length(strsplit(down_go$geneID[i],"/")[[1]])))
  down_go_cytoscape<-rbind(down_go_cytoscape,down_go_cytoscape1)
}
colnames(down_go_cytoscape)<-c("Gene","Organ")
write.csv(down_go_cytoscape,"go_down.positive regulation of cytokine production.csv")

###########################Analysis on overlapping pathway for cytoscape
organ <- c("adipose","adrenal","bone_marrow","brain","colon","heart","intestine","kidney","lacrimal",
           "liver","lung","pancreas","salivary","skeletalmuscle","spleen","stomach","thymus")

go_down_filter<-data.frame()
for (i in 1:length(organ)){
  go_down_filter1<-data.frame(rownames(down_all_go_data_union_count_order_filter),down_all_go_data_union_count_order_filter[,i],rep(organ[i],nrow(down_all_go_data_union_count_order_filter)))
  go_down_filter<-rbind(go_down_filter,go_down_filter1)
  
}
head(go_down_filter)
colnames(go_down_filter)<-c("Pathway","P_value","Organ")
go_down_filter<-na.omit(go_down_filter)
dim(go_down_filter)
table(go_down_filter$Organ)
write.csv(go_down_filter,"go_down_filter.csv")
  
###########################Analysis on all the pathways for cytoscape
go_down_all<-data.frame()
for (i in 1:length(organ)){
  go_down_all1<-data.frame(rownames(down_all_go_data_union_log),down_all_go_data_union_log[,i],rep(organ[i],nrow(down_all_go_data_union_log)))
  go_down_all<-rbind(go_down_all,go_down_all1)
  
}
head(go_down_all)
dim(go_down_all)
colnames(go_down_all)<-c("Pathway","P_value","Organ")
go_down_all<-na.omit(go_down_all)
dim(go_down_all)
table(go_down_all$Organ)
write.csv(go_down_all,"go_down_all.csv")

###########################Visualization on the enriched pathways of each organ
organ <- c("adipose","adrenal","bone_marrow","brain","colon","heart","intestine","kidney","lacrimal",
           "liver","lung","pancreas","salivary","skeletalmuscle","spleen","stomach","thymus")

down_go_list_all<-list(down_adipose_go,down_adrenal_go,down_bone_marrow_go,down_brain_go,down_colon_go,down_heart_go,down_intestine_go,
                     down_kidney_go,down_lacrimal_go,down_liver_go,down_lung_go,down_pancreas_go,down_salivary_go,down_skeletalmuscle_go,down_spleen_go,
                     down_stomach_go,down_thymus_go)
head(down_go_list_all[[1]])

down_go_all<-data.frame()
for (i in (1:length(down_go_list_all))) {
  if (nrow(down_go_list_all[[i]])>10) 
    down_go_all_1<-as.data.frame(down_go_list_all[[i]][c(1:10),]) else
      down_go_all_1<-as.data.frame(down_go_list_all[[i]])
    down_go_all_1$Organ<-rep(organ[i],nrow(down_go_all_1))
    down_go_all<-rbind(down_go_all,down_go_all_1)
}
table(down_go_all$Organ)
down_go_all$Pathway_organ<-paste0(down_go_all$Description,"_",down_go_all$Organ,sep="")
down_go_all$LogP<-round((-log(down_go_all$pvalue,10)),2)
head(down_go_all)

ggdotchart(down_go_all, x = "Pathway_organ", y = "LogP",
           color = "Organ",group = "Organ",title="Top10 enriched pathways of each organ (Down)",
           sorting = "descending", 
           palette = (mycolors[c(1:17)]),
           add = "segments",                             
           xlab="", 
           rotate = FALSE,
           dot.size = 3 )+geom_vline(xintercept =c(10,20,30,40,50,60,70,80,90,100,110,120,130,140,145,155,165),linetype="dashed")

ggdotchart(down_go_all, x = "Pathway_organ", y = "LogP",group = "Organ",color = "#3C5488FF",fill = "#3C5488FF",
           sorting = "descending", title="Top10 enriched pathways of each organ (Down)",
           add = "segments",                             
           xlab="", 
           rotate = FALSE,
           dot.size = 3)+geom_vline(xintercept =c(10,20,30,40,50,60,70,80,90,100,110,120,130,140,145,155,165),linetype="dashed")



