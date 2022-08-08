##########################Definition of significantly enriched biological pathways based on DEGs (as examplified by male-biased DEGs here)
setwd("../Pathway_enrichment/")
deg<-read.csv('../DEG.csv')
head(deg)
##########################library required packages and define color panels for downstream analysis  
##########################
##########################
##########################
##########################
##########################
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

##########################Define significantly changed pathways for each organ 
##########################
##########################
##########################
##########################
##########################

#################up_adipose
{up_adipose <- bitr(deg[which(deg$Organ=="adipose"&deg$Change=="up"&deg$Compare=="MSVSFS"),]$Gene_symbol, fromType = "SYMBOL",
                     toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                     OrgDb = org.Mm.eg.db)
head(up_adipose)
dim(up_adipose)
up_adipose_go <- enrichGO(gene      = up_adipose$ENTREZID,
                            OrgDb         = org.Mm.eg.db,
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.01,
                            qvalueCutoff  = 0.01,
                            readable      = TRUE)
dim(up_adipose_go)}
################visualization examples
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

#################save data
save(up_adipose_go,up_adrenal_go,up_brain_go,up_bone_marrow_go,up_colon_go,up_heart_go,up_intestine_go,
     up_kidney_go,up_lacrimal_go,up_liver_go,up_lung_go,up_pancreas_go,up_salivary_go,up_skeletalmuscle_go,up_spleen_go,
     up_stomach_go,up_thymus_go,file="up_go_filter.Rdata")

