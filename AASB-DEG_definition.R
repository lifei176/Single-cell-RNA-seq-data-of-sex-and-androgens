setwd("")
DEG<-read.csv("DEGs.csv")
tissue<-names(table(DEG$Tissue))

data.frame.AASB<-data.frame()
AASB_count<-data.frame()
AASB_tissue<-data.frame()
for (i in 1:17){
DEG_1<-DEG[which(DEG$Tissue==tissue[i]),]
DEG_1_MSvsFS<-droplevels(DEG_1[DEG_1$Comparison=="MSvsFS",])
DEG_1_MSvsFS_up<-droplevels(DEG_1_MSvsFS[DEG_1_MSvsFS$Change=="up",])
DEG_1_MSvsFS_down<-droplevels(DEG_1_MSvsFS[DEG_1_MSvsFS$Change=="down",])

DEG_1_MCvsMS<-droplevels(DEG_1[DEG_1$Comparison=="MCvsMS",])
DEG_1_MCvsMS_up<-droplevels(DEG_1_MCvsMS[DEG_1_MCvsMS$Change=="up",])
DEG_1_MCvsMS_down<-droplevels(DEG_1_MCvsMS[DEG_1_MCvsMS$Change=="down",])

DEG_1_FDvsFS<-droplevels(DEG_1[DEG_1$Comparison=="FDvsFS",])
DEG_1_FDvsFS_up<-droplevels(DEG_1_FDvsFS[DEG_1_FDvsFS$Change=="up",])
DEG_1_FDvsFS_down<-droplevels(DEG_1_FDvsFS[DEG_1_FDvsFS$Change=="down",])
#intersect cell types
cell_up<-intersect(names(table(DEG_1_MSvsFS_up$Cell_type)),intersect(names(table(DEG_1_MCvsMS_down$Cell_type)),
                                                           names(table(DEG_1_FDvsFS_up$Cell_type))))

cell_down<-intersect(names(table(DEG_1_MSvsFS_down$Cell_type)),intersect(names(table(DEG_1_MCvsMS_up$Cell_type)),
                                                                     names(table(DEG_1_FDvsFS_down$Cell_type))))
#positive AASB-DEGs
data.frame.up<-data.frame()
for (j in 1:length(cell_up)){
cell_up_MSvsFS<-as.character(droplevels(DEG_1_MSvsFS_up[which(DEG_1_MSvsFS_up$Cell_type==cell_up[j]),])$Gene_symbol)
cell_down_MCvsMS<-as.character(droplevels(DEG_1_MCvsMS_down[which(DEG_1_MCvsMS_down$Cell_type==cell_up[j]),])$Gene_symbol)
cell_up_FDvsFS<-as.character(droplevels(DEG_1_FDvsFS_up[which(DEG_1_FDvsFS_up$Cell_type==cell_up[j]),])$Gene_symbol)
positive_AASB_DEG<-intersect(cell_up_MSvsFS,intersect(cell_down_MCvsMS,cell_up_FDvsFS))
data.frame.up1<-data.frame(Tissue=rep(tissue[i],length(positive_AASB_DEG)),
                             Cell_type=rep(cell_up[j],length(positive_AASB_DEG)),
                             Change=rep("Positive",length(positive_AASB_DEG)),
                             Gene=positive_AASB_DEG)
data.frame.up<-rbind(data.frame.up,data.frame.up1)
}
#negative AASB-DEGs
data.frame.down<-data.frame()
for (h in 1:length(cell_down)){
  cell_down_MSvsFS<-as.character(droplevels(DEG_1_MSvsFS_down[which(DEG_1_MSvsFS_down$Cell_type==cell_down[h]),])$Gene_symbol)
  cell_up_MCvsMS<-as.character(droplevels(DEG_1_MCvsMS_up[which(DEG_1_MCvsMS_up$Cell_type==cell_down[h]),])$Gene_symbol)
  cell_down_FDvsFS<-as.character(droplevels(DEG_1_FDvsFS_down[which(DEG_1_FDvsFS_down$Cell_type==cell_down[h]),])$Gene_symbol)
  negative_AASB_DEG<-intersect(cell_down_MSvsFS,intersect(cell_up_MCvsMS,cell_down_FDvsFS))
  data.frame.down1<-data.frame(Tissue=rep(tissue[i],length(negative_AASB_DEG)),
                               Cell_type=rep(cell_down[h],length(negative_AASB_DEG)),
                               Change=rep("Negative",length(negative_AASB_DEG)),
                               Gene=negative_AASB_DEG)
  data.frame.down<-rbind(data.frame.down,data.frame.down1)
}
data.frame.AASB1<-rbind(data.frame.up,data.frame.down)
data.frame.AASB<-rbind(data.frame.AASB,data.frame.AASB1)#########all gene

AASB_count1<-data.frame(Tissue=tissue[i],Positive=length(unique(data.frame.up$Gene)),Negative=length(unique(data.frame.down$Gene)))
AASB_count<-rbind(AASB_count,AASB_count1)########unique count

AASB_tissue_positive<-data.frame(Tissue=rep(tissue[i],length(unique(data.frame.up$Gene))),
                                 Change= rep("Positive",length(unique(data.frame.up$Gene))),          
                                Gene= unique(data.frame.up$Gene))
AASB_tissue_negative<-data.frame(Tissue=rep(tissue[i],length(unique(data.frame.down$Gene))),
                                 Change= rep("Negative",length(unique(data.frame.down$Gene))),          
                                 Gene=unique(data.frame.down$Gene))                               
AASB_tissue1<-rbind(AASB_tissue_positive,AASB_tissue_negative)                                

AASB_tissue<-rbind(AASB_tissue,AASB_tissue1)#unique gene
}
dim(data.frame.AASB)

write.csv(data.frame.AASB,"../data.frame.AASB.csv")
write.csv(AASB_count,"../AASB_count.csv")
write.csv(AASB_tissue,"../AASB_tissue.csv")

#################Visualization
library(ggplot2)
library(ggrepel)
{library(RColorBrewer)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  color70 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[-c(6,17,18,19)]
  #37ç§?
  color37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")
  #pie(rep(1,37), col=sample(color37, 37))
  #20ç§?
  color20<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
  #pie(rep(1,20), col=sample(color20, 20))
  mycolors<-unique(c(color70,color37,color20))}

#########
deg_number$Tissue_number<-paste0(deg_number$Tissue,"_",deg_number$AASB_total,sep="")
p = ggplot(deg_number, aes(AASB_percentage, AASB_total)) +
  geom_point(aes(col=Tissue),size=4) +
  scale_color_manual(values=mycolors)+
  labs(title = "AASB-DEGs")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),panel.background=element_blank(),axis.line=element_line(colour='black'))
p+geom_text_repel(data=deg_number,size=4,aes(label=Tissue_number))

#########
AASB_tissue_positive<-droplevels(AASB_tissue[which(AASB_tissue$Change=="Positive"),])
AASB_tissue_positive_tissue<-data.frame(table(table(AASB_tissue_positive$Gene)))
ggplot(AASB_tissue_positive_tissue, aes(x=Var1, y=Freq)) + 
  geom_bar(stat = "identity")+theme_bw()+ggtitle("positive AASB-DEGs")+
  geom_text_repel(data=AASB_tissue_positive_tissue,size=4,aes(label=Freq))

AASB_tissue_negative<-droplevels(AASB_tissue[which(AASB_tissue$Change=="Negative"),])
AASB_tissue_negative_tissue<-data.frame(table(table(AASB_tissue_negative$Gene)))
ggplot(AASB_tissue_negative_tissue, aes(x=Var1, y=Freq)) + 
  geom_bar(stat = "identity")+theme_bw()+ggtitle("negative AASB-DEGs")+
  geom_text_repel(data=AASB_tissue_negative_tissue,size=4,aes(label=Freq))

################positive
head(AASB_tissue_positive)
table(AASB_tissue_positive$Tissue)
AASB_tissue_positive.heatmap<-data.frame(data1=c(1:length(unique(AASB_tissue_positive$Gene))),data2=c(1:length(unique(AASB_tissue_positive$Gene))))
for (i in (1:length(table(AASB_tissue_positive$Tissue)))){
  AASB_tissue_positive1<-AASB_tissue_positive[which(AASB_tissue_positive$Tissue==names(table(AASB_tissue_positive$Tissue))[i]),]
  row.names(AASB_tissue_positive1)<-AASB_tissue_positive1$Gene
  dim(AASB_tissue_positive1)
  AASB_tissue_positive1.select1<-AASB_tissue_positive1[as.character(unique(AASB_tissue_positive$Gene)),]
  AASB_tissue_positive.heatmap<-cbind(AASB_tissue_positive.heatmap,AASB_tissue_positive1.select1$Gene)
}

for (i in (1:17)){
  AASB_tissue_positive1<-AASB_tissue_positive[which(AASB_tissue_positive$Tissue==tissue[i]),]
  row.names(AASB_tissue_positive1)<-AASB_tissue_positive1$Gene
  dim(AASB_tissue_positive1)
  AASB_tissue_positive1.select1<-AASB_tissue_positive1[as.character(unique(AASB_tissue_positive$Gene)),]
  AASB_tissue_positive.heatmap<-cbind(AASB_tissue_positive.heatmap,AASB_tissue_positive1.select1$Gene)
}
AASB_tissue_positive.heatmap<-AASB_tissue_positive.heatmap[,-c(1,2)]
row.names(AASB_tissue_positive.heatmap)<-as.character(unique(AASB_tissue_positive$Gene))
colnames(AASB_tissue_positive.heatmap)<-tissue
head(AASB_tissue_positive.heatmap)
write.csv(AASB_tissue_positive.heatmap,"AASB_tissue_positive.heatmap.csv")

#test
{i=17
  AASB_tissue_positive1<-AASB_tissue_positive[which(AASB_tissue_positive$Tissue==tissue[i]),]
  row.names(AASB_tissue_positive1)<-AASB_tissue_positive1$Gene
  AASB_tissue_positive1
  AASB_tissue_positive1.select1<-AASB_tissue_positive1[as.character(unique(AASB_tissue_positive$Gene)),]
  448-nrow(AASB_tissue_positive1)-length(row.names(AASB_tissue_positive1.select1)[grep(pattern="NA",fixed = TRUE,row.names(AASB_tissue_positive1.select1))])
}
setdiff(row.names(AASB_tissue_positive1.select1)[grep(pattern=".1",fixed = TRUE,row.names(AASB_tissue_positive1.select1))],
        row.names(AASB_tissue_positive1.select1)[grep(pattern="NA",fixed = TRUE,row.names(AASB_tissue_positive1.select1))])


AASB_tissue_positive.heatmap<-read.csv("AASB_tissue_positive.heatmap.csv")
head(AASB_tissue_positive.heatmap)
AASB_tissue_positive.heatmap<-as.matrix(AASB_tissue_positive.heatmap)
head(AASB_tissue_positive.heatmap)
AASB_tissue_positive.heatmap[is.na(AASB_tissue_positive.heatmap)]<-0
AASB_tissue_positive.heatmap[AASB_tissue_positive.heatmap!="0"]=1
AASB_tissue_positive.heatmap<-apply(AASB_tissue_positive.heatmap,2,as.numeric)
head(AASB_tissue_positive.heatmap)
AASB_tissue_positive.heatmap<-as.data.frame(AASB_tissue_positive.heatmap)
row.names(AASB_tissue_positive.heatmap)<-as.character(unique(AASB_tissue_positive$Gene))
AASB_tissue_positive.heatmap<-AASB_tissue_positive.heatmap[,-1]
pheatmap::pheatmap(AASB_tissue_positive.heatmap,cluster_rows =FALSE,cluster_cols = FALSE,border=FALSE,color = colorRampPalette(c("white", "firebrick3"))(50),show_rownames=FALSE,na_col = "white")

AASB_tissue_positive.matrix<-AASB_tissue_positive.heatmap
AASB_tissue_positive.matrix$Merge_counts<-rowSums(AASB_tissue_positive.heatmap)
write.csv(AASB_tissue_positive.matrix,"AASB_tissue_positive.matrix.csv")

################negative
head(AASB_tissue_negative)
table(AASB_tissue_negative$Tissue)
AASB_tissue_negative.heatmap<-data.frame(data1=c(1:length(unique(AASB_tissue_negative$Gene))),data2=c(1:length(unique(AASB_tissue_negative$Gene))))
for (i in (1:17)){
  AASB_tissue_negative1<-AASB_tissue_negative[which(AASB_tissue_negative$Tissue==tissue[i]),]
  row.names(AASB_tissue_negative1)<-AASB_tissue_negative1$Gene
  dim(AASB_tissue_negative1)
  AASB_tissue_negative1.select1<-AASB_tissue_negative1[as.character(unique(AASB_tissue_negative$Gene)),]
  AASB_tissue_negative.heatmap<-cbind(AASB_tissue_negative.heatmap,AASB_tissue_negative1.select1$Gene)
}
AASB_tissue_negative.heatmap<-AASB_tissue_negative.heatmap[,-c(1,2)]
row.names(AASB_tissue_negative.heatmap)<-as.character(unique(AASB_tissue_negative$Gene))
colnames(AASB_tissue_negative.heatmap)<-tissue
head(AASB_tissue_negative.heatmap)
write.csv(AASB_tissue_negative.heatmap,"AASB_tissue_negative.heatmap.csv")


#test
{i=16
  AASB_tissue_negative1<-AASB_tissue_negative[which(AASB_tissue_negative$Tissue==tissue[i]),]
  row.names(AASB_tissue_negative1)<-AASB_tissue_negative1$Gene
  AASB_tissue_negative1
  AASB_tissue_negative1.select1<-AASB_tissue_negative1[as.character(unique(AASB_tissue_negative$Gene)),]
  634-nrow(AASB_tissue_negative1)-length(row.names(AASB_tissue_negative1.select1)[grep(pattern="NA",fixed = TRUE,row.names(AASB_tissue_negative1.select1))])
}
setdiff(row.names(AASB_tissue_negative1.select1)[grep(pattern=".1",fixed = TRUE,row.names(AASB_tissue_negative1.select1))],
        row.names(AASB_tissue_negative1.select1)[grep(pattern="NA",fixed = TRUE,row.names(AASB_tissue_negative1.select1))])


AASB_tissue_negative.heatmap<-read.csv("AASB_tissue_negative.heatmap.csv")
head(AASB_tissue_negative.heatmap)
AASB_tissue_negative.heatmap<-as.matrix(AASB_tissue_negative.heatmap)
head(AASB_tissue_negative.heatmap)
AASB_tissue_negative.heatmap[is.na(AASB_tissue_negative.heatmap)]<-0
AASB_tissue_negative.heatmap[AASB_tissue_negative.heatmap!="0"]=1
AASB_tissue_negative.heatmap<-apply(AASB_tissue_negative.heatmap,2,as.numeric)
head(AASB_tissue_negative.heatmap)
AASB_tissue_negative.heatmap<-as.data.frame(AASB_tissue_negative.heatmap)
row.names(AASB_tissue_negative.heatmap)<-as.character(unique(AASB_tissue_negative$Gene))
AASB_tissue_negative.heatmap<-AASB_tissue_negative.heatmap[,-1]
pheatmap::pheatmap(AASB_tissue_negative.heatmap,cluster_rows =FALSE,cluster_cols = FALSE,border=FALSE,color = colorRampPalette(c("white", "#3C5488FF"))(50),show_rownames=FALSE,na_col = "white")

AASB_tissue_negative.matrix<-AASB_tissue_negative.heatmap
AASB_tissue_negative.matrix$Merge_counts<-rowSums(AASB_tissue_negative.heatmap)
write.csv(AASB_tissue_negative.matrix,"AASB_tissue_negative.matrix.csv")
table(AASB_tissue_negative.matrix$Merge_counts)
table(AASB_tissue_positive.matrix$Merge_counts)
#

setdiff(row.names(AASB_tissue_negative1.select1)[grep(pattern=".1",row.names(AASB_tissue_positive1.select1))],
        row.names(AASB_tissue_negative1.select1)[grep(pattern="NA",row.names(AASB_tissue_positive1.select1))])



AASB_tissue_negative.heatmap<-read.csv("AASB_tissue_negative.heatmap.csv")
head(AASB_tissue_negative.heatmap)
AASB_tissue_negative.heatmap<-as.matrix(AASB_tissue_negative.heatmap)
head(AASB_tissue_negative.heatmap)
AASB_tissue_negative.heatmap[is.na(AASB_tissue_negative.heatmap)]<-0
AASB_tissue_negative.heatmap[AASB_tissue_negative.heatmap!="0"]=1
AASB_tissue_negative.heatmap<-apply(AASB_tissue_negative.heatmap,2,as.numeric)
head(AASB_tissue_negative.heatmap)
AASB_tissue_negative.heatmap<-as.data.frame(AASB_tissue_negative.heatmap)
row.names(AASB_tissue_negative.heatmap)<-as.character(unique(AASB_tissue_negative$Gene))
AASB_tissue_negative.heatmap<-AASB_tissue_negative.heatmap[,-1]
pheatmap::pheatmap(AASB_tissue_negative.heatmap,cluster_rows =FALSE,cluster_cols = FALSE,border=FALSE,color = colorRampPalette(c("white", "firebrick3"))(50),show_rownames=FALSE,na_col = "white")

