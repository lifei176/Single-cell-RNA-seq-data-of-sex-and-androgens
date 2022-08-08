my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', "#b20000",'#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')

##################################Meta.data for each tissue 
setwd("../Meta.data")
csvfile.meta.data = list.files(pattern="meta.data")
list2env(
  lapply(setNames(csvfile.meta.data, make.names(gsub("*.csv$", "", csvfile.meta.data))),
         read.csv,header=TRUE,check.names=FALSE), envir = .GlobalEnv)  
list.meta.data=list(meta.data.adipose,meta.data.adrenal,meta.data.bone_marrow,meta.data.brain,meta.data.colon,meta.data.heart,meta.data.intestine,meta.data.kidney,
                    meta.data.lacrimal,meta.data.liver,meta.data.lung,meta.data.pancreas,meta.data.salivary,meta.data.skeletalmuscle,meta.data.spleen,meta.data.stomach,meta.data.thymus)

setwd("../Cell_proportion")
organ <- c("adipose","adrenal","bone_marrow","brain","colon","heart","intestine","kidney","lacrimal",
           "liver","lung","pancreas","salivary","skeletalmuscle","spleen","stomach","thymus")


##################################Visualization for cell composition of each tissue
library(ggplot2)
i=2 #tissue selection
{cell_number_data<-data.frame(table((list.meta.data[[i]]$cell_name)))
cell_number_data[order(cell_number_data$Freq,decreasing = T),]$Var1
cell.prop<-as.data.frame(prop.table(table((list.meta.data[[i]]$cell_name),list.meta.data[[i]]$organ_condition)))
colnames(cell.prop)<-c("cluster","condition","proportion")
cell.prop$cluster<-factor(cell.prop$cluster,levels = cell_number_data[order(cell_number_data$Freq,decreasing = T),]$Var1)}

###########Barplot visualization of cell composition without cell cluster annotation
ggplot(cell.prop,aes(condition,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle(organ[i])+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values=my36colors)+theme(legend.position = 'none')

###########Barplot visualization of cell composition with cell cluster annotation
ggplot(cell.prop,aes(condition,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle(organ[i])+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values=my36colors)

###########Point visualization of cell composition
sum.FD<-sum(cell.prop[which(cell.prop$condition==names(table(cell.prop$condition))[1]),]$proportion)
sum.FS<-sum(cell.prop[which(cell.prop$condition==names(table(cell.prop$condition))[2]),]$proportion)
sum.MC<-sum(cell.prop[which(cell.prop$condition==names(table(cell.prop$condition))[3]),]$proportion)
sum.MS<-sum(cell.prop[which(cell.prop$condition==names(table(cell.prop$condition))[4]),]$proportion)
cell.prop$proportion.all<-c(rep(sum.FD,length(table(cell.prop$cluster))),rep(sum.FS,length(table(cell.prop$cluster))),rep(sum.MC,length(table(cell.prop$cluster))),rep(sum.MS,length(table(cell.prop$cluster))))
cell.prop$proportion.final<-cell.prop$proportion/cell.prop$proportion.all
ggplot(cell.prop,aes(x=condition,y=cluster,size=proportion.final)) + geom_point(aes(colour =cluster ))+
  scale_color_manual(values=my36colors)+theme_bw()

################Integrate cell proportion data of 17 tissues
cell.prop.all<-data.frame()
for (i in (1:17)){
cell_number_data<-data.frame(table((list.meta.data[[i]]$cell_name)))
  cell_number_data[order(cell_number_data$Freq,decreasing = T),]$Var1
  cell.prop<-as.data.frame(prop.table(table((list.meta.data[[i]]$cell_name),list.meta.data[[i]]$organ_condition)))
  colnames(cell.prop)<-c("cluster","condition","proportion")
  sum.FD<-sum(cell.prop[which(cell.prop$condition==names(table(cell.prop$condition))[1]),]$proportion)
  sum.FS<-sum(cell.prop[which(cell.prop$condition==names(table(cell.prop$condition))[2]),]$proportion)
  sum.MC<-sum(cell.prop[which(cell.prop$condition==names(table(cell.prop$condition))[3]),]$proportion)
  sum.MS<-sum(cell.prop[which(cell.prop$condition==names(table(cell.prop$condition))[4]),]$proportion)
  cell.prop$proportion.all<-c(rep(sum.FD,length(table(cell.prop$cluster))),rep(sum.FS,length(table(cell.prop$cluster))),rep(sum.MC,length(table(cell.prop$cluster))),rep(sum.MS,length(table(cell.prop$cluster))))
  cell.prop$proportion.final<-cell.prop$proportion/cell.prop$proportion.all  
  cell.prop.all<-rbind(cell.prop.all,cell.prop)
}
dim(cell.prop.all)
write.csv(cell.prop.all,"cell.prop.all.csv")
cell.prop.all$tissue=c(rep("adipose",128),rep("adrenal",84),rep("bone_marrow",64),rep("brain",84),rep("colon",52),
                       rep("heart",108),rep("intestine",60),rep("kidney",108),rep("lacrimal",60),rep("liver",40),
                       rep("lung",96),rep("pancreas",96),rep("salivary",80),rep("skeletalmuscle",112),rep("spleen",80),
                       rep("stomach",68),rep("thymus",44))


cell.prop.all$treatment<-c(rep("FD",table(cell.prop.all$condition)[1]),rep("FS",table(cell.prop.all$condition)[1]),rep("MC",table(cell.prop.all$condition)[1]),rep("MS",table(cell.prop.all$condition)[1]),
                           rep("FD",table(cell.prop.all$condition)[1*4+1]),rep("FS",table(cell.prop.all$condition)[1*4+1]),rep("MC",table(cell.prop.all$condition)[1*4+1]),rep("MS",table(cell.prop.all$condition)[1*4+1]),
                           rep("FD",table(cell.prop.all$condition)[2*4+1]),rep("FS",table(cell.prop.all$condition)[2*4+1]),rep("MC",table(cell.prop.all$condition)[2*4+1]),rep("MS",table(cell.prop.all$condition)[2*4+1]),
                           rep("FD",table(cell.prop.all$condition)[3*4+1]),rep("FS",table(cell.prop.all$condition)[3*4+1]),rep("MC",table(cell.prop.all$condition)[3*4+1]),rep("MS",table(cell.prop.all$condition)[3*4+1]),
                           rep("FD",table(cell.prop.all$condition)[4*4+1]),rep("FS",table(cell.prop.all$condition)[4*4+1]),rep("MC",table(cell.prop.all$condition)[4*4+1]),rep("MS",table(cell.prop.all$condition)[4*4+1]),
                           rep("FD",table(cell.prop.all$condition)[5*4+1]),rep("FS",table(cell.prop.all$condition)[5*4+1]),rep("MC",table(cell.prop.all$condition)[5*4+1]),rep("MS",table(cell.prop.all$condition)[5*4+1]),
                           rep("FD",table(cell.prop.all$condition)[6*4+1]),rep("FS",table(cell.prop.all$condition)[6*4+1]),rep("MC",table(cell.prop.all$condition)[6*4+1]),rep("MS",table(cell.prop.all$condition)[6*4+1]),
                           rep("FD",table(cell.prop.all$condition)[7*4+1]),rep("FS",table(cell.prop.all$condition)[7*4+1]),rep("MC",table(cell.prop.all$condition)[7*4+1]),rep("MS",table(cell.prop.all$condition)[7*4+1]),
                           rep("FD",table(cell.prop.all$condition)[8*4+1]),rep("FS",table(cell.prop.all$condition)[8*4+1]),rep("MC",table(cell.prop.all$condition)[8*4+1]),rep("MS",table(cell.prop.all$condition)[8*4+1]),
                           rep("FD",table(cell.prop.all$condition)[9*4+1]),rep("FS",table(cell.prop.all$condition)[9*4+1]),rep("MC",table(cell.prop.all$condition)[9*4+1]),rep("MS",table(cell.prop.all$condition)[9*4+1]),
                           rep("FD",table(cell.prop.all$condition)[10*4+1]),rep("FS",table(cell.prop.all$condition)[10*4+1]),rep("MC",table(cell.prop.all$condition)[10*4+1]),rep("MS",table(cell.prop.all$condition)[10*4+1]),
                           rep("FD",table(cell.prop.all$condition)[11*4+1]),rep("FS",table(cell.prop.all$condition)[11*4+1]),rep("MC",table(cell.prop.all$condition)[11*4+1]),rep("MS",table(cell.prop.all$condition)[11*4+1]),
                           rep("FD",table(cell.prop.all$condition)[12*4+1]),rep("FS",table(cell.prop.all$condition)[12*4+1]),rep("MC",table(cell.prop.all$condition)[12*4+1]),rep("MS",table(cell.prop.all$condition)[12*4+1]),
                           rep("FD",table(cell.prop.all$condition)[13*4+1]),rep("FS",table(cell.prop.all$condition)[13*4+1]),rep("MC",table(cell.prop.all$condition)[13*4+1]),rep("MS",table(cell.prop.all$condition)[13*4+1]),
                           rep("FD",table(cell.prop.all$condition)[14*4+1]),rep("FS",table(cell.prop.all$condition)[14*4+1]),rep("MC",table(cell.prop.all$condition)[14*4+1]),rep("MS",table(cell.prop.all$condition)[14*4+1]),
                           rep("FD",table(cell.prop.all$condition)[15*4+1]),rep("FS",table(cell.prop.all$condition)[15*4+1]),rep("MC",table(cell.prop.all$condition)[15*4+1]),rep("MS",table(cell.prop.all$condition)[15*4+1]),
                           rep("FD",table(cell.prop.all$condition)[16*4+1]),rep("FS",table(cell.prop.all$condition)[16*4+1]),rep("MC",table(cell.prop.all$condition)[16*4+1]),rep("MS",table(cell.prop.all$condition)[16*4+1]))

cell.number.tissue<-data.frame(Tissue=organ, Cell_number=as.numeric(54286, 33859,40253,58435,30414,85500,39871,39475,37354,131656,80436,30813,50887,75334,59156,28006,36964))
head(cell.prop.all)
dim(cell.prop.all)
length(table(cell.prop.all$tissue))

organ <- c("adipose","adrenal","bone_marrow","brain","colon","heart","intestine","kidney","lacrimal",
           "liver","lung","pancreas","salivary","skeletalmuscle","spleen","stomach","thymus")
organ[1]
tissue.cell.prop.all<-data.frame()
for (i in 1:17){
  tissue.cell.prop.all1<-cell.prop.all[cell.prop.all$tissue==organ[i]&cell.prop.all$treatment=="FD",]
  tissue.cell.prop.all<-rbind(tissue.cell.prop.all,tissue.cell.prop.all1)
}
dim(tissue.cell.prop.all)#341
head(tissue.cell.prop.all)

tissue.cell.prop.all.number<-c()
for (i in 1:17){
  tissue.cell.prop.all.number1<-rep(cell.number.tissue$Cell_number[i],table(tissue.cell.prop.all$tissue)[i])
  tissue.cell.prop.all.number<-c(tissue.cell.prop.all.number,tissue.cell.prop.all.number1)
}
length(tissue.cell.prop.all.number)
tissue.cell.prop.all$cell_number_tissue<-tissue.cell.prop.all.number

MSvsFS<-(cell.prop.all[cell.prop.all$treatment=="MS",]$proportion.final)/(cell.prop.all[cell.prop.all$treatment=="FS",]$proportion.final)
MSvsFS_FC<-log(MSvsFS,2)
MSvsFS_FC<-gsub(-Inf,-100,MSvsFS_FC)
MSvsFS_FC<-gsub(Inf,100,MSvsFS_FC)
MSvsFS_FC<-as.numeric(gsub(NaN,0,MSvsFS_FC))
length(MSvsFS_FC)

MCvsMS<-(cell.prop.all[cell.prop.all$treatment=="MC",]$proportion.final)/(cell.prop.all[cell.prop.all$treatment=="MS",]$proportion.final)
MCvsMS_FC<-log(MCvsMS,2)
MCvsMS_FC<-gsub(-Inf,-100,MCvsMS_FC)
MCvsMS_FC<-gsub(Inf,100,MCvsMS_FC)
MCvsMS_FC<-as.numeric(gsub(NaN,0,MCvsMS_FC))

FDvsFS<-(cell.prop.all[cell.prop.all$treatment=="FD",]$proportion.final)/(cell.prop.all[cell.prop.all$treatment=="FS",]$proportion.final)
FDvsFS_FC<-log(FDvsFS,2)
FDvsFS_FC<-gsub(-Inf,-100,FDvsFS_FC)
FDvsFS_FC<-gsub(Inf,100,FDvsFS_FC)
FDvsFS_FC<-as.numeric(gsub(NaN,0,FDvsFS_FC))

all_FC<-data.frame(MSvsFS_FC=MSvsFS_FC,MCvsMS_FC=MCvsMS_FC,FDvsFS_FC=FDvsFS_FC)#341

MS_cell_number<-cell.prop.all[cell.prop.all$treatment=="MS",]$proportion*tissue.cell.prop.all$cell_number_tissue
FS_cell_number<-cell.prop.all[cell.prop.all$treatment=="FS",]$proportion*tissue.cell.prop.all$cell_number_tissue
FD_cell_number<-cell.prop.all[cell.prop.all$treatment=="FD",]$proportion*tissue.cell.prop.all$cell_number_tissue
MC_cell_number<-cell.prop.all[cell.prop.all$treatment=="MC",]$proportion*tissue.cell.prop.all$cell_number_tissue

MS_cell_proportion<-cell.prop.all[cell.prop.all$treatment=="MS",]$proportion.final
FS_cell_proportion<-cell.prop.all[cell.prop.all$treatment=="FS",]$proportion.final
FD_cell_proportion<-cell.prop.all[cell.prop.all$treatment=="FD",]$proportion.final
MC_cell_proportion<-cell.prop.all[cell.prop.all$treatment=="MC",]$proportion.final

all_FC$cluster<-tissue.cell.prop.all$cluster
all_FC$cell_type<-tissue.cell.prop.all$cell_type
all_FC$tissue<-tissue.cell.prop.all$tissue

head(all_FC)
all_FC$cell_number_MS<-round(MS_cell_number,0)
all_FC$cell_number_FS<-round(FS_cell_number,0)
all_FC$cell_number_MC<-round(MC_cell_number,0)
all_FC$cell_number_FD<-round(FD_cell_number,0)

all_FC$cell_proportion_MS<-MS_cell_proportion
all_FC$cell_proportion_FS<-FS_cell_proportion
all_FC$cell_proportion_MC<-MC_cell_proportion
all_FC$cell_proportion_FD<-FD_cell_proportion

head(all_FC)
dim((all_FC))#341  13
write.csv(all_FC,"all_FC.csv")

all_FC_filter<-all_FC[(all_FC$MSvsFS_FC>0&all_FC$MSvsFS_FC<100)|(all_FC$MSvsFS_FC<0&all_FC$MSvsFS_FC>(-100)),]
all_FC_filter<-all_FC_filter[(all_FC_filter$MCvsMS_FC>0&all_FC_filter$MCvsMS_FC<100)|(all_FC_filter$MCvsMS_FC<0&all_FC_filter$MCvsMS_FC>(-100)),]
all_FC_filter<-all_FC_filter[(all_FC_filter$FDvsFS_FC>0&all_FC_filter$FDvsFS_FC<100)|(all_FC_filter$FDvsFS_FC<0&all_FC_filter$FDvsFS_FC>(-100)),]
dim(all_FC_filter)#337 13
head(all_FC_filter)

library(ggplot2)
library(viridis)
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

####################Pearson correlation
ggplot(all_FC_filter,aes(MSvsFS_FC,MCvsMS_FC))+geom_point(aes(colour = tissue),size=4)+
  geom_smooth(method="lm")+theme_bw()+ggtitle("pearson cor=0.5532 p-value<2.2e-16")+scale_color_manual(values=mycolors)
cor.test(all_FC_filter$MSvsFS_FC,all_FC_filter$MCvsMS_FC,method = "pearson")#pearson cor=-0.5532133  p-value < 2.2e-16

ggplot(all_FC_filter,aes(MSvsFS_FC,FDvsFS_FC))+geom_point(aes(colour = tissue),size=4)+
  geom_smooth(method="lm")+theme_bw()+ggtitle("pearson cor=0.3022 p-value=1.508e-08")+scale_color_manual(values=mycolors)
cor.test(all_FC_filter$MSvsFS_FC,all_FC_filter$FDvsFS_FC,method = "pearson")#pearson cor=0.3022377 p-value=1.508e-08



