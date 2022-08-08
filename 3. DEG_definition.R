library(Seurat)
library(ggplot2)
setwd("..")
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', "#b20000",'#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')

brain.combined.filter.integrated <- readRDS("../brain.combined.filter.integrated.rds")
lung.combined.filter.integrated <- readRDS("../lung.combined.filter.integrated.rds")
pancreas.combined.filter.integrated <- readRDS("../pancreas.combined.filter.integrated.rds")
spleen.combined.filter.integrated <- readRDS("../spleen.combined.filter.integrated.rds")
stomach.combined.filter.integrated <- readRDS("../stomach.combined.filter.integrated.rds")

adipose.combined.filter.integrated<-readRDS("../adipose.combined.filter.integrated.rds")
heart.combined.filter.integrated<-readRDS("../heart.combined.filter.integrated.rds")
kidney.combined.filter.integrated<-readRDS("../kidney.combined.filter.integrated.rds")
skeletalmuscle.combined.filter.integrated<-readRDS("../skeletalmuscle.combined.filter.integrated.rds")

adrenal.combined.filter.integrated<-readRDS("../adrenal.combined.filter.integrated.rds")
bone_marrow.combined.filter.integrated<-readRDS("../bone_marrow.combined.filter.integrated.rds")
liver.combined.filter.integrated<-readRDS("../liver.combined.filter.integrated.rds")
salivary.combined.filter.integrated<-readRDS("../salivary.combined.filter.integrated.rds")
thymus.combined.filter.integrated<-readRDS("../thymus.combined.filter.integrated.rds")

colon.combined.filter.integrated<-readRDS("../colon.combined.filter.integrated.rds")
intestine.combined.filter.integrated<-readRDS("../intestine.combined.filter.integrated.rds")
lacrimal.combined.filter.integrated<-readRDS("../lacrimal.combined.filter.integrated.rds")

library(patchwork)
library(dplyr)
library(tidyr)
###################1. adipose: DEGs definition across multiple cell types and multiple conditions
###################
###################
###################
###################
###################
DefaultAssay(adipose.combined.filter.integrated) <- "RNA"
adipose.combined.filter.integrated$celltype.stim <- paste(Idents(adipose.combined.filter.integrated), adipose.combined.filter.integrated$condition, sep = "_")
adipose.combined.filter.integrated$celltype <- Idents(adipose.combined.filter.integrated)
Idents(adipose.combined.filter.integrated) <- "celltype.stim"
adipose<-adipose.combined.filter.integrated
table(adipose$celltype.stim)
length(table(adipose$celltype.stim))#128
length(table(adipose$cell_name))#32
#############Cell_type filtering for downstream DEG analysis
cell.use<-c()
for (i in 1:length(table(adipose$cell_name))){
FD_cell_number<-nrow(adipose@meta.data[which(adipose@meta.data$cell_name==names(table(adipose$cell_name))[i]&adipose@meta.data$condition=="FD"),])
FS_cell_number<-nrow(adipose@meta.data[which(adipose@meta.data$cell_name==names(table(adipose$cell_name))[i]&adipose@meta.data$condition=="FS"),])
MC_cell_number<-nrow(adipose@meta.data[which(adipose@meta.data$cell_name==names(table(adipose$cell_name))[i]&adipose@meta.data$condition=="MC"),])
MS_cell_number<-nrow(adipose@meta.data[which(adipose@meta.data$cell_name==names(table(adipose$cell_name))[i]&adipose@meta.data$condition=="MS"),])
cell.use<-c(cell.use,FD_cell_number>2&FS_cell_number>2&MC_cell_number>2&MS_cell_number>2)
}
cell.use.data<-data.frame(table(adipose$cell_name),cell.use)
celltypes_adipose_filter<-as.character(unique(cell.use.data[which(cell.use.data$cell.use==TRUE),]$Var1))

VS <- data.frame(compare=c('MSVSFS','FDVSFS','MCVSMS'),
                 ident1=c('_MS','_FD','_MC'),
                 ident2=c('_FS','_FS','_MS'))

###################Define DEGs across multiple cell types
deg_adipose <- data.frame()
for (j in (1:dim(VS)[1])){
  eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
  for (i in (1:length(celltypes_adipose_filter))){
    assign(paste0("deg_markers_",celltypes_adipose_filter[i],"_",format(VS[j,1]),"_adipose",sep=""),FindMarkers(adipose, ident.1 = paste0(celltypes_adipose_filter[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_adipose_filter[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
    eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_adipose_filter[i]),"_",format(VS[j,1]),"_adipose)[1]>0){deg_markers_",format(celltypes_adipose_filter[i]),"_",format(VS[j,1]),"_adipose <- subset(deg_markers_",format(celltypes_adipose_filter[i]),"_",format(VS[j,1]),"_adipose,p_val_adj<0.05)};
						   if(dim(deg_markers_",format(celltypes_adipose_filter[i]),"_",format(VS[j,1]),"_adipose)[1]>0){deg_markers_",format(celltypes_adipose_filter[i]),"_",format(VS[j,1]),"_adipose$Gene_symbol <- row.names(deg_markers_",format(celltypes_adipose_filter[i]),"_",format(VS[j,1]),"_adipose);
						   deg_markers_",format(celltypes_adipose_filter[i]),"_",format(VS[j,1]),"_adipose$Cell_type <- '",format(celltypes_adipose_filter[i]),"';
						   deg_markers_",format(celltypes_adipose_filter[i]),"_",format(VS[j,1]),"_adipose$Change <- ifelse(deg_markers_",format(celltypes_adipose_filter[i]),"_",format(VS[j,1]),"_adipose$avg_log2FC > 0 ,'up','down');
						   deg_markers_",format(celltypes_adipose_filter[i]),"_",format(VS[j,1]),"_adipose$Organ <- 'adipose';
						   deg_markers_",format(celltypes_adipose_filter[i]),"_",format(VS[j,1]),"_adipose$Compare <- '",format(VS[j,1]),"';",
                             format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_adipose_filter[i]),"_",format(VS[j,1]),"_adipose)};
						   rm(deg_markers_",format(celltypes_adipose_filter[i]),"_",format(VS[j,1]),"_adipose)")))
  }
  gc()
  eval(parse(text = paste0("deg_adipose <- bind_rows(deg_adipose,",format(VS[j,1]),");
                         rm(",format(VS[j,1]),")")))
}
dim(deg_adipose)
write.csv(deg_adipose,'E:/lifei/androgen_project/20211216_Reanalyses/DEG/DEG_list/DEGs.adipose.avg_log2FC05.p_val_adj005.csv')

###################2. adrenal: DEGs definition across multiple cell types and multiple conditions
###################
###################
###################
###################
###################

 DefaultAssay(adrenal.combined.filter.integrated) <- "RNA"
  adrenal.combined.filter.integrated$celltype.stim <- paste(Idents(adrenal.combined.filter.integrated), adrenal.combined.filter.integrated$condition, sep = "_")
  adrenal.combined.filter.integrated$celltype <- Idents(adrenal.combined.filter.integrated)
  Idents(adrenal.combined.filter.integrated) <- "celltype.stim"
  adrenal<-adrenal.combined.filter.integrated
  table(adrenal$celltype.stim)
  length(table(adrenal$celltype.stim))#89
  length(table(adrenal$cell_name))#23
  #############Cell_type filtering for downstream DEG analysis
  cell.use<-c()
  for (i in 1:length(table(adrenal$cell_name))){
    FD_cell_number<-nrow(adrenal@meta.data[which(adrenal@meta.data$cell_name==names(table(adrenal$cell_name))[i]&adrenal@meta.data$condition=="FD"),])
    FS_cell_number<-nrow(adrenal@meta.data[which(adrenal@meta.data$cell_name==names(table(adrenal$cell_name))[i]&adrenal@meta.data$condition=="FS"),])
    MC_cell_number<-nrow(adrenal@meta.data[which(adrenal@meta.data$cell_name==names(table(adrenal$cell_name))[i]&adrenal@meta.data$condition=="MC"),])
    MS_cell_number<-nrow(adrenal@meta.data[which(adrenal@meta.data$cell_name==names(table(adrenal$cell_name))[i]&adrenal@meta.data$condition=="MS"),])
    cell.use<-c(cell.use,FD_cell_number>2&FS_cell_number>2&MC_cell_number>2&MS_cell_number>2)
  }
  cell.use.data<-data.frame(table(adrenal$cell_name),cell.use)
  celltypes_adrenal_filter<-as.character(unique(cell.use.data[which(cell.use.data$cell.use==TRUE),]$Var1))
  
  VS <- data.frame(compare=c('MSVSFS','FDVSFS','MCVSMS'),
                   ident1=c('_MS','_FD','_MC'),
                   ident2=c('_FS','_FS','_MS'))
  
  ###################Define DEGs across multiple cell types
  deg_adrenal <- data.frame()
  for (j in (1:dim(VS)[1])){
    eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
    for (i in (1:length(celltypes_adrenal_filter))){
      assign(paste0("deg_markers_",celltypes_adrenal_filter[i],"_",format(VS[j,1]),"_adrenal",sep=""),FindMarkers(adrenal, ident.1 = paste0(celltypes_adrenal_filter[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_adrenal_filter[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
      eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_adrenal_filter[i]),"_",format(VS[j,1]),"_adrenal)[1]>0){deg_markers_",format(celltypes_adrenal_filter[i]),"_",format(VS[j,1]),"_adrenal <- subset(deg_markers_",format(celltypes_adrenal_filter[i]),"_",format(VS[j,1]),"_adrenal,p_val_adj<0.05)};
						   if(dim(deg_markers_",format(celltypes_adrenal_filter[i]),"_",format(VS[j,1]),"_adrenal)[1]>0){deg_markers_",format(celltypes_adrenal_filter[i]),"_",format(VS[j,1]),"_adrenal$Gene_symbol <- row.names(deg_markers_",format(celltypes_adrenal_filter[i]),"_",format(VS[j,1]),"_adrenal);
						   deg_markers_",format(celltypes_adrenal_filter[i]),"_",format(VS[j,1]),"_adrenal$Cell_type <- '",format(celltypes_adrenal_filter[i]),"';
						   deg_markers_",format(celltypes_adrenal_filter[i]),"_",format(VS[j,1]),"_adrenal$Change <- ifelse(deg_markers_",format(celltypes_adrenal_filter[i]),"_",format(VS[j,1]),"_adrenal$avg_log2FC > 0 ,'up','down');
						   deg_markers_",format(celltypes_adrenal_filter[i]),"_",format(VS[j,1]),"_adrenal$Organ <- 'adrenal';
						   deg_markers_",format(celltypes_adrenal_filter[i]),"_",format(VS[j,1]),"_adrenal$Compare <- '",format(VS[j,1]),"';",
                               format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_adrenal_filter[i]),"_",format(VS[j,1]),"_adrenal)};
						   rm(deg_markers_",format(celltypes_adrenal_filter[i]),"_",format(VS[j,1]),"_adrenal)")))
    }
    gc()
    eval(parse(text = paste0("deg_adrenal <- bind_rows(deg_adrenal,",format(VS[j,1]),");
                         rm(",format(VS[j,1]),")")))
  }
  dim(deg_adrenal)
  write.csv(deg_adrenal,'E:/lifei/androgen_project/20211216_Reanalyses/DEG/DEG_list/DEGs.adrenal.avg_log2FC05.p_val_adj005.csv')

###################3. bone_marrow: DEGs definition across multiple cell types and multiple conditions
###################
###################
###################
###################
###################
  DefaultAssay(bone_marrow.combined.filter.integrated) <- "RNA"
  bone_marrow.combined.filter.integrated$celltype.stim <- paste(Idents(bone_marrow.combined.filter.integrated), bone_marrow.combined.filter.integrated$condition, sep = "_")
  bone_marrow.combined.filter.integrated$celltype <- Idents(bone_marrow.combined.filter.integrated)
  Idents(bone_marrow.combined.filter.integrated) <- "celltype.stim"
  bone_marrow<-bone_marrow.combined.filter.integrated
  table(bone_marrow$celltype.stim)
  length(table(bone_marrow$celltype.stim))#89
  length(table(bone_marrow$cell_name))#23
  #############Cell_type filtering for downstream DEG analysis
  cell.use<-c()
  for (i in 1:length(table(bone_marrow$cell_name))){
    FD_cell_number<-nrow(bone_marrow@meta.data[which(bone_marrow@meta.data$cell_name==names(table(bone_marrow$cell_name))[i]&bone_marrow@meta.data$condition=="FD"),])
    FS_cell_number<-nrow(bone_marrow@meta.data[which(bone_marrow@meta.data$cell_name==names(table(bone_marrow$cell_name))[i]&bone_marrow@meta.data$condition=="FS"),])
    MC_cell_number<-nrow(bone_marrow@meta.data[which(bone_marrow@meta.data$cell_name==names(table(bone_marrow$cell_name))[i]&bone_marrow@meta.data$condition=="MC"),])
    MS_cell_number<-nrow(bone_marrow@meta.data[which(bone_marrow@meta.data$cell_name==names(table(bone_marrow$cell_name))[i]&bone_marrow@meta.data$condition=="MS"),])
    cell.use<-c(cell.use,FD_cell_number>2&FS_cell_number>2&MC_cell_number>2&MS_cell_number>2)
  }
  cell.use.data<-data.frame(table(bone_marrow$cell_name),cell.use)
  celltypes_bone_marrow_filter<-as.character(unique(cell.use.data[which(cell.use.data$cell.use==TRUE),]$Var1))
  
  VS <- data.frame(compare=c('MSVSFS','FDVSFS','MCVSMS'),
                   ident1=c('_MS','_FD','_MC'),
                   ident2=c('_FS','_FS','_MS'))
  
  ###################Define DEGs across multiple cell types
  deg_bone_marrow <- data.frame()
  for (j in (1:dim(VS)[1])){
    eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
    for (i in (1:length(celltypes_bone_marrow_filter))){
      assign(paste0("deg_markers_",celltypes_bone_marrow_filter[i],"_",format(VS[j,1]),"_bone_marrow",sep=""),FindMarkers(bone_marrow, ident.1 = paste0(celltypes_bone_marrow_filter[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_bone_marrow_filter[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
      eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_bone_marrow_filter[i]),"_",format(VS[j,1]),"_bone_marrow)[1]>0){deg_markers_",format(celltypes_bone_marrow_filter[i]),"_",format(VS[j,1]),"_bone_marrow <- subset(deg_markers_",format(celltypes_bone_marrow_filter[i]),"_",format(VS[j,1]),"_bone_marrow,p_val_adj<0.05)};
						   if(dim(deg_markers_",format(celltypes_bone_marrow_filter[i]),"_",format(VS[j,1]),"_bone_marrow)[1]>0){deg_markers_",format(celltypes_bone_marrow_filter[i]),"_",format(VS[j,1]),"_bone_marrow$Gene_symbol <- row.names(deg_markers_",format(celltypes_bone_marrow_filter[i]),"_",format(VS[j,1]),"_bone_marrow);
						   deg_markers_",format(celltypes_bone_marrow_filter[i]),"_",format(VS[j,1]),"_bone_marrow$Cell_type <- '",format(celltypes_bone_marrow_filter[i]),"';
						   deg_markers_",format(celltypes_bone_marrow_filter[i]),"_",format(VS[j,1]),"_bone_marrow$Change <- ifelse(deg_markers_",format(celltypes_bone_marrow_filter[i]),"_",format(VS[j,1]),"_bone_marrow$avg_log2FC > 0 ,'up','down');
						   deg_markers_",format(celltypes_bone_marrow_filter[i]),"_",format(VS[j,1]),"_bone_marrow$Organ <- 'bone_marrow';
						   deg_markers_",format(celltypes_bone_marrow_filter[i]),"_",format(VS[j,1]),"_bone_marrow$Compare <- '",format(VS[j,1]),"';",
                               format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_bone_marrow_filter[i]),"_",format(VS[j,1]),"_bone_marrow)};
						   rm(deg_markers_",format(celltypes_bone_marrow_filter[i]),"_",format(VS[j,1]),"_bone_marrow)")))
    }
    gc()
    eval(parse(text = paste0("deg_bone_marrow <- bind_rows(deg_bone_marrow,",format(VS[j,1]),");
                         rm(",format(VS[j,1]),")")))
  }
  dim(deg_bone_marrow)
  write.csv(deg_bone_marrow,'E:/lifei/androgen_project/20211216_Reanalyses/DEG/DEG_list/DEGs.bone_marrow.avg_log2FC05.p_val_adj005.csv')

###################4. brain: DEGs definition across multiple cell types and multiple conditions
###################
###################
###################
###################
###################

  DefaultAssay(brain.combined.filter.integrated) <- "RNA"
  brain.combined.filter.integrated$celltype.stim <- paste(Idents(brain.combined.filter.integrated), brain.combined.filter.integrated$condition, sep = "_")
  brain.combined.filter.integrated$celltype <- Idents(brain.combined.filter.integrated)
  Idents(brain.combined.filter.integrated) <- "celltype.stim"
  brain<-brain.combined.filter.integrated
  table(brain$celltype.stim)
  length(table(brain$celltype.stim))#89
  length(table(brain$cell_name))#23
  #############Cell_type filtering for downstream DEG analysis
  cell.use<-c()
  for (i in 1:length(table(brain$cell_name))){
    FD_cell_number<-nrow(brain@meta.data[which(brain@meta.data$cell_name==names(table(brain$cell_name))[i]&brain@meta.data$condition=="FD"),])
    FS_cell_number<-nrow(brain@meta.data[which(brain@meta.data$cell_name==names(table(brain$cell_name))[i]&brain@meta.data$condition=="FS"),])
    MC_cell_number<-nrow(brain@meta.data[which(brain@meta.data$cell_name==names(table(brain$cell_name))[i]&brain@meta.data$condition=="MC"),])
    MS_cell_number<-nrow(brain@meta.data[which(brain@meta.data$cell_name==names(table(brain$cell_name))[i]&brain@meta.data$condition=="MS"),])
    cell.use<-c(cell.use,FD_cell_number>2&FS_cell_number>2&MC_cell_number>2&MS_cell_number>2)
  }
  cell.use.data<-data.frame(table(brain$cell_name),cell.use)
  celltypes_brain_filter<-as.character(unique(cell.use.data[which(cell.use.data$cell.use==TRUE),]$Var1))
  
  VS <- data.frame(compare=c('MSVSFS','FDVSFS','MCVSMS'),
                   ident1=c('_MS','_FD','_MC'),
                   ident2=c('_FS','_FS','_MS'))
  
  ###################Define DEGs across multiple cell types
  deg_brain <- data.frame()
  for (j in (1:dim(VS)[1])){
    eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
    for (i in (1:length(celltypes_brain_filter))){
      assign(paste0("deg_markers_",celltypes_brain_filter[i],"_",format(VS[j,1]),"_brain",sep=""),FindMarkers(brain, ident.1 = paste0(celltypes_brain_filter[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_brain_filter[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
      eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_brain_filter[i]),"_",format(VS[j,1]),"_brain)[1]>0){deg_markers_",format(celltypes_brain_filter[i]),"_",format(VS[j,1]),"_brain <- subset(deg_markers_",format(celltypes_brain_filter[i]),"_",format(VS[j,1]),"_brain,p_val_adj<0.05)};
						   if(dim(deg_markers_",format(celltypes_brain_filter[i]),"_",format(VS[j,1]),"_brain)[1]>0){deg_markers_",format(celltypes_brain_filter[i]),"_",format(VS[j,1]),"_brain$Gene_symbol <- row.names(deg_markers_",format(celltypes_brain_filter[i]),"_",format(VS[j,1]),"_brain);
						   deg_markers_",format(celltypes_brain_filter[i]),"_",format(VS[j,1]),"_brain$Cell_type <- '",format(celltypes_brain_filter[i]),"';
						   deg_markers_",format(celltypes_brain_filter[i]),"_",format(VS[j,1]),"_brain$Change <- ifelse(deg_markers_",format(celltypes_brain_filter[i]),"_",format(VS[j,1]),"_brain$avg_log2FC > 0 ,'up','down');
						   deg_markers_",format(celltypes_brain_filter[i]),"_",format(VS[j,1]),"_brain$Organ <- 'brain';
						   deg_markers_",format(celltypes_brain_filter[i]),"_",format(VS[j,1]),"_brain$Compare <- '",format(VS[j,1]),"';",
                               format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_brain_filter[i]),"_",format(VS[j,1]),"_brain)};
						   rm(deg_markers_",format(celltypes_brain_filter[i]),"_",format(VS[j,1]),"_brain)")))
    }
    gc()
    eval(parse(text = paste0("deg_brain <- bind_rows(deg_brain,",format(VS[j,1]),");
                         rm(",format(VS[j,1]),")")))
  }
  dim(deg_brain)
  write.csv(deg_brain,'E:/lifei/androgen_project/20211216_Reanalyses/DEG/DEG_list/DEGs.brain.avg_log2FC05.p_val_adj005.csv')

###################5. colon: DEGs definition across multiple cell types and multiple conditions
###################
###################
###################
###################
###################

  DefaultAssay(colon.combined.filter.integrated) <- "RNA"
  colon.combined.filter.integrated$celltype.stim <- paste(Idents(colon.combined.filter.integrated), colon.combined.filter.integrated$condition, sep = "_")
  colon.combined.filter.integrated$celltype <- Idents(colon.combined.filter.integrated)
  Idents(colon.combined.filter.integrated) <- "celltype.stim"
  colon<-colon.combined.filter.integrated
  table(colon$celltype.stim)
  length(table(colon$celltype.stim))#89
  length(table(colon$cell_name))#23
  #############Cell_type filtering for downstream DEG analysis
  cell.use<-c()
  for (i in 1:length(table(colon$cell_name))){
    FD_cell_number<-nrow(colon@meta.data[which(colon@meta.data$cell_name==names(table(colon$cell_name))[i]&colon@meta.data$condition=="FD"),])
    FS_cell_number<-nrow(colon@meta.data[which(colon@meta.data$cell_name==names(table(colon$cell_name))[i]&colon@meta.data$condition=="FS"),])
    MC_cell_number<-nrow(colon@meta.data[which(colon@meta.data$cell_name==names(table(colon$cell_name))[i]&colon@meta.data$condition=="MC"),])
    MS_cell_number<-nrow(colon@meta.data[which(colon@meta.data$cell_name==names(table(colon$cell_name))[i]&colon@meta.data$condition=="MS"),])
    cell.use<-c(cell.use,FD_cell_number>2&FS_cell_number>2&MC_cell_number>2&MS_cell_number>2)
  }
  cell.use.data<-data.frame(table(colon$cell_name),cell.use)
  celltypes_colon_filter<-as.character(unique(cell.use.data[which(cell.use.data$cell.use==TRUE),]$Var1))
  
  VS <- data.frame(compare=c('MSVSFS','FDVSFS','MCVSMS'),
                   ident1=c('_MS','_FD','_MC'),
                   ident2=c('_FS','_FS','_MS'))
  
  ###################Define DEGs across multiple cell types
  deg_colon <- data.frame()
  for (j in (1:dim(VS)[1])){
    eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
    for (i in (1:length(celltypes_colon_filter))){
      assign(paste0("deg_markers_",celltypes_colon_filter[i],"_",format(VS[j,1]),"_colon",sep=""),FindMarkers(colon, ident.1 = paste0(celltypes_colon_filter[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_colon_filter[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
      eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_colon_filter[i]),"_",format(VS[j,1]),"_colon)[1]>0){deg_markers_",format(celltypes_colon_filter[i]),"_",format(VS[j,1]),"_colon <- subset(deg_markers_",format(celltypes_colon_filter[i]),"_",format(VS[j,1]),"_colon,p_val_adj<0.05)};
						   if(dim(deg_markers_",format(celltypes_colon_filter[i]),"_",format(VS[j,1]),"_colon)[1]>0){deg_markers_",format(celltypes_colon_filter[i]),"_",format(VS[j,1]),"_colon$Gene_symbol <- row.names(deg_markers_",format(celltypes_colon_filter[i]),"_",format(VS[j,1]),"_colon);
						   deg_markers_",format(celltypes_colon_filter[i]),"_",format(VS[j,1]),"_colon$Cell_type <- '",format(celltypes_colon_filter[i]),"';
						   deg_markers_",format(celltypes_colon_filter[i]),"_",format(VS[j,1]),"_colon$Change <- ifelse(deg_markers_",format(celltypes_colon_filter[i]),"_",format(VS[j,1]),"_colon$avg_log2FC > 0 ,'up','down');
						   deg_markers_",format(celltypes_colon_filter[i]),"_",format(VS[j,1]),"_colon$Organ <- 'colon';
						   deg_markers_",format(celltypes_colon_filter[i]),"_",format(VS[j,1]),"_colon$Compare <- '",format(VS[j,1]),"';",
                               format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_colon_filter[i]),"_",format(VS[j,1]),"_colon)};
						   rm(deg_markers_",format(celltypes_colon_filter[i]),"_",format(VS[j,1]),"_colon)")))
    }
    gc()
    eval(parse(text = paste0("deg_colon <- bind_rows(deg_colon,",format(VS[j,1]),");
                         rm(",format(VS[j,1]),")")))
  }
  dim(deg_colon)
  write.csv(deg_colon,'E:/lifei/androgen_project/20211216_Reanalyses/DEG/DEG_list/DEGs.colon.avg_log2FC05.p_val_adj005.csv')

###################6. heart: DEGs definition across multiple cell types and multiple conditions
###################
###################
###################
###################
###################

  DefaultAssay(heart.combined.filter.integrated) <- "RNA"
  heart.combined.filter.integrated$celltype.stim <- paste(Idents(heart.combined.filter.integrated), heart.combined.filter.integrated$condition, sep = "_")
  heart.combined.filter.integrated$celltype <- Idents(heart.combined.filter.integrated)
  Idents(heart.combined.filter.integrated) <- "celltype.stim"
  heart<-heart.combined.filter.integrated
  table(heart$celltype.stim)
  length(table(heart$celltype.stim))#89
  length(table(heart$cell_name))#23
  #############Cell_type filtering for downstream DEG analysis
  cell.use<-c()
  for (i in 1:length(table(heart$cell_name))){
    FD_cell_number<-nrow(heart@meta.data[which(heart@meta.data$cell_name==names(table(heart$cell_name))[i]&heart@meta.data$condition=="FD"),])
    FS_cell_number<-nrow(heart@meta.data[which(heart@meta.data$cell_name==names(table(heart$cell_name))[i]&heart@meta.data$condition=="FS"),])
    MC_cell_number<-nrow(heart@meta.data[which(heart@meta.data$cell_name==names(table(heart$cell_name))[i]&heart@meta.data$condition=="MC"),])
    MS_cell_number<-nrow(heart@meta.data[which(heart@meta.data$cell_name==names(table(heart$cell_name))[i]&heart@meta.data$condition=="MS"),])
    cell.use<-c(cell.use,FD_cell_number>2&FS_cell_number>2&MC_cell_number>2&MS_cell_number>2)
  }
  cell.use.data<-data.frame(table(heart$cell_name),cell.use)
  celltypes_heart_filter<-as.character(unique(cell.use.data[which(cell.use.data$cell.use==TRUE),]$Var1))
  
  VS <- data.frame(compare=c('MSVSFS','FDVSFS','MCVSMS'),
                   ident1=c('_MS','_FD','_MC'),
                   ident2=c('_FS','_FS','_MS'))
  
  ###################Define DEGs across multiple cell types
  deg_heart <- data.frame()
  for (j in (1:dim(VS)[1])){
    eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
    for (i in (1:length(celltypes_heart_filter))){
      assign(paste0("deg_markers_",celltypes_heart_filter[i],"_",format(VS[j,1]),"_heart",sep=""),FindMarkers(heart, ident.1 = paste0(celltypes_heart_filter[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_heart_filter[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
      eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_heart_filter[i]),"_",format(VS[j,1]),"_heart)[1]>0){deg_markers_",format(celltypes_heart_filter[i]),"_",format(VS[j,1]),"_heart <- subset(deg_markers_",format(celltypes_heart_filter[i]),"_",format(VS[j,1]),"_heart,p_val_adj<0.05)};
						   if(dim(deg_markers_",format(celltypes_heart_filter[i]),"_",format(VS[j,1]),"_heart)[1]>0){deg_markers_",format(celltypes_heart_filter[i]),"_",format(VS[j,1]),"_heart$Gene_symbol <- row.names(deg_markers_",format(celltypes_heart_filter[i]),"_",format(VS[j,1]),"_heart);
						   deg_markers_",format(celltypes_heart_filter[i]),"_",format(VS[j,1]),"_heart$Cell_type <- '",format(celltypes_heart_filter[i]),"';
						   deg_markers_",format(celltypes_heart_filter[i]),"_",format(VS[j,1]),"_heart$Change <- ifelse(deg_markers_",format(celltypes_heart_filter[i]),"_",format(VS[j,1]),"_heart$avg_log2FC > 0 ,'up','down');
						   deg_markers_",format(celltypes_heart_filter[i]),"_",format(VS[j,1]),"_heart$Organ <- 'heart';
						   deg_markers_",format(celltypes_heart_filter[i]),"_",format(VS[j,1]),"_heart$Compare <- '",format(VS[j,1]),"';",
                               format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_heart_filter[i]),"_",format(VS[j,1]),"_heart)};
						   rm(deg_markers_",format(celltypes_heart_filter[i]),"_",format(VS[j,1]),"_heart)")))
    }
    gc()
    eval(parse(text = paste0("deg_heart <- bind_rows(deg_heart,",format(VS[j,1]),");
                         rm(",format(VS[j,1]),")")))
  }
  dim(deg_heart)
  write.csv(deg_heart,'E:/lifei/androgen_project/20211216_Reanalyses/DEG/DEG_list/DEGs.heart.avg_log2FC05.p_val_adj005.csv')

###################7. intestine: DEGs definition across multiple cell types and multiple conditions
###################
###################
###################
###################
###################

  DefaultAssay(intestine.combined.filter.integrated) <- "RNA"
  intestine.combined.filter.integrated$celltype.stim <- paste(Idents(intestine.combined.filter.integrated), intestine.combined.filter.integrated$condition, sep = "_")
  intestine.combined.filter.integrated$celltype <- Idents(intestine.combined.filter.integrated)
  Idents(intestine.combined.filter.integrated) <- "celltype.stim"
  intestine<-intestine.combined.filter.integrated
  table(intestine$celltype.stim)
  length(table(intestine$celltype.stim))#89
  length(table(intestine$cell_name))#23
  #############Cell_type filtering for downstream DEG analysis
  cell.use<-c()
  for (i in 1:length(table(intestine$cell_name))){
    FD_cell_number<-nrow(intestine@meta.data[which(intestine@meta.data$cell_name==names(table(intestine$cell_name))[i]&intestine@meta.data$condition=="FD"),])
    FS_cell_number<-nrow(intestine@meta.data[which(intestine@meta.data$cell_name==names(table(intestine$cell_name))[i]&intestine@meta.data$condition=="FS"),])
    MC_cell_number<-nrow(intestine@meta.data[which(intestine@meta.data$cell_name==names(table(intestine$cell_name))[i]&intestine@meta.data$condition=="MC"),])
    MS_cell_number<-nrow(intestine@meta.data[which(intestine@meta.data$cell_name==names(table(intestine$cell_name))[i]&intestine@meta.data$condition=="MS"),])
    cell.use<-c(cell.use,FD_cell_number>2&FS_cell_number>2&MC_cell_number>2&MS_cell_number>2)
  }
  cell.use.data<-data.frame(table(intestine$cell_name),cell.use)
  celltypes_intestine_filter<-as.character(unique(cell.use.data[which(cell.use.data$cell.use==TRUE),]$Var1))
  
  VS <- data.frame(compare=c('MSVSFS','FDVSFS','MCVSMS'),
                   ident1=c('_MS','_FD','_MC'),
                   ident2=c('_FS','_FS','_MS'))
  
  ###################Define DEGs across multiple cell types
  deg_intestine <- data.frame()
  for (j in (1:dim(VS)[1])){
    eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
    for (i in (1:length(celltypes_intestine_filter))){
      assign(paste0("deg_markers_",celltypes_intestine_filter[i],"_",format(VS[j,1]),"_intestine",sep=""),FindMarkers(intestine, ident.1 = paste0(celltypes_intestine_filter[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_intestine_filter[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
      eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_intestine_filter[i]),"_",format(VS[j,1]),"_intestine)[1]>0){deg_markers_",format(celltypes_intestine_filter[i]),"_",format(VS[j,1]),"_intestine <- subset(deg_markers_",format(celltypes_intestine_filter[i]),"_",format(VS[j,1]),"_intestine,p_val_adj<0.05)};
						   if(dim(deg_markers_",format(celltypes_intestine_filter[i]),"_",format(VS[j,1]),"_intestine)[1]>0){deg_markers_",format(celltypes_intestine_filter[i]),"_",format(VS[j,1]),"_intestine$Gene_symbol <- row.names(deg_markers_",format(celltypes_intestine_filter[i]),"_",format(VS[j,1]),"_intestine);
						   deg_markers_",format(celltypes_intestine_filter[i]),"_",format(VS[j,1]),"_intestine$Cell_type <- '",format(celltypes_intestine_filter[i]),"';
						   deg_markers_",format(celltypes_intestine_filter[i]),"_",format(VS[j,1]),"_intestine$Change <- ifelse(deg_markers_",format(celltypes_intestine_filter[i]),"_",format(VS[j,1]),"_intestine$avg_log2FC > 0 ,'up','down');
						   deg_markers_",format(celltypes_intestine_filter[i]),"_",format(VS[j,1]),"_intestine$Organ <- 'intestine';
						   deg_markers_",format(celltypes_intestine_filter[i]),"_",format(VS[j,1]),"_intestine$Compare <- '",format(VS[j,1]),"';",
                               format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_intestine_filter[i]),"_",format(VS[j,1]),"_intestine)};
						   rm(deg_markers_",format(celltypes_intestine_filter[i]),"_",format(VS[j,1]),"_intestine)")))
    }
    gc()
    eval(parse(text = paste0("deg_intestine <- bind_rows(deg_intestine,",format(VS[j,1]),");
                         rm(",format(VS[j,1]),")")))
  }
  dim(deg_intestine)
  write.csv(deg_intestine,'E:/lifei/androgen_project/20211216_Reanalyses/DEG/DEG_list/DEGs.intestine.avg_log2FC05.p_val_adj005.csv')

###################8. kidney: DEGs definition across multiple cell types and multiple conditions
###################
###################
###################
###################
###################

  DefaultAssay(kidney.combined.filter.integrated) <- "RNA"
  kidney.combined.filter.integrated$celltype.stim <- paste(Idents(kidney.combined.filter.integrated), kidney.combined.filter.integrated$condition, sep = "_")
  kidney.combined.filter.integrated$celltype <- Idents(kidney.combined.filter.integrated)
  Idents(kidney.combined.filter.integrated) <- "celltype.stim"
  kidney<-kidney.combined.filter.integrated
  table(kidney$celltype.stim)
  length(table(kidney$celltype.stim))#89
  length(table(kidney$cell_name))#23
  #############Cell_type filtering for downstream DEG analysis
  cell.use<-c()
  for (i in 1:length(table(kidney$cell_name))){
    FD_cell_number<-nrow(kidney@meta.data[which(kidney@meta.data$cell_name==names(table(kidney$cell_name))[i]&kidney@meta.data$condition=="FD"),])
    FS_cell_number<-nrow(kidney@meta.data[which(kidney@meta.data$cell_name==names(table(kidney$cell_name))[i]&kidney@meta.data$condition=="FS"),])
    MC_cell_number<-nrow(kidney@meta.data[which(kidney@meta.data$cell_name==names(table(kidney$cell_name))[i]&kidney@meta.data$condition=="MC"),])
    MS_cell_number<-nrow(kidney@meta.data[which(kidney@meta.data$cell_name==names(table(kidney$cell_name))[i]&kidney@meta.data$condition=="MS"),])
    cell.use<-c(cell.use,FD_cell_number>2&FS_cell_number>2&MC_cell_number>2&MS_cell_number>2)
  }
  cell.use.data<-data.frame(table(kidney$cell_name),cell.use)
  celltypes_kidney_filter<-as.character(unique(cell.use.data[which(cell.use.data$cell.use==TRUE),]$Var1))
  
  VS <- data.frame(compare=c('MSVSFS','FDVSFS','MCVSMS'),
                   ident1=c('_MS','_FD','_MC'),
                   ident2=c('_FS','_FS','_MS'))
  
  ###################Define DEGs across multiple cell types
  deg_kidney <- data.frame()
  for (j in (1:dim(VS)[1])){
    eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
    for (i in (1:length(celltypes_kidney_filter))){
      assign(paste0("deg_markers_",celltypes_kidney_filter[i],"_",format(VS[j,1]),"_kidney",sep=""),FindMarkers(kidney, ident.1 = paste0(celltypes_kidney_filter[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_kidney_filter[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
      eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_kidney_filter[i]),"_",format(VS[j,1]),"_kidney)[1]>0){deg_markers_",format(celltypes_kidney_filter[i]),"_",format(VS[j,1]),"_kidney <- subset(deg_markers_",format(celltypes_kidney_filter[i]),"_",format(VS[j,1]),"_kidney,p_val_adj<0.05)};
						   if(dim(deg_markers_",format(celltypes_kidney_filter[i]),"_",format(VS[j,1]),"_kidney)[1]>0){deg_markers_",format(celltypes_kidney_filter[i]),"_",format(VS[j,1]),"_kidney$Gene_symbol <- row.names(deg_markers_",format(celltypes_kidney_filter[i]),"_",format(VS[j,1]),"_kidney);
						   deg_markers_",format(celltypes_kidney_filter[i]),"_",format(VS[j,1]),"_kidney$Cell_type <- '",format(celltypes_kidney_filter[i]),"';
						   deg_markers_",format(celltypes_kidney_filter[i]),"_",format(VS[j,1]),"_kidney$Change <- ifelse(deg_markers_",format(celltypes_kidney_filter[i]),"_",format(VS[j,1]),"_kidney$avg_log2FC > 0 ,'up','down');
						   deg_markers_",format(celltypes_kidney_filter[i]),"_",format(VS[j,1]),"_kidney$Organ <- 'kidney';
						   deg_markers_",format(celltypes_kidney_filter[i]),"_",format(VS[j,1]),"_kidney$Compare <- '",format(VS[j,1]),"';",
                               format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_kidney_filter[i]),"_",format(VS[j,1]),"_kidney)};
						   rm(deg_markers_",format(celltypes_kidney_filter[i]),"_",format(VS[j,1]),"_kidney)")))
    }
    gc()
    eval(parse(text = paste0("deg_kidney <- bind_rows(deg_kidney,",format(VS[j,1]),");
                         rm(",format(VS[j,1]),")")))
  }
  dim(deg_kidney)
  write.csv(deg_kidney,'E:/lifei/androgen_project/20211216_Reanalyses/DEG/DEG_list/DEGs.kidney.avg_log2FC05.p_val_adj005.csv')

###################9. lacrimal: DEGs definition across multiple cell types and multiple conditions
###################
###################
###################
###################
###################

  DefaultAssay(lacrimal.combined.filter.integrated) <- "RNA"
  lacrimal.combined.filter.integrated$celltype.stim <- paste(Idents(lacrimal.combined.filter.integrated), lacrimal.combined.filter.integrated$condition, sep = "_")
  lacrimal.combined.filter.integrated$celltype <- Idents(lacrimal.combined.filter.integrated)
  Idents(lacrimal.combined.filter.integrated) <- "celltype.stim"
  lacrimal<-lacrimal.combined.filter.integrated
  table(lacrimal$celltype.stim)
  length(table(lacrimal$celltype.stim))#89
  length(table(lacrimal$cell_name))#23
  #############Cell_type filtering for downstream DEG analysis
  cell.use<-c()
  for (i in 1:length(table(lacrimal$cell_name))){
    FD_cell_number<-nrow(lacrimal@meta.data[which(lacrimal@meta.data$cell_name==names(table(lacrimal$cell_name))[i]&lacrimal@meta.data$condition=="FD"),])
    FS_cell_number<-nrow(lacrimal@meta.data[which(lacrimal@meta.data$cell_name==names(table(lacrimal$cell_name))[i]&lacrimal@meta.data$condition=="FS"),])
    MC_cell_number<-nrow(lacrimal@meta.data[which(lacrimal@meta.data$cell_name==names(table(lacrimal$cell_name))[i]&lacrimal@meta.data$condition=="MC"),])
    MS_cell_number<-nrow(lacrimal@meta.data[which(lacrimal@meta.data$cell_name==names(table(lacrimal$cell_name))[i]&lacrimal@meta.data$condition=="MS"),])
    cell.use<-c(cell.use,FD_cell_number>2&FS_cell_number>2&MC_cell_number>2&MS_cell_number>2)
  }
  cell.use.data<-data.frame(table(lacrimal$cell_name),cell.use)
  celltypes_lacrimal_filter<-as.character(unique(cell.use.data[which(cell.use.data$cell.use==TRUE),]$Var1))
  
  VS <- data.frame(compare=c('MSVSFS','FDVSFS','MCVSMS'),
                   ident1=c('_MS','_FD','_MC'),
                   ident2=c('_FS','_FS','_MS'))
  
  ###################Define DEGs across multiple cell types
  deg_lacrimal <- data.frame()
  for (j in (1:dim(VS)[1])){
    eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
    for (i in (1:length(celltypes_lacrimal_filter))){
      assign(paste0("deg_markers_",celltypes_lacrimal_filter[i],"_",format(VS[j,1]),"_lacrimal",sep=""),FindMarkers(lacrimal, ident.1 = paste0(celltypes_lacrimal_filter[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_lacrimal_filter[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
      eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_lacrimal_filter[i]),"_",format(VS[j,1]),"_lacrimal)[1]>0){deg_markers_",format(celltypes_lacrimal_filter[i]),"_",format(VS[j,1]),"_lacrimal <- subset(deg_markers_",format(celltypes_lacrimal_filter[i]),"_",format(VS[j,1]),"_lacrimal,p_val_adj<0.05)};
						   if(dim(deg_markers_",format(celltypes_lacrimal_filter[i]),"_",format(VS[j,1]),"_lacrimal)[1]>0){deg_markers_",format(celltypes_lacrimal_filter[i]),"_",format(VS[j,1]),"_lacrimal$Gene_symbol <- row.names(deg_markers_",format(celltypes_lacrimal_filter[i]),"_",format(VS[j,1]),"_lacrimal);
						   deg_markers_",format(celltypes_lacrimal_filter[i]),"_",format(VS[j,1]),"_lacrimal$Cell_type <- '",format(celltypes_lacrimal_filter[i]),"';
						   deg_markers_",format(celltypes_lacrimal_filter[i]),"_",format(VS[j,1]),"_lacrimal$Change <- ifelse(deg_markers_",format(celltypes_lacrimal_filter[i]),"_",format(VS[j,1]),"_lacrimal$avg_log2FC > 0 ,'up','down');
						   deg_markers_",format(celltypes_lacrimal_filter[i]),"_",format(VS[j,1]),"_lacrimal$Organ <- 'lacrimal';
						   deg_markers_",format(celltypes_lacrimal_filter[i]),"_",format(VS[j,1]),"_lacrimal$Compare <- '",format(VS[j,1]),"';",
                               format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_lacrimal_filter[i]),"_",format(VS[j,1]),"_lacrimal)};
						   rm(deg_markers_",format(celltypes_lacrimal_filter[i]),"_",format(VS[j,1]),"_lacrimal)")))
    }
    gc()
    eval(parse(text = paste0("deg_lacrimal <- bind_rows(deg_lacrimal,",format(VS[j,1]),");
                         rm(",format(VS[j,1]),")")))
  }
  dim(deg_lacrimal)
  write.csv(deg_lacrimal,'E:/lifei/androgen_project/20211216_Reanalyses/DEG/DEG_list/DEGs.lacrimal.avg_log2FC05.p_val_adj005.csv')

###################10. liver: DEGs definition across multiple cell types and multiple conditions
###################
###################
###################
###################
###################

  DefaultAssay(liver.combined.filter.integrated) <- "RNA"
  liver.combined.filter.integrated$celltype.stim <- paste(Idents(liver.combined.filter.integrated), liver.combined.filter.integrated$condition, sep = "_")
  liver.combined.filter.integrated$celltype <- Idents(liver.combined.filter.integrated)
  Idents(liver.combined.filter.integrated) <- "celltype.stim"
  liver<-liver.combined.filter.integrated
  table(liver$celltype.stim)
  length(table(liver$celltype.stim))#89
  length(table(liver$cell_name))#23
  #############Cell_type filtering for downstream DEG analysis
  cell.use<-c()
  for (i in 1:length(table(liver$cell_name))){
    FD_cell_number<-nrow(liver@meta.data[which(liver@meta.data$cell_name==names(table(liver$cell_name))[i]&liver@meta.data$condition=="FD"),])
    FS_cell_number<-nrow(liver@meta.data[which(liver@meta.data$cell_name==names(table(liver$cell_name))[i]&liver@meta.data$condition=="FS"),])
    MC_cell_number<-nrow(liver@meta.data[which(liver@meta.data$cell_name==names(table(liver$cell_name))[i]&liver@meta.data$condition=="MC"),])
    MS_cell_number<-nrow(liver@meta.data[which(liver@meta.data$cell_name==names(table(liver$cell_name))[i]&liver@meta.data$condition=="MS"),])
    cell.use<-c(cell.use,FD_cell_number>2&FS_cell_number>2&MC_cell_number>2&MS_cell_number>2)
  }
  cell.use.data<-data.frame(table(liver$cell_name),cell.use)
  celltypes_liver_filter<-as.character(unique(cell.use.data[which(cell.use.data$cell.use==TRUE),]$Var1))
  
  VS <- data.frame(compare=c('MSVSFS','FDVSFS','MCVSMS'),
                   ident1=c('_MS','_FD','_MC'),
                   ident2=c('_FS','_FS','_MS'))
  
  ###################Define DEGs across multiple cell types
  deg_liver <- data.frame()
  for (j in (1:dim(VS)[1])){
    eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
    for (i in (1:length(celltypes_liver_filter))){
      assign(paste0("deg_markers_",celltypes_liver_filter[i],"_",format(VS[j,1]),"_liver",sep=""),FindMarkers(liver, ident.1 = paste0(celltypes_liver_filter[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_liver_filter[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
      eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_liver_filter[i]),"_",format(VS[j,1]),"_liver)[1]>0){deg_markers_",format(celltypes_liver_filter[i]),"_",format(VS[j,1]),"_liver <- subset(deg_markers_",format(celltypes_liver_filter[i]),"_",format(VS[j,1]),"_liver,p_val_adj<0.05)};
						   if(dim(deg_markers_",format(celltypes_liver_filter[i]),"_",format(VS[j,1]),"_liver)[1]>0){deg_markers_",format(celltypes_liver_filter[i]),"_",format(VS[j,1]),"_liver$Gene_symbol <- row.names(deg_markers_",format(celltypes_liver_filter[i]),"_",format(VS[j,1]),"_liver);
						   deg_markers_",format(celltypes_liver_filter[i]),"_",format(VS[j,1]),"_liver$Cell_type <- '",format(celltypes_liver_filter[i]),"';
						   deg_markers_",format(celltypes_liver_filter[i]),"_",format(VS[j,1]),"_liver$Change <- ifelse(deg_markers_",format(celltypes_liver_filter[i]),"_",format(VS[j,1]),"_liver$avg_log2FC > 0 ,'up','down');
						   deg_markers_",format(celltypes_liver_filter[i]),"_",format(VS[j,1]),"_liver$Organ <- 'liver';
						   deg_markers_",format(celltypes_liver_filter[i]),"_",format(VS[j,1]),"_liver$Compare <- '",format(VS[j,1]),"';",
                               format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_liver_filter[i]),"_",format(VS[j,1]),"_liver)};
						   rm(deg_markers_",format(celltypes_liver_filter[i]),"_",format(VS[j,1]),"_liver)")))
    }
    gc()
    eval(parse(text = paste0("deg_liver <- bind_rows(deg_liver,",format(VS[j,1]),");
                         rm(",format(VS[j,1]),")")))
  }
  dim(deg_liver)
  write.csv(deg_liver,'E:/lifei/androgen_project/20211216_Reanalyses/DEG/DEG_list/DEGs.liver.avg_log2FC05.p_val_adj005.csv')

###################11. lung: DEGs definition across multiple cell types and multiple conditions
###################
###################
###################
###################
###################

  DefaultAssay(lung.combined.filter.integrated) <- "RNA"
  lung.combined.filter.integrated$celltype.stim <- paste(Idents(lung.combined.filter.integrated), lung.combined.filter.integrated$condition, sep = "_")
  lung.combined.filter.integrated$celltype <- Idents(lung.combined.filter.integrated)
  Idents(lung.combined.filter.integrated) <- "celltype.stim"
  lung<-lung.combined.filter.integrated
  table(lung$celltype.stim)
  length(table(lung$celltype.stim))#89
  length(table(lung$cell_name))#23
  #############Cell_type filtering for downstream DEG analysis
  cell.use<-c()
  for (i in 1:length(table(lung$cell_name))){
    FD_cell_number<-nrow(lung@meta.data[which(lung@meta.data$cell_name==names(table(lung$cell_name))[i]&lung@meta.data$condition=="FD"),])
    FS_cell_number<-nrow(lung@meta.data[which(lung@meta.data$cell_name==names(table(lung$cell_name))[i]&lung@meta.data$condition=="FS"),])
    MC_cell_number<-nrow(lung@meta.data[which(lung@meta.data$cell_name==names(table(lung$cell_name))[i]&lung@meta.data$condition=="MC"),])
    MS_cell_number<-nrow(lung@meta.data[which(lung@meta.data$cell_name==names(table(lung$cell_name))[i]&lung@meta.data$condition=="MS"),])
    cell.use<-c(cell.use,FD_cell_number>2&FS_cell_number>2&MC_cell_number>2&MS_cell_number>2)
  }
  cell.use.data<-data.frame(table(lung$cell_name),cell.use)
  celltypes_lung_filter<-as.character(unique(cell.use.data[which(cell.use.data$cell.use==TRUE),]$Var1))
  
  VS <- data.frame(compare=c('MSVSFS','FDVSFS','MCVSMS'),
                   ident1=c('_MS','_FD','_MC'),
                   ident2=c('_FS','_FS','_MS'))
  
  ###################Define DEGs across multiple cell types
  deg_lung <- data.frame()
  for (j in (1:dim(VS)[1])){
    eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
    for (i in (1:length(celltypes_lung_filter))){
      assign(paste0("deg_markers_",celltypes_lung_filter[i],"_",format(VS[j,1]),"_lung",sep=""),FindMarkers(lung, ident.1 = paste0(celltypes_lung_filter[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_lung_filter[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
      eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_lung_filter[i]),"_",format(VS[j,1]),"_lung)[1]>0){deg_markers_",format(celltypes_lung_filter[i]),"_",format(VS[j,1]),"_lung <- subset(deg_markers_",format(celltypes_lung_filter[i]),"_",format(VS[j,1]),"_lung,p_val_adj<0.05)};
						   if(dim(deg_markers_",format(celltypes_lung_filter[i]),"_",format(VS[j,1]),"_lung)[1]>0){deg_markers_",format(celltypes_lung_filter[i]),"_",format(VS[j,1]),"_lung$Gene_symbol <- row.names(deg_markers_",format(celltypes_lung_filter[i]),"_",format(VS[j,1]),"_lung);
						   deg_markers_",format(celltypes_lung_filter[i]),"_",format(VS[j,1]),"_lung$Cell_type <- '",format(celltypes_lung_filter[i]),"';
						   deg_markers_",format(celltypes_lung_filter[i]),"_",format(VS[j,1]),"_lung$Change <- ifelse(deg_markers_",format(celltypes_lung_filter[i]),"_",format(VS[j,1]),"_lung$avg_log2FC > 0 ,'up','down');
						   deg_markers_",format(celltypes_lung_filter[i]),"_",format(VS[j,1]),"_lung$Organ <- 'lung';
						   deg_markers_",format(celltypes_lung_filter[i]),"_",format(VS[j,1]),"_lung$Compare <- '",format(VS[j,1]),"';",
                               format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_lung_filter[i]),"_",format(VS[j,1]),"_lung)};
						   rm(deg_markers_",format(celltypes_lung_filter[i]),"_",format(VS[j,1]),"_lung)")))
    }
    gc()
    eval(parse(text = paste0("deg_lung <- bind_rows(deg_lung,",format(VS[j,1]),");
                         rm(",format(VS[j,1]),")")))
  }
  dim(deg_lung)
  write.csv(deg_lung,'E:/lifei/androgen_project/20211216_Reanalyses/DEG/DEG_list/DEGs.lung.avg_log2FC05.p_val_adj005.csv')

###################12. pancreas: DEGs definition across multiple cell types and multiple conditions
###################
###################
###################
###################
###################

  DefaultAssay(pancreas.combined.filter.integrated) <- "RNA"
  pancreas.combined.filter.integrated$celltype.stim <- paste(Idents(pancreas.combined.filter.integrated), pancreas.combined.filter.integrated$condition, sep = "_")
  pancreas.combined.filter.integrated$celltype <- Idents(pancreas.combined.filter.integrated)
  Idents(pancreas.combined.filter.integrated) <- "celltype.stim"
  pancreas<-pancreas.combined.filter.integrated
  table(pancreas$celltype.stim)
  length(table(pancreas$celltype.stim))#89
  length(table(pancreas$cell_name))#23
  #############Cell_type filtering for downstream DEG analysis
  cell.use<-c()
  for (i in 1:length(table(pancreas$cell_name))){
    FD_cell_number<-nrow(pancreas@meta.data[which(pancreas@meta.data$cell_name==names(table(pancreas$cell_name))[i]&pancreas@meta.data$condition=="FD"),])
    FS_cell_number<-nrow(pancreas@meta.data[which(pancreas@meta.data$cell_name==names(table(pancreas$cell_name))[i]&pancreas@meta.data$condition=="FS"),])
    MC_cell_number<-nrow(pancreas@meta.data[which(pancreas@meta.data$cell_name==names(table(pancreas$cell_name))[i]&pancreas@meta.data$condition=="MC"),])
    MS_cell_number<-nrow(pancreas@meta.data[which(pancreas@meta.data$cell_name==names(table(pancreas$cell_name))[i]&pancreas@meta.data$condition=="MS"),])
    cell.use<-c(cell.use,FD_cell_number>2&FS_cell_number>2&MC_cell_number>2&MS_cell_number>2)
  }
  cell.use.data<-data.frame(table(pancreas$cell_name),cell.use)
  celltypes_pancreas_filter<-as.character(unique(cell.use.data[which(cell.use.data$cell.use==TRUE),]$Var1))
  
  VS <- data.frame(compare=c('MSVSFS','FDVSFS','MCVSMS'),
                   ident1=c('_MS','_FD','_MC'),
                   ident2=c('_FS','_FS','_MS'))
  
  ###################Define DEGs across multiple cell types
  deg_pancreas <- data.frame()
  for (j in (1:dim(VS)[1])){
    eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
    for (i in (1:length(celltypes_pancreas_filter))){
      assign(paste0("deg_markers_",celltypes_pancreas_filter[i],"_",format(VS[j,1]),"_pancreas",sep=""),FindMarkers(pancreas, ident.1 = paste0(celltypes_pancreas_filter[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_pancreas_filter[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
      eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_pancreas_filter[i]),"_",format(VS[j,1]),"_pancreas)[1]>0){deg_markers_",format(celltypes_pancreas_filter[i]),"_",format(VS[j,1]),"_pancreas <- subset(deg_markers_",format(celltypes_pancreas_filter[i]),"_",format(VS[j,1]),"_pancreas,p_val_adj<0.05)};
						   if(dim(deg_markers_",format(celltypes_pancreas_filter[i]),"_",format(VS[j,1]),"_pancreas)[1]>0){deg_markers_",format(celltypes_pancreas_filter[i]),"_",format(VS[j,1]),"_pancreas$Gene_symbol <- row.names(deg_markers_",format(celltypes_pancreas_filter[i]),"_",format(VS[j,1]),"_pancreas);
						   deg_markers_",format(celltypes_pancreas_filter[i]),"_",format(VS[j,1]),"_pancreas$Cell_type <- '",format(celltypes_pancreas_filter[i]),"';
						   deg_markers_",format(celltypes_pancreas_filter[i]),"_",format(VS[j,1]),"_pancreas$Change <- ifelse(deg_markers_",format(celltypes_pancreas_filter[i]),"_",format(VS[j,1]),"_pancreas$avg_log2FC > 0 ,'up','down');
						   deg_markers_",format(celltypes_pancreas_filter[i]),"_",format(VS[j,1]),"_pancreas$Organ <- 'pancreas';
						   deg_markers_",format(celltypes_pancreas_filter[i]),"_",format(VS[j,1]),"_pancreas$Compare <- '",format(VS[j,1]),"';",
                               format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_pancreas_filter[i]),"_",format(VS[j,1]),"_pancreas)};
						   rm(deg_markers_",format(celltypes_pancreas_filter[i]),"_",format(VS[j,1]),"_pancreas)")))
    }
    gc()
    eval(parse(text = paste0("deg_pancreas <- bind_rows(deg_pancreas,",format(VS[j,1]),");
                         rm(",format(VS[j,1]),")")))
  }
  dim(deg_pancreas)
  write.csv(deg_pancreas,'E:/lifei/androgen_project/20211216_Reanalyses/DEG/DEG_list/DEGs.pancreas.avg_log2FC05.p_val_adj005.csv')

###################13. salivary: DEGs definition across multiple cell types and multiple conditions
###################
###################
###################
###################
###################

  DefaultAssay(salivary.combined.filter.integrated) <- "RNA"
  salivary.combined.filter.integrated$celltype.stim <- paste(Idents(salivary.combined.filter.integrated), salivary.combined.filter.integrated$condition, sep = "_")
  salivary.combined.filter.integrated$celltype <- Idents(salivary.combined.filter.integrated)
  Idents(salivary.combined.filter.integrated) <- "celltype.stim"
  salivary<-salivary.combined.filter.integrated
  table(salivary$celltype.stim)
  length(table(salivary$celltype.stim))#89
  length(table(salivary$cell_name))#23
  #############Cell_type filtering for downstream DEG analysis
  cell.use<-c()
  for (i in 1:length(table(salivary$cell_name))){
    FD_cell_number<-nrow(salivary@meta.data[which(salivary@meta.data$cell_name==names(table(salivary$cell_name))[i]&salivary@meta.data$condition=="FD"),])
    FS_cell_number<-nrow(salivary@meta.data[which(salivary@meta.data$cell_name==names(table(salivary$cell_name))[i]&salivary@meta.data$condition=="FS"),])
    MC_cell_number<-nrow(salivary@meta.data[which(salivary@meta.data$cell_name==names(table(salivary$cell_name))[i]&salivary@meta.data$condition=="MC"),])
    MS_cell_number<-nrow(salivary@meta.data[which(salivary@meta.data$cell_name==names(table(salivary$cell_name))[i]&salivary@meta.data$condition=="MS"),])
    cell.use<-c(cell.use,FD_cell_number>2&FS_cell_number>2&MC_cell_number>2&MS_cell_number>2)
  }
  cell.use.data<-data.frame(table(salivary$cell_name),cell.use)
  celltypes_salivary_filter<-as.character(unique(cell.use.data[which(cell.use.data$cell.use==TRUE),]$Var1))
  
  VS <- data.frame(compare=c('MSVSFS','FDVSFS','MCVSMS'),
                   ident1=c('_MS','_FD','_MC'),
                   ident2=c('_FS','_FS','_MS'))
  
  ###################Define DEGs across multiple cell types
  deg_salivary <- data.frame()
  for (j in (1:dim(VS)[1])){
    eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
    for (i in (1:length(celltypes_salivary_filter))){
      assign(paste0("deg_markers_",celltypes_salivary_filter[i],"_",format(VS[j,1]),"_salivary",sep=""),FindMarkers(salivary, ident.1 = paste0(celltypes_salivary_filter[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_salivary_filter[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
      eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_salivary_filter[i]),"_",format(VS[j,1]),"_salivary)[1]>0){deg_markers_",format(celltypes_salivary_filter[i]),"_",format(VS[j,1]),"_salivary <- subset(deg_markers_",format(celltypes_salivary_filter[i]),"_",format(VS[j,1]),"_salivary,p_val_adj<0.05)};
						   if(dim(deg_markers_",format(celltypes_salivary_filter[i]),"_",format(VS[j,1]),"_salivary)[1]>0){deg_markers_",format(celltypes_salivary_filter[i]),"_",format(VS[j,1]),"_salivary$Gene_symbol <- row.names(deg_markers_",format(celltypes_salivary_filter[i]),"_",format(VS[j,1]),"_salivary);
						   deg_markers_",format(celltypes_salivary_filter[i]),"_",format(VS[j,1]),"_salivary$Cell_type <- '",format(celltypes_salivary_filter[i]),"';
						   deg_markers_",format(celltypes_salivary_filter[i]),"_",format(VS[j,1]),"_salivary$Change <- ifelse(deg_markers_",format(celltypes_salivary_filter[i]),"_",format(VS[j,1]),"_salivary$avg_log2FC > 0 ,'up','down');
						   deg_markers_",format(celltypes_salivary_filter[i]),"_",format(VS[j,1]),"_salivary$Organ <- 'salivary';
						   deg_markers_",format(celltypes_salivary_filter[i]),"_",format(VS[j,1]),"_salivary$Compare <- '",format(VS[j,1]),"';",
                               format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_salivary_filter[i]),"_",format(VS[j,1]),"_salivary)};
						   rm(deg_markers_",format(celltypes_salivary_filter[i]),"_",format(VS[j,1]),"_salivary)")))
    }
    gc()
    eval(parse(text = paste0("deg_salivary <- bind_rows(deg_salivary,",format(VS[j,1]),");
                         rm(",format(VS[j,1]),")")))
  }
  dim(deg_salivary)
  write.csv(deg_salivary,'E:/lifei/androgen_project/20211216_Reanalyses/DEG/DEG_list/DEGs.salivary.avg_log2FC05.p_val_adj005.csv')

###################14. skeletalmuscle: DEGs definition across multiple cell types and multiple conditions
###################
###################
###################
###################
###################

  DefaultAssay(skeletalmuscle.combined.filter.integrated) <- "RNA"
  skeletalmuscle.combined.filter.integrated$celltype.stim <- paste(Idents(skeletalmuscle.combined.filter.integrated), skeletalmuscle.combined.filter.integrated$condition, sep = "_")
  skeletalmuscle.combined.filter.integrated$celltype <- Idents(skeletalmuscle.combined.filter.integrated)
  Idents(skeletalmuscle.combined.filter.integrated) <- "celltype.stim"
  skeletalmuscle<-skeletalmuscle.combined.filter.integrated
  table(skeletalmuscle$celltype.stim)
  length(table(skeletalmuscle$celltype.stim))#89
  length(table(skeletalmuscle$cell_name))#23
  #############Cell_type filtering for downstream DEG analysis
  cell.use<-c()
  for (i in 1:length(table(skeletalmuscle$cell_name))){
    FD_cell_number<-nrow(skeletalmuscle@meta.data[which(skeletalmuscle@meta.data$cell_name==names(table(skeletalmuscle$cell_name))[i]&skeletalmuscle@meta.data$condition=="FD"),])
    FS_cell_number<-nrow(skeletalmuscle@meta.data[which(skeletalmuscle@meta.data$cell_name==names(table(skeletalmuscle$cell_name))[i]&skeletalmuscle@meta.data$condition=="FS"),])
    MC_cell_number<-nrow(skeletalmuscle@meta.data[which(skeletalmuscle@meta.data$cell_name==names(table(skeletalmuscle$cell_name))[i]&skeletalmuscle@meta.data$condition=="MC"),])
    MS_cell_number<-nrow(skeletalmuscle@meta.data[which(skeletalmuscle@meta.data$cell_name==names(table(skeletalmuscle$cell_name))[i]&skeletalmuscle@meta.data$condition=="MS"),])
    cell.use<-c(cell.use,FD_cell_number>2&FS_cell_number>2&MC_cell_number>2&MS_cell_number>2)
  }
  cell.use.data<-data.frame(table(skeletalmuscle$cell_name),cell.use)
  celltypes_skeletalmuscle_filter<-as.character(unique(cell.use.data[which(cell.use.data$cell.use==TRUE),]$Var1))
  
  VS <- data.frame(compare=c('MSVSFS','FDVSFS','MCVSMS'),
                   ident1=c('_MS','_FD','_MC'),
                   ident2=c('_FS','_FS','_MS'))
  
  ###################Define DEGs across multiple cell types
  deg_skeletalmuscle <- data.frame()
  for (j in (1:dim(VS)[1])){
    eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
    for (i in (1:length(celltypes_skeletalmuscle_filter))){
      assign(paste0("deg_markers_",celltypes_skeletalmuscle_filter[i],"_",format(VS[j,1]),"_skeletalmuscle",sep=""),FindMarkers(skeletalmuscle, ident.1 = paste0(celltypes_skeletalmuscle_filter[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_skeletalmuscle_filter[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
      eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_skeletalmuscle_filter[i]),"_",format(VS[j,1]),"_skeletalmuscle)[1]>0){deg_markers_",format(celltypes_skeletalmuscle_filter[i]),"_",format(VS[j,1]),"_skeletalmuscle <- subset(deg_markers_",format(celltypes_skeletalmuscle_filter[i]),"_",format(VS[j,1]),"_skeletalmuscle,p_val_adj<0.05)};
						   if(dim(deg_markers_",format(celltypes_skeletalmuscle_filter[i]),"_",format(VS[j,1]),"_skeletalmuscle)[1]>0){deg_markers_",format(celltypes_skeletalmuscle_filter[i]),"_",format(VS[j,1]),"_skeletalmuscle$Gene_symbol <- row.names(deg_markers_",format(celltypes_skeletalmuscle_filter[i]),"_",format(VS[j,1]),"_skeletalmuscle);
						   deg_markers_",format(celltypes_skeletalmuscle_filter[i]),"_",format(VS[j,1]),"_skeletalmuscle$Cell_type <- '",format(celltypes_skeletalmuscle_filter[i]),"';
						   deg_markers_",format(celltypes_skeletalmuscle_filter[i]),"_",format(VS[j,1]),"_skeletalmuscle$Change <- ifelse(deg_markers_",format(celltypes_skeletalmuscle_filter[i]),"_",format(VS[j,1]),"_skeletalmuscle$avg_log2FC > 0 ,'up','down');
						   deg_markers_",format(celltypes_skeletalmuscle_filter[i]),"_",format(VS[j,1]),"_skeletalmuscle$Organ <- 'skeletalmuscle';
						   deg_markers_",format(celltypes_skeletalmuscle_filter[i]),"_",format(VS[j,1]),"_skeletalmuscle$Compare <- '",format(VS[j,1]),"';",
                               format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_skeletalmuscle_filter[i]),"_",format(VS[j,1]),"_skeletalmuscle)};
						   rm(deg_markers_",format(celltypes_skeletalmuscle_filter[i]),"_",format(VS[j,1]),"_skeletalmuscle)")))
    }
    gc()
    eval(parse(text = paste0("deg_skeletalmuscle <- bind_rows(deg_skeletalmuscle,",format(VS[j,1]),");
                         rm(",format(VS[j,1]),")")))
  }
  dim(deg_skeletalmuscle)
  write.csv(deg_skeletalmuscle,'E:/lifei/androgen_project/20211216_Reanalyses/DEG/DEG_list/DEGs.skeletalmuscle.avg_log2FC05.p_val_adj005.csv')

###################15. spleen: DEGs definition across multiple cell types and multiple conditions
###################
###################
###################
###################
###################

  DefaultAssay(spleen.combined.filter.integrated) <- "RNA"
  spleen.combined.filter.integrated$celltype.stim <- paste(Idents(spleen.combined.filter.integrated), spleen.combined.filter.integrated$condition, sep = "_")
  spleen.combined.filter.integrated$celltype <- Idents(spleen.combined.filter.integrated)
  Idents(spleen.combined.filter.integrated) <- "celltype.stim"
  spleen<-spleen.combined.filter.integrated
  table(spleen$celltype.stim)
  length(table(spleen$celltype.stim))#89
  length(table(spleen$cell_name))#23
  #############Cell_type filtering for downstream DEG analysis
  cell.use<-c()
  for (i in 1:length(table(spleen$cell_name))){
    FD_cell_number<-nrow(spleen@meta.data[which(spleen@meta.data$cell_name==names(table(spleen$cell_name))[i]&spleen@meta.data$condition=="FD"),])
    FS_cell_number<-nrow(spleen@meta.data[which(spleen@meta.data$cell_name==names(table(spleen$cell_name))[i]&spleen@meta.data$condition=="FS"),])
    MC_cell_number<-nrow(spleen@meta.data[which(spleen@meta.data$cell_name==names(table(spleen$cell_name))[i]&spleen@meta.data$condition=="MC"),])
    MS_cell_number<-nrow(spleen@meta.data[which(spleen@meta.data$cell_name==names(table(spleen$cell_name))[i]&spleen@meta.data$condition=="MS"),])
    cell.use<-c(cell.use,FD_cell_number>2&FS_cell_number>2&MC_cell_number>2&MS_cell_number>2)
  }
  cell.use.data<-data.frame(table(spleen$cell_name),cell.use)
  celltypes_spleen_filter<-as.character(unique(cell.use.data[which(cell.use.data$cell.use==TRUE),]$Var1))
  
  VS <- data.frame(compare=c('MSVSFS','FDVSFS','MCVSMS'),
                   ident1=c('_MS','_FD','_MC'),
                   ident2=c('_FS','_FS','_MS'))
  
  ###################Define DEGs across multiple cell types
  deg_spleen <- data.frame()
  for (j in (1:dim(VS)[1])){
    eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
    for (i in (1:length(celltypes_spleen_filter))){
      assign(paste0("deg_markers_",celltypes_spleen_filter[i],"_",format(VS[j,1]),"_spleen",sep=""),FindMarkers(spleen, ident.1 = paste0(celltypes_spleen_filter[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_spleen_filter[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
      eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_spleen_filter[i]),"_",format(VS[j,1]),"_spleen)[1]>0){deg_markers_",format(celltypes_spleen_filter[i]),"_",format(VS[j,1]),"_spleen <- subset(deg_markers_",format(celltypes_spleen_filter[i]),"_",format(VS[j,1]),"_spleen,p_val_adj<0.05)};
						   if(dim(deg_markers_",format(celltypes_spleen_filter[i]),"_",format(VS[j,1]),"_spleen)[1]>0){deg_markers_",format(celltypes_spleen_filter[i]),"_",format(VS[j,1]),"_spleen$Gene_symbol <- row.names(deg_markers_",format(celltypes_spleen_filter[i]),"_",format(VS[j,1]),"_spleen);
						   deg_markers_",format(celltypes_spleen_filter[i]),"_",format(VS[j,1]),"_spleen$Cell_type <- '",format(celltypes_spleen_filter[i]),"';
						   deg_markers_",format(celltypes_spleen_filter[i]),"_",format(VS[j,1]),"_spleen$Change <- ifelse(deg_markers_",format(celltypes_spleen_filter[i]),"_",format(VS[j,1]),"_spleen$avg_log2FC > 0 ,'up','down');
						   deg_markers_",format(celltypes_spleen_filter[i]),"_",format(VS[j,1]),"_spleen$Organ <- 'spleen';
						   deg_markers_",format(celltypes_spleen_filter[i]),"_",format(VS[j,1]),"_spleen$Compare <- '",format(VS[j,1]),"';",
                               format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_spleen_filter[i]),"_",format(VS[j,1]),"_spleen)};
						   rm(deg_markers_",format(celltypes_spleen_filter[i]),"_",format(VS[j,1]),"_spleen)")))
    }
    gc()
    eval(parse(text = paste0("deg_spleen <- bind_rows(deg_spleen,",format(VS[j,1]),");
                         rm(",format(VS[j,1]),")")))
  }
  dim(deg_spleen)
  write.csv(deg_spleen,'E:/lifei/androgen_project/20211216_Reanalyses/DEG/DEG_list/DEGs.spleen.avg_log2FC05.p_val_adj005.csv')

###################16. stomach: DEGs definition across multiple cell types and multiple conditions
###################
###################
###################
###################
###################

  DefaultAssay(stomach.combined.filter.integrated) <- "RNA"
  stomach.combined.filter.integrated$celltype.stim <- paste(Idents(stomach.combined.filter.integrated), stomach.combined.filter.integrated$condition, sep = "_")
  stomach.combined.filter.integrated$celltype <- Idents(stomach.combined.filter.integrated)
  Idents(stomach.combined.filter.integrated) <- "celltype.stim"
  stomach<-stomach.combined.filter.integrated
  table(stomach$celltype.stim)
  length(table(stomach$celltype.stim))#89
  length(table(stomach$cell_name))#23
  #############Cell_type filtering for downstream DEG analysis
  cell.use<-c()
  for (i in 1:length(table(stomach$cell_name))){
    FD_cell_number<-nrow(stomach@meta.data[which(stomach@meta.data$cell_name==names(table(stomach$cell_name))[i]&stomach@meta.data$condition=="FD"),])
    FS_cell_number<-nrow(stomach@meta.data[which(stomach@meta.data$cell_name==names(table(stomach$cell_name))[i]&stomach@meta.data$condition=="FS"),])
    MC_cell_number<-nrow(stomach@meta.data[which(stomach@meta.data$cell_name==names(table(stomach$cell_name))[i]&stomach@meta.data$condition=="MC"),])
    MS_cell_number<-nrow(stomach@meta.data[which(stomach@meta.data$cell_name==names(table(stomach$cell_name))[i]&stomach@meta.data$condition=="MS"),])
    cell.use<-c(cell.use,FD_cell_number>2&FS_cell_number>2&MC_cell_number>2&MS_cell_number>2)
  }
  cell.use.data<-data.frame(table(stomach$cell_name),cell.use)
  celltypes_stomach_filter<-as.character(unique(cell.use.data[which(cell.use.data$cell.use==TRUE),]$Var1))
  
  VS <- data.frame(compare=c('MSVSFS','FDVSFS','MCVSMS'),
                   ident1=c('_MS','_FD','_MC'),
                   ident2=c('_FS','_FS','_MS'))
  
  ###################Define DEGs across multiple cell types
  deg_stomach <- data.frame()
  for (j in (1:dim(VS)[1])){
    eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
    for (i in (1:length(celltypes_stomach_filter))){
      assign(paste0("deg_markers_",celltypes_stomach_filter[i],"_",format(VS[j,1]),"_stomach",sep=""),FindMarkers(stomach, ident.1 = paste0(celltypes_stomach_filter[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_stomach_filter[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
      eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_stomach_filter[i]),"_",format(VS[j,1]),"_stomach)[1]>0){deg_markers_",format(celltypes_stomach_filter[i]),"_",format(VS[j,1]),"_stomach <- subset(deg_markers_",format(celltypes_stomach_filter[i]),"_",format(VS[j,1]),"_stomach,p_val_adj<0.05)};
						   if(dim(deg_markers_",format(celltypes_stomach_filter[i]),"_",format(VS[j,1]),"_stomach)[1]>0){deg_markers_",format(celltypes_stomach_filter[i]),"_",format(VS[j,1]),"_stomach$Gene_symbol <- row.names(deg_markers_",format(celltypes_stomach_filter[i]),"_",format(VS[j,1]),"_stomach);
						   deg_markers_",format(celltypes_stomach_filter[i]),"_",format(VS[j,1]),"_stomach$Cell_type <- '",format(celltypes_stomach_filter[i]),"';
						   deg_markers_",format(celltypes_stomach_filter[i]),"_",format(VS[j,1]),"_stomach$Change <- ifelse(deg_markers_",format(celltypes_stomach_filter[i]),"_",format(VS[j,1]),"_stomach$avg_log2FC > 0 ,'up','down');
						   deg_markers_",format(celltypes_stomach_filter[i]),"_",format(VS[j,1]),"_stomach$Organ <- 'stomach';
						   deg_markers_",format(celltypes_stomach_filter[i]),"_",format(VS[j,1]),"_stomach$Compare <- '",format(VS[j,1]),"';",
                               format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_stomach_filter[i]),"_",format(VS[j,1]),"_stomach)};
						   rm(deg_markers_",format(celltypes_stomach_filter[i]),"_",format(VS[j,1]),"_stomach)")))
    }
    gc()
    eval(parse(text = paste0("deg_stomach <- bind_rows(deg_stomach,",format(VS[j,1]),");
                         rm(",format(VS[j,1]),")")))
  }
  dim(deg_stomach)
  write.csv(deg_stomach,'E:/lifei/androgen_project/20211216_Reanalyses/DEG/DEG_list/DEGs.stomach.avg_log2FC05.p_val_adj005.csv')

###################17. thymus: DEGs definition across multiple cell types and multiple conditions
###################
###################
###################
###################
###################

  DefaultAssay(thymus.combined.filter.integrated) <- "RNA"
  thymus.combined.filter.integrated$celltype.stim <- paste(Idents(thymus.combined.filter.integrated), thymus.combined.filter.integrated$condition, sep = "_")
  thymus.combined.filter.integrated$celltype <- Idents(thymus.combined.filter.integrated)
  Idents(thymus.combined.filter.integrated) <- "celltype.stim"
  thymus<-thymus.combined.filter.integrated
  table(thymus$celltype.stim)
  length(table(thymus$celltype.stim))#89
  length(table(thymus$cell_name))#23
  #############Cell_type filtering for downstream DEG analysis
  cell.use<-c()
  for (i in 1:length(table(thymus$cell_name))){
    FD_cell_number<-nrow(thymus@meta.data[which(thymus@meta.data$cell_name==names(table(thymus$cell_name))[i]&thymus@meta.data$condition=="FD"),])
    FS_cell_number<-nrow(thymus@meta.data[which(thymus@meta.data$cell_name==names(table(thymus$cell_name))[i]&thymus@meta.data$condition=="FS"),])
    MC_cell_number<-nrow(thymus@meta.data[which(thymus@meta.data$cell_name==names(table(thymus$cell_name))[i]&thymus@meta.data$condition=="MC"),])
    MS_cell_number<-nrow(thymus@meta.data[which(thymus@meta.data$cell_name==names(table(thymus$cell_name))[i]&thymus@meta.data$condition=="MS"),])
    cell.use<-c(cell.use,FD_cell_number>2&FS_cell_number>2&MC_cell_number>2&MS_cell_number>2)
  }
  cell.use.data<-data.frame(table(thymus$cell_name),cell.use)
  celltypes_thymus_filter<-as.character(unique(cell.use.data[which(cell.use.data$cell.use==TRUE),]$Var1))
  
  VS <- data.frame(compare=c('MSVSFS','FDVSFS','MCVSMS'),
                   ident1=c('_MS','_FD','_MC'),
                   ident2=c('_FS','_FS','_MS'))
  
  ###################Define DEGs across multiple cell types
  deg_thymus <- data.frame()
  for (j in (1:dim(VS)[1])){
    eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
    for (i in (1:length(celltypes_thymus_filter))){
      assign(paste0("deg_markers_",celltypes_thymus_filter[i],"_",format(VS[j,1]),"_thymus",sep=""),FindMarkers(thymus, ident.1 = paste0(celltypes_thymus_filter[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_thymus_filter[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
      eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_thymus_filter[i]),"_",format(VS[j,1]),"_thymus)[1]>0){deg_markers_",format(celltypes_thymus_filter[i]),"_",format(VS[j,1]),"_thymus <- subset(deg_markers_",format(celltypes_thymus_filter[i]),"_",format(VS[j,1]),"_thymus,p_val_adj<0.05)};
						   if(dim(deg_markers_",format(celltypes_thymus_filter[i]),"_",format(VS[j,1]),"_thymus)[1]>0){deg_markers_",format(celltypes_thymus_filter[i]),"_",format(VS[j,1]),"_thymus$Gene_symbol <- row.names(deg_markers_",format(celltypes_thymus_filter[i]),"_",format(VS[j,1]),"_thymus);
						   deg_markers_",format(celltypes_thymus_filter[i]),"_",format(VS[j,1]),"_thymus$Cell_type <- '",format(celltypes_thymus_filter[i]),"';
						   deg_markers_",format(celltypes_thymus_filter[i]),"_",format(VS[j,1]),"_thymus$Change <- ifelse(deg_markers_",format(celltypes_thymus_filter[i]),"_",format(VS[j,1]),"_thymus$avg_log2FC > 0 ,'up','down');
						   deg_markers_",format(celltypes_thymus_filter[i]),"_",format(VS[j,1]),"_thymus$Organ <- 'thymus';
						   deg_markers_",format(celltypes_thymus_filter[i]),"_",format(VS[j,1]),"_thymus$Compare <- '",format(VS[j,1]),"';",
                               format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_thymus_filter[i]),"_",format(VS[j,1]),"_thymus)};
						   rm(deg_markers_",format(celltypes_thymus_filter[i]),"_",format(VS[j,1]),"_thymus)")))
    }
    gc()
    eval(parse(text = paste0("deg_thymus <- bind_rows(deg_thymus,",format(VS[j,1]),");
                         rm(",format(VS[j,1]),")")))
  }
  dim(deg_thymus)
  write.csv(deg_thymus,'E:/lifei/androgen_project/20211216_Reanalyses/DEG/DEG_list/DEGs.thymus.avg_log2FC05.p_val_adj005.csv')

  
###################Combine all the DEGs for 17 organs
###################
###################
###################
###################
###################
organ <- c("adipose","adrenal","bone_marrow","brain","colon","heart","intestine","kidney","lacrimal",
           "liver","lung","pancreas","salivary","skeletalmuscle","spleen","stomach","thymus")
deg <- data.frame()
for (i in 1:length(organ)) {
  eval(parse(text = paste0("deg_",format(organ[i])," <- read.csv('","E:/lifei/androgen_project/20211216_Reanalyses/DEG/DEG_list/DEGs.",format(organ[i]),".avg_log2FC05.p_val_adj005.csv',header=T,row.names=1);
                           deg <- bind_rows(deg,deg_",format(organ[i]),");
  rm(deg_",format(organ[i]),")")))
}
write.csv(deg,'DEG.csv')

