#Check out Jack Humphrey cell deconvolution data from the same cohort

#Written by: Jarrett Eshima
#Date: 11-3-23

#Reference: https://www.nature.com/articles/s41593-022-01205-3
#MuSiC Estimates for ~400 samples are provided in Supplemental Table 4

setwd("G:/SpinalCord/Publication/CellDecon")
CellProps = read.csv("MuSiC_Estimates.csv") #Available from: https://www.nature.com/articles/s41593-022-01205-3#Sec27 (Supplemental Table 3)
colnames(CellProps) = CellProps[1,]
CellProps = CellProps[-1,]
CellProps = CellProps[,1:9]
CellProps$sample = gsub("-",".",CellProps$sample)

SubtypeLabels = read.csv("Supplemental_Dataset2_Page2.csv") #Use the second page of Supplemental Dataset 2

CellProps$Subtype = NA

for(i in 1:nrow(CellProps)){
  ind = which(SubtypeLabels$ExternalSampleId == CellProps$sample[i])
  if(length(ind)>0){
    CellProps$Subtype[i] = SubtypeLabels$Subtype[ind]
  }
  
}

CleanCellProps = CellProps[-which(is.na(CellProps$Subtype)),]
rownames(CleanCellProps) = CleanCellProps$sample; CleanCellProps = CleanCellProps[,-1]

Subtype = CleanCellProps$Subtype
FinalCellProps = CleanCellProps[,-9]

ggdat = data.frame(matrix(NA,nrow=(nrow(FinalCellProps)*ncol(FinalCellProps)),ncol = 3))
colnames(ggdat) = c("Subtype","Cell","Percentage")

ggdat$Subtype = rep(Subtype,ncol(FinalCellProps))
ggdat$Cell = c(rep("Astrocytes",nrow(FinalCellProps)),rep("Endothelial Cells",nrow(FinalCellProps)),rep("Excitatory Neurons",nrow(FinalCellProps)),rep("Inhibitory Neurons",nrow(FinalCellProps)),rep("Microglia",nrow(FinalCellProps)),rep("Oligodendrocytes",nrow(FinalCellProps)),rep("Oligodendrocytes Precursor Cells",nrow(FinalCellProps)),rep("Pericytes",nrow(FinalCellProps)))
ggdat$Percentage = c(FinalCellProps[,1],FinalCellProps[,2],FinalCellProps[,3],FinalCellProps[,4],FinalCellProps[,5],FinalCellProps[,6],FinalCellProps[,7],FinalCellProps[,8])

ggdat2 = ggdat[-which(ggdat$Cell == "Oligodendrocytes Precursor Cells"),]

ggdat3 = ggdat2
ggdat3$Percentage = as.numeric(ggdat3$Percentage)
ggdat3$Percentage[which(ggdat3$Cell == "Excitatory Neurons")] = ggdat3$Percentage[which(ggdat3$Cell == "Excitatory Neurons")] + ggdat3$Percentage[which(ggdat3$Cell == "Inhibitory Neurons")]
ggdat3$Cell[which(ggdat3$Cell == "Excitatory Neurons")] = "Neurons"
ggdat3 = ggdat3[-which(ggdat3$Cell == "Inhibitory Neurons"),]

ggdat3$Subtype = as.factor(ggdat3$Subtype)
ggdat3$Cell = as.factor(ggdat3$Cell)
ggdat3$Percentage = as.numeric(ggdat3$Percentage)

customcol = rep(NA,nrow(ggdat3))
customshape = rep(NA,nrow(ggdat3))
for(i in 1:nrow(ggdat3)){
  if(ggdat3$Subtype[i] == "GLIA"){
    customcol[i] = "#d7ed87"
    customshape[i] = 21
  }else if(ggdat3$Subtype[i] == "OX"){
    customcol[i] = "#3ca2de"
    customshape[i] = 22
  }else if(ggdat3$Subtype[i] == "TD"){
    customcol[i] = "#e8844a"
    customshape[i] = 24
  }
}


#Plot
library(ggplot2)
p = ggplot(ggdat3,aes(x=Cell,y=Percentage,fill=Subtype)) + geom_boxplot(outlier.shape = NA)
p = p+scale_fill_manual(values = c("goldenrod1","navy","firebrick"))
p = p+theme(panel.background = element_rect(fill = 'white', color = 'white'),panel.grid.major = element_line(color = 'gray75'),panel.grid.minor = element_line(color = 'gray75'))
p = p+ylim(0,1)
p = p+ theme(plot.margin = unit(c(1,1,1,1), "cm"))
p = p+theme(axis.title.x=element_text(vjust=-2))
p = p+theme(axis.title.y=element_text(angle=90, vjust=6))
p = p + geom_point(aes(x=Cell,y=Percentage,fill=Subtype),shape=21,size=1,position = position_jitterdodge(),col=customcol)
#p = p + geom_point(aes(x=Cell,y=Percentage,fill=Subtype),shape=customshape,size=1,position = position_jitterdodge(),col=customcol)
p

######## p-values

pval = data.frame(matrix(NA,nrow=3,ncol=6))
colnames(pval) = c("Astrocytes","Endothelial Cells","Microglia","Neurons","Oligodendrocytes","Pericytes")
rownames(pval) = c("GvO","GvT","OvT")

for(i in 1:ncol(pval)){
  
  FiltPercentages = ggdat3[which(ggdat3$Cell == colnames(pval)[i]),]
  
  Glia = FiltPercentages$Percentage[which(FiltPercentages$Subtype == "GLIA")]
  Ox = FiltPercentages$Percentage[which(FiltPercentages$Subtype == "OX")]
  TD = FiltPercentages$Percentage[which(FiltPercentages$Subtype == "TD")]
  
  GvO = wilcox.test(Glia,Ox,alternative = "two.sided",paired = F)
  GvT = wilcox.test(Glia,TD,alternative = "two.sided",paired = F)
  OvT = wilcox.test(Ox,TD,alternative = "two.sided",paired = F)
  
  pval[1,i] = GvO$p.value
  pval[2,i] = GvT$p.value
  pval[3,i] = OvT$p.value
  
}

pval2 = as.numeric(unlist(pval))

adjp = p.adjust(pval2,method="bonferroni")

adjpmat = data.frame(matrix(adjp,nrow=3,ncol=6))
colnames(adjpmat) = colnames(pval)
rownames(adjpmat) = rownames(pval)


metadata = read.csv("GEOCollaborator.csv")
metadata$ExternalSampleId = gsub("-",".",metadata$ExternalSampleId)
nmd = metadata[metadata$ExternalSampleId %in% rownames(CleanCellProps),]

table(nmd$Sample.Source)

# table(metadata$Sample.Source)
# 
# check = metadata[which(metadata$Sample.Source == "Cortex_Frontal"),]
# check$SITE = NA
# 
# sradata = read.csv("SraRunTable.txt")
# sradata$sample_id_alt = gsub("-",".",sradata$sample_id_alt)
# 
# for(i in 1:nrow(check)){
#   ind = which(sradata$sample_id_alt == check$ExternalSampleId[i])
#   if(length(ind)>0){
#     check$SITE[i] = sradata$Project[ind]
#   }
# }
# 
# table(check$SITE)

#Unspecified motor cortex:
#NYGC: 63     TargetALS: 11

#Lateral Motor Cortex:
#NYGC: 16     TargetALS: 121

#Medial Motor Cortex:
#NYGC: 16     TargetALS: 121

#Frontal Cortex:
#NYGC: 183     TargetALS: 144



##################################### Tissue level - scRNA-seq way better for this
TissueProps = CleanCellProps
TissueProps$Tissue = NA

for(i in 1:nrow(TissueProps)){
  ind = which(SubtypeLabels$ExternalSampleId == rownames(TissueProps)[i])
  if(length(ind)>0){
    TissueProps$Tissue[i] = SubtypeLabels$Sample.Source[ind]
  }
}

Subtype = TissueProps[,9]
Tissue = TissueProps[,10]
TissueProps = TissueProps[,-9:-10]

ggdat = data.frame(matrix(NA,nrow=(nrow(TissueProps)*ncol(TissueProps)),ncol = 3))
colnames(ggdat) = c("Tissue","Cell","Percentage")

ggdat$Tissue = rep(Tissue,ncol(TissueProps))
ggdat$Cell = c(rep("Astrocytes",nrow(TissueProps)),rep("Endothelial Cells",nrow(TissueProps)),rep("Excitatory Neurons",nrow(TissueProps)),rep("Inhibitory Neurons",nrow(TissueProps)),rep("Microglia",nrow(TissueProps)),rep("Oligodendrocytes",nrow(TissueProps)),rep("Oligodendrocytes Precursor Cells",nrow(TissueProps)),rep("Pericytes",nrow(TissueProps)))
ggdat$Percentage = c(TissueProps[,1],TissueProps[,2],TissueProps[,3],TissueProps[,4],TissueProps[,5],TissueProps[,6],TissueProps[,7],TissueProps[,8])

ggdat2 = ggdat[-which(ggdat$Cell == "Oligodendrocytes Precursor Cells"),]

ggdat3 = ggdat2
ggdat3$Percentage = as.numeric(ggdat3$Percentage)
ggdat3$Percentage[which(ggdat3$Cell == "Excitatory Neurons")] = ggdat3$Percentage[which(ggdat3$Cell == "Excitatory Neurons")] + ggdat3$Percentage[which(ggdat3$Cell == "Inhibitory Neurons")]
ggdat3$Cell[which(ggdat3$Cell == "Excitatory Neurons")] = "Neurons"
ggdat3 = ggdat3[-which(ggdat3$Cell == "Inhibitory Neurons"),]

ggdat3$Tissue = as.factor(ggdat3$Tissue)
ggdat3$Cell = as.factor(ggdat3$Cell)
ggdat3$Percentage = as.numeric(ggdat3$Percentage)

customcol = rep(NA,length(Subtype))
customshape = rep(NA,length(Subtype))
for(i in 1:length(Subtype)){
  if(Subtype[i] == "GLIA"){
    customcol[i] = "goldenrod1"
    customshape[i] = 21
  }else if(Subtype[i] == "OX"){
    customcol[i] = "navy"
    customshape[i] = 22
  }else if(Subtype[i] == "TD"){
    customcol[i] = "firebrick"
    customshape[i] = 24
  }
}

customcol = rep(customcol,length(names(table(ggdat3$Cell))))
customshape = rep(customshape,length(names(table(ggdat3$Cell))))

#Plot
library(ggplot2)
p = ggplot(ggdat3,aes(x=Cell,y=Percentage,fill=Tissue)) + geom_boxplot(outlier.shape = NA)
p = p+scale_fill_manual(values = c("#7cd6c7","#3caac7","#738aff"))
p = p+theme(panel.background = element_rect(fill = 'white', color = 'white'),panel.grid.major = element_line(color = 'gray75'),panel.grid.minor = element_line(color = 'gray75'))
p = p+ylim(0,1)
p = p+ theme(plot.margin = unit(c(1,1,1,1), "cm"))
p = p+theme(axis.title.x=element_text(vjust=-2))
p = p+theme(axis.title.y=element_text(angle=90, vjust=6))
p = p + geom_point(aes(x=Cell,y=Percentage,fill=Tissue),shape=21,size=2.5,position = position_jitterdodge(),col=customcol)
p