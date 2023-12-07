
library(ggplot2)
#Written: November 6th, 2023

setwd("G:/SpinalCord/Publication/Manuscript/Tables")

Cortex = read.csv("Cortex_HC_Pheno.csv") #Provided in Supplementary Dataset 2 from previous paper (https://www.nature.com/articles/s41467-022-35494-w)
Spinal = read.csv("SampleLevelPhenotype.csv") #Provided in Supplementary Dataset 1 from this paper

# Cortex_HC = Cortex$Subject_ID
# Spinal_HC = Spinal$ExternalSubjectId[which(Spinal$Subtype == "CONTROL")]
# Clean_HC = names(table(c(Cortex_HC,Spinal_HC)))
#n=88 biologically independent non-neurological control patients


## Truncated Stathmin 2 - Nice plot

ggdat = data.frame(matrix(NA,nrow=nrow(Spinal)*3,ncol=3))
colnames(ggdat) = c("Subtype","CountScale","Value")
ggdat$Subtype = rep(Spinal$Subtype,3)
ggdat$CountScale = c(rep("STMN2_TPM",nrow(Spinal)),rep("tSTMN2_counts",nrow(Spinal)),rep("tSTMN2_TPM",nrow(Spinal)))
ggdat$Value = c(Spinal$STMN2_TPM,Spinal$tSTMN2_counts,Spinal$tSTMN2_TPM)

ggdat$Subtype = as.factor(ggdat$Subtype)
ggdat$CountScale = as.factor(ggdat$CountScale)
ggdat$Value = as.numeric(ggdat$Value)

customcol = rep(NA,nrow(ggdat))
customshape = rep(NA,nrow(ggdat))
for(i in 1:nrow(ggdat)){
  if(ggdat$Subtype[i] == "GLIA"){
    customcol[i] = "#faee66"
    customshape[i] = 21
  }else if(ggdat$Subtype[i] == "OX"){
    customcol[i] = "#3ca2de"
    customshape[i] = 22
  }else if(ggdat$Subtype[i] == "TD"){
    customcol[i] = "#e8844a"
    customshape[i] = 24
  }else{
    customcol[i] = "gray80"
  }
}


#Plots

#With all count scales
p = ggplot(ggdat,aes(x=CountScale,y=Value,fill=Subtype)) + geom_boxplot(outlier.shape = NA)
p = p+scale_fill_manual(values = c("gray40","goldenrod1","navy","firebrick"))
p = p+theme(panel.background = element_rect(fill = 'white', color = 'white'),panel.grid.major = element_line(color = 'gray75'),panel.grid.minor = element_line(color = 'gray75'))
p = p+ theme(plot.margin = unit(c(1,1,1,1), "cm"))
p = p+ylim(0,150)
p = p+theme(axis.title.x=element_text(vjust=-2))
p = p+theme(axis.title.y=element_text(angle=90, vjust=6))
p = p + geom_point(aes(x=CountScale,y=Value,fill=Subtype),shape=21,size=1,position = position_jitterdodge(),col=customcol)
p

#Zoom-in on the tSTMN2_TPM count scale
ggdat2 = ggdat[-which(ggdat$CountScale == "STMN2_TPM"),]
ggdat2 = ggdat2[-which(ggdat2$CountScale == "tSTMN2_counts"),]
p = ggplot(ggdat2,aes(x=CountScale,y=Value,fill=Subtype)) + geom_boxplot(outlier.shape = NA)
p = p+scale_fill_manual(values = c("gray40","goldenrod1","navy","firebrick"))
p = p+theme(panel.background = element_rect(fill = 'white', color = 'white'),panel.grid.major = element_line(color = 'gray75'),panel.grid.minor = element_line(color = 'gray75'))
p = p+ylim(0,2)
p = p+ theme(plot.margin = unit(c(1,1,1,1), "cm"))
p = p+theme(axis.title.x=element_text(vjust=-2))
p = p+theme(axis.title.y=element_text(angle=90, vjust=6))
p = p + geom_point(aes(x=CountScale,y=Value,fill=Subtype),shape=21,size=1,position = position_jitterdodge(),col=customcol[nrow(Spinal)*2+1:nrow(Spinal)])
p

########## STATS

pval = data.frame(matrix(NA,nrow=6,ncol=3))
colnames(pval) = c("STMN2_TPM","tSTMN2_counts","tSTMN2_TPM")
rownames(pval) = c("GvO","GvT","OvT","GvC","OvC","TvC")

for(i in 1:ncol(pval)){
  
  FiltCounts = ggdat[which(ggdat$CountScale == colnames(pval)[i]),]
  
  Glia = FiltCounts$Value[which(FiltCounts$Subtype == "GLIA")]
  Ox = FiltCounts$Value[which(FiltCounts$Subtype == "OX")]
  TD = FiltCounts$Value[which(FiltCounts$Subtype == "TD")]
  CTR = FiltCounts$Value[which(FiltCounts$Subtype == "CONTROL")]
  
  GvO = wilcox.test(Glia,Ox,alternative = "two.sided",paired = F)
  GvT = wilcox.test(Glia,TD,alternative = "two.sided",paired = F)
  OvT = wilcox.test(Ox,TD,alternative = "two.sided",paired = F)
  GvC = wilcox.test(Glia,CTR,alternative = "two.sided",paired = F)
  OvC = wilcox.test(Ox,CTR,alternative = "two.sided",paired = F)
  TvC = wilcox.test(TD,CTR,alternative = "two.sided",paired = F)
  
  pval[1,i] = GvO$p.value
  pval[2,i] = GvT$p.value
  pval[3,i] = OvT$p.value
  pval[4,i] = GvC$p.value
  pval[5,i] = OvC$p.value
  pval[6,i] = TvC$p.value
  
}

pval2 = as.numeric(unlist(pval))

adjp = p.adjust(pval2,method="bonferroni")

adjpmat = data.frame(matrix(adjp,nrow=6,ncol=3))
colnames(adjpmat) = colnames(pval)
rownames(adjpmat) = rownames(pval)