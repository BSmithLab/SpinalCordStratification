
#Short script to plot RPKM normalized expression of ALS-Ox marker genes, grouped by % ALS-Ox

#Written by: Jarrett Eshima
#Date: 1-7-24

library(ggplot2)

## Read in the data

#Files found at: https://figshare.com/authors/Jarrett_Eshima/13813720

setwd("G:/SpinalCord/Publication/SummaryPlot")

SpinalPheno = read.csv("SpinalCord_Phenotype.csv")
CortexPheno = read.csv("Cortex_Phenotype.csv")
PatientPheno = read.csv("Supplemental_Dataset2.csv")

SpinalExpression = read.csv("Spinal_RPKM_OX_7GeneClassifier.csv")
rownames(SpinalExpression) = SpinalExpression$X; SpinalExpression = SpinalExpression[,-1]
CortexExpression = read.csv("Cortex_RPKM_OX_7GeneClassifier.csv")
rownames(CortexExpression) = CortexExpression$X; CortexExpression = CortexExpression[,-1]

for(i in 1:nrow(CortexPheno)){
  if(CortexPheno$Subtype[i] == "TE"){
    CortexPheno$Subtype[i] = "TD"
  }else if(CortexPheno$Subtype[i] == "Control"){
    CortexPheno$Subtype[i] = "CTR"
  }
}

Patients = names(table(SpinalPheno$ExternalSubjectId))[names(table(SpinalPheno$ExternalSubjectId)) %in% names(table(CortexPheno$Patient))] #unique patients with obs from both cortex and spinal cord
Controls = names(table(SpinalPheno$ExternalSubjectId[which(SpinalPheno$Subtype == "CTR")]))[names(table(SpinalPheno$ExternalSubjectId[which(SpinalPheno$Subtype == "CTR")])) %in% names(table(CortexPheno$Patient[which(CortexPheno$Subtype == "CTR")]))]
ALS = Patients[!Patients %in% Controls]

PatientPheno$OxPercent = NA
PatientPheno$OxCort = NA
PatientPheno$OxSpinal = NA

for(i in 1:nrow(PatientPheno)){
  
  subtypes = PatientPheno[i,2:8]
  subtypes = subtypes[!is.na(subtypes)]
  
  nOx = as.numeric(table(subtypes)[which(names(table(subtypes)) == "OX")])
  if(length(nOx)>0){
    PatientPheno$OxPercent[i] = nOx/sum(as.numeric(table(subtypes)))
  }else{
    PatientPheno$OxPercent[i] = 0/sum(as.numeric(table(subtypes)))
  }
  
  cortsubtypes = PatientPheno[i,2:5]
  cortsubtypes = cortsubtypes[!is.na(cortsubtypes)]
  
  if(length(cortsubtypes)>0){
    nOxc = as.numeric(table(cortsubtypes)[which(names(table(cortsubtypes)) == "OX")])
    if(length(nOxc)>0){
      PatientPheno$OxCort[i] = nOxc/sum(as.numeric(table(cortsubtypes)))
    }else{
      PatientPheno$OxCort[i] = 0/sum(as.numeric(table(cortsubtypes)))
    }
  }
  
  spinalsubtypes = PatientPheno[i,6:8]
  spinalsubtypes = spinalsubtypes[!is.na(spinalsubtypes)]
  
  if(length(spinalsubtypes)>0){
    nOxs = as.numeric(table(spinalsubtypes)[which(names(table(spinalsubtypes)) == "OX")])
    if(length(nOxs)>0){
      PatientPheno$OxSpinal[i] = nOxs/sum(as.numeric(table(spinalsubtypes)))
    }else{
      PatientPheno$OxSpinal[i] = 0/sum(as.numeric(table(spinalsubtypes)))
    }
  }
}

CortexPheno$OxPercent = NA
CortexPheno$OxCort = NA
CortexPheno$OxSpinal = NA

for(i in 1:nrow(CortexPheno)){
  ind = which(PatientPheno$Patient == CortexPheno$Patient[i])
  if(length(ind)>0){
    CortexPheno$OxPercent[i] = PatientPheno$OxPercent[ind]
    CortexPheno$OxCort[i] = PatientPheno$OxCort[ind]
    CortexPheno$OxSpinal[i] = PatientPheno$OxSpinal[ind]
  }
}

SpinalPheno$OxPercent = NA
SpinalPheno$OxCort = NA
SpinalPheno$OxSpinal = NA

for(i in 1:nrow(SpinalPheno)){
  ind = which(PatientPheno$Patient == SpinalPheno$ExternalSubjectId[i])
  if(length(ind)>0){
    SpinalPheno$OxPercent[i] = PatientPheno$OxPercent[ind]
    SpinalPheno$OxCort[i] = PatientPheno$OxCort[ind]
    SpinalPheno$OxSpinal[i] = PatientPheno$OxSpinal[ind]
  }
}


tmpind = NA
for(i in 1:nrow(CortexPheno)){
  if(CortexPheno$Subtype[i] != "CTR" && is.na(CortexPheno$OxPercent[i])){
    tmpind = c(tmpind,i)
  }
}
tmpind = tmpind[!is.na(tmpind)]
FiltCortexPheno = CortexPheno[-tmpind,] #no corresponding spinal cord samples


tmpind = NA
for(i in 1:nrow(SpinalPheno)){
  if(SpinalPheno$Subtype[i] != "CTR" && is.na(SpinalPheno$OxPercent[i])){
    tmpind = c(tmpind,i)
  }
}
tmpind = tmpind[!is.na(tmpind)]
FiltSpinalPheno = SpinalPheno #All spinal samples have corresponding cortex

## Plot 1 - Individual samples plotted and grouped by % ALS-Ox using factors

FiltCortexPheno$OxFactor = NA

for(i in 1:nrow(FiltCortexPheno)){
  if(FiltCortexPheno$Subtype[i] == "CTR"){
    FiltCortexPheno$OxFactor[i] = 4
  }else{
    if(FiltCortexPheno$OxPercent[i] == 1){
      FiltCortexPheno$OxFactor[i] = 1
    }else if(FiltCortexPheno$OxPercent[i] >= 0.5 && FiltCortexPheno$OxPercent[i] < 1){
      FiltCortexPheno$OxFactor[i] = 2
    }else if(FiltCortexPheno$OxPercent[i] < 0.5){
      FiltCortexPheno$OxFactor[i] = 3
    }else{
      FiltCortexPheno$OxFactor[i] = NA
    }
  }
}

FiltSpinalPheno$OxFactor = NA

for(i in 1:nrow(FiltSpinalPheno)){
  if(FiltSpinalPheno$Subtype[i] == "CTR"){
    FiltSpinalPheno$OxFactor[i] = 4
  }else{
    if(FiltSpinalPheno$OxPercent[i] == 1){
      FiltSpinalPheno$OxFactor[i] = 1
    }else if(FiltSpinalPheno$OxPercent[i] >= 0.5 && FiltSpinalPheno$OxPercent[i] < 1){
      FiltSpinalPheno$OxFactor[i] = 2
    }else if(FiltSpinalPheno$OxPercent[i] < 0.5){
      FiltSpinalPheno$OxFactor[i] = 3
    }else{
      FiltSpinalPheno$OxFactor[i] = NA
    }
  }
}

CortexExpression = t(CortexExpression)
FiltCortexExpression = CortexExpression[rownames(CortexExpression) %in% FiltCortexPheno$Subject,]
SpinalExpression = t(SpinalExpression)

FullExpression = rbind(FiltCortexExpression,SpinalExpression)

#Reorder to match sample order in expression
FactorTable = data.frame(matrix(NA,nrow=sum(nrow(FiltCortexPheno),nrow(FiltSpinalPheno)),ncol = 2))
colnames(FactorTable) = c("Sample","Factor")
FactorTable$Sample = rownames(FullExpression)

for(i in 1:nrow(FactorTable)){
  ind1 = which(FiltCortexPheno$Subject == FactorTable$Sample[i])
  ind2 = which(FiltSpinalPheno$ExternalSampleId == FactorTable$Sample[i])
  if(length(ind1)>0){
    FactorTable$Factor[i] = FiltCortexPheno$OxFactor[ind1]
  }else if(length(ind2)>0){
    FactorTable$Factor[i] = FiltSpinalPheno$OxFactor[ind2]
  }else{
    cat(paste("Problem with sample: ",FactorTable$Sample[i],"\n",sep=""))
  }
}

ggdat = cbind(FullExpression,FactorTable$Factor)
colnames(ggdat)[ncol(ggdat)] = "Factor"
ggdat2 = data.frame(ggdat)
for(i in 1:nrow(ggdat2)){
  if(ggdat2$Factor[i] == 1){
    ggdat2$Factor[i] = "Concordant ALS-Ox"
  }else if(ggdat2$Factor[i] == 2){
    ggdat2$Factor[i] = "At least 50% ALS-Ox"
  }else if(ggdat2$Factor[i] == 3){
    ggdat2$Factor[i] = "Generally not ALS-Ox"
  }else if(ggdat2$Factor[i] == 4){
    ggdat2$Factor[i] = "Control"
  }
}

#Now prep for ggplot 
ngene = 7
ggdat3 = data.frame(matrix(NA,nrow=nrow(ggdat2)*(ncol(ggdat2)-1),ncol=4))
colnames(ggdat3) = c("Sample","Group","Gene","Expression")
ggdat3$Sample = rep(rownames(ggdat2),ngene)
ggdat3$Gene = c(rep(colnames(ggdat2)[1],nrow(ggdat2)),rep(colnames(ggdat2)[2],nrow(ggdat2)),rep(colnames(ggdat2)[3],nrow(ggdat2)),rep(colnames(ggdat2)[4],nrow(ggdat2)),rep(colnames(ggdat2)[5],nrow(ggdat2)),rep(colnames(ggdat2)[6],nrow(ggdat2)),rep(colnames(ggdat2)[7],nrow(ggdat2)))
ggdat3$Expression = c(ggdat2$GABRA1,ggdat2$SLC17A6,ggdat2$HTR2A,ggdat2$B4GALT6,ggdat2$GAD2,ggdat2$GLRA3,ggdat2$PCSK1)
ggdat3$Group = rep(ggdat2$Factor,ngene)
ggdat3$Group = factor(ggdat3$Group,levels = c("Concordant ALS-Ox","At least 50% ALS-Ox","Generally not ALS-Ox","Control"))

customcol = c(rep("#bf0808",nrow(ggdat2)),rep("#ebbc57",nrow(ggdat2)),rep("#46d44d",nrow(ggdat2)),rep("#21b570",nrow(ggdat2)),rep("#52bce3",nrow(ggdat2)),rep("#a078eb",nrow(ggdat2)),rep("#e872e8",nrow(ggdat2)))

p = ggplot(ggdat3,aes(x=Group,y=Expression,fill=Gene)) + geom_boxplot(outlier.shape = NA)
p = p+scale_fill_manual(values = c("#bf0808","#ebbc57","#46d44d","#409180","#73adeb","#a078eb","#e872e8"))
p = p+theme(panel.background = element_rect(fill = 'white', color = 'white'),panel.grid.major = element_line(color = 'gray75'),panel.grid.minor = element_line(color = 'gray75'))
#p = p+ylim(0,1)
p = p+ theme(plot.margin = unit(c(1,1,1,1), "cm"))
p = p+theme(axis.title.x=element_text(vjust=-2))
p = p+theme(axis.title.y=element_text(angle=90, vjust=6))
p = p + geom_point(aes(x=Group,y=Expression,fill=Gene,colour=Group),shape=19,size=0.5,position = position_jitterdodge()) + facet_wrap(~Gene,scales = "free")
p


table(FiltCortexPheno$Patient[which(FiltCortexPheno$Subtype == "CTR")]) %in% table(FiltSpinalPheno$ExternalSubjectId[which(FiltSpinalPheno$Subtype == "CTR")])
