#High Concordance Heatmap
setwd("G:/SpinalCord/Publication/Manuscript/Tables/Concordance")
TD_Pheno_Final = read.csv("ALS-TD_HighlyConcordantCNS_4Covar_9-21-23_MatchingCNS.csv") #Supplemental Dataset 4, rows 24-46
TD_Pheno_Final = TD_Pheno_Final[-which(TD_Pheno_Final$SubjectID == ""),]
Ox_Pheno_Final = read.csv("ALS-Ox_HighlyConcordantCNS_4Covar_9-21-23_MatchingCNS.csv") #Supplemental Dataset 4, rows 10-23
Ox_Pheno_Final = Ox_Pheno_Final[-which(Ox_Pheno_Final$SubjectID == ""),]
Glia_Pheno_Final = read.csv("ALS-Glia_HighlyConcordantCNS_4Covar_9-21-23_MatchingCNS.csv") #Supplemental Dataset 4, rows 2-9
Glia_Pheno_Final = Glia_Pheno_Final[-which(Glia_Pheno_Final$SubjectID == ""),]

### ALS-Glia

TissueLoc = c("FrontalCortexSubtype","MedialMotorCortexSubtype","LateralMotorCortexSubtype","UnspecifiedMotorCortexSubtype","CervicalSubtype","ThoracicSubtype","LumbarSubtype")
tmpheat = data.frame(matrix(NA,nrow=nrow(Glia_Pheno_Final)*length(TissueLoc),ncol = 3))
colnames(tmpheat) = c("Module","Response","Value")

for(i in 1:nrow(Glia_Pheno_Final)){
  
  tmp = rep(paste(Glia_Pheno_Final$SubjectID[i],sep=""),length(TissueLoc))
  if(i == 1){
    ModLabs = tmp
  }else{
    ModLabs = c(ModLabs,tmp)
  }
  
}

tmpheat$Module = ModLabs
tmpheat$Response = rep(TissueLoc,nrow(Glia_Pheno_Final))

for(i in 1:nrow(tmpheat)){
  
  rind = which(Glia_Pheno_Final$SubjectID == tmpheat$Module[i])
  cind = which(colnames(Glia_Pheno_Final) == tmpheat$Response[i])
  
  tmpheat$Value[i] = Glia_Pheno_Final[rind,cind]
  
}

GLIA_Heat = tmpheat

### ALS-Ox

TissueLoc = c("FrontalCortexSubtype","MedialMotorCortexSubtype","LateralMotorCortexSubtype","UnspecifiedMotorCortexSubtype","CervicalSubtype","ThoracicSubtype","LumbarSubtype")
tmpheat = data.frame(matrix(NA,nrow=nrow(Ox_Pheno_Final)*length(TissueLoc),ncol = 3))
colnames(tmpheat) = c("Module","Response","Value")

for(i in 1:nrow(Ox_Pheno_Final)){
  
  tmp = rep(paste(Ox_Pheno_Final$SubjectID[i],sep=""),length(TissueLoc))
  if(i == 1){
    ModLabs = tmp
  }else{
    ModLabs = c(ModLabs,tmp)
  }
  
}

tmpheat$Module = ModLabs
tmpheat$Response = rep(TissueLoc,nrow(Ox_Pheno_Final))

for(i in 1:nrow(tmpheat)){
  
  rind = which(Ox_Pheno_Final$SubjectID == tmpheat$Module[i])
  cind = which(colnames(Ox_Pheno_Final) == tmpheat$Response[i])
  
  tmpheat$Value[i] = Ox_Pheno_Final[rind,cind]
  
}

OX_Heat = tmpheat


### ALS-TD

TissueLoc = c("FrontalCortexSubtype","MedialMotorCortexSubtype","LateralMotorCortexSubtype","UnspecifiedMotorCortexSubtype","CervicalSubtype","ThoracicSubtype","LumbarSubtype")
tmpheat = data.frame(matrix(NA,nrow=nrow(TD_Pheno_Final)*length(TissueLoc),ncol = 3))
colnames(tmpheat) = c("Module","Response","Value")

for(i in 1:nrow(TD_Pheno_Final)){
  
  tmp = rep(paste(TD_Pheno_Final$SubjectID[i],sep=""),length(TissueLoc))
  if(i == 1){
    ModLabs = tmp
  }else{
    ModLabs = c(ModLabs,tmp)
  }

}

tmpheat$Module = ModLabs
tmpheat$Response = rep(TissueLoc,nrow(TD_Pheno_Final))

for(i in 1:nrow(tmpheat)){
  
  rind = which(TD_Pheno_Final$SubjectID == tmpheat$Module[i])
  cind = which(colnames(TD_Pheno_Final) == tmpheat$Response[i])
  
  tmpheat$Value[i] = TD_Pheno_Final[rind,cind]
  
}

TD_Heat = tmpheat

ConHeat = rbind(GLIA_Heat,OX_Heat,TD_Heat)

for(i in 1:length(ConHeat$Value)){
  if(is.na(ConHeat$Value[i])){
    ConHeat$Value[i] = ""
  }
}

# for(i in 1:length(ConHeat$Value)){
#   if(ConHeat$Value[i] == "TD"){
#     ConHeat$Value[i]
#   }
# }

ConHeat$Value = as.factor(ConHeat$Value)
ConHeat$Response = factor(ConHeat$Response,levels = TissueLoc)
colors = c("gray80","goldenrod1","navy","firebrick")

p = ggplot(ConHeat,aes(Response,Module)) + geom_tile(aes(fill=Value),color="black")
p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "white"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray50"))
p

p = ggplot(ConHeat,aes(Response,Module)) + geom_tile(aes(fill=Value),color="black") + scale_fill_manual(values = colors)
p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "white"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray50"))
p


###################################################################

#### Plot Clinical Variables

#Age at Onset
boxplot(Glia_Pheno_Final$AgeofOnset,Ox_Pheno_Final$AgeofOnset,TD_Pheno_Final$AgeofOnset,col=c("goldenrod1","darkblue","firebrick"),xaxt="n",cex.axis=1.5,cex=1.5,cex.lab=1.5,outline=F,ylim=c(25,90),ylab="Age at Onset (yrs)")
axis(1,at=c(1,2,3),labels = c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis=1.5)
points(rep(1,length(Glia_Pheno_Final$AgeofOnset)),Glia_Pheno_Final$AgeofOnset,pch=19,col="gray50")
points(rep(2,length(Ox_Pheno_Final$AgeofOnset)),Ox_Pheno_Final$AgeofOnset,pch=19,col="gray50")
points(rep(3,length(TD_Pheno_Final$AgeofOnset)),TD_Pheno_Final$AgeofOnset,pch=19,col="gray50")

#Age at Death
boxplot(Glia_Pheno_Final$AgeofDeath,Ox_Pheno_Final$AgeofDeath,TD_Pheno_Final$AgeofDeath,col=c("goldenrod1","darkblue","firebrick"),xaxt="n",cex.axis=1.5,cex=1.5,cex.lab=1.5,outline=F,ylim=c(25,90),ylab="Age at Death (yrs)")
axis(1,at=c(1,2,3),labels = c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis=1.5)
points(rep(1,length(Glia_Pheno_Final$AgeofDeath)),Glia_Pheno_Final$AgeofDeath,pch=19,col="gray50")
points(rep(2,length(Ox_Pheno_Final$AgeofDeath)),Ox_Pheno_Final$AgeofDeath,pch=19,col="gray50")
points(rep(3,length(TD_Pheno_Final$AgeofDeath)),TD_Pheno_Final$AgeofDeath,pch=19,col="gray50")

#Duration
boxplot(Glia_Pheno_Final$Duration,Ox_Pheno_Final$Duration,TD_Pheno_Final$Duration,col=c("goldenrod1","darkblue","firebrick"),xaxt="n",cex.axis=1.5,cex=1.5,cex.lab=1.5,outline=F,ylim=c(0,140),ylab="Disease Duration (months)")
axis(1,at=c(1,2,3),labels = c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis=1.5)
points(rep(1,length(Glia_Pheno_Final$Duration)),Glia_Pheno_Final$Duration,pch=19,col="gray50")
points(rep(2,length(Ox_Pheno_Final$Duration)),Ox_Pheno_Final$Duration,pch=19,col="gray50")
points(rep(3,length(TD_Pheno_Final$Duration)),TD_Pheno_Final$Duration,pch=19,col="gray50")

#Site of Symptom Onset
#In excel


