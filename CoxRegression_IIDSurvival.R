#Cox Proportional Hazard Regression Model for Survival (Frailty model)
#Written by: Jarrett Eshima
#Date: December 27th, 2023

####################################################################################################################################
#Load Packages
library(survival)
library(survminer)
library(beeswarm)
library(ggplot2)

#Load Data
setwd("G:/SpinalCord/Publication/FricksSurvival")
CNS_MEM = read.csv("Supplemental_Dataset3.csv") #Submitted with manuscript

#You can skip to line 325 if you don't care how the PC loadings are calculated.

####################################################################################################################################

#Files found at: https://figshare.com/authors/Jarrett_Eshima/13813720

########### VERIFY THE HEATMAP GIVES THE PCA -- YES

HiSeq_Spinal = read.csv("SpinalHeatmap_sake_v3_HiSeq.csv")
rownames(HiSeq_Spinal) = HiSeq_Spinal$Gene; HiSeq_Spinal = HiSeq_Spinal[,-1]
tmp = prcomp(t(HiSeq_Spinal))
plot(tmp$x[,1],tmp$x[,2])

customcols = data.frame(matrix(NA,nrow=ncol(HiSeq_Spinal),ncol = 3))
colnames(customcols) = c("Sample","Group","Color")
customcols$Sample = gsub("-",".",colnames(HiSeq_Spinal))

for(i in 1:nrow(customcols)){
  for(j in 1:nrow(CNS_MEM)){
    if(customcols$Sample[i] == CNS_MEM$Sample[j]){
      customcols$Group[i] = CNS_MEM$Subtype[j]
    }
  }
}

for(i in 1:nrow(customcols)){
  if(customcols$Group[i] == "TD"){
    customcols$Color[i] = "firebrick"
  }else if(customcols$Group[i] == "OX"){
    customcols$Color[i] = "navy"
  }else{
    customcols$Color[i] = "goldenrod1"
  }
}

#Inversing components gives SAKE output - not sure this matters much
plot(-1*tmp$x[,1],-1*tmp$x[,2],col = customcols$Color,pch=20) #The original manuscript version from SAKE
#Without inverse 
plot(tmp$x[,1],tmp$x[,2],col = customcols$Color,pch=20) 

##########

NovaSeq_Spinal = read.csv("SpinalHeatmap_sake_v3_NovaSeq.csv")
rownames(NovaSeq_Spinal) = NovaSeq_Spinal$Gene; NovaSeq_Spinal = NovaSeq_Spinal[,-1]
NovaSeq_Spinal = NovaSeq_Spinal[,-which(colnames(NovaSeq_Spinal) == "CGND.HRA.00290")] #Duplicate patient sample 
tmp = prcomp(t(NovaSeq_Spinal))
plot(tmp$x[,1],tmp$x[,2])

customcols2 = data.frame(matrix(NA,nrow=ncol(NovaSeq_Spinal),ncol = 3))
colnames(customcols2) = c("Sample","Group","Color")
customcols2$Sample = gsub("-",".",colnames(NovaSeq_Spinal))

for(i in 1:nrow(customcols2)){
  for(j in 1:nrow(CNS_MEM)){
    if(customcols2$Sample[i] == CNS_MEM$Sample[j]){
      customcols2$Group[i] = CNS_MEM$Subtype[j]
    }
  }
}

for(i in 1:nrow(customcols2)){
  if(customcols2$Group[i] == "TD"){
    customcols2$Color[i] = "firebrick"
  }else if(customcols2$Group[i] == "OX"){
    customcols2$Color[i] = "navy"
  }else{
    customcols2$Color[i] = "goldenrod1"
  }
}

#Inversing component 1 gives SAKE output - not sure this matters much
plot(-1*tmp$x[,1],tmp$x[,2],col = customcols2$Color,pch=20) #The original manuscript version from SAKE

#Without inverse 
plot(tmp$x[,1],tmp$x[,2],col = customcols2$Color,pch=20) 


HiSeq_Cortex = read.csv("CortexHeatmap_sake_v3_HiSeq.csv")
rownames(HiSeq_Cortex) = HiSeq_Cortex$Gene; HiSeq_Cortex = HiSeq_Cortex[,-1]
for(i in 1:ncol(HiSeq_Cortex)){
  colnames(HiSeq_Cortex)[i] = substr(colnames(HiSeq_Cortex)[i],6,nchar(colnames(HiSeq_Cortex)[i]))
}
HiSeq_Cortex = HiSeq_Cortex[,-which(colnames(HiSeq_Cortex) == "CGND.HRA.00300")] #this is the only pre-fALS patient and causes problems with hazard ratio estimation of disease group covariate

tmp = prcomp(t(HiSeq_Cortex))
plot(tmp$x[,1],tmp$x[,2])

customcols3 = data.frame(matrix(NA,nrow=ncol(HiSeq_Cortex),ncol = 3))
colnames(customcols3) = c("Sample","Group","Color")
customcols3$Sample = gsub("-",".",colnames(HiSeq_Cortex))

for(i in 1:nrow(customcols3)){
  for(j in 1:nrow(CNS_MEM)){
    if(customcols3$Sample[i] == CNS_MEM$Sample[j]){
      customcols3$Group[i] = CNS_MEM$Subtype[j]
    }
  }
}

for(i in 1:nrow(customcols3)){
  if(customcols3$Group[i] == "TD"){
    customcols3$Color[i] = "firebrick"
  }else if(customcols3$Group[i] == "OX"){
    customcols3$Color[i] = "navy"
  }else if(customcols3$Group[i] == "GLIA"){
    customcols3$Color[i] = "goldenrod1"
  }
}

plot(tmp$x[,1],tmp$x[,2],col = customcols3$Color,pch=20) 


NovaSeq_Cortex = read.csv("CortexHeatmap_sake_v3_NovaSeq.csv")
rownames(NovaSeq_Cortex) = NovaSeq_Cortex$Gene; NovaSeq_Cortex = NovaSeq_Cortex[,-1]
for(i in 1:ncol(NovaSeq_Cortex)){
  colnames(NovaSeq_Cortex)[i] = substr(colnames(NovaSeq_Cortex)[i],6,nchar(colnames(NovaSeq_Cortex)[i]))
}
#HiSeq_Cortex = HiSeq_Cortex[,-which(colnames(HiSeq_Cortex) == "CGND.HRA.00300")] #this is the only pre-fALS patient and causes problems with hazard ratio estimation of disease group covariate

tmp = prcomp(t(NovaSeq_Cortex))
plot(tmp$x[,1],tmp$x[,2])

customcols4 = data.frame(matrix(NA,nrow=ncol(NovaSeq_Cortex),ncol = 3))
colnames(customcols4) = c("Sample","Group","Color")
customcols4$Sample = gsub("-",".",colnames(NovaSeq_Cortex))

for(i in 1:nrow(customcols4)){
  for(j in 1:nrow(CNS_MEM)){
    if(customcols4$Sample[i] == CNS_MEM$Sample[j]){
      customcols4$Group[i] = CNS_MEM$Subtype[j]
    }
  }
}

for(i in 1:nrow(customcols4)){
  if(customcols4$Group[i] == "TD"){
    customcols4$Color[i] = "firebrick"
  }else if(customcols4$Group[i] == "OX"){
    customcols4$Color[i] = "navy"
  }else if(customcols4$Group[i] == "GLIA"){
    customcols4$Color[i] = "goldenrod1"
  }
}

plot(tmp$x[,1],tmp$x[,2],col = customcols4$Color,pch=20) 

#Check HISEQ cutoff = ~0.1

Spinal_HiSeq_Scores = read.csv("HiSeq_SpinalCord_nsNMF_AvgFeatureScores_9-27-23.csv")
HiSeq_Filt_ScoresSpinal = Spinal_HiSeq_Scores[Spinal_HiSeq_Scores$X %in% rownames(HiSeq_Spinal),]
table(HiSeq_Filt_ScoresSpinal$Rep1Group)
#Check NOVASEQ cutoff = ~0.05

Spinal_NovaSeq_Scores = read.csv("NovaSeq_SpinalCord_nsNMF_AvgFeatureScores_9-27-23.csv")
NovaSeq_Filt_ScoresSpinal = Spinal_NovaSeq_Scores[Spinal_NovaSeq_Scores$X %in% rownames(NovaSeq_Spinal),]
table(NovaSeq_Filt_ScoresSpinal$Rep1Group)


Cortex_HiSeq_Scores = read.csv("Cortex_HiSeq_SpinalCord_nsNMF_AvgFeatureScores_2021.csv")
HiSeq_Filt_ScoresCort = Cortex_HiSeq_Scores[Cortex_HiSeq_Scores$X %in% rownames(HiSeq_Cortex),]
table(HiSeq_Filt_ScoresCort$Rep1Group)
#Check NOVASEQ cutoff = ~0.05

Cortex_NovaSeq_Scores = read.csv("Cortex_NovaSeq_SpinalCord_nsNMF_AvgFeatureScores_2021.csv")
NovaSeq_Filt_ScoresCort = Cortex_NovaSeq_Scores[Cortex_NovaSeq_Scores$X %in% rownames(NovaSeq_Cortex),]
table(NovaSeq_Filt_ScoresCort$Rep1Group)

#Using fewer features still gives good separation
# test = NovaSeq_Spinal[rownames(NovaSeq_Spinal) %in% NovaSeq_Filt_Scores$X[1:500],]
# tmp = prcomp(t(test))
# plot(-1*tmp$x[,1],tmp$x[,2],col = customcols2$Color,pch=20) 

############################ Option 1 - Top 500 genes ##############################

NovaSeq_Filt_ScoresCort = NovaSeq_Filt_ScoresCort[order(NovaSeq_Filt_ScoresCort$AverageScore,decreasing = T),]
HiSeq_Filt_ScoresCort = HiSeq_Filt_ScoresCort[order(HiSeq_Filt_ScoresCort$AverageScore,decreasing = T),]
NovaSeq_Filt_ScoresSpinal = NovaSeq_Filt_ScoresSpinal[order(NovaSeq_Filt_ScoresSpinal$AverageScore,decreasing = T),]
HiSeq_Filt_ScoresSpinal = HiSeq_Filt_ScoresSpinal[order(HiSeq_Filt_ScoresSpinal$AverageScore,decreasing = T),]

Nova_genes_Cortex = NovaSeq_Filt_ScoresCort$X[1:500]
Nova_genes_Spinal = NovaSeq_Filt_ScoresSpinal$X[1:500]
Hi_genes_Cortex = HiSeq_Filt_ScoresCort$X[1:500]
Hi_genes_Spinal = HiSeq_Filt_ScoresSpinal$X[1:500]

V1_Nova_Cortex = NovaSeq_Cortex[rownames(NovaSeq_Cortex) %in% Nova_genes_Cortex,]
V1_Nova_Spinal = NovaSeq_Spinal[rownames(NovaSeq_Spinal) %in% Nova_genes_Spinal,]
V1_Hi_Cortex = HiSeq_Cortex[rownames(HiSeq_Cortex) %in% Hi_genes_Cortex,]
V1_Hi_Spinal = HiSeq_Spinal[rownames(HiSeq_Spinal) %in% Hi_genes_Spinal,]

#check overlap
table(rownames(V1_Nova_Cortex) %in% rownames(V1_Nova_Spinal))
table(rownames(V1_Hi_Cortex) %in% rownames(V1_Hi_Spinal)) 
table(rownames(V1_Nova_Cortex) %in% rownames(V1_Hi_Cortex))
table(rownames(V1_Nova_Spinal) %in% rownames(V1_Hi_Spinal))

#GET PCA COMPS
Nova_Cortex_PCA = prcomp(t(V1_Nova_Cortex))
plot(Nova_Cortex_PCA$x[,1],Nova_Cortex_PCA$x[,2],col=customcols4$Color,pch=20,cex=1.4,main="NovaSeq - Cortex PCA")
Nova_Spinal_PCA = prcomp(t(V1_Nova_Spinal))
plot(Nova_Spinal_PCA$x[,1],Nova_Spinal_PCA$x[,2],col=customcols2$Color,pch=20,cex=1.4,main="NovaSeq - Spinal Cord PCA")
Hi_Cortex_PCA = prcomp(t(V1_Hi_Cortex))
plot(Hi_Cortex_PCA$x[,1],Hi_Cortex_PCA$x[,2],col=customcols3$Color,pch=20,cex=1.4,main="HiSeq - Cortex PCA")
Hi_Spinal_PCA = prcomp(t(V1_Hi_Spinal))
plot(Hi_Spinal_PCA$x[,1],Hi_Spinal_PCA$x[,2],col=customcols$Color,pch=20,cex=1.4,main="HiSeq - Spinal Cord PCA")



#Use varimax to produce standardized scores 
ncomp = 2

#NovaSeq Cortex
pca_varimax_NC = prcomp(t(V1_Nova_Cortex), center=T, scale=T)
rawLoadings = pca_varimax_NC$rotation[,1:ncomp] %*% diag(pca_varimax_NC$sdev, ncomp, ncomp)
scores_NC = scale(pca_varimax_NC$x[,1:2]) %*% varimax(rawLoadings)$rotmat
plot(scores_NC,col=customcols4$Color,pch=20,cex=1.5,main="Varimax rotated PCA - Novaseq Cortex")

#NovaSeq Spinal Cord
pca_varimax_NS = prcomp(t(V1_Nova_Spinal), center=T, scale=T)
rawLoadings = pca_varimax_NS$rotation[,1:ncomp] %*% diag(pca_varimax_NS$sdev, ncomp, ncomp)
scores_NS = scale(pca_varimax_NS$x[,1:2]) %*% varimax(rawLoadings)$rotmat
plot(scores_NS,col=customcols2$Color,pch=20,cex=1.5,main="Varimax rotated PCA - Novaseq Spinal Cord")

#HiSeq Cortex
pca_varimax_HC = prcomp(t(V1_Hi_Cortex), center=T, scale=T)
rawLoadings     <- pca_varimax_HC$rotation[,1:ncomp] %*% diag(pca_varimax_HC$sdev, ncomp, ncomp)
scores_HC <- scale(pca_varimax_HC$x[,1:2]) %*% varimax(rawLoadings)$rotmat
plot(scores_HC,col=customcols3$Color,pch=20,cex=1.5,main="Varimax rotated PCA - HiSeq Cortex")

#HiSeq Spinal Cord
pca_varimax_HS = prcomp(t(V1_Hi_Spinal), center=T, scale=T)
rawLoadings     <- pca_varimax_HS$rotation[,1:ncomp] %*% diag(pca_varimax_HS$sdev, ncomp, ncomp)
scores_HS <- scale(pca_varimax_HS$x[,1:2]) %*% varimax(rawLoadings)$rotmat
plot(scores_HS,col=customcols$Color,pch=20,cex=1.5,main="Varimax rotated PCA - HiSeq Spinal Cord")

#rotate PCA spaces so that the spatial region corresponds to the same subtype across PCAs

#NovaSeq Cortex
angle = 25
rads = angle*(pi/180)
zshift = matrix(c(cos(rads),-sin(rads),sin(rads),cos(rads)),nrow=2)
RotStandardScores_NC = scores_NC %*% zshift
plot(RotStandardScores_NC,col=customcols4$Color,pch=20,cex=1.5,main="Varimax and 25 degree Z rotated PCA - NovaSeq Cortex")

#NovaSeq Spinal Cord
angle = -25
rads = angle*(pi/180)
zshift = matrix(c(cos(rads),-sin(rads),sin(rads),cos(rads)),nrow=2)
RotStandardScores_NS = scores_NS %*% zshift
plot(RotStandardScores_NS,col=customcols2$Color,pch=20,cex=1.5,main="Varimax and -25 degree Z rotated PCA - NovaSeq Spinal Cord")

#HiSeq Cortex
angle = 90
rads = angle*(pi/180)
zshift = matrix(c(cos(rads),-sin(rads),sin(rads),cos(rads)),nrow=2)
RotStandardScores_HC = scores_HC %*% zshift
plot(RotStandardScores_HC,col=customcols3$Color,pch=20,cex=1.5,main="Varimax and 90 degree Z rotated PCA - HiSeq Cortex")

#HiSeq Spinal Cord
angle = -30
rads = angle*(pi/180)
zshift = matrix(c(cos(rads),-sin(rads),sin(rads),cos(rads)),nrow=2)
RotStandardScores_HS = scores_HS %*% zshift
plot(RotStandardScores_HS,col=customcols$Color,pch=20,cex=1.5,main="Varimax and -30 degree Z rotated PCA - HiSeq Spinal Cord")

#Add to CNS_MEM

CNS_MEM$PC1 = NA
CNS_MEM$PC2 = NA

for(i in 1:nrow(CNS_MEM)){
  counter = 0
  NC_ind = which(rownames(RotStandardScores_NC) == CNS_MEM$Sample[i])
  NS_ind = which(rownames(RotStandardScores_NS) == CNS_MEM$Sample[i])
  HC_ind = which(rownames(RotStandardScores_HC) == CNS_MEM$Sample[i])
  HS_ind = which(rownames(RotStandardScores_HS) == CNS_MEM$Sample[i])

  if(length(NC_ind)>0){
    CNS_MEM$PC1[i] = RotStandardScores_NC[NC_ind,1]
    CNS_MEM$PC2[i] = RotStandardScores_NC[NC_ind,2]
    counter = counter+1
  }

  if(length(NS_ind)>0){
    CNS_MEM$PC1[i] = RotStandardScores_NS[NS_ind,1]
    CNS_MEM$PC2[i] = RotStandardScores_NS[NS_ind,2]
    counter = counter+1
  }

  if(length(HC_ind)>0){
    CNS_MEM$PC1[i] = RotStandardScores_HC[HC_ind,1]
    CNS_MEM$PC2[i] = RotStandardScores_HC[HC_ind,2]
    counter = counter+1
  }

  if(length(HS_ind)>0){
    CNS_MEM$PC1[i] = RotStandardScores_HS[HS_ind,1]
    CNS_MEM$PC2[i] = RotStandardScores_HS[HS_ind,2]
    counter = counter+1
  }

  if(counter>1){
    cat(paste("warning multiple matching row: ",i,", corresponding to sample:",CNS_MEM$Sample[i],sep = ""))
  }
}

#write.csv(CNS_MEM,"CNS_MEM_wPCAcomponents.csv")
#write.csv(CNS_MEM,"Supplemental_Dataset3.csv")

########################################### COX PH MODEL ##################################################################

#Prep Data
CNS_MEM$SubjectSex = as.factor(CNS_MEM$SubjectSex)
CNS_MEM$SubjectSex = relevel(CNS_MEM$SubjectSex,"Male")
CNS_MEM$disease_group = as.factor(CNS_MEM$disease_group)
CNS_MEM$disease_group = relevel(CNS_MEM$disease_group,"ALS-TDP")
CNS_MEM$Subtype = as.factor(CNS_MEM$Subtype)
CNS_MEM$Subtype = relevel(CNS_MEM$Subtype,"TD")
CNS_MEM$Site = as.factor(CNS_MEM$Site)
CNS_MEM$AgeOnset = as.numeric(CNS_MEM$AgeOnset)


#Cox prop hazard model
coxfit = coxph(Surv(time) ~ SubjectSex + AgeOnset + disease_group + PC1 + PC2 + cluster(PatientID), data = CNS_MEM)
summary(coxfit)

library(coxme)
#with 1|PatientID
coxfit = coxme(Surv(time) ~ SubjectSex + AgeOnset + disease_group + PC1 + PC2 + (1|PatientID) , data = CNS_MEM)
summary(coxfit)

#Similar results

###################################################### LEVEL GROUPING APPROACH #######################################################
#Too many combinations when examined at the sample level
#6 tissues, 3 groups, permuted = >700 unique combinations

###Must take Meta assignment approach - which falls into similar repeat measure limitations as the original survival

#Switch phenotype file to patient level (n=221)
CNS_MEM_Patient = data.frame(matrix(NA,nrow = length(table(CNS_MEM$PatientID)),ncol=7))
colnames(CNS_MEM_Patient) = c("Sex","AgeOnset","disease_group","time","CortexMeta","SpinalMeta","SubtypeComb")
rownames(CNS_MEM_Patient) = names(table(CNS_MEM$PatientID))

for(i in 1:nrow(CNS_MEM_Patient)){
  ind = which(CNS_MEM$PatientID == rownames(CNS_MEM_Patient)[i])
  
  if(length(ind)>1){
    
    keepind_sex = which(table(CNS_MEM$SubjectSex[ind])!= 0)
    keepind_age = which(table(CNS_MEM$AgeOnset[ind])!= 0)
    keepind_disease = which(table(CNS_MEM$disease_group[ind])!= 0)
    keepind_time = which(table(CNS_MEM$time[ind]) != 0)
    
    if(length(keepind_sex)==1){
      CNS_MEM_Patient$Sex[i] = names(keepind_sex)
    }else if(length(keepind_sex)<1){
      CNS_MEM_Patient$Sex[i] = NA
    }else{
      cat(paste("Warning: multiple sex levels matching patient",rownames(CNS_MEM_Patient)[i],"\n",sep=""))
    }
    
    if(length(keepind_age)==1){
      CNS_MEM_Patient$AgeOnset[i] = names(keepind_age)
    }else if(length(keepind_age)<1){
      CNS_MEM_Patient$AgeOnset[i] = NA 
    }else{
      cat(paste("Warning: multiple ages matching patient",rownames(CNS_MEM_Patient)[i],"\n",sep=""))
    }
    
    if(length(keepind_disease)==1){
      CNS_MEM_Patient$disease_group[i] = names(keepind_disease)
    }else if(length(keepind_disease)<1){
      CNS_MEM_Patient$disease_group[i] = NA 
    }else{
      cat(paste("Warning: multiple disease group levels matching patient",rownames(CNS_MEM_Patient)[i],"\n",sep=""))
    }
    
    if(length(keepind_time)==1){
      CNS_MEM_Patient$time[i] = names(keepind_time)
    }else if(length(keepind_time)<1){
      CNS_MEM_Patient$time[i] = NA 
    }else{
      cat(paste("Warning: multiple disease group levels matching patient",rownames(CNS_MEM_Patient)[i],"\n",sep=""))
    }
    
    
  }else{
    CNS_MEM_Patient$Sex[i] = CNS_MEM$SubjectSex[ind]
    CNS_MEM_Patient$AgeOnset[i] = CNS_MEM$AgeOnset[ind]
    CNS_MEM_Patient$disease_group[i] = CNS_MEM$disease_group[ind]
    CNS_MEM_Patient$time[i] = CNS_MEM$time[ind]
  }
}


#Majority Agreement for CNS-level subtype
cortextissues = c("Frontal Cortex","Lateral Motor Cortex","Medial Motor Cortex","Other Motor Cortex")
spinaltissues = c("Spinal_Cord_Cervical","Spinal_Cord_Thoracic","Spinal_Cord_Lumbar")

Cortex_MEM = CNS_MEM[CNS_MEM$tissue %in% cortextissues,]
Spinal_MEM = CNS_MEM[CNS_MEM$tissue %in% spinaltissues,]

for(i in 1:nrow(CNS_MEM_Patient)){
  
  ind_cort = which(Cortex_MEM$PatientID == rownames(CNS_MEM_Patient)[i])
  if(length(ind_cort)>0){
    stypes_cort =  Cortex_MEM$Subtype[ind_cort]
    tmp = names(table(stypes_cort)[which(table(stypes_cort) == max(table(stypes_cort)))])
    if(length(tmp) == 1){
      CNS_MEM_Patient$CortexMeta[i] = tmp
    }else{
      CNS_MEM_Patient$CortexMeta[i] = "Discordant"
    }
  }
  
  ind_spin = which(Spinal_MEM$PatientID == rownames(CNS_MEM_Patient)[i])
  if(length(ind_spin)>0){
    stypes_spin = Spinal_MEM$Subtype[ind_spin]
    tmp2 = names(table(stypes_spin)[which(table(stypes_spin) == max(table(stypes_spin)))])
    if(length(tmp2) == 1){
      CNS_MEM_Patient$SpinalMeta[i] = tmp2
    }else{
      CNS_MEM_Patient$SpinalMeta[i] = "Discordant"
    }
  }
  
}

Final_MEM = CNS_MEM_Patient[-which(is.na(CNS_MEM_Patient$CortexMeta)),]
Final_MEM = Final_MEM[-which(is.na(Final_MEM$SpinalMeta)),]

Final_MEM$Subtype = paste(Final_MEM$CortexMeta,Final_MEM$SpinalMeta,sep="")
Final_MEM$PatientID = rownames(Final_MEM)

#Prep Data
Final_MEM$Sex = as.factor(Final_MEM$Sex)
Final_MEM$Sex = relevel(Final_MEM$Sex,"Male")
Final_MEM$disease_group = as.factor(Final_MEM$disease_group)
Final_MEM$disease_group = relevel(Final_MEM$disease_group,"ALS-TDP")
Final_MEM$AgeOnset = as.numeric(Final_MEM$AgeOnset)
Final_MEM$Subtype = as.factor(Final_MEM$Subtype)
Final_MEM$PatientID = as.factor(Final_MEM$PatientID)
Final_MEM$time = as.numeric(Final_MEM$time)


#Cox prop hazard model

coxfit = coxph(Surv(time) ~ Sex + AgeOnset + disease_group + Subtype + cluster(PatientID), data = Final_MEM)
summary(coxfit)

library(coxme)
#with 1|PatientID
coxfit = coxme(Surv(time) ~ Sex + AgeOnset + disease_group + Subtype + (1|PatientID), data = Final_MEM)
summary(coxfit)
