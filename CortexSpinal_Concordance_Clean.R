library(stringr)

wd = "G:/SpinalCord/Publication/SystemicAnalysis"
setwd(wd)

##############################################################################################################################################################################################################
#PATIENT IDs

SC_Nova = read.csv("SpinalCord_NovaSeq_RobustSubtypeAssignment_4Covariate_9-18-23_majority.csv") #file derived from previous script

#colnames(SC_Nova) = str_sub(colnames(SC_Nova),6,nchar(colnames(SC_Nova))) 
rownames(SC_Nova) = SC_Nova[,1]; SC_Nova = SC_Nova[,-1]

#CGND.HRA.00290 and CGND.HRA.00042 are genuine duplicates (same patient, same tissue, different sequencing platforms)
#CGND.HRA.00290 is dropped so that all patient samples are analyzed on the same sequencing platform
which(colnames(SC_Nova) == "CGND.HRA.00290")
SC_Nova = SC_Nova[,-49]

#Spinal Cord Pie Chart
tmp = t(SC_Nova)
pie(table(tmp[,12]),col=c("goldenrod1","navy","firebrick"),main="NovaSeq Spinal Cord Labels")

SC_Hi = read.csv("SpinalCord_HiSeq_RobustSubtypeAssignment_4Covariate_9-18-23_majority.csv") #file derived from previous script

#colnames(SC_Hi) = str_sub(colnames(SC_Hi),6,nchar(colnames(SC_Hi)))
rownames(SC_Hi) = SC_Hi[,1]; SC_Hi = SC_Hi[,-1]

#Spinal Cord Pie Chart
tmp2 = t(SC_Hi)
pie(table(tmp2[,12]),col=c("goldenrod1","navy","firebrick"),main="HiSeq Spinal Cord Labels")

#Spinal Cord Pie Chart
pie(table(c(tmp[,12],tmp2[,12])),col=c("goldenrod1","navy","firebrick"),main="Combined Platform Spinal Cord Labels")


Cortex = read.csv("Table_S10.csv") #Available at: https://www.nature.com/articles/s41467-022-35494-w#Sec33
Pheno = read.csv("CLINICAL_DATA_PRUDENCIO.csv") #By request from NYGC

Pheno$ExternalSampleId = gsub("-",".",Pheno$ExternalSampleId)

#Get Novaseq PATIENT IDs
NovaPatients = rep(NA,ncol(SC_Nova))
passcount = 0
for(i in 1:ncol(SC_Nova)){
  
  for(j in 1:nrow(Pheno)){
    
    if(colnames(SC_Nova)[i] == Pheno$ExternalSampleId[j]){
      NovaPatients[i] = Pheno$ExternalSubjectId[j]
      passcount = passcount+1
    }
    
    
  }
  
  
}

CleanNovaPatients = names(table(NovaPatients))

#Get HiSeq PATIENT IDs
HiPatients = rep(NA,ncol(SC_Hi))
passcount = 0
for(i in 1:ncol(SC_Hi)){
  
  for(j in 1:nrow(Pheno)){
    
    if(colnames(SC_Hi)[i] == Pheno$ExternalSampleId[j]){
      HiPatients[i] = Pheno$ExternalSubjectId[j]
      passcount = passcount+1
    }
    
    
  }
  
  
}

CleanHiPatients = names(table(HiPatients))

##############################################################################################################################################################################################################
wd = "G:/SpinalCord/Publication/SystemicAnalysis/4Covar"
setwd(wd)

#HiSeq

#Filter phenotype file

SpinalCordLabels = Pheno[Pheno$ExternalSubjectId %in% CleanHiPatients,]

table(SpinalCordLabels$Sample.Source)

SpinalString = c("Spinal_Cord_Cervical","Spinal_Cord_Lumbar","Spinal_Cord_Thoracic")

CleanSpinalCordLabels = SpinalCordLabels[SpinalCordLabels$Sample.Source %in% SpinalString,]

CleanSpinalCordLabels2 = CleanSpinalCordLabels[CleanSpinalCordLabels$ExternalSampleId %in% colnames(SC_Hi),]

FinalLabels = CleanSpinalCordLabels2


Cerv = FinalLabels[which(FinalLabels$Sample.Source == "Spinal_Cord_Cervical"),]
max(table(Cerv$ExternalSubjectId))
Lumb = FinalLabels[which(FinalLabels$Sample.Source == "Spinal_Cord_Lumbar"),]
max(table(Lumb$ExternalSubjectId))
Thor = FinalLabels[which(FinalLabels$Sample.Source == "Spinal_Cord_Thoracic"),]
max(table(Thor$ExternalSubjectId))

FinalLabels$Subtype = NA

for(i in 1:ncol(SC_Hi)){
  for(j in 1:nrow(FinalLabels)){
    if(colnames(SC_Hi)[i] == FinalLabels$ExternalSampleId[j]){
      FinalLabels$Subtype[j] = SC_Hi[12,i]
    }
  }
}


ref = data.frame(matrix(rbind(FinalLabels$ExternalSubjectId,FinalLabels$ExternalSampleId,FinalLabels$Subtype,FinalLabels$Sample.Source),nrow=4))



FinalPheno = data.frame(matrix(NA,nrow = length(CleanHiPatients),ncol = 7))
colnames(FinalPheno) = c("ExternalSubjectId","Sample1","Sample2","Sample3","CervicalSubtype","LumbarSubtype","ThoracicSubtype")
FinalPheno$ExternalSubjectId = CleanHiPatients

for(i in 1:nrow(FinalPheno)){
  for(j in 1:ncol(ref)){
    
    filldat = data.frame(ref[,which(ref[1,] == FinalPheno$ExternalSubjectId[i])])
    
    for(k in 1:ncol(filldat)){
      
      if(ncol(filldat) == 1){
        FinalPheno$Sample1[i] = filldat[2,1]
        
        if(filldat[4,1] == "Spinal_Cord_Cervical"){
          FinalPheno$CervicalSubtype[i] = filldat[3,1]
        }else if(filldat[4,1] == "Spinal_Cord_Lumbar"){
          FinalPheno$LumbarSubtype[i] = filldat[3,1]
        }else if(filldat[4,1] == "Spinal_Cord_Thoracic"){
          FinalPheno$ThoracicSubtype[i] = filldat[3,1]
        }
        
      }else if(ncol(filldat) == 2){
        FinalPheno$Sample1[i] = filldat[2,1]
        FinalPheno$Sample2[i] = filldat[2,2]
        
        if(filldat[4,1] == "Spinal_Cord_Cervical"){
          FinalPheno$CervicalSubtype[i] = filldat[3,1]
        }else if(filldat[4,1] == "Spinal_Cord_Lumbar"){
          FinalPheno$LumbarSubtype[i] = filldat[3,1]
        }else if(filldat[4,1] == "Spinal_Cord_Thoracic"){
          FinalPheno$ThoracicSubtype[i] = filldat[3,1]
        }
        
        if(filldat[4,2] == "Spinal_Cord_Cervical"){
          FinalPheno$CervicalSubtype[i] = filldat[3,2]
        }else if(filldat[4,2] == "Spinal_Cord_Lumbar"){
          FinalPheno$LumbarSubtype[i] = filldat[3,2]
        }else if(filldat[4,2] == "Spinal_Cord_Thoracic"){
          FinalPheno$ThoracicSubtype[i] = filldat[3,2]
        }
        
        
      }else if(ncol(filldat) == 3){
        FinalPheno$Sample1[i] = filldat[2,1]
        FinalPheno$Sample2[i] = filldat[2,2]
        FinalPheno$Sample3[i] = filldat[2,3]
        
        if(filldat[4,1] == "Spinal_Cord_Cervical"){
          FinalPheno$CervicalSubtype[i] = filldat[3,1]
        }else if(filldat[4,1] == "Spinal_Cord_Lumbar"){
          FinalPheno$LumbarSubtype[i] = filldat[3,1]
        }else if(filldat[4,1] == "Spinal_Cord_Thoracic"){
          FinalPheno$ThoracicSubtype[i] = filldat[3,1]
        }
        
        if(filldat[4,2] == "Spinal_Cord_Cervical"){
          FinalPheno$CervicalSubtype[i] = filldat[3,2]
        }else if(filldat[4,2] == "Spinal_Cord_Lumbar"){
          FinalPheno$LumbarSubtype[i] = filldat[3,2]
        }else if(filldat[4,2] == "Spinal_Cord_Thoracic"){
          FinalPheno$ThoracicSubtype[i] = filldat[3,2]
        }
        
        if(filldat[4,3] == "Spinal_Cord_Cervical"){
          FinalPheno$CervicalSubtype[i] = filldat[3,3]
        }else if(filldat[4,3] == "Spinal_Cord_Lumbar"){
          FinalPheno$LumbarSubtype[i] = filldat[3,3]
        }else if(filldat[4,3] == "Spinal_Cord_Thoracic"){
          FinalPheno$ThoracicSubtype[i] = filldat[3,3]
        }
      }
      
    }
    
  }
  
  if((i %% 1) == 0) cat("Patient Completed: ",i,"\n")
}



FinalPheno$Sex = NA
FinalPheno$Subject.Group = NA
FinalPheno$Subject.Group.Subcategory = NA
FinalPheno$GEO = NA
FinalPheno$SiteMotorOnset = NA
FinalPheno$SiteMotorOnsetDetail = NA
FinalPheno$AgeOnset = NA
FinalPheno$AgeDeath = NA
FinalPheno$SampleSource = NA
#FinalPheno$RIN = NA #Cant have patient-wise RINs
FinalPheno$Duration = NA
FinalPheno$Genotype = NA
FinalPheno$Prep = NA
FinalPheno$DiseaseGroup = NA


for(i in 1:nrow(FinalPheno)){
  
  for(j in 1:nrow(Pheno)){
    
    if(FinalPheno$ExternalSubjectId[i] == Pheno$ExternalSubjectId[j]){
      FinalPheno$Sex[i] = Pheno$Sex[j]
      FinalPheno$Subject.Group[i] = Pheno$Subject.Group[j]
      FinalPheno$Subject.Group.Subcategory[i] = Pheno$Subject.Group.Subcategory[j]
      FinalPheno$GEO[i] = Pheno$GEO[j]
      FinalPheno$SiteMotorOnset[i] = Pheno$Site.of.Motor.Onset[j]
      FinalPheno$SiteMotorOnsetDetail[i] = Pheno$Site.of.Motor.Onset.Detail[j]
      FinalPheno$AgeOnset[i] = Pheno$Age.at.Symptom.Onset[j]
      FinalPheno$AgeDeath[i] = Pheno$Age.at.Death[j]
      FinalPheno$SampleSource[i] = Pheno$Sample.Source[j]
      #FinalPheno$RIN[i] = Pheno$RIN[j]
      FinalPheno$Duration[i] = Pheno$Disease.Duration.in.Months[j]
      FinalPheno$Genotype[i] = Pheno$Reported.Genomic.Mutations..from.sites...NOT.in.any.way.associated.with.WGS.data.from.NYGC.[j]
      FinalPheno$Prep[i] = Pheno$Prep[j]
      #FinalPheno$Platform[i] = Pheno$Platform[j]
      FinalPheno$DiseaseGroup[i] = Pheno$disease_group[j]
      #FinalPheno$LibrarySize[i] = Pheno$library_size[j]
    }
    
    
  }
  if((i %% 1) == 0) cat("Patient Completed: ",i,"\n")
}

#Check that all patient samples are in the phenotype file
sum(!is.na(c(FinalPheno$Sample1,FinalPheno$Sample2,FinalPheno$Sample3))) == length(HiPatients)

write.csv(FinalPheno,"HiSeq_SubtypePheno_4Covar_9-18-23.csv")


#Majority Subtype assignment (not IID as determined from discussion with Drs. Plaisier and Fricks)
Dat = FinalPheno
Dat$Majority = NA

#This is the majority agreement approach
for(i in 1:nrow(Dat)){
  if(length( which(table(c(Dat$CervicalSubtype[i],Dat$LumbarSubtype[i],Dat$ThoracicSubtype[i])) == max(table(c(Dat$CervicalSubtype[i],Dat$LumbarSubtype[i],Dat$ThoracicSubtype[i]))))) == 1){
    Dat$Majority[i] = names(which(table(c(Dat$CervicalSubtype[i],Dat$LumbarSubtype[i],Dat$ThoracicSubtype[i])) == max(table(c(Dat$CervicalSubtype[i],Dat$LumbarSubtype[i],Dat$ThoracicSubtype[i])))))
  }else{
    Dat$Majority[i] = "Discordant"
  }
}

write.csv(Dat,"HiSeq_SubtypePheno_4Covar_9-18-23_MajorityAgreement.csv")

##############################################################################################################################################################################################################
#NovaSeq

#Filter phenotype file

SpinalCordLabels = Pheno[Pheno$ExternalSubjectId %in% CleanNovaPatients,]

table(SpinalCordLabels$Sample.Source)

SpinalString = c("Spinal_Cord_Cervical","Spinal_Cord_Lumbar","Spinal_Cord_Thoracic")

CleanSpinalCordLabels = SpinalCordLabels[SpinalCordLabels$Sample.Source %in% SpinalString,]

CleanSpinalCordLabels2 = CleanSpinalCordLabels[CleanSpinalCordLabels$ExternalSampleId %in% colnames(SC_Nova),]

FinalLabels = CleanSpinalCordLabels2


Cerv = FinalLabels[which(FinalLabels$Sample.Source == "Spinal_Cord_Cervical"),]
max(table(Cerv$ExternalSubjectId))
Lumb = FinalLabels[which(FinalLabels$Sample.Source == "Spinal_Cord_Lumbar"),]
max(table(Lumb$ExternalSubjectId))
Thor = FinalLabels[which(FinalLabels$Sample.Source == "Spinal_Cord_Thoracic"),]
max(table(Thor$ExternalSubjectId))

FinalLabels$Subtype = NA

for(i in 1:ncol(SC_Nova)){
  for(j in 1:nrow(FinalLabels)){
    if(colnames(SC_Nova)[i] == FinalLabels$ExternalSampleId[j]){
      FinalLabels$Subtype[j] = SC_Nova[12,i]
    }
  }
}


ref = data.frame(matrix(rbind(FinalLabels$ExternalSubjectId,FinalLabels$ExternalSampleId,FinalLabels$Subtype,FinalLabels$Sample.Source),nrow=4))



FinalPheno = data.frame(matrix(NA,nrow = length(CleanNovaPatients),ncol = 7))
colnames(FinalPheno) = c("ExternalSubjectId","Sample1","Sample2","Sample3","CervicalSubtype","LumbarSubtype","ThoracicSubtype")
FinalPheno$ExternalSubjectId = CleanNovaPatients

for(i in 1:nrow(FinalPheno)){
  for(j in 1:ncol(ref)){
    
    filldat = data.frame(ref[,which(ref[1,] == FinalPheno$ExternalSubjectId[i])])
    
    for(k in 1:ncol(filldat)){
      
      if(ncol(filldat) == 1){
        FinalPheno$Sample1[i] = filldat[2,1]
        
        if(filldat[4,1] == "Spinal_Cord_Cervical"){
          FinalPheno$CervicalSubtype[i] = filldat[3,1]
        }else if(filldat[4,1] == "Spinal_Cord_Lumbar"){
          FinalPheno$LumbarSubtype[i] = filldat[3,1]
        }else if(filldat[4,1] == "Spinal_Cord_Thoracic"){
          FinalPheno$ThoracicSubtype[i] = filldat[3,1]
        }
        
      }else if(ncol(filldat) == 2){
        FinalPheno$Sample1[i] = filldat[2,1]
        FinalPheno$Sample2[i] = filldat[2,2]
        
        if(filldat[4,1] == "Spinal_Cord_Cervical"){
          FinalPheno$CervicalSubtype[i] = filldat[3,1]
        }else if(filldat[4,1] == "Spinal_Cord_Lumbar"){
          FinalPheno$LumbarSubtype[i] = filldat[3,1]
        }else if(filldat[4,1] == "Spinal_Cord_Thoracic"){
          FinalPheno$ThoracicSubtype[i] = filldat[3,1]
        }
        
        if(filldat[4,2] == "Spinal_Cord_Cervical"){
          FinalPheno$CervicalSubtype[i] = filldat[3,2]
        }else if(filldat[4,2] == "Spinal_Cord_Lumbar"){
          FinalPheno$LumbarSubtype[i] = filldat[3,2]
        }else if(filldat[4,2] == "Spinal_Cord_Thoracic"){
          FinalPheno$ThoracicSubtype[i] = filldat[3,2]
        }
        
        
      }else if(ncol(filldat) == 3){
        FinalPheno$Sample1[i] = filldat[2,1]
        FinalPheno$Sample2[i] = filldat[2,2]
        FinalPheno$Sample3[i] = filldat[2,3]
        
        if(filldat[4,1] == "Spinal_Cord_Cervical"){
          FinalPheno$CervicalSubtype[i] = filldat[3,1]
        }else if(filldat[4,1] == "Spinal_Cord_Lumbar"){
          FinalPheno$LumbarSubtype[i] = filldat[3,1]
        }else if(filldat[4,1] == "Spinal_Cord_Thoracic"){
          FinalPheno$ThoracicSubtype[i] = filldat[3,1]
        }
        
        if(filldat[4,2] == "Spinal_Cord_Cervical"){
          FinalPheno$CervicalSubtype[i] = filldat[3,2]
        }else if(filldat[4,2] == "Spinal_Cord_Lumbar"){
          FinalPheno$LumbarSubtype[i] = filldat[3,2]
        }else if(filldat[4,2] == "Spinal_Cord_Thoracic"){
          FinalPheno$ThoracicSubtype[i] = filldat[3,2]
        }
        
        if(filldat[4,3] == "Spinal_Cord_Cervical"){
          FinalPheno$CervicalSubtype[i] = filldat[3,3]
        }else if(filldat[4,3] == "Spinal_Cord_Lumbar"){
          FinalPheno$LumbarSubtype[i] = filldat[3,3]
        }else if(filldat[4,3] == "Spinal_Cord_Thoracic"){
          FinalPheno$ThoracicSubtype[i] = filldat[3,3]
        }
      }
      
    }
    
  }
  
  if((i %% 1) == 0) cat("Patient Completed: ",i,"\n")
}



FinalPheno$Sex = NA
FinalPheno$Subject.Group = NA
FinalPheno$Subject.Group.Subcategory = NA
FinalPheno$GEO = NA
FinalPheno$SiteMotorOnset = NA
FinalPheno$SiteMotorOnsetDetail = NA
FinalPheno$AgeOnset = NA
FinalPheno$AgeDeath = NA
FinalPheno$SampleSource = NA
#FinalPheno$RIN = NA
FinalPheno$Duration = NA
FinalPheno$Genotype = NA
FinalPheno$Prep = NA
FinalPheno$DiseaseGroup = NA


for(i in 1:nrow(FinalPheno)){
  
  for(j in 1:nrow(Pheno)){
    
    if(FinalPheno$ExternalSubjectId[i] == Pheno$ExternalSubjectId[j]){
      FinalPheno$Sex[i] = Pheno$Sex[j]
      FinalPheno$Subject.Group[i] = Pheno$Subject.Group[j]
      FinalPheno$Subject.Group.Subcategory[i] = Pheno$Subject.Group.Subcategory[j]
      FinalPheno$GEO[i] = Pheno$GEO[j]
      FinalPheno$SiteMotorOnset[i] = Pheno$Site.of.Motor.Onset[j]
      FinalPheno$SiteMotorOnsetDetail[i] = Pheno$Site.of.Motor.Onset.Detail[j]
      FinalPheno$AgeOnset[i] = Pheno$Age.at.Symptom.Onset[j]
      FinalPheno$AgeDeath[i] = Pheno$Age.at.Death[j]
      FinalPheno$SampleSource[i] = Pheno$Sample.Source[j]
      #FinalPheno$RIN[i] = Pheno$RIN[j]
      FinalPheno$Duration[i] = Pheno$Disease.Duration.in.Months[j]
      FinalPheno$Genotype[i] = Pheno$Reported.Genomic.Mutations..from.sites...NOT.in.any.way.associated.with.WGS.data.from.NYGC.[j]
      FinalPheno$Prep[i] = Pheno$Prep[j]
      #FinalPheno$Platform[i] = Pheno$Platform[j]
      FinalPheno$DiseaseGroup[i] = Pheno$disease_group[j]
      #FinalPheno$LibrarySize[i] = Pheno$library_size[j]
    }
    
    
  }
  if((i %% 1) == 0) cat("Patient Completed: ",i,"\n")
}

#Check that all patient samples are in the phenotype file
sum(!is.na(c(FinalPheno$Sample1,FinalPheno$Sample2,FinalPheno$Sample3))) == length(NovaPatients)

write.csv(FinalPheno,"NovaSeq_SubtypePheno_4Covar_9-18-23.csv")


#Majority Subtype assignment (not IID as determined from discussion with Drs. Plaisier and Fricks)
Dat = FinalPheno
Dat$Majority = NA

#This is the majority agreement approach
for(i in 1:nrow(Dat)){
  if(length( which(table(c(Dat$CervicalSubtype[i],Dat$LumbarSubtype[i],Dat$ThoracicSubtype[i])) == max(table(c(Dat$CervicalSubtype[i],Dat$LumbarSubtype[i],Dat$ThoracicSubtype[i]))))) == 1){
    Dat$Majority[i] = names(which(table(c(Dat$CervicalSubtype[i],Dat$LumbarSubtype[i],Dat$ThoracicSubtype[i])) == max(table(c(Dat$CervicalSubtype[i],Dat$LumbarSubtype[i],Dat$ThoracicSubtype[i])))))
  }else{
    Dat$Majority[i] = "Discordant"
  }
}

write.csv(Dat,"NovaSeq_SubtypePheno_4Covar_9-18-23_MajorityAgreement.csv")


#################################################################################################################################

#More informative sample labels
setwd("G:/SpinalCord/Publication/SystemicAnalysis")
Pheno = read.csv("CLINICAL_DATA_PRUDENCIO.csv")
setwd(wd)
Pheno$ExternalSampleId = gsub("-",".",Pheno$ExternalSampleId)

##NovaSeq
wd = "G:/SpinalCord/Publication/SystemicAnalysis/4Covar"
setwd(wd)
Dat = read.csv("NovaSeq_SubtypePheno_4Covar_9-18-23_MajorityAgreement.csv")

Dat$CervicalSample = NA
Dat$LumbarSample = NA
Dat$ThoracicSample = NA

for(i in 1:nrow(Dat)){
  
  psamps = c(Dat$Sample1[i],Dat$Sample2[i],Dat$Sample3[i])[which(nchar(c(Dat$Sample1[i],Dat$Sample2[i],Dat$Sample3[i]))>1)]
  
  for(j in 1:length(psamps)){
    
    tmpind = which(Pheno$ExternalSampleId == psamps[j])
    
    if(Pheno$Sample.Source[tmpind] == "Spinal_Cord_Cervical"){
      Dat$CervicalSample[i] = psamps[j]
    }else if(Pheno$Sample.Source[tmpind] == "Spinal_Cord_Lumbar"){
      Dat$LumbarSample[i] = psamps[j]
    }else if(Pheno$Sample.Source[tmpind] == "Spinal_Cord_Thoracic"){
      Dat$ThoracicSample[i] = psamps[j]
    }
    
  }
  
}

write.csv(Dat,"NovaSeq_PatientPheno_4Covar_SpinalCord_SuppData_9-18-23.csv")

NovaSeq_Pheno = Dat

##HiSeq

Dat = read.csv("HiSeq_SubtypePheno_4Covar_9-18-23_MajorityAgreement.csv")

Dat$CervicalSample = NA
Dat$LumbarSample = NA
Dat$ThoracicSample = NA

for(i in 1:nrow(Dat)){
  
  psamps = c(Dat$Sample1[i],Dat$Sample2[i],Dat$Sample3[i])[which(nchar(c(Dat$Sample1[i],Dat$Sample2[i],Dat$Sample3[i]))>1)]
  
  for(j in 1:length(psamps)){
    
    tmpind = which(Pheno$ExternalSampleId == psamps[j])
    
    if(Pheno$Sample.Source[tmpind] == "Spinal_Cord_Cervical"){
      Dat$CervicalSample[i] = psamps[j]
    }else if(Pheno$Sample.Source[tmpind] == "Spinal_Cord_Lumbar"){
      Dat$LumbarSample[i] = psamps[j]
    }else if(Pheno$Sample.Source[tmpind] == "Spinal_Cord_Thoracic"){
      Dat$ThoracicSample[i] = psamps[j]
    }
    
  }
  
}

write.csv(Dat,"HiSeq_PatientPheno_4Covar_SpinalCord_SuppData_9-18-23.csv")
HiSeq_Pheno = Dat

###################################################################################################################################################
#Clean up previous study subtype label file (Table S10)
setwd("G:/SpinalCord/Publication/SystemicAnalysis")
Cortex = read.csv("Table_S10.csv")
Pheno = read.csv("CLINICAL_DATA_PRUDENCIO.csv")
#Pheno$ExternalSampleId = gsub("-",".",Pheno$ExternalSampleId)


Dat = Cortex

Dat$FrontalSubtype = NA
Dat$FrontalSample = NA
Dat$MMSubtype = NA
Dat$MMSample = NA
Dat$LMSubtype = NA
Dat$LMSample = NA
Dat$UMSubtype = NA
Dat$UMSample = NA

for(i in 1:nrow(Dat)){
  
  psamps = c(Dat$RNAseqSample1[i],Dat$RNAseqSample2[i],Dat$RNAseqSample3[i])[which(nchar(c(Dat$RNAseqSample1[i],Dat$RNAseqSample2[i],Dat$RNAseqSample3[i]))>1)]
  
  for(j in 1:length(psamps)){
    
    tmpind = which(Pheno$ExternalSampleId == psamps[j])
    
    if(Pheno$Sample.Source[tmpind] == "Cortex_Frontal"){
      Dat$FrontalSample[i] = psamps[j]
      if(j == 1){
        Dat$FrontalSubtype[i] = Dat$SubtypeLabel1[i]
      }else if(j == 2){
        Dat$FrontalSubtype[i] = Dat$SubtypeLabel2[i]
      }else if(j == 3){
        Dat$FrontalSubtype[i] = Dat$SubtypeLabel3[i]
      }
      
    }else if(Pheno$Sample.Source[tmpind] == "Cortex_Motor_Lateral"){
      Dat$MMSample[i] = psamps[j]
      if(j == 1){
        Dat$MMSubtype[i] = Dat$SubtypeLabel1[i]
      }else if(j == 2){
        Dat$MMSubtype[i] = Dat$SubtypeLabel2[i]
      }else if(j == 3){
        Dat$MMSubtype[i] = Dat$SubtypeLabel3[i]
      }
    }else if(Pheno$Sample.Source[tmpind] == "Cortex_Motor_Medial"){
      Dat$LMSample[i] = psamps[j]
      if(j == 1){
        Dat$LMSubtype[i] = Dat$SubtypeLabel1[i]
      }else if(j == 2){
        Dat$LMSubtype[i] = Dat$SubtypeLabel2[i]
      }else if(j == 3){
        Dat$LMSubtype[i] = Dat$SubtypeLabel3[i]
      }
    }else if(Pheno$Sample.Source[tmpind] == "Cortex_Motor_Unspecified"){
      Dat$UMSample[i] = psamps[j]
      if(j == 1){
        Dat$UMSubtype[i] = Dat$SubtypeLabel1[i]
      }else if(j == 2){
        Dat$UMSubtype[i] = Dat$SubtypeLabel2[i]
      }else if(j == 3){
        Dat$UMSubtype[i] = Dat$SubtypeLabel3[i]
      }
    }
    
  }
  
}


nocortex = NovaSeq_Pheno$ExternalSubjectId[!NovaSeq_Pheno$ExternalSubjectId %in% Dat$SubjectID]
nocortexmat = data.frame(matrix(NA,nrow=length(nocortex),ncol = ncol(Dat)))
colnames(nocortexmat) = colnames(Dat)
nocortexmat$SubjectID = nocortex

for(i in 1:nrow(nocortexmat)){
  for(j in 1:nrow(Pheno)){
    if(nocortexmat$SubjectID[i] == Pheno$ExternalSubjectId[j]){
      nocortexmat$Sex[i] = Pheno$Sex[j]
      nocortexmat$Group[i] = Pheno$Subject.Group[j]
      nocortexmat$Subcategory[i] = Pheno$Subject.Group.Subcategory[j]
      nocortexmat$SiteofOnset[i] = Pheno$Site.of.Motor.Onset[j]
      nocortexmat$SiteExtra[i] = Pheno$Site.of.Motor.Onset.Detail[j]
      nocortexmat$AgeofOnset[i] = Pheno$Age.at.Symptom.Onset[j]
      nocortexmat$AgeofDeath[i] = Pheno$Age.at.Death[j]
      nocortexmat$Duration[i] = Pheno$Disease.Duration.in.Months[j]
      nocortexmat$Mutation[i] = Pheno$Reported.Genomic.Mutations..from.sites...NOT.in.any.way.associated.with.WGS.data.from.NYGC.[j]
      nocortexmat$DiseaseGroup[i] = Pheno$disease_group[j]
    }
  }
}

Dat2 = rbind(Dat,nocortexmat)

write.csv(Dat2,"ALS_Part1_CleanCortexSubtypes_withSpinalPatients_9-19-23.csv")
Cortex = Dat2

##########################################################################################################################
wd = "G:/SpinalCord/Publication/SystemicAnalysis/4Covar"
setwd(wd)

HiSeq = HiSeq_Pheno
NovaSeq = NovaSeq_Pheno

#Check that all patient samples are in the phenotype file
sum(!is.na(c(NovaSeq_Pheno$Sample1,NovaSeq_Pheno$Sample2,NovaSeq_Pheno$Sample3))) == length(NovaPatients)
sum(!is.na(c(HiSeq_Pheno$Sample1,HiSeq_Pheno$Sample2,HiSeq_Pheno$Sample3))) == length(HiPatients)

HiOverlap = Cortex$SubjectID[Cortex$SubjectID %in% HiSeq$ExternalSubjectId]
NovaOverlap = Cortex$SubjectID[Cortex$SubjectID %in% NovaSeq$ExternalSubjectId]

Overlap = c(HiOverlap,NovaOverlap)

CleanOverlap = names(table(Overlap))

Concordance = data.frame(matrix(NA,nrow = length(CleanOverlap),ncol=8))
colnames(Concordance) = c("Patient","FrontalCortex","MedialMotorCortex","LateralMotorCortex","UnspecMotorCortex","CervicalSubtype","LumbarSubtype","ThoracicSubtype")
Concordance$Patient = CleanOverlap

for(i in 1:nrow(Concordance)){
  
  for(j in 1:nrow(Cortex)){
    
    if(Concordance$Patient[i] == Cortex$SubjectID[j]){
      Concordance$FrontalCortex[i] = Cortex$FrontalSubtype[j]
      Concordance$MedialMotorCortex[i] = Cortex$MMSubtype[j]
      Concordance$LateralMotorCortex[i] = Cortex$LMSubtype[j]
      Concordance$UnspecMotorCortex[i] = Cortex$UMSubtype[j]
    }
    
  }
  
}

table(colnames(HiSeq) == colnames(NovaSeq))

#Cant directly row bind... Some patients are analyzed on both sequencing platforms. Sequencing-platform concordance is worth looking into later. 
#SpinalCord = rbind(NovaSeq,HiSeq)

Both = NovaSeq$ExternalSubjectId[NovaSeq$ExternalSubjectId %in% HiSeq$ExternalSubjectId]
tmpfix = HiSeq[HiSeq$ExternalSubjectId %in% Both,]
CleanNovaSeq = NovaSeq

for(i in 1:nrow(CleanNovaSeq)){
  for(j in 1:nrow(tmpfix)){
    if(CleanNovaSeq$ExternalSubjectId[i] == tmpfix$ExternalSubjectId[j]){
      if(!is.na(tmpfix$CervicalSample[j])){
        CleanNovaSeq$CervicalSample[i] = tmpfix$CervicalSample[j]
      }
      if(!is.na(tmpfix$CervicalSubtype[j])){
        CleanNovaSeq$CervicalSubtype[i] = tmpfix$CervicalSubtype[j]
      }
      if(!is.na(tmpfix$ThoracicSample[j])){
        CleanNovaSeq$ThoracicSample[i] = tmpfix$ThoracicSample[j]
      }
      if(!is.na(tmpfix$ThoracicSubtype[j])){
        CleanNovaSeq$ThoracicSubtype[i] = tmpfix$ThoracicSubtype[j]
      }
      if(!is.na(tmpfix$LumbarSample[j])){
        CleanNovaSeq$LumbarSample[i] = tmpfix$LumbarSample[j]
      }
      if(!is.na(tmpfix$LumbarSubtype[j])){
        CleanNovaSeq$LumbarSubtype[i] = tmpfix$LumbarSubtype[j]
      }
    }
  }
}
CleanHiSeq = HiSeq[! HiSeq$ExternalSubjectId %in% Both,]
SpinalCord = rbind(CleanNovaSeq,CleanHiSeq)
table(table(SpinalCord$ExternalSubjectId)>1) #Fixed :)

for(i in 1:nrow(Concordance)){
  for(j in 1:nrow(SpinalCord)){
    if(Concordance$Patient[i] == SpinalCord$ExternalSubjectId[j]){
      Concordance$CervicalSubtype[i] = SpinalCord$CervicalSubtype[j]
      Concordance$LumbarSubtype[i] = SpinalCord$LumbarSubtype[j]
      Concordance$ThoracicSubtype[i] = SpinalCord$ThoracicSubtype[j]
      
    }
  }
}

#Check that all patient samples are in the phenotype file
sum(!is.na(c(Concordance$CervicalSubtype,Concordance$ThoracicSubtype,Concordance$LumbarSubtype)))
sum(c(length(NovaPatients),length(HiPatients)))

#Single missing file:
# test = c(SpinalCord$CervicalSample,SpinalCord$ThoracicSample,SpinalCord$LumbarSample)
# test = test[!is.na(test)]
# allsamples = c(NovaSeq_Pheno$CervicalSample,NovaSeq_Pheno$ThoracicSample,NovaSeq_Pheno$LumbarSample,HiSeq_Pheno$CervicalSample,HiSeq_Pheno$ThoracicSample,HiSeq_Pheno$LumbarSample)
# allsamples = allsamples[!is.na(allsamples)]
# 
# table(table(allsamples)>1)
# allsamples[!allsamples %in% test]
# 
# #MISSING: CGND.HRA.00290
# which(NovaSeq_Pheno$CervicalSample == "CGND.HRA.00290")
# NovaSeq_Pheno[115,]
# NovaSeq_Pheno$ExternalSubjectId[115]
# 
# which(HiSeq_Pheno$ExternalSubjectId == "NEUYV496XLP")
# HiSeq_Pheno[64,]
# #Duplicate? CGND.HRA.00042
# 
# Pheno$ExternalSampleId = gsub("-",".",Pheno$ExternalSampleId)
# which(Pheno$ExternalSampleId == "CGND.HRA.00042")
# which(Pheno$ExternalSampleId == "CGND.HRA.00290")
# 
# Pheno[c(1481,1516),]
# 
# which(NovaSeq_Pheno$CervicalSample == "CGND.HRA.00290")
# NovaSeq_Pheno[115,]
# which(HiSeq_Pheno$CervicalSample == "CGND.HRA.00042")
# HiSeq_Pheno[64,]

#Genuine duplicate: delete one of the two samples from the clustering matrix - CHECK CONCORDANCE BETWEEN DUPLICATES


sum(!is.na(c(Concordance$CervicalSubtype,Concordance$ThoracicSubtype,Concordance$LumbarSubtype))) == sum(c(length(NovaPatients),length(HiPatients)))

write.csv(Concordance,"ALS_Cortex_SpinalCord_4Covar_SubtypeConcordance_9-18-23.csv")

FinalHiSeq = Concordance[Concordance$Patient %in% CleanHiSeq$ExternalSubjectId,]
write.csv(FinalHiSeq,"SpinalCord_Concordance_HiSeq_9-18-23.csv")
FinalNovaSeq = Concordance[Concordance$Patient %in% CleanNovaSeq$ExternalSubjectId,]
write.csv(FinalNovaSeq,"SpinalCord_Concordance_NovaSeq_9-18-23.csv")

#Proceed to ConcordanceMeta_TISSUE_and_PLATFORM_LEVEL.R

