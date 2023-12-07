

wd = "G:/SpinalCord/Publication/Manuscript/Tables"
setwd(wd)

###################################################################################################################################################################################################################################################################################################################################
#Run this section once.

Cortex = read.csv("G:/SpinalCord/Publication/SystemicAnalysis/Table_S10.csv")

CNSSubtype = read.csv("ALS_Cortex_SpinalCord_4Covar_SubtypeConcordance_9-18-23.csv")
#CNSSubtype = read.csv("ALS_Cortex_SpinalCord_NoGliaNoRIN_SubtypeConcordance_8-31-23.csv")

#Fill Assigned Subtype

CNSSubtype$CortexSubtype = NA

for(i in 1:nrow(CNSSubtype)){
  for(j in 1:nrow(Cortex)){
    if(CNSSubtype$Patient[i] == Cortex$SubjectID[j]){
      CNSSubtype$CortexSubtype[i] = Cortex$FinalSubtype[j]
    }
  }
}

Dat = CNSSubtype
Dat$SpinalSubtype = NA

for(i in 1:nrow(Dat)){
  if(length( which(table(c(Dat$CervicalSubtype[i],Dat$LumbarSubtype[i],Dat$ThoracicSubtype[i])) == max(table(c(Dat$CervicalSubtype[i],Dat$LumbarSubtype[i],Dat$ThoracicSubtype[i]))))) == 1){
    Dat$SpinalSubtype[i] = names(which(table(c(Dat$CervicalSubtype[i],Dat$LumbarSubtype[i],Dat$ThoracicSubtype[i])) == max(table(c(Dat$CervicalSubtype[i],Dat$LumbarSubtype[i],Dat$ThoracicSubtype[i])))))
  }else{
    Dat$SpinalSubtype[i] = "Discordant"
  }
}

for(i in 1:nrow(Dat)){
  for(j in 1:ncol(Dat)){
    if(!is.na(Dat[i,j])){
      if(Dat[i,j] == "TE"){
        Dat[i,j] = "TD"
      }
    }
  }
}

write.csv(Dat,"ALS_Cortex_SpinalCord_4Covar_SubtypeConcordance_AssignedSubtype_9-20-23.csv")

######################################################################################################################################
#Load in Output

#Without oligodendrocyte markers
#CNSSubtype = read.csv("G:/SpinalCord/Publication/SystemicAnalysis/NoOligo/ALS_Cortex_SpinalCordNoOligo_SubtypeConcordance_AssignedSubtype_8-16-23.csv")

#Without glial markers (Jack Humphrey paper)
#CNSSubtype = read.csv("G:/SpinalCord/Publication/Manuscript/Tables/ALS_Cortex_SpinalCord_NoGlia_SubtypeConcordance_AssignedSubtype_8-22-23.csv")

#Without glial markers (Jack Humphrey paper) and RIN-dependent genes
#CNSSubtype = read.csv("G:/SpinalCord/Publication/Manuscript/Tables/ALS_Cortex_SpinalCord_NoGliaNoRIN05_SubtypeConcordance_AssignedSubtype_8-31-23.csv")

#WITH glial markers (Jack Humphrey paper) and WITHOUT RIN-dependent genes (alpha = 0.05)
#CNSSubtype = read.csv("G:/SpinalCord/Publication/Manuscript/Tables/ALS_Cortex_SpinalCord_NoGliaNoRIN05_SubtypeConcordance_AssignedSubtype_9-20-23.csv")

#Sex, Site, Tissue, and RIN dependent genes expression removal
CNSSubtype = read.csv("G:/SpinalCord/Publication/Manuscript/Tables/ALS_Cortex_SpinalCord_4Covar_SubtypeConcordance_AssignedSubtype_9-20-23.csv")

CNSSubtype = CNSSubtype[-which(CNSSubtype$CortexSubtype == "Discordant"),]
CNSSubtype = CNSSubtype[-which(CNSSubtype$SpinalSubtype == "Discordant"),]
#Remove the patients without cortex samples 
CNSSubtype = CNSSubtype[-which(CNSSubtype$CortexSubtype == ""),]

###############################################  META SUBTYPE without discordant  #######################################################################################

#Concordance by spinal cord section - USING THE "META" SUBTYPE ASSIGNED BY CORTEX AGREEMENT (PREVIOUS PAPER)
concount = discount = 0
miscount = 0 #Add in a counter to keep track of unavailable tissue samples
GtoT = GtoO = TtoO = TtoG = OtoT = OtoG = 0
GoodG = GoodO = GoodT = 0
for(i in 1:nrow(CNSSubtype)){
  if(CNSSubtype$SpinalSubtype[i] == "" || CNSSubtype$CortexSubtype == ""){
    miscount = miscount+1
  }else if(CNSSubtype$CortexSubtype[i] == CNSSubtype$SpinalSubtype[i]){
    concount = concount+1
    if(CNSSubtype$SpinalSubtype[i] =="GLIA"){
      GoodG = GoodG+1
    }else if(CNSSubtype$SpinalSubtype[i] == "OX"){
      GoodO = GoodO+1
    }else if(CNSSubtype$SpinalSubtype[i] == "TD"){
      GoodT = GoodT+1
    }
  }else if(CNSSubtype$CortexSubtype[i] != CNSSubtype$SpinalSubtype[i]){
    discount = discount+1
    
    if(CNSSubtype$CortexSubtype[i] == "GLIA" && CNSSubtype$SpinalSubtype[i] == "TD"){
      GtoT = GtoT+1
    }else if(CNSSubtype$CortexSubtype[i] == "GLIA" && CNSSubtype$SpinalSubtype[i] == "OX"){
      GtoO = GtoO+1
    }
    
    if(CNSSubtype$CortexSubtype[i] == "TD" && CNSSubtype$SpinalSubtype[i] == "GLIA"){
      TtoG = TtoG+1
    }else if(CNSSubtype$CortexSubtype[i] == "TD" && CNSSubtype$SpinalSubtype[i] == "OX"){
      TtoO = TtoO+1
    }
    
    if(CNSSubtype$CortexSubtype[i] == "OX" && CNSSubtype$SpinalSubtype[i] == "TD"){
      OtoT = OtoT+1
    }else if(CNSSubtype$CortexSubtype[i] == "OX" && CNSSubtype$SpinalSubtype[i] == "GLIA"){
      OtoG = OtoG+1
    }
    
  }
  
}

#NoGlia Angle: 140
#NoGliaNoRin Angle: 110
ConcordantPercent = concount/sum(concount,discount)
DiscordantPercent = 1-ConcordantPercent
pie(c(ConcordantPercent*100,DiscordantPercent*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 95,main="Subtype Concordance - Cortex vs Spinal Cord - META",cex=2)
text(-0.5,0,"54",cex=2.5)
text(0.5,0,"61",cex=2.5)
table(CNSSubtype$CortexSubtype == CNSSubtype$SpinalSubtype)

##################################### Clustering Reassignment: Cortex vs Spinal Cord ###################################################################################################

################################# Cortex vs Cervical Spinal Cord

#ALS-Glia
a = GoodG/sum(GoodG,GtoO,GtoT)
b = GtoO/sum(GoodG,GtoO,GtoT)
c = GtoT/sum(GoodG,GtoO,GtoT)
pie(c(a,b,c),main="ALS-Glia Clustering: Cortex vs Spinal Cord - META",col = c("goldenrod1","chartreuse3","darkorange1"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 90,cex=1.25)
text(-0.5,0,"9",cex=2.5)
text(0.4,0.35,"5",cex=2.5)
text(0.4,-0.4,"4",cex=2.5)

#ALS-Ox
a = GoodO/sum(GoodO,OtoG,OtoT)
b = OtoG/sum(GoodO,OtoG,OtoT)
c = OtoT/sum(GoodO,OtoG,OtoT)
pie(c(a,b,c),main="ALS-Ox Clustering: Cortex vs Spinal Cord - META",col = c("navy","chartreuse3","#8d2ca3"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 110,cex=1.25)
text(-0.5,0,"22",cex=2.5,col="white")
text(0.45,0.25,"27",cex=2.5,col="white")
text(0.05,-0.6,"8",cex=2.5,col="white")

#ALS-TD
a = GoodT/sum(GoodT,TtoG,TtoO)
b = TtoG/sum(GoodT,TtoG,TtoO)
c = TtoO/sum(GoodT,TtoG,TtoO)
pie(c(a,b,c),main="ALS-TD Clustering: Cortex vs Spinal Cord - META",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 77,cex=1.25)
text(-0.5,0,"23",cex=2.5,col="white")
text(0.5,0.35,"9",cex=2.5,col="white")
text(0.475,-0.375,"8",cex=2.5,col="white")

######################################################################################################################################


###############################################  META SUBTYPE WITH discordant  #######################################################################################

#Sex, Site, Tissue, and RIN dependent genes expression removal
CNSSubtype = read.csv("G:/SpinalCord/Publication/Manuscript/Tables/ALS_Cortex_SpinalCord_4Covar_SubtypeConcordance_AssignedSubtype_9-20-23.csv")
CNSSubtype = CNSSubtype[-which(CNSSubtype$CortexSubtype == ""),]

#Concordance by spinal cord section - USING THE "META" SUBTYPE ASSIGNED BY CORTEX AGREEMENT (PREVIOUS PAPER)
concount = discount = 0
miscount = 0 #Add in a counter to keep track of unavailable tissue samples
GtoT = GtoO = TtoO = TtoG = OtoT = OtoG = DtoG = DtoO = DtoT = GtoD = OtoD = TtoD = 0
GoodG = GoodO = GoodT = GoodDis = 0
for(i in 1:nrow(CNSSubtype)){
  if(CNSSubtype$SpinalSubtype[i] == "" || CNSSubtype$CortexSubtype == ""){
    miscount = miscount+1
  }else if(CNSSubtype$CortexSubtype[i] == CNSSubtype$SpinalSubtype[i]){
    concount = concount+1
    if(CNSSubtype$SpinalSubtype[i] =="GLIA"){
      GoodG = GoodG+1
    }else if(CNSSubtype$SpinalSubtype[i] == "OX"){
      GoodO = GoodO+1
    }else if(CNSSubtype$SpinalSubtype[i] == "TD"){
      GoodT = GoodT+1
    }else if(CNSSubtype$SpinalSubtype[i] == "Discordant"){
      GoodDis = GoodDis+1
    }
  }else if(CNSSubtype$CortexSubtype[i] != CNSSubtype$SpinalSubtype[i]){
    discount = discount+1
    
    if(CNSSubtype$CortexSubtype[i] == "GLIA" && CNSSubtype$SpinalSubtype[i] == "TD"){
      GtoT = GtoT+1
    }else if(CNSSubtype$CortexSubtype[i] == "GLIA" && CNSSubtype$SpinalSubtype[i] == "OX"){
      GtoO = GtoO+1
    }else if(CNSSubtype$CortexSubtype[i] == "GLIA" && CNSSubtype$SpinalSubtype[i] == "Discordant"){
      GtoD = GtoD+1
    }
    
    if(CNSSubtype$CortexSubtype[i] == "TD" && CNSSubtype$SpinalSubtype[i] == "GLIA"){
      TtoG = TtoG+1
    }else if(CNSSubtype$CortexSubtype[i] == "TD" && CNSSubtype$SpinalSubtype[i] == "OX"){
      TtoO = TtoO+1
    }else if(CNSSubtype$CortexSubtype[i] == "TD" && CNSSubtype$SpinalSubtype[i] == "Discordant"){
      TtoD = TtoD+1
    }
    
    if(CNSSubtype$CortexSubtype[i] == "OX" && CNSSubtype$SpinalSubtype[i] == "TD"){
      OtoT = OtoT+1
    }else if(CNSSubtype$CortexSubtype[i] == "OX" && CNSSubtype$SpinalSubtype[i] == "GLIA"){
      OtoG = OtoG+1
    }else if(CNSSubtype$CortexSubtype[i] == "OX" && CNSSubtype$SpinalSubtype[i] == "Discordant"){
      OtoD = OtoD+1
    }
    
    if(CNSSubtype$CortexSubtype[i] == "Discordant" && CNSSubtype$SpinalSubtype[i] == "GLIA"){
      DtoG = DtoG+1
    }else if(CNSSubtype$CortexSubtype[i] == "Discordant" && CNSSubtype$SpinalSubtype[i] == "OX"){
      DtoO = DtoO+1
    }else if(CNSSubtype$CortexSubtype[i] == "Discordant" && CNSSubtype$SpinalSubtype[i] == "TD"){
      DtoT = DtoT+1
    }
    
  }
  
}

#NoGlia Angle: 140
#NoGliaNoRin Angle: 110
ConcordantPercent = concount/sum(concount,discount)
DiscordantPercent = 1-ConcordantPercent
pie(c(ConcordantPercent*100,DiscordantPercent*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 122,main="Subtype Concordance - Cortex vs Spinal Cord - META with Discordant")
text(-0.5,0,"61",cex=2.5)
text(0.5,0,"131",cex=2.5)
table(CNSSubtype$CortexSubtype == CNSSubtype$SpinalSubtype)

##################################### Clustering Reassignment: Cortex vs Spinal Cord ###################################################################################################

################################# Cortex vs Cervical Spinal Cord

#ALS-Glia
a = GoodG/sum(GoodG,GtoO,GtoT,GtoD)
b = GtoO/sum(GoodG,GtoO,GtoT,GtoD)
c = GtoT/sum(GoodG,GtoO,GtoT,GtoD)
d = GtoD/sum(GoodG,GtoO,GtoT,GtoD)
pie(c(a,b,c,d),main="ALS-Glia Clustering: Cortex vs Spinal Cord",col = c("goldenrod1","chartreuse3","darkorange1","gray50"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)","Glia (Cortex) to Discordant (Spinal)"),init.angle = 125)
text(-0.6,0,"9",cex=2.5)
text(-0.05,-0.6,"4",cex=2.5)
text(0.5,-0.4,"5",cex=2.5)
text(0.3,0.5,"10",cex=2.5)

#ALS-Ox
a = GoodO/sum(GoodO,OtoG,OtoT,OtoD)
b = OtoG/sum(GoodO,OtoG,OtoT,OtoD)
c = OtoT/sum(GoodO,OtoG,OtoT,OtoD)
d = OtoD/sum(GoodO,OtoG,OtoT,OtoD)
pie(c(a,b,c,d),main="ALS-Ox Clustering: Cortex vs Spinal Cord",col = c("navy","chartreuse3","#8d2ca3","gray50"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)","Ox (Cortex) to Discordant (Spinal)"),init.angle = 137)
text(-0.6,0,"22",cex=2.5,col="white")
text(-0.25,-0.6,"8",cex=2.5,col="white")
text(0.45,-0.375,"27",cex=2.5,col="white")
text(0.175,0.6,"29",cex=2.5,col="white")

#ALS-TD
a = GoodT/sum(GoodT,TtoG,TtoO,TtoD)
b = TtoG/sum(GoodT,TtoG,TtoO,TtoD)
c = TtoO/sum(GoodT,TtoG,TtoO,TtoD)
d = TtoD/sum(GoodT,TtoG,TtoO,TtoD)
pie(c(a,b,c,d),main="ALS-TD Clustering: Cortex vs Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3","gray50"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)","TD (Cortex) to Discordant (Spinal)"),init.angle = 95)
text(-0.6,0,"23",cex=2.5,col="white")
text(0.3,-0.55,"8",cex=2.5,col="white")
text(0.65,0,"9",cex=2.5,col="white")
text(0.275,0.575,"8",cex=2.5,col="white")

#Discordant
a = GoodDis/sum(GoodDis,DtoG,DtoO,DtoT)
b = DtoG/sum(GoodDis,DtoG,DtoO,DtoT)
c = DtoO/sum(GoodDis,DtoG,DtoO,DtoT)
d = DtoT/sum(GoodDis,DtoG,DtoO,DtoT)
pie(c(a,b,c,d),main="ALS-Discordant Clustering: Cortex vs Spinal Cord",col = c("gray50","#f2ed8d","#a5e3f2","#f7929e"),labels=c("Concordant","Discordant (Cortex) to Glia (Spinal)","Discordant (Cortex) to Ox (Spinal)","Discordant (Cortex) to TD (Spinal)"),init.angle = 140)
text(-0.6,0,"7",cex=2.5)
text(-0.1,-0.6,"6",cex=2.5)
text(0.6,-0.15,"8",cex=2.5)
text(0.05,0.625,"9",cex=2.5)

######################################################################################################################################


######################################################################################################################################
#Highly concordant patients

#Without glial markers (Jack Humphrey paper) and RIN-dependent genes
CNSSubtype = read.csv("G:/SpinalCord/Publication/Manuscript/Tables/ALS_Cortex_SpinalCord_4Covar_SubtypeConcordance_AssignedSubtype_9-20-23.csv")
#CNSSubtype$CortexSubtype[which(CNSSubtype$CortexSubtype == "TE")] = "TD"
#CNSSubtype = CNSSubtype[,-1]
rownames(CNSSubtype) = CNSSubtype$Patient
CNSSubtype = CNSSubtype[,-1]

#Full Agreement across all available tissues:
GLIAind = OXind = TDind = rep(NA,nrow(CNSSubtype))
for(i in 1:nrow(CNSSubtype)){
  if(length(table(unlist(unname(CNSSubtype[i,])))) == 2){
    if(names(table(unlist(unname(CNSSubtype[i,])))[2]) == "GLIA"){
      GLIAind[i] = i
    }else if(names(table(unlist(unname(CNSSubtype[i,])))[2]) == "OX"){
      OXind[i] = i
    }else if(names(table(unlist(unname(CNSSubtype[i,])))[2]) == "TD"){
      TDind[i] = i
    }
  }
}

GLIAind = GLIAind[!is.na(GLIAind)]
OXind = OXind[!is.na(OXind)]
TDind = TDind[!is.na(TDind)]

Glia_Concordant = rownames(CNSSubtype)[GLIAind]
OX_Concordant = rownames(CNSSubtype)[OXind]
TD_Concordant = rownames(CNSSubtype)[TDind]

setwd("G:/SpinalCord/Publication/SystemicAnalysis")
Cortex = read.csv("ALS_Part1_CleanCortexSubtypes_withSpinalPatients_9-19-23.csv")

#Concordant Patient Supplemental Tables:

Glia_Pheno = Cortex[Cortex$SubjectID %in% Glia_Concordant,]
Ox_Pheno = Cortex[Cortex$SubjectID %in% OX_Concordant,]
TD_Pheno = Cortex[Cortex$SubjectID %in% TD_Concordant,]


Glia_Spinal = CNSSubtype[GLIAind,]
Ox_Spinal = CNSSubtype[OXind,]
TD_Spinal = CNSSubtype[TDind,]

rownames(TD_Spinal) == TD_Pheno$SubjectID

Glia_Cortex = data.frame(matrix(NA,nrow(Glia_Pheno),ncol(Glia_Pheno)))
colnames(Glia_Cortex) = colnames(Glia_Pheno)
Ox_Cortex = data.frame(matrix(NA,nrow(Ox_Pheno),ncol(Ox_Pheno)))
colnames(Ox_Cortex) = colnames(Ox_Pheno)
TD_Cortex = data.frame(matrix(NA,nrow(TD_Pheno),ncol(TD_Pheno)))
colnames(TD_Cortex) = colnames(TD_Pheno)

Glia_Cortex$SubjectID = rownames(Glia_Spinal)
Ox_Cortex$SubjectID = rownames(Ox_Spinal)
TD_Cortex$SubjectID = rownames(TD_Spinal)

for(i in 1:nrow(Glia_Cortex)){
  for(j in 1:nrow(Glia_Pheno)){
    if(Glia_Cortex$SubjectID[i] == Glia_Pheno$SubjectID[j]){
      Glia_Cortex[i,] = Glia_Pheno[j,]
    }
  }
}

for(i in 1:nrow(Ox_Cortex)){
  for(j in 1:nrow(Ox_Pheno)){
    if(Ox_Cortex$SubjectID[i] == Ox_Pheno$SubjectID[j]){
      Ox_Cortex[i,] = Ox_Pheno[j,]
    }
  }
}

for(i in 1:nrow(TD_Cortex)){
  for(j in 1:nrow(TD_Pheno)){
    if(TD_Cortex$SubjectID[i] == TD_Pheno$SubjectID[j]){
      TD_Cortex[i,] = TD_Pheno[j,]
    }
  }
}


Glia_Pheno_Final = cbind(Glia_Cortex,Glia_Spinal)
Ox_Pheno_Final = cbind(Ox_Cortex,Ox_Spinal)
TD_Pheno_Final = cbind(TD_Cortex,TD_Spinal)


#Add sample IDs
setwd("G:/SpinalCord/Publication/SystemicAnalysis/4Covar")

HiSeq = read.csv("HiSeq_PatientPheno_4Covar_SpinalCord_SuppData_9-18-23.csv")
NovaSeq = read.csv("NovaSeq_PatientPheno_4Covar_SpinalCord_SuppData_9-18-23.csv")

Glia_Pheno_Final$CervicalSample = NA
Glia_Pheno_Final$ThoracicSample = NA
Glia_Pheno_Final$LumbarSample = NA

for(i in 1:nrow(Glia_Pheno_Final)){
  for(j in 1:nrow(HiSeq)){
    if(Glia_Pheno_Final$SubjectID[i] == HiSeq$ExternalSubjectId[j]){
      if(!is.na(HiSeq$CervicalSample[j])){
        Glia_Pheno_Final$CervicalSample[i] = HiSeq$CervicalSample[j]
      }
      if(!is.na(HiSeq$ThoracicSample[j])){
        Glia_Pheno_Final$ThoracicSample[i] = HiSeq$ThoracicSample[j]
      }
      
      if(!is.na(HiSeq$LumbarSample[j])){
        Glia_Pheno_Final$LumbarSample[i] = HiSeq$LumbarSample[j]
      }
      
    }
  }
  for(j in 1:nrow(NovaSeq)){
    if(Glia_Pheno_Final$SubjectID[i] == NovaSeq$ExternalSubjectId[j]){
      if(!is.na(NovaSeq$CervicalSample[j])){
        Glia_Pheno_Final$CervicalSample[i] = NovaSeq$CervicalSample[j]
      }
      if(!is.na(NovaSeq$ThoracicSample[j])){
        Glia_Pheno_Final$ThoracicSample[i] = NovaSeq$ThoracicSample[j]
      }
      if(!is.na(NovaSeq$LumbarSample[j])){
        Glia_Pheno_Final$LumbarSample[i] = NovaSeq$LumbarSample[j]
      }
      
    }
  }
}

write.csv(Glia_Pheno_Final,"ALS-Glia_HighlyConcordantCNS_4Covar_9-21-23.csv")



Ox_Pheno_Final$CervicalSample = NA
Ox_Pheno_Final$ThoracicSample = NA
Ox_Pheno_Final$LumbarSample = NA

for(i in 1:nrow(Ox_Pheno_Final)){
  for(j in 1:nrow(HiSeq)){
    if(Ox_Pheno_Final$SubjectID[i] == HiSeq$ExternalSubjectId[j]){
      if(!is.na(NovaSeq$CervicalSample[j])){
        Ox_Pheno_Final$CervicalSample[i] = HiSeq$CervicalSample[j]
      }
      
      if(!is.na(NovaSeq$ThoracicSample[j])){
        Ox_Pheno_Final$ThoracicSample[i] = HiSeq$ThoracicSample[j]
      }
      
      if(!is.na(NovaSeq$LumbarSample[j])){
        Ox_Pheno_Final$LumbarSample[i] = HiSeq$LumbarSample[j]
      }
      
    }
  }
  for(j in 1:nrow(NovaSeq)){
    if(Ox_Pheno_Final$SubjectID[i] == NovaSeq$ExternalSubjectId[j]){
      if(!is.na(NovaSeq$CervicalSample[j])){
        Ox_Pheno_Final$CervicalSample[i] = NovaSeq$CervicalSample[j]
      }
      if(!is.na(NovaSeq$ThoracicSample[j])){
        Ox_Pheno_Final$ThoracicSample[i] = NovaSeq$ThoracicSample[j]
      }
      if(!is.na(NovaSeq$LumbarSample[j])){
        Ox_Pheno_Final$LumbarSample[i] = NovaSeq$LumbarSample[j]
      }
      
    }
  }
}

write.csv(Ox_Pheno_Final,"ALS-OX_HighlyConcordantCNS_4Covar_9-21-23.csv")



TD_Pheno_Final$CervicalSample = NA
TD_Pheno_Final$ThoracicSample = NA
TD_Pheno_Final$LumbarSample = NA

for(i in 1:nrow(TD_Pheno_Final)){
  for(j in 1:nrow(HiSeq)){
    if(TD_Pheno_Final$SubjectID[i] == HiSeq$ExternalSubjectId[j]){
      if(!is.na(HiSeq$CervicalSample[j])){
        TD_Pheno_Final$CervicalSample[i] = HiSeq$CervicalSample[j]
      }
      
      if(!is.na(HiSeq$ThoracicSample[j])){
        TD_Pheno_Final$ThoracicSample[i] = HiSeq$ThoracicSample[j]
      }
      
      if(!is.na(HiSeq$LumbarSample[j])){
        TD_Pheno_Final$LumbarSample[i] = HiSeq$LumbarSample[j]
      }
      
    }
  }
  for(j in 1:nrow(NovaSeq)){
    if(TD_Pheno_Final$SubjectID[i] == NovaSeq$ExternalSubjectId[j]){
      if(!is.na(NovaSeq$CervicalSample[j])){
        TD_Pheno_Final$CervicalSample[i] = NovaSeq$CervicalSample[j]
      }
      
      if(!is.na(NovaSeq$ThoracicSample[j])){
        TD_Pheno_Final$ThoracicSample[i] = NovaSeq$ThoracicSample[j]
      }
      
      if(!is.na(NovaSeq$LumbarSample[j])){
        TD_Pheno_Final$LumbarSample[i] = NovaSeq$LumbarSample[j]
      }
      
    }
  }
}

write.csv(TD_Pheno_Final,"ALS-TD_HighlyConcordantCNS_4Covar_9-21-23.csv")



####### AT THE TISSUE LEVEL
#Use ConcordanceMeta_TISSUE_LEVEL.R
