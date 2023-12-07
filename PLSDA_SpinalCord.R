
## Spinal Cord classification using PLS-DA
#Written by: Jarrett Eshima
#Date: October 23, 2023
#Smith Research Lab

library(ggplot2)
library(rgl)
library(pls)
library(mixOmics)
library(MLmetrics)
library(edgeR)

#Split point plots - Ref: https://stackoverflow.com/questions/31538534/plotting-half-circles-in-r
upper.half.circle <- function(x,y,r,nsteps=100,...){  
  rs <- seq(0,pi,len=nsteps) 
  xc <- x+r*cos(rs) 
  yc <- y+r*sin(rs) 
  polygon(xc,yc,...) 
} 

lower.half.circle <- function(x,y,r,nsteps=100,...){ 
  rs <- seq(0,pi,len=nsteps) 
  xc <- x-r*cos(rs) 
  yc <- y-r*sin(rs) 
  polygon(xc,yc,...) 
} 

###########################################  LOAD DATA  ########################################################################
#Cortex Expression
load("G:/SpinalCord/Publication/DifferentialExpression/CNSPlots/ALSPatientStratification_UnivariateDatasets_RINSite_PeerReview.RData") #Available at: https://figshare.com/authors/Jarrett_Eshima/13813720

#Save Cortex Level Expression
Cortex_GO = filt.GO.sig
Cortex_GT = filt.GT.sig
Cortex_GC = filt.glia.sig
Cortex_TO = filt.TO.sig
Cortex_OC = filt.ox.sig
Cortex_TC = filt.TE.sig
Cortex_Gond = filt.glia.sig.ond
Cortex_Oond = filt.ox.sig.ond
Cortex_Tond = filt.TE.sig.ond
Cortex_Cond = filt.COND.sig

CortexExpr = NormCounts
CortexRaw = rCountData_rinsite
CortexPheno = FullPheno_sr

#Add patient ID to cortex Phenotype
clinicaldata = read.csv("G:/SpinalCord/Publication/CLINICAL_DATA_PRUDENCIO.csv") #Provided by NYGC by request
clinicaldata$ExternalSampleId = gsub("-",".",clinicaldata$ExternalSampleId)
ExtraCortex = clinicaldata[clinicaldata$ExternalSampleId %in% CortexPheno$Subject,]

CortexPheno$Patient = NA

for(i in 1:nrow(CortexPheno)){
  ind = which(ExtraCortex$ExternalSampleId == CortexPheno$Subject[i])
  CortexPheno$Patient[i] = ExtraCortex$ExternalSubjectId[ind]
}


#Spinal Cord Expression
load("G:/SpinalCord/Publication/DifferentialExpression/4Covar/SpinalCord_DifferentialExpression_FULLCohort_Top5000_10-18-23.RData") #Available at: https://figshare.com/authors/Jarrett_Eshima/13813720

Spinal_GO = filt.GO.sig2
Spinal_GT = filt.GT.sig2
Spinal_GC = filt.glia.sig2
Spinal_TO = filt.TO.sig2
Spinal_OC = filt.ox.sig2
Spinal_TC = filt.TE.sig2

SpinalExpr = NormCounts
SpinalRaw = rCountData_rinsite
SpinalPheno = FullPheno_sr

CleanPatients = names(table(c(CortexPheno$Patient,SpinalPheno$ExternalSubjectId)))

##################################### Three-Gene Classifier for ALS-Ox #############################################################################################################

#Candidates: SLC17A6, GLRA3, GAD2, GABRA1, B4GALT6, PCSK1, HTR2A

################## RPKM COUNTS (More robust against instrument/user/site variability - not good for cross sample comparisons as discussed in many recent papers but probably optimal for a cohort-blind classifier)

#PUBLICATION Fig. S12

CortexHoldoutGenes = c("HTR2A","SLC17A6","B4GALT6")
SpinalHoldout = c("GAD2","GLRA3","SLC17A6")
HiSeqHoldoutGenes = c("B4GALT6","GLRA3","SLC17A6")

#Use gold transcript colors (stable) where possible: https://useast.ensembl.org/Help/View?id=151
#GRCh38.p12 

#Estimate lengths from single transcript variant 
GABRA1Len = 4238 #NCBI Reference Sequence: NM_001127644.2; ENST00000393943.10
SLCLen = 3949 #UniProt: Q9P2U8; ENST00000263160.3
GAD2Len = 2462 #ENST00000259271.7
GLRALen = 8697 #ENST00000274093.8; NCBI Reference Sequence: NM_006529.4
PCSKLen = 5136 #ENST00000311106.8; NCBI Reference Sequence: NM_000439.5
HTRLen = 4835 #ENST00000378688.8
MYLLen = 2786 #ENST00000279022.7; NCBI Reference Sequence: NM_006097.5
GALNACLen = 2170 # ENST00000225276.9
TAGLen = 1477 #ENST00000278968.10
B4GALTLen = 2128 #ENST00000237019.11

HiSeqHoldout = F
CortexHoldout = T
iter = 100 #Number of cross-validation rounds
train = 0.8
test = 1-train
set.seed(123456)
testauc_c1 = testauc_c2 = testf1_c1 = testf1_c2 = rep(NA,iter)
testf1_c1_OX = testf1_c2_OX = testf1_c1_NO = testf1_c2_NO = rep(NA,iter)

############################################# Top 3 PUB ################################################################################

Candidates = SpinalHoldout

Cortex_Class_Raw = CortexRaw[rownames(CortexRaw) %in% Candidates,] #These are RAW (and not normalized expression counts)
Spinal_Class_Raw = SpinalRaw[rownames(SpinalRaw) %in% Candidates,] #These are RAW (and not normalized expression counts)

Cortex_Class_RPKM = Cortex_Class_Raw
Spinal_Class_RPKM = Spinal_Class_Raw #In the same order as Cortex

CortexPheno$library_size = as.numeric(CortexPheno$library_size)
SpinalPheno$library_size = as.numeric(SpinalPheno$library_size)

#Adjust first gene to RPKM scale
if(rownames(Cortex_Class_RPKM)[1] == "GLRA3"){
  Cortex_Class_RPKM[1,] = (Cortex_Class_Raw[1,]/(CortexPheno$library_size * GLRALen))*10^9
  Spinal_Class_RPKM[1,] = (Spinal_Class_Raw[1,]/(SpinalPheno$library_size * GLRALen))*10^9
}else if(rownames(Cortex_Class_RPKM)[1] == "PCSK1"){
  Cortex_Class_RPKM[1,] = (Cortex_Class_Raw[1,]/(CortexPheno$library_size * PCSKLen))*10^9
  Spinal_Class_RPKM[1,] = (Spinal_Class_Raw[1,]/(SpinalPheno$library_size * PCSKLen))*10^9
}else if(rownames(Cortex_Class_RPKM)[1] == "HTR2A"){
  Cortex_Class_RPKM[1,] = (Cortex_Class_Raw[1,]/(CortexPheno$library_size * HTRLen))*10^9
  Spinal_Class_RPKM[1,] = (Spinal_Class_Raw[1,]/(SpinalPheno$library_size * HTRLen))*10^9
}else if(rownames(Cortex_Class_RPKM)[1] == "GABRA1"){
  Cortex_Class_RPKM[1,] = (Cortex_Class_Raw[1,]/(CortexPheno$library_size * GABRA1Len))*10^9
  Spinal_Class_RPKM[1,] = (Spinal_Class_Raw[1,]/(SpinalPheno$library_size * GABRA1Len))*10^9
}else if(rownames(Cortex_Class_RPKM)[1] == "GAD2"){
  Cortex_Class_RPKM[1,] = (Cortex_Class_Raw[1,]/(CortexPheno$library_size * GAD2Len))*10^9
  Spinal_Class_RPKM[1,] = (Spinal_Class_Raw[1,]/(SpinalPheno$library_size * GAD2Len))*10^9
}else if(rownames(Cortex_Class_RPKM)[1] == "SLC17A6"){
  Cortex_Class_RPKM[1,] = (Cortex_Class_Raw[1,]/(CortexPheno$library_size * SLCLen))*10^9
  Spinal_Class_RPKM[1,] = (Spinal_Class_Raw[1,]/(SpinalPheno$library_size * SLCLen))*10^9
}else if(rownames(Cortex_Class_RPKM)[1] == "B4GALT6"){
  Cortex_Class_RPKM[1,] = (Cortex_Class_Raw[1,]/(CortexPheno$library_size * B4GALTLen))*10^9
  Spinal_Class_RPKM[1,] = (Spinal_Class_Raw[1,]/(SpinalPheno$library_size * B4GALTLen))*10^9
}

#Adjust second gene to RPKM scale
if(rownames(Cortex_Class_RPKM)[2] == "GLRA3"){
  Cortex_Class_RPKM[2,] = (Cortex_Class_Raw[2,]/(CortexPheno$library_size * GLRALen))*10^9
  Spinal_Class_RPKM[2,] = (Spinal_Class_Raw[2,]/(SpinalPheno$library_size * GLRALen))*10^9
}else if(rownames(Cortex_Class_RPKM)[2] == "PCSK1"){
  Cortex_Class_RPKM[2,] = (Cortex_Class_Raw[2,]/(CortexPheno$library_size * PCSKLen))*10^9
  Spinal_Class_RPKM[2,] = (Spinal_Class_Raw[2,]/(SpinalPheno$library_size * PCSKLen))*10^9
}else if(rownames(Cortex_Class_RPKM)[2] == "HTR2A"){
  Cortex_Class_RPKM[2,] = (Cortex_Class_Raw[2,]/(CortexPheno$library_size * HTRLen))*10^9
  Spinal_Class_RPKM[2,] = (Spinal_Class_Raw[2,]/(SpinalPheno$library_size * HTRLen))*10^9
}else if(rownames(Cortex_Class_RPKM)[2] == "GABRA1"){
  Cortex_Class_RPKM[2,] = (Cortex_Class_Raw[2,]/(CortexPheno$library_size * GABRA1Len))*10^9
  Spinal_Class_RPKM[2,] = (Spinal_Class_Raw[2,]/(SpinalPheno$library_size * GABRA1Len))*10^9
}else if(rownames(Cortex_Class_RPKM)[2] == "GAD2"){
  Cortex_Class_RPKM[2,] = (Cortex_Class_Raw[2,]/(CortexPheno$library_size * GAD2Len))*10^9
  Spinal_Class_RPKM[2,] = (Spinal_Class_Raw[2,]/(SpinalPheno$library_size * GAD2Len))*10^9
}else if(rownames(Cortex_Class_RPKM)[2] == "SLC17A6"){
  Cortex_Class_RPKM[2,] = (Cortex_Class_Raw[2,]/(CortexPheno$library_size * SLCLen))*10^9
  Spinal_Class_RPKM[2,] = (Spinal_Class_Raw[2,]/(SpinalPheno$library_size * SLCLen))*10^9
}else if(rownames(Cortex_Class_RPKM)[2] == "B4GALT6"){
  Cortex_Class_RPKM[2,] = (Cortex_Class_Raw[2,]/(CortexPheno$library_size * B4GALTLen))*10^9
  Spinal_Class_RPKM[2,] = (Spinal_Class_Raw[2,]/(SpinalPheno$library_size * B4GALTLen))*10^9
}

#Adjust third gene to RPKM scale
if(rownames(Cortex_Class_RPKM)[3] == "GLRA3"){
  Cortex_Class_RPKM[3,] = (Cortex_Class_Raw[3,]/(CortexPheno$library_size * GLRALen))*10^9
  Spinal_Class_RPKM[3,] = (Spinal_Class_Raw[3,]/(SpinalPheno$library_size * GLRALen))*10^9
}else if(rownames(Cortex_Class_RPKM)[3] == "PCSK1"){
  Cortex_Class_RPKM[3,] = (Cortex_Class_Raw[3,]/(CortexPheno$library_size * PCSKLen))*10^9
  Spinal_Class_RPKM[3,] = (Spinal_Class_Raw[3,]/(SpinalPheno$library_size * PCSKLen))*10^9
}else if(rownames(Cortex_Class_RPKM)[3] == "HTR2A"){
  Cortex_Class_RPKM[3,] = (Cortex_Class_Raw[3,]/(CortexPheno$library_size * HTRLen))*10^9
  Spinal_Class_RPKM[3,] = (Spinal_Class_Raw[3,]/(SpinalPheno$library_size * HTRLen))*10^9
}else if(rownames(Cortex_Class_RPKM)[3] == "GABRA1"){
  Cortex_Class_RPKM[3,] = (Cortex_Class_Raw[3,]/(CortexPheno$library_size * GABRA1Len))*10^9
  Spinal_Class_RPKM[3,] = (Spinal_Class_Raw[3,]/(SpinalPheno$library_size * GABRA1Len))*10^9
}else if(rownames(Cortex_Class_RPKM)[3] == "GAD2"){
  Cortex_Class_RPKM[3,] = (Cortex_Class_Raw[3,]/(CortexPheno$library_size * GAD2Len))*10^9
  Spinal_Class_RPKM[3,] = (Spinal_Class_Raw[3,]/(SpinalPheno$library_size * GAD2Len))*10^9
}else if(rownames(Cortex_Class_RPKM)[3] == "SLC17A6"){
  Cortex_Class_RPKM[3,] = (Cortex_Class_Raw[3,]/(CortexPheno$library_size * SLCLen))*10^9
  Spinal_Class_RPKM[3,] = (Spinal_Class_Raw[3,]/(SpinalPheno$library_size * SLCLen))*10^9
}else if(rownames(Cortex_Class_RPKM)[3] == "B4GALT6"){
  Cortex_Class_RPKM[3,] = (Cortex_Class_Raw[3,]/(CortexPheno$library_size * B4GALTLen))*10^9
  Spinal_Class_RPKM[3,] = (Spinal_Class_Raw[3,]/(SpinalPheno$library_size * B4GALTLen))*10^9
}

############################################# Top 7 PUB ################################################################################
#7-gene PLSDA
#PUBLICATION Fig. S13

Candidates7 = c("SLC17A6","GAD2","GABRA1","GLRA3","PCSK1","HTR2A","B4GALT6")
Cortex_Class_Raw = CortexRaw[rownames(CortexRaw) %in% Candidates7,] #These are RAW (and not normalized expression counts)
Spinal_Class_Raw = SpinalRaw[rownames(SpinalRaw) %in% Candidates7,] #These are RAW (and not normalized expression counts)

Cortex_Class_RPKM = Cortex_Class_Raw
Spinal_Class_RPKM = Spinal_Class_Raw #In the same order as Cortex

Cortex_Class_RPKM[1,] = (Cortex_Class_Raw[1,]/(CortexPheno$library_size * GABRA1Len))*10^9
Cortex_Class_RPKM[2,] = (Cortex_Class_Raw[2,]/(CortexPheno$library_size * SLCLen))*10^9
Cortex_Class_RPKM[3,] = (Cortex_Class_Raw[3,]/(CortexPheno$library_size * HTRLen))*10^9
Cortex_Class_RPKM[4,] = (Cortex_Class_Raw[4,]/(CortexPheno$library_size * B4GALTLen))*10^9
Cortex_Class_RPKM[5,] = (Cortex_Class_Raw[5,]/(CortexPheno$library_size * GAD2Len))*10^9
Cortex_Class_RPKM[6,] = (Cortex_Class_Raw[6,]/(CortexPheno$library_size * GLRALen))*10^9
Cortex_Class_RPKM[7,] = (Cortex_Class_Raw[7,]/(CortexPheno$library_size * PCSKLen))*10^9

Spinal_Class_RPKM[1,] = (Spinal_Class_Raw[1,]/(SpinalPheno$library_size * GABRA1Len))*10^9
Spinal_Class_RPKM[2,] = (Spinal_Class_Raw[2,]/(SpinalPheno$library_size * SLCLen))*10^9
Spinal_Class_RPKM[3,] = (Spinal_Class_Raw[3,]/(SpinalPheno$library_size * HTRLen))*10^9
Spinal_Class_RPKM[4,] = (Spinal_Class_Raw[4,]/(SpinalPheno$library_size * B4GALTLen))*10^9
Spinal_Class_RPKM[5,] = (Spinal_Class_Raw[5,]/(SpinalPheno$library_size * GAD2Len))*10^9
Spinal_Class_RPKM[6,] = (Spinal_Class_Raw[6,]/(SpinalPheno$library_size * GLRALen))*10^9
Spinal_Class_RPKM[7,] = (Spinal_Class_Raw[7,]/(SpinalPheno$library_size * PCSKLen))*10^9

Candidates = Candidates7
######################################################################################################################################

ncand = length(Candidates)
CortexMagRPKM = rep(NA,ncol(Cortex_Class_RPKM))
SpinalMagRPKM = rep(NA,ncol(Spinal_Class_RPKM))
for(i in 1:ncol(Cortex_Class_RPKM)){
  CortexMagRPKM[i] = sqrt(sum(Cortex_Class_RPKM[seq(1:ncand),i]^2)) #Verified manually
}

for(i in 1:ncol(Spinal_Class_RPKM)){
  SpinalMagRPKM[i] = sqrt(sum(Spinal_Class_RPKM[seq(1:ncand),i]^2)) #Verified manually
}

#CortexMagRPKM = sqrt(Cortex_Class_RPKM[1,]^2 + Cortex_Class_RPKM[2,]^2 + Cortex_Class_RPKM[3,]^2)
#SpinalMagRPKM = sqrt(Spinal_Class_RPKM[1,]^2 + Spinal_Class_RPKM[2,]^2 + Spinal_Class_RPKM[3,]^2)

names(CortexMagRPKM) = colnames(Cortex_Class_RPKM)
names(SpinalMagRPKM) = colnames(Spinal_Class_RPKM)

CNSClassifier = data.frame(matrix(NA,nrow=length(CleanPatients),ncol=2))
colnames(CNSClassifier) = c("CortexMagnitude","SpinalMagnitude")
rownames(CNSClassifier) = CleanPatients

#For patients with >1 sample per CNS region (Cortex, Spinal Cord), take the average of the magnitude

for(j in 1:nrow(CNSClassifier)){
  
  #Cortex
  sampind = which(CortexPheno$Patient == rownames(CNSClassifier)[j])
  
  if(length(sampind)>0){
    cortexsamples = CortexPheno$Subject[sampind]
    
    Mags = as.numeric(CortexMagRPKM[names(CortexMagRPKM) %in% cortexsamples])
    if(length(Mags)>1){
      fillCortexMag = mean(Mags)
    }else{
      fillCortexMag = Mags
    }
    
    CNSClassifier$CortexMagnitude[j] = fillCortexMag
  }
  
  #Spinal Cord
  sampind2 = which(SpinalPheno$ExternalSubjectId == rownames(CNSClassifier)[j])
  
  if(length(sampind2)>0){
    spinalsamples = SpinalPheno$ExternalSampleId[sampind2]
    
    Mags2 = as.numeric(SpinalMagRPKM[names(SpinalMagRPKM) %in% spinalsamples])
    if(length(Mags2)>1){
      fillSpinalMag = mean(Mags2)
    }else{
      fillSpinalMag = Mags2
    }
    
    CNSClassifier$SpinalMagnitude[j] = fillSpinalMag
  }
}

CNSClassifier_Clean = CNSClassifier[-which(is.na(apply(CNSClassifier,1,mean))),]

#Read in fully concordant patient data
ConcordantPatients = read.csv("G:/SpinalCord/Publication/Manuscript/Tables/Publication/Supplemental_Dataset2.csv")
OxPatients = ConcordantPatients$SubjectID[which(ConcordantPatients$MetaSpinalSubtype == "OX")]

AllLabels = read.csv("G:/SpinalCord/Publication/Manuscript/Tables/Publication/Supplemental_Dataset1.csv")


#Split Points - Fancy Plot

SplitCols = data.frame(matrix(NA,nrow = nrow(CNSClassifier_Clean),ncol = 2))
colnames(SplitCols) = c("CortexColor","SpinalColor")
rownames(SplitCols) = rownames(CNSClassifier_Clean)

for(j in 1:nrow(SplitCols)){
  
  ind = which(AllLabels$Patient == rownames(SplitCols)[j])
  ind2 = which(clinicaldata$ExternalSubjectId == rownames(SplitCols)[j])[1]
  
  if(length(ind)>0){
    if(AllLabels$MetaCortex[ind] == "GLIA"){
      SplitCols$CortexColor[j] = "goldenrod1"
    }else if(AllLabels$MetaCortex[ind] == "OX"){
      SplitCols$CortexColor[j] = "navy"
    }else if(AllLabels$MetaCortex[ind] == "TD"){
      SplitCols$CortexColor[j] = "firebrick"
    }else if(AllLabels$MetaCortex[ind] == "Discordant"){
      SplitCols$CortexColor[j] = "gray50"
    }
  }
  
  if(length(ind)>0){
    if(AllLabels$MetaSpinal[ind] == "GLIA"){
      SplitCols$SpinalColor[j] = "goldenrod1"
    }else if(AllLabels$MetaSpinal[ind] == "OX"){
      SplitCols$SpinalColor[j] = "navy"
    }else if(AllLabels$MetaSpinal[ind] == "TD"){
      SplitCols$SpinalColor[j] = "firebrick"
    }else if(AllLabels$MetaSpinal[ind] == "Discordant"){
      SplitCols$SpinalColor[j] = "gray50"
    }
  }
  
  if(length(ind2)>0){
    if(clinicaldata$Subject.Group[ind2] == "Non-Neurological Control"){
      SplitCols$CortexColor[j] = "white"
      SplitCols$SpinalColor[j] = "white"
    }else if(clinicaldata$Subject.Group[ind2] == "Other Neurological Disorders"){
      SplitCols$CortexColor[j] = "black"
      SplitCols$SpinalColor[j] = "black"
    }
  }
  
}

SplitCols[rownames(SplitCols) %in% OxPatients,] = c("#429ac9","#429ac9")


plot(1, type="n",xlab="CortexMagnitude", ylab="SpinalMagnitude",main = paste("RPKM Magnitude Classifier: ",Candidates,sep=""),xlim=c(0,max(CNSClassifier_Clean$CortexMagnitude)),ylim=c(0,max(CNSClassifier_Clean$SpinalMagnitude)),asp = 1)
circsize = 1
for(k in 1:nrow(CNSClassifier_Clean)){
  uppercol = SplitCols$CortexColor[k]
  lowercol = SplitCols$SpinalColor[k]
  
  upper.half.circle(CNSClassifier_Clean$CortexMagnitude[k],CNSClassifier_Clean$SpinalMagnitude[k],circsize,nsteps=1000,col=uppercol)
  lower.half.circle(CNSClassifier_Clean$CortexMagnitude[k],CNSClassifier_Clean$SpinalMagnitude[k],circsize,nsteps=1000,col=lowercol)
  
}

########## PLS-DA

Cortex_pls = Cortex_Class_RPKM
Cortex_pls = data.frame(Cortex_pls)

CortexPheno$PLSSubtype = NA
customcols = rep(NA,nrow(CortexPheno))

for(m in 1:nrow(CortexPheno)){
  if(CortexPheno$Subtype[m] == "OX"){
    CortexPheno$PLSSubtype[m] = "OX"
    customcols[m] = "navy"
  }else{
    CortexPheno$PLSSubtype[m] = "NotOX"
    customcols[m] = "gray50"
  }
}

Spinal_pls = Spinal_Class_RPKM
Spinal_pls = data.frame(Spinal_pls)

SpinalPheno$PLSSubtype = NA
customcols2 = rep(NA,nrow(SpinalPheno))

for(n in 1:nrow(SpinalPheno)){
  if(SpinalPheno$Subtype[n] == "OX"){
    SpinalPheno$PLSSubtype[n] = "OX"
    customcols2[n] = "navy"
  }else{
    SpinalPheno$PLSSubtype[n] = "NotOX"
    customcols2[n] = "gray50"
  }
}


if(HiSeqHoldout == T){
  
  C_NS_ind = which(CortexPheno$plaform == "NovaSeq")
  C_HS_ind = which(CortexPheno$platform == "HiSeq")
  
  S_NS_ind = which(SpinalPheno$Platform == "NovaSeq")
  S_HS_ind = which(SpinalPheno$Platform == "HiSeq")
  
  Cortex_NS = Cortex_pls[,C_NS_ind]
  Cortex_HS = Cortex_pls[,C_HS_ind]
  Spinal_NS = Spinal_pls[,S_NS_ind]
  Spinal_HS = Spinal_pls[,S_HS_ind]
  
  Nova = cbind(Cortex_NS,Spinal_NS)
  Hi = cbind(Cortex_HS,Spinal_HS)
  
  NovaPheno = data.frame(matrix(NA,nrow=sum(ncol(Cortex_NS),ncol(Spinal_NS)),ncol=5))
  colnames(NovaPheno) = c("sample","tissue","group","platform","PLSSubtype")
  NovaPheno$sample = c(CortexPheno$Subject[C_NS_ind],SpinalPheno$ExternalSampleId[S_NS_ind])
  NovaPheno$tissue = c(CortexPheno$tissue[C_NS_ind],SpinalPheno$Tissue[S_NS_ind])
  NovaPheno$group = c(CortexPheno$disease_group[C_NS_ind],SpinalPheno$disease_group[S_NS_ind])
  NovaPheno$platform = c(CortexPheno$platform[C_NS_ind],SpinalPheno$Platform[S_NS_ind])
  NovaPheno$PLSSubtype = c(CortexPheno$PLSSubtype[C_NS_ind],SpinalPheno$PLSSubtype[S_NS_ind])
  
  
  HiPheno = data.frame(matrix(NA,nrow=sum(ncol(Cortex_HS),ncol(Spinal_HS)),ncol=5))
  colnames(HiPheno) = c("sample","tissue","group","platform","PLSSubtype")
  HiPheno$sample = c(CortexPheno$Subject[C_HS_ind],SpinalPheno$ExternalSampleId[S_HS_ind])
  HiPheno$tissue = c(CortexPheno$tissue[C_HS_ind],SpinalPheno$Tissue[S_HS_ind])
  HiPheno$group = c(CortexPheno$disease_group[C_HS_ind],SpinalPheno$disease_group[S_HS_ind])
  HiPheno$platform = c(CortexPheno$platform[C_HS_ind],SpinalPheno$Platform[S_HS_ind])
  HiPheno$PLSSubtype = c(CortexPheno$PLSSubtype[C_HS_ind],SpinalPheno$PLSSubtype[S_HS_ind])
  
  for(i in 1:iter){
    #80/20 split
    trainind = sample(seq(1,ncol(Nova)),size=round(ncol(Nova)*train,0),replace = F)
    Nova_pls_train = Nova[,trainind]
    Nova_pls_test = Nova[,-trainind]
    
    #Build model
    NovaPLS = plsda(t(Nova_pls_train),NovaPheno$PLSSubtype[trainind])
    
    #Predict test dataset
    pred = predict(NovaPLS, newdata = t(Nova_pls_test))
    predc1 = pred$class$max.dist[,1]
    predc2 = pred$class$max.dist[,2]
    
    testPheno = NovaPheno[-trainind,]
    OXind = which(testPheno$PLSSubtype == "OX") #this works bc prediction and phenotype samples in the same order
    Nind = which(testPheno$PLSSubtype == "NotOX") #this works bc prediction and phenotype samples in the same order
    
    #Overall
    testf1_c1[i] = MLmetrics::F1_Score(y_true = as.factor(NovaPheno$PLSSubtype[-trainind]),y_pred = as.factor(predc1))
    testf1_c2[i] = MLmetrics::F1_Score(y_true = as.factor(NovaPheno$PLSSubtype[-trainind]),y_pred = as.factor(predc2))
    
    #ALS-Ox
    testf1_c1_OX[i] = MLmetrics::F1_Score(y_true = factor(testPheno$PLSSubtype[OXind],levels = c("OX","NotOX")),y_pred = factor(predc1[OXind],levels = c("OX","NotOX")))
    testf1_c2_OX[i] = MLmetrics::F1_Score(y_true = factor(testPheno$PLSSubtype[OXind],levels = c("OX","NotOX")),y_pred = factor(predc2[OXind],levels = c("OX","NotOX")))
    
    #Not Ox
    testf1_c1_NO[i] = MLmetrics::F1_Score(y_true = factor(testPheno$PLSSubtype[Nind],levels = c("NotOX","OX")),y_pred = factor(predc1[Nind],levels = c("NotOX","OX")))
    testf1_c2_NO[i] = MLmetrics::F1_Score(y_true = factor(testPheno$PLSSubtype[Nind],levels = c("NotOX","OX")),y_pred = factor(predc2[Nind],levels = c("NotOX","OX")))
    
    #Predict Test Data
    TESTres = auroc(NovaPLS,newdata = t(Nova_pls_test),outcome.test = NovaPheno$PLSSubtype[-trainind],roc.comp = 1,title="NovaSeq Test Data",plot=F)
    testauc_c1[i] = TESTres$Comp1[[1]]
    testauc_c2[i] = TESTres$Comp2[[1]]
  }
  
  #Rebuild the model one last time using all training data
  NovaPLS = plsda(t(Nova),NovaPheno$PLSSubtype)
  
  #Get performance and F1
  NovaModRes = perf(NovaPLS)
  
  #Predict "Validation" Data - ROC PLOT
  auroc(NovaPLS,newdata = t(Hi),outcome.test = HiPheno$PLSSubtype,roc.comp = 1,title="HiSeq Validation Data")
  Sys.sleep(1)
  SenSpec = auroc(NovaPLS,newdata = t(Hi),outcome.test = HiPheno$PLSSubtype,roc.comp = 1,title="HiSeq Validation Data")
  
  #View correlation of variables with response (subtype)
  network(NovaPLS)
  Sys.sleep(2)
  
  #View PLSDA Space
  plotIndiv(NovaPLS,legend = T,main = paste("PLS-DA: ",paste(Candidates,sep=","),sep=""))
  #As points...
  plotIndiv(NovaPLS,legend = T,main = paste("PLS-DA: ",paste(Candidates,sep=","),sep=""),ind.names = F)
  
  #Performance metrics in the test cohort - Component 1
  boxplot(testauc_c1,testf1_c1,col=c("#bd265d","#7b2cbf"),ylim=c(0,1),xaxt='n',main="AUC and F1 metrics: NovaSeq Test Cohort Component 1",outline=FALSE)
  axis(1,at=c(1,2),labels = c("AUC","F1"),cex.axis=1.5)
  points(rep(1,length(testauc_c1)),testauc_c1,col="#800e37")
  points(rep(2,length(testf1_c1)),testf1_c1,col="#400661")
  
  #Performance metrics in the test cohort - Component 2
  boxplot(testauc_c2,testf1_c2,col=c("#bd265d","#7b2cbf"),ylim=c(0,1),xaxt='n',main="AUC and F1 metrics: NovaSeq Test Cohort Component 2",outline=FALSE)
  axis(1,at=c(1,2),labels = c("AUC","F1"),cex.axis=1.5)
  points(rep(1,length(testauc_c2)),testauc_c2,col="#800e37")
  points(rep(2,length(testf1_c2)),testf1_c2,col="#400661")
  
  #Subtype level F1 performance metrics - Component 1
  boxplot(testf1_c1_NO,testf1_c1_OX,col=c("gray50","#4975cc"),ylim=c(0,1),ylab="F1 Score",xaxt='n',main="Subtype F1 NovaSeq Test Cohort C1",outline=FALSE)
  axis(1,at=c(1,2),labels = c("Not ALS-Ox","ALS-Ox"),cex.axis=1.5)
  points(rep(1,length(testf1_c1_NO)),testf1_c1_NO,col="gray20",pch=19)
  points(rep(2,length(testf1_c1_OX)),testf1_c1_OX,col="navy",pch=19)
  
  #Subtype level F1 performance metrics - Component 2
  boxplot(testf1_c2_NO,testf1_c2_OX,col=c("gray50","#4975cc"),ylim=c(0,1),ylab="F1 Score",xaxt='n',main="Subtype F1 NovaSeq Test Cohort C2",outline=FALSE)
  axis(1,at=c(1,2),labels = c("Not ALS-Ox","ALS-Ox"),cex.axis=1.5)
  points(rep(1,length(testf1_c2_NO)),testf1_c2_NO,col="gray20",pch=19)
  points(rep(2,length(testf1_c2_OX)),testf1_c2_OX,col="navy",pch=19)
  
}else if(CortexHoldout == T){
  #use Cortex as validation cohort
  for(i in 1:iter){
    trainind = sample(seq(1,ncol(Spinal_pls)),size=round(ncol(Spinal_pls)*train,0),replace = F)
    
    Spinal_pls_train = Spinal_pls[,trainind]
    Spinal_pls_test = Spinal_pls[,-trainind]
    
    SpinalPLS = plsda(t(Spinal_pls_train),SpinalPheno$PLSSubtype[trainind])
    
    #Predict test dataset
    pred = predict(SpinalPLS, newdata = t(Spinal_pls_test))
    predc1 = pred$class$max.dist[,1]
    predc2 = pred$class$max.dist[,2]
    
    testPheno = SpinalPheno[-trainind,]
    OXind = which(testPheno$PLSSubtype == "OX") #this works bc prediction and phenotype samples in the same order
    Nind = which(testPheno$PLSSubtype == "NotOX") #this works bc prediction and phenotype samples in the same order
    
    #Overall
    testf1_c1[i] = MLmetrics::F1_Score(y_true = as.factor(SpinalPheno$PLSSubtype[-trainind]),y_pred = as.factor(predc1))
    testf1_c2[i] = MLmetrics::F1_Score(y_true = as.factor(SpinalPheno$PLSSubtype[-trainind]),y_pred = as.factor(predc2))
    
    #ALS-Ox
    testf1_c1_OX[i] = MLmetrics::F1_Score(y_true = factor(testPheno$PLSSubtype[OXind],levels = c("OX","NotOX")),y_pred = factor(predc1[OXind],levels = c("OX","NotOX")))
    testf1_c2_OX[i] = MLmetrics::F1_Score(y_true = factor(testPheno$PLSSubtype[OXind],levels = c("OX","NotOX")),y_pred = factor(predc2[OXind],levels = c("OX","NotOX")))
    
    #Not Ox
    testf1_c1_NO[i] = MLmetrics::F1_Score(y_true = factor(testPheno$PLSSubtype[Nind],levels = c("NotOX","OX")),y_pred = factor(predc1[Nind],levels = c("NotOX","OX")))
    testf1_c2_NO[i] = MLmetrics::F1_Score(y_true = factor(testPheno$PLSSubtype[Nind],levels = c("NotOX","OX")),y_pred = factor(predc2[Nind],levels = c("NotOX","OX")))
    
    
    #Predict Test Data
    TESTres = auroc(SpinalPLS,newdata = t(Spinal_pls_test),outcome.test = SpinalPheno$PLSSubtype[-trainind],roc.comp = 1,title="Spinal Test Data",plot=F)
    testauc_c1[i] = TESTres$Comp1[[1]]
    testauc_c2[i] = TESTres$Comp2[[1]]
    
  }
  
  #Predict Test Data - ROC PLOT
  # auroc(SpinalPLS,newdata = t(Cortex_pls_test),outcome.test = CortexPheno$PLSSubtype[-trainind],roc.comp = 1,title="Cortex Test Data")
  # Sys.sleep(1)
  
  #Rebuild the model one last time using all training data
  SpinalPLS = plsda(t(Spinal_pls),SpinalPheno$PLSSubtype)
  
  #Get performance and F1
  CortexModRes = perf(SpinalPLS)
  
  #Predict "Validation" Data - ROC PLOT
  auroc(SpinalPLS,newdata = t(Cortex_pls),outcome.test = CortexPheno$PLSSubtype,roc.comp = 1,title="Cortex Validation Data")
  Sys.sleep(1)
  SenSpec = auroc(SpinalPLS,newdata = t(Cortex_pls),outcome.test = CortexPheno$PLSSubtype,roc.comp = 1,title="Cortex Validation Data")
  
  #View correlation of variables with response (subtype)
  network(SpinalPLS)
  Sys.sleep(2)
  
  #View PLSDA Space
  plotIndiv(SpinalPLS,legend = T,main = paste("PLS-DA: ",paste(Candidates,sep=","),sep=""))
  #As points...
  plotIndiv(SpinalPLS,legend = T,main = paste("PLS-DA: ",paste(Candidates,sep=","),sep=""),ind.names = F)
  
  #Performance metrics in the test cohort - Component 1
  boxplot(testauc_c1,testf1_c1,col=c("#bd265d","#7b2cbf"),ylim=c(0,1),xaxt='n',main="AUC and F1 metrics: Spinal Test Cohort Component 1",outline=FALSE)
  axis(1,at=c(1,2),labels = c("AUC","F1"),cex.axis=1.5)
  points(rep(1,length(testauc_c1)),testauc_c1,col="#800e37")
  points(rep(2,length(testf1_c1)),testf1_c1,col="#400661")
  
  #Performance metrics in the test cohort - Component 2
  boxplot(testauc_c2,testf1_c2,col=c("#bd265d","#7b2cbf"),ylim=c(0,1),xaxt='n',main="AUC and F1 metrics: Spinal Test Cohort Component 2",outline=FALSE)
  axis(1,at=c(1,2),labels = c("AUC","F1"),cex.axis=1.5)
  points(rep(1,length(testauc_c2)),testauc_c2,col="#800e37")
  points(rep(2,length(testf1_c2)),testf1_c2,col="#400661")
  
  #Subtype level F1 performance metrics - Component 1
  boxplot(testf1_c1_NO,testf1_c1_OX,col=c("gray50","#4975cc"),ylim=c(0,1),ylab="F1 Score",xaxt='n',main="Subtype F1 Spinal Test Cohort C1",outline=FALSE)
  axis(1,at=c(1,2),labels = c("Not ALS-Ox","ALS-Ox"),cex.axis=1.5)
  points(rep(1,length(testf1_c1_NO)),testf1_c1_NO,col="gray20",pch=19)
  points(rep(2,length(testf1_c1_OX)),testf1_c1_OX,col="navy",pch=19)
  
  #Subtype level F1 performance metrics - Component 2
  boxplot(testf1_c2_NO,testf1_c2_OX,col=c("gray50","#4975cc"),ylim=c(0,1),ylab="F1 Score",xaxt='n',main="Subtype F1 Spinal Test Cohort C2",outline=FALSE)
  axis(1,at=c(1,2),labels = c("Not ALS-Ox","ALS-Ox"),cex.axis=1.5)
  points(rep(1,length(testf1_c2_NO)),testf1_c2_NO,col="gray20",pch=19)
  points(rep(2,length(testf1_c2_OX)),testf1_c2_OX,col="navy",pch=19)
  
}else{
  #use spinal cord as validation cohort
  for(i in 1:iter){
    trainind = sample(seq(1,ncol(Cortex_pls)),size=round(ncol(Cortex_pls)*train,0),replace = F)
    
    Cortex_pls_train = Cortex_pls[,trainind]
    Cortex_pls_test = Cortex_pls[,-trainind]
    
    CortexPLS = plsda(t(Cortex_pls_train),CortexPheno$PLSSubtype[trainind])
    
    #Predict test dataset
    pred = predict(CortexPLS, newdata = t(Cortex_pls_test))
    predc1 = pred$class$max.dist[,1]
    predc2 = pred$class$max.dist[,2]
    
    testPheno = CortexPheno[-trainind,]
    OXind = which(testPheno$PLSSubtype == "OX") #this works bc prediction and phenotype samples in the same order
    Nind = which(testPheno$PLSSubtype == "NotOX") #this works bc prediction and phenotype samples in the same order
    
    #Overall
    testf1_c1[i] = MLmetrics::F1_Score(y_true = as.factor(CortexPheno$PLSSubtype[-trainind]),y_pred = as.factor(predc1))
    testf1_c2[i] = MLmetrics::F1_Score(y_true = as.factor(CortexPheno$PLSSubtype[-trainind]),y_pred = as.factor(predc2))
    
    #ALS-Ox
    testf1_c1_OX[i] = MLmetrics::F1_Score(y_true = factor(testPheno$PLSSubtype[OXind],levels = c("OX","NotOX")),y_pred = factor(predc1[OXind],levels = c("OX","NotOX")))
    testf1_c2_OX[i] = MLmetrics::F1_Score(y_true = factor(testPheno$PLSSubtype[OXind],levels = c("OX","NotOX")),y_pred = factor(predc2[OXind],levels = c("OX","NotOX")))
    
    #Not Ox
    testf1_c1_NO[i] = MLmetrics::F1_Score(y_true = factor(testPheno$PLSSubtype[Nind],levels = c("NotOX","OX")),y_pred = factor(predc1[Nind],levels = c("NotOX","OX")))
    testf1_c2_NO[i] = MLmetrics::F1_Score(y_true = factor(testPheno$PLSSubtype[Nind],levels = c("NotOX","OX")),y_pred = factor(predc2[Nind],levels = c("NotOX","OX")))
    
    #Predict Test Data
    TESTres = auroc(CortexPLS,newdata = t(Cortex_pls_test),outcome.test = CortexPheno$PLSSubtype[-trainind],roc.comp = 1,title="Cortex Test Data",plot=F)
    testauc_c1[i] = TESTres$Comp1[[1]]
    testauc_c2[i] = TESTres$Comp2[[1]]
    
  }
  
  #Predict Test Data - ROC PLOT
  # auroc(CortexPLS,newdata = t(Cortex_pls_test),outcome.test = CortexPheno$PLSSubtype[-trainind],roc.comp = 1,title="Cortex Test Data")
  # Sys.sleep(1)
  
  #Rebuild the model one last time using all training data
  CortexPLS = plsda(t(Cortex_pls),CortexPheno$PLSSubtype)
  
  #Get performance and F1
  CortexModRes = perf(CortexPLS)
  
  #Predict "Validation" Data - ROC PLOT
  auroc(CortexPLS,newdata = t(Spinal_pls),outcome.test = SpinalPheno$PLSSubtype,roc.comp = 1,title="Spinal Cord Validation Data")
  Sys.sleep(1)
  SenSpec = auroc(CortexPLS,newdata = t(Spinal_pls),outcome.test = SpinalPheno$PLSSubtype,roc.comp = 1,title="Spinal Cord Validation Data")
  
  #View correlation of variables with response (subtype)
  network(CortexPLS)
  Sys.sleep(2)
  
  #View PLSDA Space
  plotIndiv(CortexPLS,legend = T,main = paste("PLS-DA: ",paste(Candidates,sep=","),sep=""))
  #As points...
  plotIndiv(CortexPLS,legend = T,main = paste("PLS-DA: ",paste(Candidates,sep=","),sep=""),ind.names = F)
  
  #Performance metrics in the test cohort - Component 1
  boxplot(testauc_c1,testf1_c1,col=c("#bd265d","#7b2cbf"),ylim=c(0,1),xaxt='n',main="AUC and F1 metrics: Cortex Test Cohort Component 1",outline=FALSE)
  axis(1,at=c(1,2),labels = c("AUC","F1"),cex.axis=1.5)
  points(rep(1,length(testauc_c1)),testauc_c1,col="#800e37")
  points(rep(2,length(testf1_c1)),testf1_c1,col="#400661")
  
  #Performance metrics in the test cohort - Component 2
  boxplot(testauc_c2,testf1_c2,col=c("#bd265d","#7b2cbf"),ylim=c(0,1),xaxt='n',main="AUC and F1 metrics: Cortex Test Cohort Component 2",outline=FALSE)
  axis(1,at=c(1,2),labels = c("AUC","F1"),cex.axis=1.5)
  points(rep(1,length(testauc_c2)),testauc_c2,col="#800e37")
  points(rep(2,length(testf1_c2)),testf1_c2,col="#400661")
  
  #Subtype level F1 performance metrics - Component 1
  boxplot(testf1_c1_NO,testf1_c1_OX,col=c("gray50","#4975cc"),ylim=c(0,1),ylab="F1 Score",xaxt='n',main="Subtype F1 Cortex Test Cohort C1",outline=FALSE)
  axis(1,at=c(1,2),labels = c("Not ALS-Ox","ALS-Ox"),cex.axis=1.5)
  points(rep(1,length(testf1_c1_NO)),testf1_c1_NO,col="gray20",pch=19)
  points(rep(2,length(testf1_c1_OX)),testf1_c1_OX,col="navy",pch=19)
  
  #Subtype level F1 performance metrics - Component 2
  boxplot(testf1_c2_NO,testf1_c2_OX,col=c("gray50","#4975cc"),ylim=c(0,1),ylab="F1 Score",xaxt='n',main="Subtype F1 Cortex Test Cohort C2",outline=FALSE)
  axis(1,at=c(1,2),labels = c("Not ALS-Ox","ALS-Ox"),cex.axis=1.5)
  points(rep(1,length(testf1_c2_NO)),testf1_c2_NO,col="gray20",pch=19)
  points(rep(2,length(testf1_c2_OX)),testf1_c2_OX,col="navy",pch=19)
  
}

#Get Sensitivity and Specificity Metrics:
SenSpecDat = SenSpec$graph.Comp1$data
SenSpecDat$Specificity = (100-SenSpecDat$Specificity)

############################################################################################################################################################
################################# CYCLE COMBINATIONS TO IDENTIFY BEST GENES FOR THIS COHORT ################################################################
############################################################################################################################################################

library(mixOmics)
library(MLmetrics)

#Use gold transcript colors (stable) where possible: https://useast.ensembl.org/Help/View?id=151

#Estimate lengths from single transcript variant 
GABRA1Len = 4238 #NCBI Reference Sequence: NM_001127644.2; ENST00000393943.10
SLCLen = 3949 #UniProt: Q9P2U8; ENST00000263160.3
GAD2Len = 2462 #ENST00000259271.7
GLRALen = 8697 #ENST00000274093.8; NCBI Reference Sequence: NM_006529.4
PCSKLen = 5136 #ENST00000311106.8; NCBI Reference Sequence: NM_000439.5
HTRLen = 4835 #ENST00000378688.8
B4GALTLen = 2128 #ENST00000237019.11

SpinalPheno$library_size = as.numeric(SpinalPheno$library_size)
CortexPheno$library_size = as.numeric(CortexPheno$library_size)

CandList = c("GLRA3","PCSK1","HTR2A","GABRA1","GAD2","SLC17A6","B4GALT6")
Combinations = combn(CandList,3,simplify = F)
#keep=c(5,7,9)
#Combinations = Combinations[keep]

#Parameters
set.seed(123456)
HiSeqHoldout = T
CortexHoldout = F
iter = 100 #Number of cross-validation rounds
train = 0.8
test = 1-train

#Result containers
testauc_c1 = testauc_c2 = testf1_c1 = testf1_c2 = rep(NA,iter)
testf1_c1_cent = testf1_c2_cent = testf1_c1_maha = testf1_c2_maha = rep(NA,iter)
testauc_c1_mean = testauc_c2_mean = testf1_c1_mean = testf1_c2_mean = rep(NA,length(Combinations))
testf1_c1_cent_mean = testf1_c2_cent_mean = testf1_c1_maha_mean = testf1_c2_maha_mean = rep(NA,length(Combinations))
val_c1_auc = rep(NA,length(Combinations))

for(z in 1:length(Combinations)){
  
  #Cycle 3-gene combinations
  Candidates = Combinations[[z]]
  #Test 9-gene 
  #Candidates = FullCandList
  cat("Current Genes:",Candidates,"\n")
  
  Cortex_Class_Raw = CortexRaw[rownames(CortexRaw) %in% Candidates,] #These are RAW (and not normalized expression counts)
  Spinal_Class_Raw = SpinalRaw[rownames(SpinalRaw) %in% Candidates,] #These are RAW (and not normalized expression counts)
  
  Cortex_Class_RPKM = Cortex_Class_Raw
  Spinal_Class_RPKM = Spinal_Class_Raw #In the same order as Cortex
  
  #Adjust first gene to RPKM scale
  if(rownames(Cortex_Class_RPKM)[1] == "GLRA3"){
    Cortex_Class_RPKM[1,] = (Cortex_Class_Raw[1,]/(CortexPheno$library_size * GLRALen))*10^9
    Spinal_Class_RPKM[1,] = (Spinal_Class_Raw[1,]/(SpinalPheno$library_size * GLRALen))*10^9
  }else if(rownames(Cortex_Class_RPKM)[1] == "PCSK1"){
    Cortex_Class_RPKM[1,] = (Cortex_Class_Raw[1,]/(CortexPheno$library_size * PCSKLen))*10^9
    Spinal_Class_RPKM[1,] = (Spinal_Class_Raw[1,]/(SpinalPheno$library_size * PCSKLen))*10^9
  }else if(rownames(Cortex_Class_RPKM)[1] == "HTR2A"){
    Cortex_Class_RPKM[1,] = (Cortex_Class_Raw[1,]/(CortexPheno$library_size * HTRLen))*10^9
    Spinal_Class_RPKM[1,] = (Spinal_Class_Raw[1,]/(SpinalPheno$library_size * HTRLen))*10^9
  }else if(rownames(Cortex_Class_RPKM)[1] == "GABRA1"){
    Cortex_Class_RPKM[1,] = (Cortex_Class_Raw[1,]/(CortexPheno$library_size * GABRA1Len))*10^9
    Spinal_Class_RPKM[1,] = (Spinal_Class_Raw[1,]/(SpinalPheno$library_size * GABRA1Len))*10^9
  }else if(rownames(Cortex_Class_RPKM)[1] == "GAD2"){
    Cortex_Class_RPKM[1,] = (Cortex_Class_Raw[1,]/(CortexPheno$library_size * GAD2Len))*10^9
    Spinal_Class_RPKM[1,] = (Spinal_Class_Raw[1,]/(SpinalPheno$library_size * GAD2Len))*10^9
  }else if(rownames(Cortex_Class_RPKM)[1] == "SLC17A6"){
    Cortex_Class_RPKM[1,] = (Cortex_Class_Raw[1,]/(CortexPheno$library_size * SLCLen))*10^9
    Spinal_Class_RPKM[1,] = (Spinal_Class_Raw[1,]/(SpinalPheno$library_size * SLCLen))*10^9
  }else if(rownames(Cortex_Class_RPKM)[1] == "B4GALT6"){
    Cortex_Class_RPKM[1,] = (Cortex_Class_Raw[1,]/(CortexPheno$library_size * B4GALTLen))*10^9
    Spinal_Class_RPKM[1,] = (Spinal_Class_Raw[1,]/(SpinalPheno$library_size * B4GALTLen))*10^9
  }
  
  #Adjust second gene to RPKM scale
  if(rownames(Cortex_Class_RPKM)[2] == "GLRA3"){
    Cortex_Class_RPKM[2,] = (Cortex_Class_Raw[2,]/(CortexPheno$library_size * GLRALen))*10^9
    Spinal_Class_RPKM[2,] = (Spinal_Class_Raw[2,]/(SpinalPheno$library_size * GLRALen))*10^9
  }else if(rownames(Cortex_Class_RPKM)[2] == "PCSK1"){
    Cortex_Class_RPKM[2,] = (Cortex_Class_Raw[2,]/(CortexPheno$library_size * PCSKLen))*10^9
    Spinal_Class_RPKM[2,] = (Spinal_Class_Raw[2,]/(SpinalPheno$library_size * PCSKLen))*10^9
  }else if(rownames(Cortex_Class_RPKM)[2] == "HTR2A"){
    Cortex_Class_RPKM[2,] = (Cortex_Class_Raw[2,]/(CortexPheno$library_size * HTRLen))*10^9
    Spinal_Class_RPKM[2,] = (Spinal_Class_Raw[2,]/(SpinalPheno$library_size * HTRLen))*10^9
  }else if(rownames(Cortex_Class_RPKM)[2] == "GABRA1"){
    Cortex_Class_RPKM[2,] = (Cortex_Class_Raw[2,]/(CortexPheno$library_size * GABRA1Len))*10^9
    Spinal_Class_RPKM[2,] = (Spinal_Class_Raw[2,]/(SpinalPheno$library_size * GABRA1Len))*10^9
  }else if(rownames(Cortex_Class_RPKM)[2] == "GAD2"){
    Cortex_Class_RPKM[2,] = (Cortex_Class_Raw[2,]/(CortexPheno$library_size * GAD2Len))*10^9
    Spinal_Class_RPKM[2,] = (Spinal_Class_Raw[2,]/(SpinalPheno$library_size * GAD2Len))*10^9
  }else if(rownames(Cortex_Class_RPKM)[2] == "SLC17A6"){
    Cortex_Class_RPKM[2,] = (Cortex_Class_Raw[2,]/(CortexPheno$library_size * SLCLen))*10^9
    Spinal_Class_RPKM[2,] = (Spinal_Class_Raw[2,]/(SpinalPheno$library_size * SLCLen))*10^9
  }else if(rownames(Cortex_Class_RPKM)[2] == "B4GALT6"){
    Cortex_Class_RPKM[2,] = (Cortex_Class_Raw[2,]/(CortexPheno$library_size * B4GALTLen))*10^9
    Spinal_Class_RPKM[2,] = (Spinal_Class_Raw[2,]/(SpinalPheno$library_size * B4GALTLen))*10^9
  }
  
  #Adjust third gene to RPKM scale
  if(rownames(Cortex_Class_RPKM)[3] == "GLRA3"){
    Cortex_Class_RPKM[3,] = (Cortex_Class_Raw[3,]/(CortexPheno$library_size * GLRALen))*10^9
    Spinal_Class_RPKM[3,] = (Spinal_Class_Raw[3,]/(SpinalPheno$library_size * GLRALen))*10^9
  }else if(rownames(Cortex_Class_RPKM)[3] == "PCSK1"){
    Cortex_Class_RPKM[3,] = (Cortex_Class_Raw[3,]/(CortexPheno$library_size * PCSKLen))*10^9
    Spinal_Class_RPKM[3,] = (Spinal_Class_Raw[3,]/(SpinalPheno$library_size * PCSKLen))*10^9
  }else if(rownames(Cortex_Class_RPKM)[3] == "HTR2A"){
    Cortex_Class_RPKM[3,] = (Cortex_Class_Raw[3,]/(CortexPheno$library_size * HTRLen))*10^9
    Spinal_Class_RPKM[3,] = (Spinal_Class_Raw[3,]/(SpinalPheno$library_size * HTRLen))*10^9
  }else if(rownames(Cortex_Class_RPKM)[3] == "GABRA1"){
    Cortex_Class_RPKM[3,] = (Cortex_Class_Raw[3,]/(CortexPheno$library_size * GABRA1Len))*10^9
    Spinal_Class_RPKM[3,] = (Spinal_Class_Raw[3,]/(SpinalPheno$library_size * GABRA1Len))*10^9
  }else if(rownames(Cortex_Class_RPKM)[3] == "GAD2"){
    Cortex_Class_RPKM[3,] = (Cortex_Class_Raw[3,]/(CortexPheno$library_size * GAD2Len))*10^9
    Spinal_Class_RPKM[3,] = (Spinal_Class_Raw[3,]/(SpinalPheno$library_size * GAD2Len))*10^9
  }else if(rownames(Cortex_Class_RPKM)[3] == "SLC17A6"){
    Cortex_Class_RPKM[3,] = (Cortex_Class_Raw[3,]/(CortexPheno$library_size * SLCLen))*10^9
    Spinal_Class_RPKM[3,] = (Spinal_Class_Raw[3,]/(SpinalPheno$library_size * SLCLen))*10^9
  }else if(rownames(Cortex_Class_RPKM)[3] == "B4GALT6"){
    Cortex_Class_RPKM[3,] = (Cortex_Class_Raw[3,]/(CortexPheno$library_size * B4GALTLen))*10^9
    Spinal_Class_RPKM[3,] = (Spinal_Class_Raw[3,]/(SpinalPheno$library_size * B4GALTLen))*10^9
  }
  
  
  CortexMagRPKM = sqrt(Cortex_Class_RPKM[1,]^2 + Cortex_Class_RPKM[2,]^2 + Cortex_Class_RPKM[3,]^2)
  SpinalMagRPKM = sqrt(Spinal_Class_RPKM[1,]^2 + Spinal_Class_RPKM[2,]^2 + Spinal_Class_RPKM[3,]^2)
  
  CNSClassifier = data.frame(matrix(NA,nrow=length(CleanPatients),ncol=2))
  colnames(CNSClassifier) = c("CortexMagnitude","SpinalMagnitude")
  rownames(CNSClassifier) = CleanPatients
  
  #For patients with >1 sample per CNS region (Cortex, Spinal Cord), take the average of the magnitude
  
  for(j in 1:nrow(CNSClassifier)){
    
    #Cortex
    sampind = which(CortexPheno$Patient == rownames(CNSClassifier)[j])
    
    if(length(sampind)>0){
      cortexsamples = CortexPheno$Subject[sampind]
      
      Mags = as.numeric(CortexMagRPKM[names(CortexMagRPKM) %in% cortexsamples])
      if(length(Mags)>1){
        fillCortexMag = mean(Mags)
      }else{
        fillCortexMag = Mags
      }
      
      CNSClassifier$CortexMagnitude[j] = fillCortexMag
    }
    
    #Spinal Cord
    sampind2 = which(SpinalPheno$ExternalSubjectId == rownames(CNSClassifier)[j])
    
    if(length(sampind2)>0){
      spinalsamples = SpinalPheno$ExternalSampleId[sampind2]
      
      Mags2 = as.numeric(SpinalMagRPKM[names(SpinalMagRPKM) %in% spinalsamples])
      if(length(Mags2)>1){
        fillSpinalMag = mean(Mags2)
      }else{
        fillSpinalMag = Mags2
      }
      
      CNSClassifier$SpinalMagnitude[j] = fillSpinalMag
    }
  }
  
  CNSClassifier_Clean = CNSClassifier[-which(is.na(apply(CNSClassifier,1,mean))),]
  
  #Read in fully concordant patient data
  ConcordantPatients = read.csv("G:/SpinalCord/Publication/Manuscript/Tables/Publication/Supplemental_Dataset2.csv")
  OxPatients = ConcordantPatients$SubjectID[which(ConcordantPatients$MetaSpinalSubtype == "OX")]
  
  AllLabels = read.csv("G:/SpinalCord/Publication/Manuscript/Tables/Publication/Supplemental_Dataset1.csv")
  
  
  #Split Points - Fancy Plot
  
  SplitCols = data.frame(matrix(NA,nrow = nrow(CNSClassifier_Clean),ncol = 2))
  colnames(SplitCols) = c("CortexColor","SpinalColor")
  rownames(SplitCols) = rownames(CNSClassifier_Clean)
  
  for(j in 1:nrow(SplitCols)){
    
    ind = which(AllLabels$Patient == rownames(SplitCols)[j])
    ind2 = which(clinicaldata$ExternalSubjectId == rownames(SplitCols)[j])[1]
    
    if(length(ind)>0){
      if(AllLabels$MetaCortex[ind] == "GLIA"){
        SplitCols$CortexColor[j] = "goldenrod1"
      }else if(AllLabels$MetaCortex[ind] == "OX"){
        SplitCols$CortexColor[j] = "navy"
      }else if(AllLabels$MetaCortex[ind] == "TD"){
        SplitCols$CortexColor[j] = "firebrick"
      }else if(AllLabels$MetaCortex[ind] == "Discordant"){
        SplitCols$CortexColor[j] = "gray50"
      }
    }
    
    if(length(ind)>0){
      if(AllLabels$MetaSpinal[ind] == "GLIA"){
        SplitCols$SpinalColor[j] = "goldenrod1"
      }else if(AllLabels$MetaSpinal[ind] == "OX"){
        SplitCols$SpinalColor[j] = "navy"
      }else if(AllLabels$MetaSpinal[ind] == "TD"){
        SplitCols$SpinalColor[j] = "firebrick"
      }else if(AllLabels$MetaSpinal[ind] == "Discordant"){
        SplitCols$SpinalColor[j] = "gray50"
      }
    }
    
    if(clinicaldata$Subject.Group[ind2] == "Non-Neurological Control"){
      SplitCols$CortexColor[j] = "white"
      SplitCols$SpinalColor[j] = "white"
    }else if(clinicaldata$Subject.Group[ind2] == "Other Neurological Disorders"){
      SplitCols$CortexColor[j] = "black"
      SplitCols$SpinalColor[j] = "black"
    }
    
  }
  
  SplitCols[rownames(SplitCols) %in% OxPatients,] = c("#429ac9","#429ac9")
  
  
  plot(1, type="n",xlab="CortexMagnitude", ylab="SpinalMagnitude",main = paste("RPKM Magnitude Classifier: ",Candidates,sep=""),xlim=c(0,max(CNSClassifier_Clean$CortexMagnitude)),ylim=c(0,max(CNSClassifier_Clean$SpinalMagnitude)),asp = 1)
  for(k in 1:nrow(CNSClassifier_Clean)){
    uppercol = SplitCols$CortexColor[k]
    lowercol = SplitCols$SpinalColor[k]
    
    upper.half.circle(CNSClassifier_Clean$CortexMagnitude[k],CNSClassifier_Clean$SpinalMagnitude[k],1,nsteps=1000,col=uppercol)
    lower.half.circle(CNSClassifier_Clean$CortexMagnitude[k],CNSClassifier_Clean$SpinalMagnitude[k],1,nsteps=1000,col=lowercol)
    
  }
  
  Sys.sleep(5)
  
  ########## PLS-DA
  
  Cortex_pls = Cortex_Class_RPKM
  Cortex_pls = data.frame(Cortex_pls)
  
  CortexPheno$PLSSubtype = NA
  customcols = rep(NA,nrow(CortexPheno))
  
  for(m in 1:nrow(CortexPheno)){
    if(CortexPheno$Subtype[m] == "OX"){
      CortexPheno$PLSSubtype[m] = "OX"
      customcols[m] = "navy"
    }else{
      CortexPheno$PLSSubtype[m] = "NotOX"
      customcols[m] = "gray50"
    }
  }
  
  Spinal_pls = Spinal_Class_RPKM
  Spinal_pls = data.frame(Spinal_pls)
  
  SpinalPheno$PLSSubtype = NA
  customcols2 = rep(NA,nrow(SpinalPheno))
  
  for(n in 1:nrow(SpinalPheno)){
    if(SpinalPheno$Subtype[n] == "OX"){
      SpinalPheno$PLSSubtype[n] = "OX"
      customcols2[n] = "navy"
    }else{
      SpinalPheno$PLSSubtype[n] = "NotOX"
      customcols2[n] = "gray50"
    }
  }
  
  
  if(HiSeqHoldout == T){
    
    C_NS_ind = which(CortexPheno$plaform == "NovaSeq")
    C_HS_ind = which(CortexPheno$platform == "HiSeq")
    
    S_NS_ind = which(SpinalPheno$Platform == "NovaSeq")
    S_HS_ind = which(SpinalPheno$Platform == "HiSeq")
    
    Cortex_NS = Cortex_pls[,C_NS_ind]
    Cortex_HS = Cortex_pls[,C_HS_ind]
    Spinal_NS = Spinal_pls[,S_NS_ind]
    Spinal_HS = Spinal_pls[,S_HS_ind]
    
    Nova = cbind(Cortex_NS,Spinal_NS)
    Hi = cbind(Cortex_HS,Spinal_HS)
    
    NovaPheno = data.frame(matrix(NA,nrow=sum(ncol(Cortex_NS),ncol(Spinal_NS)),ncol=5))
    colnames(NovaPheno) = c("sample","tissue","group","platform","PLSSubtype")
    NovaPheno$sample = c(CortexPheno$Subject[C_NS_ind],SpinalPheno$ExternalSampleId[S_NS_ind])
    NovaPheno$tissue = c(CortexPheno$tissue[C_NS_ind],SpinalPheno$Tissue[S_NS_ind])
    NovaPheno$group = c(CortexPheno$disease_group[C_NS_ind],SpinalPheno$disease_group[S_NS_ind])
    NovaPheno$platform = c(CortexPheno$platform[C_NS_ind],SpinalPheno$Platform[S_NS_ind])
    NovaPheno$PLSSubtype = c(CortexPheno$PLSSubtype[C_NS_ind],SpinalPheno$PLSSubtype[S_NS_ind])
    
    
    HiPheno = data.frame(matrix(NA,nrow=sum(ncol(Cortex_HS),ncol(Spinal_HS)),ncol=5))
    colnames(HiPheno) = c("sample","tissue","group","platform","PLSSubtype")
    HiPheno$sample = c(CortexPheno$Subject[C_HS_ind],SpinalPheno$ExternalSampleId[S_HS_ind])
    HiPheno$tissue = c(CortexPheno$tissue[C_HS_ind],SpinalPheno$Tissue[S_HS_ind])
    HiPheno$group = c(CortexPheno$disease_group[C_HS_ind],SpinalPheno$disease_group[S_HS_ind])
    HiPheno$platform = c(CortexPheno$platform[C_HS_ind],SpinalPheno$Platform[S_HS_ind])
    HiPheno$PLSSubtype = c(CortexPheno$PLSSubtype[C_HS_ind],SpinalPheno$PLSSubtype[S_HS_ind])
    
    for(i in 1:iter){
      #80/20 split
      trainind = sample(seq(1,ncol(Nova)),size=round(ncol(Nova)*train,0),replace = F)
      Nova_pls_train = Nova[,trainind]
      Nova_pls_test = Nova[,-trainind]
      
      #Build model
      NovaPLS = plsda(t(Nova_pls_train),NovaPheno$PLSSubtype[trainind])
      
      #Predict test dataset
      pred = predict(NovaPLS, newdata = t(Nova_pls_test))
      predc1 = pred$class$max.dist[,1]
      predc2 = pred$class$max.dist[,2]
      
      predcen1 = pred$class$centroids.dist[,1]
      predcen2 = pred$class$centroids.dist[,2]
      predmaha1 = pred$class$mahalanobis.dist[,1]
      predmaha2 = pred$class$mahalanobis.dist[,2]
      
      #Max dist
      testf1_c1[i] = MLmetrics::F1_Score(y_true = as.factor(NovaPheno$PLSSubtype[-trainind]),y_pred = as.factor(predc1))
      testf1_c2[i] = MLmetrics::F1_Score(y_true = as.factor(NovaPheno$PLSSubtype[-trainind]),y_pred = as.factor(predc2))
      
      #Centroid and Mahalanobis distances
      testf1_c1_cent[i] = MLmetrics::F1_Score(y_true = as.factor(NovaPheno$PLSSubtype[-trainind]),y_pred = as.factor(predcen1))
      testf1_c2_cent[i] = MLmetrics::F1_Score(y_true = as.factor(NovaPheno$PLSSubtype[-trainind]),y_pred = as.factor(predcen2))
      testf1_c1_maha[i] = MLmetrics::F1_Score(y_true = as.factor(NovaPheno$PLSSubtype[-trainind]),y_pred = as.factor(predmaha1))
      testf1_c2_maha[i] = MLmetrics::F1_Score(y_true = as.factor(NovaPheno$PLSSubtype[-trainind]),y_pred = as.factor(predmaha2))
      
      #Predict Test Data
      TESTres = auroc(NovaPLS,newdata = t(Nova_pls_test),outcome.test = NovaPheno$PLSSubtype[-trainind],roc.comp = 1,title="NovaSeq Test Data",plot=F)
      testauc_c1[i] = TESTres$Comp1[[1]]
      testauc_c2[i] = TESTres$Comp2[[1]]
    }
    
    #Rebuild the model one last time using all training data
    NovaPLS = plsda(t(Nova),NovaPheno$PLSSubtype)
    
    #Get performance and F1
    NovaModRes = perf(NovaPLS)
    
    #Predict "Validation" Data - ROC PLOT
    auroc(NovaPLS,newdata = t(Hi),outcome.test = HiPheno$PLSSubtype,roc.comp = 1,title="HiSeq Validation Data")
    Sys.sleep(1)
    
    tmp = auroc(NovaPLS,newdata = t(Hi),outcome.test = HiPheno$PLSSubtype,roc.comp = 1,title="HiSeq Validation Data")
    val_c1_auc[z] = tmp$Comp1[1]
    
    #View correlation of variables with response (subtype)
    network(NovaPLS)
    Sys.sleep(2)
    
    #View PLSDA Space
    plotIndiv(NovaPLS,legend = T,main = paste("PLS-DA: ",paste(Candidates,sep=","),sep=""))
    
    #Performance metrics in the test cohort - Component 1
    boxplot(testauc_c1,testf1_c1,col=c("#bd265d","#7b2cbf"),xaxt='n',ylim=c(0,1),main="AUC and F1 metrics: NovaSeq Test Cohort Component 1",outline=FALSE)
    axis(1,at=c(1,2),labels=c("Max Distance AUC","Max Distance F1"))
    points(rep(1,length(testauc_c1)),testauc_c1,col="#800e37")
    points(rep(2,length(testf1_c1)),testf1_c1,col="#400661")
    
    #Performance metrics in the test cohort - Component 2
    boxplot(testauc_c2,testf1_c2,col=c("#bd265d","#7b2cbf"),xaxt='n',ylim=c(0,1),main="AUC and F1 metrics: NovaSeq Test Cohort Component 2",outline=FALSE)
    axis(1,at=c(1,2),labels=c("Max Distance AUC","Max Distance F1"))
    points(rep(1,length(testauc_c2)),testauc_c2,col="#800e37")
    points(rep(2,length(testf1_c2)),testf1_c2,col="#400661")
    
    #Other Distance Metrics - Component 1
    boxplot(testf1_c1_cent,testf1_c1_maha,col=c("#bd265d","#7b2cbf"),xaxt='n',ylim=c(0,1),main="AUC and F1 metrics: Cortex Test Cohort Component 1",outline=FALSE)
    axis(1,at=c(1,2),labels=c("Centroid F1","Mahalanobis F1"))
    points(rep(1,length(testf1_c1_cent)),testf1_c1_cent,col="#800e37")
    points(rep(2,length(testf1_c1_maha)),testf1_c1_maha,col="#400661")
    
    #Other Distance Metrics - Component 2
    boxplot(testf1_c2_cent,testf1_c2_maha,col=c("#bd265d","#7b2cbf"),xaxt='n',ylim=c(0,1),main="AUC and F1 metrics: Cortex Test Cohort Component 2",outline=FALSE)
    axis(1,at=c(1,2),labels=c("Centroid F1","Mahalanobis F1"))
    points(rep(1,length(testf1_c2_cent)),testf1_c2_cent,col="#800e37")
    points(rep(2,length(testf1_c2_maha)),testf1_c2_maha,col="#400661")
    
  }else if(CortexHoldout == T){
    #use Cortex as validation cohort
    for(i in 1:iter){
      trainind = sample(seq(1,ncol(Spinal_pls)),size=round(ncol(Spinal_pls)*train,0),replace = F)
      
      Spinal_pls_train = Spinal_pls[,trainind]
      Spinal_pls_test = Spinal_pls[,-trainind]
      
      SpinalPLS = plsda(t(Spinal_pls_train),SpinalPheno$PLSSubtype[trainind])
      
      #Predict test dataset
      pred = predict(SpinalPLS, newdata = t(Spinal_pls_test))
      predc1 = pred$class$max.dist[,1]
      predc2 = pred$class$max.dist[,2]
      
      predcen1 = pred$class$centroids.dist[,1]
      predcen2 = pred$class$centroids.dist[,2]
      predmaha1 = pred$class$mahalanobis.dist[,1]
      predmaha2 = pred$class$mahalanobis.dist[,2]
      
      testf1_c1[i] = MLmetrics::F1_Score(y_true = as.factor(SpinalPheno$PLSSubtype[-trainind]),y_pred = as.factor(predc1))
      testf1_c2[i] = MLmetrics::F1_Score(y_true = as.factor(SpinalPheno$PLSSubtype[-trainind]),y_pred = as.factor(predc2))
      
      #Centroid and Mahalanobis distances
      testf1_c1_cent[i] = MLmetrics::F1_Score(y_true = as.factor(SpinalPheno$PLSSubtype[-trainind]),y_pred = as.factor(predcen1))
      testf1_c2_cent[i] = MLmetrics::F1_Score(y_true = as.factor(SpinalPheno$PLSSubtype[-trainind]),y_pred = as.factor(predcen2))
      testf1_c1_maha[i] = MLmetrics::F1_Score(y_true = as.factor(SpinalPheno$PLSSubtype[-trainind]),y_pred = as.factor(predmaha1))
      testf1_c2_maha[i] = MLmetrics::F1_Score(y_true = as.factor(SpinalPheno$PLSSubtype[-trainind]),y_pred = as.factor(predmaha2))
      
      
      #Predict Test Data
      TESTres = auroc(SpinalPLS,newdata = t(Spinal_pls_test),outcome.test = SpinalPheno$PLSSubtype[-trainind],roc.comp = 1,title="Spinal Test Data",plot=F)
      testauc_c1[i] = TESTres$Comp1[[1]]
      testauc_c2[i] = TESTres$Comp2[[1]]
      
    }
    
    #Predict Test Data - ROC PLOT
    # auroc(SpinalPLS,newdata = t(Cortex_pls_test),outcome.test = CortexPheno$PLSSubtype[-trainind],roc.comp = 1,title="Cortex Test Data")
    # Sys.sleep(1)
    
    #Rebuild the model one last time using all training data
    SpinalPLS = plsda(t(Spinal_pls),SpinalPheno$PLSSubtype)
    
    #Get performance and F1
    CortexModRes = perf(SpinalPLS)
    
    #Predict "Validation" Data - ROC PLOT
    auroc(SpinalPLS,newdata = t(Cortex_pls),outcome.test = CortexPheno$PLSSubtype,roc.comp = 1,title="Cortex Validation Data")
    Sys.sleep(1)
    
    tmp = auroc(SpinalPLS,newdata = t(Cortex_pls),outcome.test = CortexPheno$PLSSubtype,roc.comp = 1,title="Cortex Validation Data")
    val_c1_auc[z] = tmp$Comp1[1]
    
    #View correlation of variables with response (subtype)
    network(SpinalPLS)
    Sys.sleep(2)
    
    #View PLSDA Space
    plotIndiv(SpinalPLS,legend = T,main = paste("PLS-DA: ",paste(Candidates,sep=","),sep=""))
    #As points...
    plotIndiv(SpinalPLS,legend = T,main = paste("PLS-DA: ",paste(Candidates,sep=","),sep=""),ind.names = F)
    
    
    #Performance metrics in the test cohort - Component 1
    boxplot(testauc_c1,testf1_c1,col=c("#bd265d","#7b2cbf"),ylim=c(0,1),xaxt='n',main="AUC and F1 metrics: Cortex Test Cohort Component 1",outline=FALSE)
    axis(1,at=c(1,2),labels=c("Max Distance AUC","Max Distance F1"))
    points(rep(1,length(testauc_c1)),testauc_c1,col="#800e37")
    points(rep(2,length(testf1_c1)),testf1_c1,col="#400661")
    
    #Performance metrics in the test cohort - Component 2
    boxplot(testauc_c2,testf1_c2,col=c("#bd265d","#7b2cbf"),ylim=c(0,1),xaxt='n',main="AUC and F1 metrics: Cortex Test Cohort Component 2",outline=FALSE)
    axis(1,at=c(1,2),labels=c("Max Distance AUC","Max Distance F1"))
    points(rep(1,length(testauc_c2)),testauc_c2,col="#800e37")
    points(rep(2,length(testf1_c2)),testf1_c2,col="#400661")
    
    #Other Distance Metrics - Component 1
    boxplot(testf1_c1_cent,testf1_c1_maha,col=c("#bd265d","#7b2cbf"),xaxt='n',ylim=c(0,1),main="AUC and F1 metrics: Cortex Test Cohort Component 1",outline=FALSE)
    axis(1,at=c(1,2),labels=c("Centroid F1","Mahalanobis F1"))
    points(rep(1,length(testf1_c1_cent)),testf1_c1_cent,col="#800e37")
    points(rep(2,length(testf1_c1_maha)),testf1_c1_maha,col="#400661")
    
    #Other Distance Metrics - Component 2
    boxplot(testf1_c2_cent,testf1_c2_maha,col=c("#bd265d","#7b2cbf"),xaxt='n',ylim=c(0,1),main="AUC and F1 metrics: Cortex Test Cohort Component 2",outline=FALSE)
    axis(1,at=c(1,2),labels=c("Centroid F1","Mahalanobis F1"))
    points(rep(1,length(testf1_c2_cent)),testf1_c2_cent,col="#800e37")
    points(rep(2,length(testf1_c2_maha)),testf1_c2_maha,col="#400661")
  }else{
    #use spinal cord as validation cohort
    for(i in 1:iter){
      trainind = sample(seq(1,ncol(Cortex_pls)),size=round(ncol(Cortex_pls)*train,0),replace = F)
      
      Cortex_pls_train = Cortex_pls[,trainind]
      Cortex_pls_test = Cortex_pls[,-trainind]
      
      CortexPLS = plsda(t(Cortex_pls_train),CortexPheno$PLSSubtype[trainind])
      
      #Predict test dataset
      pred = predict(CortexPLS, newdata = t(Cortex_pls_test))
      predc1 = pred$class$max.dist[,1]
      predc2 = pred$class$max.dist[,2]
      
      predcen1 = pred$class$centroids.dist[,1]
      predcen2 = pred$class$centroids.dist[,2]
      predmaha1 = pred$class$mahalanobis.dist[,1]
      predmaha2 = pred$class$mahalanobis.dist[,2]
      
      #Max distance
      testf1_c1[i] = MLmetrics::F1_Score(y_true = as.factor(CortexPheno$PLSSubtype[-trainind]),y_pred = as.factor(predc1))
      testf1_c2[i] = MLmetrics::F1_Score(y_true = as.factor(CortexPheno$PLSSubtype[-trainind]),y_pred = as.factor(predc2))
      
      #Centroid and Mahalanobis distances
      testf1_c1_cent[i] = MLmetrics::F1_Score(y_true = as.factor(CortexPheno$PLSSubtype[-trainind]),y_pred = as.factor(predcen1))
      testf1_c2_cent[i] = MLmetrics::F1_Score(y_true = as.factor(CortexPheno$PLSSubtype[-trainind]),y_pred = as.factor(predcen2))
      testf1_c1_maha[i] = MLmetrics::F1_Score(y_true = as.factor(CortexPheno$PLSSubtype[-trainind]),y_pred = as.factor(predmaha1))
      testf1_c2_maha[i] = MLmetrics::F1_Score(y_true = as.factor(CortexPheno$PLSSubtype[-trainind]),y_pred = as.factor(predmaha2))
      
      
      #Predict Test Data
      TESTres = auroc(CortexPLS,newdata = t(Cortex_pls_test),outcome.test = CortexPheno$PLSSubtype[-trainind],roc.comp = 1,title="Cortex Test Data",plot=F)
      testauc_c1[i] = TESTres$Comp1[[1]]
      testauc_c2[i] = TESTres$Comp2[[1]]
      
    }
    
    #Predict Test Data - ROC PLOT
    # auroc(CortexPLS,newdata = t(Cortex_pls_test),outcome.test = CortexPheno$PLSSubtype[-trainind],roc.comp = 1,title="Cortex Test Data")
    # Sys.sleep(1)
    
    #Rebuild the model one last time using all training data
    CortexPLS = plsda(t(Cortex_pls),CortexPheno$PLSSubtype)
    
    #Get performance and F1
    CortexModRes = perf(CortexPLS)
    
    #Predict "Validation" Data - ROC PLOT
    auroc(CortexPLS,newdata = t(Spinal_pls),outcome.test = SpinalPheno$PLSSubtype,roc.comp = 1,title="Spinal Cord Validation Data")
    Sys.sleep(1)
    
    tmp = auroc(CortexPLS,newdata = t(Spinal_pls),outcome.test = SpinalPheno$PLSSubtype,roc.comp = 1,title="Spinal Cord Validation Data")
    val_c1_auc[z] = tmp$Comp1[1]
    
    #View correlation of variables with response (subtype)
    network(CortexPLS)
    Sys.sleep(2)
    
    #View PLSDA Space
    plotIndiv(CortexPLS,legend = T,main = paste("PLS-DA: ",paste(Candidates,sep=","),sep=""))
    
    #Performance metrics in the test cohort - Component 1
    boxplot(testauc_c1,testf1_c1,col=c("#bd265d","#7b2cbf"),ylim=c(0,1),xaxt='n',main="AUC and F1 metrics: Cortex Test Cohort Component 1",outline=FALSE)
    axis(1,at=c(1,2),labels=c("Max Distance AUC","Max Distance F1"))
    points(rep(1,length(testauc_c1)),testauc_c1,col="#800e37")
    points(rep(2,length(testf1_c1)),testf1_c1,col="#400661")
    
    #Performance metrics in the test cohort - Component 2
    boxplot(testauc_c2,testf1_c2,col=c("#bd265d","#7b2cbf"),ylim=c(0,1),xaxt='n',main="AUC and F1 metrics: Cortex Test Cohort Component 2",outline=FALSE)
    axis(1,at=c(1,2),labels=c("Max Distance AUC","Max Distance F1"))
    points(rep(1,length(testauc_c2)),testauc_c2,col="#800e37")
    points(rep(2,length(testf1_c2)),testf1_c2,col="#400661")
    
    #Other Distance Metrics - Component 1
    boxplot(testf1_c1_cent,testf1_c1_maha,col=c("#bd265d","#7b2cbf"),xaxt='n',ylim=c(0,1),main="AUC and F1 metrics: Cortex Test Cohort Component 1",outline=FALSE)
    axis(1,at=c(1,2),labels=c("Centroid F1","Mahalanobis F1"))
    points(rep(1,length(testf1_c1_cent)),testf1_c1_cent,col="#800e37")
    points(rep(2,length(testf1_c1_maha)),testf1_c1_maha,col="#400661")
    
    #Other Distance Metrics - Component 2
    boxplot(testf1_c2_cent,testf1_c2_maha,col=c("#bd265d","#7b2cbf"),xaxt='n',ylim=c(0,1),main="AUC and F1 metrics: Cortex Test Cohort Component 2",outline=FALSE)
    axis(1,at=c(1,2),labels=c("Centroid F1","Mahalanobis F1"))
    points(rep(1,length(testf1_c2_cent)),testf1_c2_cent,col="#800e37")
    points(rep(2,length(testf1_c2_maha)),testf1_c2_maha,col="#400661")
  }
  
  testauc_c1_mean[z] = mean(testauc_c1,na.rm = T)
  testauc_c2_mean[z] = mean(testauc_c2,na.rm = T)
  testf1_c1_mean[z] = mean(testf1_c1,na.rm = T)
  testf1_c2_mean[z] = mean(testf1_c2,na.rm = T)
  
  testf1_c1_cent_mean[z] = mean(testf1_c1_cent,na.rm = T)
  testf1_c2_cent_mean[z] = mean(testf1_c2_cent,na.rm = T)
  testf1_c1_maha_mean[z] = mean(testf1_c1_maha,na.rm = T)
  testf1_c2_maha_mean[z] = mean(testf1_c2_maha,na.rm = T)
  
}

#Check out Combinations using component 1 and F1 scores on test cohort

tmp = testf1_c1_mean
names(tmp) = Combinations
tmp[order(tmp)]

names(val_c1_auc) = Combinations
val_c1_auc[order(val_c1_auc)]

####################################################################    END     #####################################################################################################################################################

######### Supervised Learning expression matrices

#Use gold transcript colors: https://useast.ensembl.org/Help/View?id=151

#Estimate lengths from single transcript variant 
# GABRA1Len = 4238 #NCBI Reference Sequence: NM_001127644.2; ENST00000393943.10
# SLCLen = 3949 #UniProt: Q9P2U8; ENST00000263160.3
# GAD2Len = 2462 #ENST00000259271.7
# GLRALen = 8697 #ENST00000274093.8; NCBI Reference Sequence: NM_006529.4
# PCSKLen = 5136 #ENST00000311106.8; NCBI Reference Sequence: NM_000439.5
# HTRLen = 4835 #ENST00000378688.8
# MYLLen = 2786 #ENST00000279022.7; NCBI Reference Sequence: NM_006097.5
# GALNACLen = 2170 # ENST00000225276.9
# TAGLen = 1477 #ENST00000278968.10
# B4GALTLen = 2128 #ENST00000237019.11
# 
# #Alternative number of features for machine learning
# 
# Candidates7 = c("SLC17A6","GAD2","GABRA1","GLRA3","PCSK1","HTR2A","B4GALT6")
# Cortex_Class_Raw = CortexRaw[rownames(CortexRaw) %in% Candidates7,] #These are RAW (and not normalized expression counts)
# Spinal_Class_Raw = SpinalRaw[rownames(SpinalRaw) %in% Candidates7,] #These are RAW (and not normalized expression counts)
# 
# Cortex_Class_RPKM = Cortex_Class_Raw
# Spinal_Class_RPKM = Spinal_Class_Raw #In the same order as Cortex
# 
# Cortex_Class_RPKM[1,] = (Cortex_Class_Raw[1,]/(CortexPheno$library_size * GABRA1Len))*10^9
# Cortex_Class_RPKM[2,] = (Cortex_Class_Raw[2,]/(CortexPheno$library_size * SLCLen))*10^9
# Cortex_Class_RPKM[3,] = (Cortex_Class_Raw[3,]/(CortexPheno$library_size * HTRLen))*10^9
# Cortex_Class_RPKM[4,] = (Cortex_Class_Raw[4,]/(CortexPheno$library_size * B4GALTLen))*10^9
# Cortex_Class_RPKM[5,] = (Cortex_Class_Raw[5,]/(CortexPheno$library_size * GAD2Len))*10^9
# Cortex_Class_RPKM[6,] = (Cortex_Class_Raw[6,]/(CortexPheno$library_size * GLRALen))*10^9
# Cortex_Class_RPKM[7,] = (Cortex_Class_Raw[7,]/(CortexPheno$library_size * PCSKLen))*10^9
# 
# Spinal_Class_RPKM[1,] = (Spinal_Class_Raw[1,]/(SpinalPheno$library_size * GABRA1Len))*10^9
# Spinal_Class_RPKM[2,] = (Spinal_Class_Raw[2,]/(SpinalPheno$library_size * SLCLen))*10^9
# Spinal_Class_RPKM[3,] = (Spinal_Class_Raw[3,]/(SpinalPheno$library_size * HTRLen))*10^9
# Spinal_Class_RPKM[4,] = (Spinal_Class_Raw[4,]/(SpinalPheno$library_size * B4GALTLen))*10^9
# Spinal_Class_RPKM[5,] = (Spinal_Class_Raw[5,]/(SpinalPheno$library_size * GAD2Len))*10^9
# Spinal_Class_RPKM[6,] = (Spinal_Class_Raw[6,]/(SpinalPheno$library_size * GLRALen))*10^9
# Spinal_Class_RPKM[7,] = (Spinal_Class_Raw[7,]/(SpinalPheno$library_size * PCSKLen))*10^9
# 
# write.csv(Cortex_Class_RPKM,"Cortex_RPKM_OX_7GeneClassifier.csv")
# write.csv(Spinal_Class_RPKM,"Spinal_RPKM_OX_7GeneClassifier.csv")
# 
# 
# Candidates10 = c("SLC17A6","GAD2","GABRA1","GLRA3","PCSK1","HTR2A","MYL9","ST6GALNAC2","TAGLN","B4GALT6")
# Cortex_Class_Raw = CortexRaw[rownames(CortexRaw) %in% Candidates10,] #These are RAW (and not normalized expression counts)
# Spinal_Class_Raw = SpinalRaw[rownames(SpinalRaw) %in% Candidates10,] #These are RAW (and not normalized expression counts)
# 
# Cortex_Class_RPKM = Cortex_Class_Raw
# Spinal_Class_RPKM = Spinal_Class_Raw #In the same order as Cortex
# 
# table(CortexPheno$Subject == colnames(Cortex_Class_Raw))
# 
# Cortex_Class_RPKM = Cortex_Class_Raw
# Spinal_Class_RPKM = Spinal_Class_Raw
# 
# Cortex_Class_RPKM[1,] = (Cortex_Class_Raw[1,]/(CortexPheno$library_size * GABRA1Len))*10^9
# Cortex_Class_RPKM[2,] = (Cortex_Class_Raw[2,]/(CortexPheno$library_size * GALNACLen))*10^9
# Cortex_Class_RPKM[3,] = (Cortex_Class_Raw[3,]/(CortexPheno$library_size * SLCLen))*10^9
# Cortex_Class_RPKM[4,] = (Cortex_Class_Raw[4,]/(CortexPheno$library_size * MYLLen))*10^9
# Cortex_Class_RPKM[5,] = (Cortex_Class_Raw[5,]/(CortexPheno$library_size * HTRLen))*10^9
# Cortex_Class_RPKM[6,] = (Cortex_Class_Raw[6,]/(CortexPheno$library_size * B4GALTLen))*10^9
# Cortex_Class_RPKM[7,] = (Cortex_Class_Raw[7,]/(CortexPheno$library_size * GAD2Len))*10^9
# Cortex_Class_RPKM[8,] = (Cortex_Class_Raw[8,]/(CortexPheno$library_size * GLRALen))*10^9
# Cortex_Class_RPKM[9,] = (Cortex_Class_Raw[9,]/(CortexPheno$library_size * TAGLen))*10^9
# Cortex_Class_RPKM[10,] = (Cortex_Class_Raw[10,]/(CortexPheno$library_size * PCSKLen))*10^9
# 
# SpinalPheno$library_size = as.numeric(SpinalPheno$library_size)
# Spinal_Class_RPKM[1,] = (Spinal_Class_Raw[1,]/(SpinalPheno$library_size * GABRA1Len))*10^9
# Spinal_Class_RPKM[2,] = (Spinal_Class_Raw[2,]/(SpinalPheno$library_size * GALNACLen))*10^9
# Spinal_Class_RPKM[3,] = (Spinal_Class_Raw[3,]/(SpinalPheno$library_size * SLCLen))*10^9
# Spinal_Class_RPKM[4,] = (Spinal_Class_Raw[4,]/(SpinalPheno$library_size * MYLLen))*10^9
# Spinal_Class_RPKM[5,] = (Spinal_Class_Raw[5,]/(SpinalPheno$library_size * HTRLen))*10^9
# Spinal_Class_RPKM[6,] = (Spinal_Class_Raw[6,]/(SpinalPheno$library_size * B4GALTLen))*10^9
# Spinal_Class_RPKM[7,] = (Spinal_Class_Raw[7,]/(SpinalPheno$library_size * GAD2Len))*10^9
# Spinal_Class_RPKM[8,] = (Spinal_Class_Raw[8,]/(SpinalPheno$library_size * GLRALen))*10^9
# Spinal_Class_RPKM[9,] = (Spinal_Class_Raw[9,]/(SpinalPheno$library_size * TAGLen))*10^9
# Spinal_Class_RPKM[10,] = (Spinal_Class_Raw[10,]/(SpinalPheno$library_size * PCSKLen))*10^9
# 
# write.csv(Cortex_Class_RPKM,"Cortex_RPKM_10GeneClassifier.csv")
# write.csv(Spinal_Class_RPKM,"Spinal_RPKM_10GeneClassifier.csv")

################################### OTHER COUNT SCALES - NOT INCLUDED IN PAPER ##################################################
################## MoR COUNTS

# Candidates = c("SLC17A6","GAD2","GABRA1")
# #Candidates = c("SLC17A6","GAD2","GABRA1","GLRA3","PCSK1","HTR2A")
# 
# Cortex_Class_Raw = CortexRaw[rownames(CortexRaw) %in% Candidates,] #These are RAW (and not normalized expression counts)
# Cortex_Class = CortexExpr[rownames(CortexExpr) %in% Candidates,] #These are DESeq2 median-of-ratio counts (genes fit to MLE using negative binomial)
# 
# Spinal_Class_Raw = SpinalRaw[rownames(SpinalRaw) %in% Candidates,] #These are RAW (and not normalized expression counts)
# Spinal_Class = SpinalExpr[rownames(SpinalExpr) %in% Candidates,] #These are DESeq2 median-of-ratio counts (genes fit to MLE using negative binomial)
# 
# #visualize 
# 
# subcolors = rep(NA,ncol(Cortex_Class))
# for(i in 1:nrow(CortexPheno)){
#   if(CortexPheno$Subtype[i] == "GLIA"){
#     subcolors[i] = "goldenrod1"
#   }else if(CortexPheno$Subtype[i] == "OX"){
#     subcolors[i] = "navy"
#   }else if(CortexPheno$Subtype[i] == "TE"){
#     subcolors[i] = "firebrick"
#   }else{
#     subcolors[i] = "gray50"
#   }
# }
# 
# plot3d(Cortex_Class[1,],Cortex_Class[2,],Cortex_Class[3,],col = subcolors,size=5)
# 
# CortexMag = sqrt(Cortex_Class[1,]^2 + Cortex_Class[2,]^2 + Cortex_Class[3,]^2)
# SpinalMag = sqrt(Spinal_Class[1,]^2 + Spinal_Class[2,]^2 + Spinal_Class[3,]^2)
# 
# 
# #Coldata matrix for DESeq2
# clinicaldata = read.csv("G:/SpinalCord/Publication/CLINICAL_DATA_PRUDENCIO.csv") 
# clinicaldata$ExternalSampleId = gsub("-",".",clinicaldata$ExternalSampleId)
# ExtraCortex = clinicaldata[clinicaldata$ExternalSampleId %in% CortexPheno$Subject,]
# 
# CortexPheno$Patient = NA
# 
# for(i in 1:nrow(CortexPheno)){
#   ind = which(ExtraCortex$ExternalSampleId == CortexPheno$Subject[i])
#   CortexPheno$Patient[i] = ExtraCortex$ExternalSubjectId[ind]
# }
# 
# CleanPatients = names(table(c(CortexPheno$Patient,SpinalPheno$ExternalSubjectId)))
# 
# #### Match Magnitude vectors to patients
# 
# CNSClassifier = data.frame(matrix(NA,nrow=length(CleanPatients),ncol=2))
# colnames(CNSClassifier) = c("CortexMagnitude","SpinalMagnitude")
# rownames(CNSClassifier) = CleanPatients
# 
# #For patients with >1 sample per CNS region (Cortex, Spinal Cord), take the average of the magnitude
# 
# for(i in 1:nrow(CNSClassifier)){
#   
#   #Cortex
#   sampind = which(CortexPheno$Patient == rownames(CNSClassifier)[i])
#   
#   if(length(sampind)>0){
#     cortexsamples = CortexPheno$Subject[sampind]
#     
#     Mags = as.numeric(CortexMag[names(CortexMag) %in% cortexsamples])
#     if(length(Mags)>1){
#       fillCortexMag = mean(Mags)
#     }else{
#       fillCortexMag = Mags
#     }
#     
#     CNSClassifier$CortexMagnitude[i] = fillCortexMag
#   }
#   
#   #Spinal Cord
#   sampind2 = which(SpinalPheno$ExternalSubjectId == rownames(CNSClassifier)[i])
#   
#   if(length(sampind2)>0){
#     spinalsamples = SpinalPheno$ExternalSampleId[sampind2]
#     
#     Mags2 = as.numeric(SpinalMag[names(SpinalMag) %in% spinalsamples])
#     if(length(Mags2)>1){
#       fillSpinalMag = mean(Mags2)
#     }else{
#       fillSpinalMag = Mags2
#     }
#     
#     CNSClassifier$SpinalMagnitude[i] = fillSpinalMag
#   }
# }
# 
# CNSClassifier_Clean = CNSClassifier[-which(is.na(apply(CNSClassifier,1,mean))),]
# 
# #Read in fully concordant patient data
# ConcordantPatients = read.csv("G:/SpinalCord/Publication/Manuscript/Tables/Publication/Supplemental_Dataset2.csv")
# OxPatients = ConcordantPatients$SubjectID[which(ConcordantPatients$MetaSpinalSubtype == "OX")]
# 
# AllLabels = read.csv("G:/SpinalCord/Publication/Manuscript/Tables/Publication/Supplemental_Dataset1.csv")
# 
# 
# #Split Points - Fancy Plot
# 
# SplitCols = data.frame(matrix(NA,nrow = nrow(CNSClassifier_Clean),ncol = 2))
# colnames(SplitCols) = c("CortexColor","SpinalColor")
# rownames(SplitCols) = rownames(CNSClassifier_Clean)
# 
# for(i in 1:nrow(SplitCols)){
#   
#   ind = which(AllLabels$Patient == rownames(SplitCols)[i])
#   ind2 = which(clinicaldata$ExternalSubjectId == rownames(SplitCols)[i])[1]
#   
#   if(length(ind)>0){
#     if(AllLabels$MetaCortex[ind] == "GLIA"){
#       SplitCols$CortexColor[i] = "goldenrod1"
#     }else if(AllLabels$MetaCortex[ind] == "OX"){
#       SplitCols$CortexColor[i] = "navy"
#     }else if(AllLabels$MetaCortex[ind] == "TD"){
#       SplitCols$CortexColor[i] = "firebrick"
#     }else if(AllLabels$MetaCortex[ind] == "Discordant"){
#       SplitCols$CortexColor[i] = "gray25"
#     }else if(is.na(AllLabels$MetaCortex[ind])){
#       SplitCols$CortexColor[i] = "gray75"
#     }
#   }
#   
#   if(length(ind)>0){
#     if(AllLabels$MetaSpinal[ind] == "GLIA"){
#       SplitCols$SpinalColor[i] = "goldenrod1"
#     }else if(AllLabels$MetaSpinal[ind] == "OX"){
#       SplitCols$SpinalColor[i] = "navy"
#     }else if(AllLabels$MetaSpinal[ind] == "TD"){
#       SplitCols$SpinalColor[i] = "firebrick"
#     }else if(AllLabels$MetaSpinal[ind] == "Discordant"){
#       SplitCols$SpinalColor[i] = "gray25"
#     }else if(is.na(AllLabels$MetaSpinal[ind])){
#       SplitCols$SpinalColor[i] = "gray75"
#     }
#   }
#   
#   if(clinicaldata$Subject.Group[ind2] == "Non-Neurological Control" || clinicaldata$Subject.Group[ind2] == "Other Neurological Disorders"){
#     SplitCols$CortexColor[i] = "white"
#     SplitCols$SpinalColor[i] = "white"
#   }
#   
# }
# 
# 
# plot(1, type="n",xlab="CortexMagnitude", ylab="SpinalMagnitude",main = "MoR Magnitude Classifier - GABRA1, GAD2, SLC17A6",xlim=c(0,max(CNSClassifier_Clean$CortexMagnitude)),ylim=c(0,max(CNSClassifier_Clean$SpinalMagnitude)), asp = 1)
# for(i in 1:nrow(CNSClassifier_Clean)){
#   uppercol = SplitCols$CortexColor[i]
#   lowercol = SplitCols$SpinalColor[i]
#   
#   upper.half.circle(CNSClassifier_Clean$CortexMagnitude[i],CNSClassifier_Clean$SpinalMagnitude[i],100,nsteps=1000,col=uppercol)
#   lower.half.circle(CNSClassifier_Clean$CortexMagnitude[i],CNSClassifier_Clean$SpinalMagnitude[i],100,nsteps=1000,col=lowercol)
#   
# }
# 
# 
# 
# #################### RAW COUNTS (highly susceptible to instrument/user/site variability)
# 
# CortexMagRaw = sqrt(Cortex_Class_Raw[1,]^2 + Cortex_Class_Raw[2,]^2 + Cortex_Class_Raw[3,]^2)
# SpinalMagRaw = sqrt(Spinal_Class_Raw[1,]^2 + Spinal_Class_Raw[2,]^2 + Spinal_Class_Raw[3,]^2)
# 
# CNSClassifier = data.frame(matrix(NA,nrow=length(CleanPatients),ncol=2))
# colnames(CNSClassifier) = c("CortexMagnitude","SpinalMagnitude")
# rownames(CNSClassifier) = CleanPatients
# 
# #For patients with >1 sample per CNS region (Cortex, Spinal Cord), take the average of the magnitude
# 
# for(i in 1:nrow(CNSClassifier)){
#   
#   #Cortex
#   sampind = which(CortexPheno$Patient == rownames(CNSClassifier)[i])
#   
#   if(length(sampind)>0){
#     cortexsamples = CortexPheno$Subject[sampind]
#     
#     Mags = as.numeric(CortexMagRaw[names(CortexMagRaw) %in% cortexsamples])
#     if(length(Mags)>1){
#       fillCortexMag = mean(Mags)
#     }else{
#       fillCortexMag = Mags
#     }
#     
#     CNSClassifier$CortexMagnitude[i] = fillCortexMag
#   }
#   
#   #Spinal Cord
#   sampind2 = which(SpinalPheno$ExternalSubjectId == rownames(CNSClassifier)[i])
#   
#   if(length(sampind2)>0){
#     spinalsamples = SpinalPheno$ExternalSampleId[sampind2]
#     
#     Mags2 = as.numeric(SpinalMagRaw[names(SpinalMagRaw) %in% spinalsamples])
#     if(length(Mags2)>1){
#       fillSpinalMag = mean(Mags2)
#     }else{
#       fillSpinalMag = Mags2
#     }
#     
#     CNSClassifier$SpinalMagnitude[i] = fillSpinalMag
#   }
# }
# 
# CNSClassifier_Clean = CNSClassifier[-which(is.na(apply(CNSClassifier,1,mean))),]
# 
# #Read in fully concordant patient data
# ConcordantPatients = read.csv("G:/SpinalCord/Publication/Manuscript/Tables/Publication/Supplemental_Dataset2.csv")
# OxPatients = ConcordantPatients$SubjectID[which(ConcordantPatients$MetaSpinalSubtype == "OX")]
# 
# AllLabels = read.csv("G:/SpinalCord/Publication/Manuscript/Tables/Publication/Supplemental_Dataset1.csv")
# 
# 
# #Split Points - Fancy Plot
# 
# SplitCols = data.frame(matrix(NA,nrow = nrow(CNSClassifier_Clean),ncol = 2))
# colnames(SplitCols) = c("CortexColor","SpinalColor")
# rownames(SplitCols) = rownames(CNSClassifier_Clean)
# 
# for(i in 1:nrow(SplitCols)){
#   
#   ind = which(AllLabels$Patient == rownames(SplitCols)[i])
#   ind2 = which(clinicaldata$ExternalSubjectId == rownames(SplitCols)[i])[1]
#   
#   if(length(ind)>0){
#     if(AllLabels$MetaCortex[ind] == "GLIA"){
#       SplitCols$CortexColor[i] = "goldenrod1"
#     }else if(AllLabels$MetaCortex[ind] == "OX"){
#       SplitCols$CortexColor[i] = "navy"
#     }else if(AllLabels$MetaCortex[ind] == "TD"){
#       SplitCols$CortexColor[i] = "firebrick"
#     }else if(AllLabels$MetaCortex[ind] == "Discordant"){
#       SplitCols$CortexColor[i] = "gray25"
#     }else if(is.na(AllLabels$MetaCortex[ind])){
#       SplitCols$CortexColor[i] = "gray75"
#     }
#   }
#   
#   if(length(ind)>0){
#     if(AllLabels$MetaSpinal[ind] == "GLIA"){
#       SplitCols$SpinalColor[i] = "goldenrod1"
#     }else if(AllLabels$MetaSpinal[ind] == "OX"){
#       SplitCols$SpinalColor[i] = "navy"
#     }else if(AllLabels$MetaSpinal[ind] == "TD"){
#       SplitCols$SpinalColor[i] = "firebrick"
#     }else if(AllLabels$MetaSpinal[ind] == "Discordant"){
#       SplitCols$SpinalColor[i] = "gray25"
#     }else if(is.na(AllLabels$MetaSpinal[ind])){
#       SplitCols$SpinalColor[i] = "gray75"
#     }
#   }
#   
#   if(clinicaldata$Subject.Group[ind2] == "Non-Neurological Control" || clinicaldata$Subject.Group[ind2] == "Other Neurological Disorders"){
#     SplitCols$CortexColor[i] = "white"
#     SplitCols$SpinalColor[i] = "white"
#   }
#   
# }
# 
# 
# plot(1, type="n",xlab="CortexMagnitude", ylab="SpinalMagnitude",main = "Raw Magnitude Classifier - GABRA1, GAD2, SLC17A6",xlim=c(0,max(CNSClassifier_Clean$CortexMagnitude)),ylim=c(0,max(CNSClassifier_Clean$SpinalMagnitude)), asp = 1)
# for(i in 1:nrow(CNSClassifier_Clean)){
#   uppercol = SplitCols$CortexColor[i]
#   lowercol = SplitCols$SpinalColor[i]
#   
#   upper.half.circle(CNSClassifier_Clean$CortexMagnitude[i],CNSClassifier_Clean$SpinalMagnitude[i],100,nsteps=1000,col=uppercol)
#   lower.half.circle(CNSClassifier_Clean$CortexMagnitude[i],CNSClassifier_Clean$SpinalMagnitude[i],100,nsteps=1000,col=lowercol)
#   
# }


############################################# Top 10 ################################################################################
#10-gene PLSDA
#Exploratory - worse than 7 gene

# Candidates10 = c("SLC17A6","GAD2","GABRA1","GLRA3","PCSK1","HTR2A","MYL9","ST6GALNAC2","TAGLN","B4GALT6")
# Cortex_Class_Raw = CortexRaw[rownames(CortexRaw) %in% Candidates10,] #These are RAW (and not normalized expression counts)
# Spinal_Class_Raw = SpinalRaw[rownames(SpinalRaw) %in% Candidates10,] #These are RAW (and not normalized expression counts)
# 
# Cortex_Class_RPKM = Cortex_Class_Raw
# Spinal_Class_RPKM = Spinal_Class_Raw #In the same order as Cortex
# 
# Cortex_Class_RPKM[1,] = (Cortex_Class_Raw[1,]/(CortexPheno$library_size * GABRA1Len))*10^9
# Cortex_Class_RPKM[2,] = (Cortex_Class_Raw[2,]/(CortexPheno$library_size * GALNACLen))*10^9
# Cortex_Class_RPKM[3,] = (Cortex_Class_Raw[3,]/(CortexPheno$library_size * SLCLen))*10^9
# Cortex_Class_RPKM[4,] = (Cortex_Class_Raw[4,]/(CortexPheno$library_size * MYLLen))*10^9
# Cortex_Class_RPKM[5,] = (Cortex_Class_Raw[5,]/(CortexPheno$library_size * HTRLen))*10^9
# Cortex_Class_RPKM[6,] = (Cortex_Class_Raw[6,]/(CortexPheno$library_size * B4GALTLen))*10^9
# Cortex_Class_RPKM[7,] = (Cortex_Class_Raw[7,]/(CortexPheno$library_size * GAD2Len))*10^9
# Cortex_Class_RPKM[8,] = (Cortex_Class_Raw[8,]/(CortexPheno$library_size * GLRALen))*10^9
# Cortex_Class_RPKM[9,] = (Cortex_Class_Raw[9,]/(CortexPheno$library_size * TAGLen))*10^9
# Cortex_Class_RPKM[10,] = (Cortex_Class_Raw[10,]/(CortexPheno$library_size * PCSKLen))*10^9
# 
# SpinalPheno$library_size = as.numeric(SpinalPheno$library_size)
# Spinal_Class_RPKM[1,] = (Spinal_Class_Raw[1,]/(SpinalPheno$library_size * GABRA1Len))*10^9
# Spinal_Class_RPKM[2,] = (Spinal_Class_Raw[2,]/(SpinalPheno$library_size * GALNACLen))*10^9
# Spinal_Class_RPKM[3,] = (Spinal_Class_Raw[3,]/(SpinalPheno$library_size * SLCLen))*10^9
# Spinal_Class_RPKM[4,] = (Spinal_Class_Raw[4,]/(SpinalPheno$library_size * MYLLen))*10^9
# Spinal_Class_RPKM[5,] = (Spinal_Class_Raw[5,]/(SpinalPheno$library_size * HTRLen))*10^9
# Spinal_Class_RPKM[6,] = (Spinal_Class_Raw[6,]/(SpinalPheno$library_size * B4GALTLen))*10^9
# Spinal_Class_RPKM[7,] = (Spinal_Class_Raw[7,]/(SpinalPheno$library_size * GAD2Len))*10^9
# Spinal_Class_RPKM[8,] = (Spinal_Class_Raw[8,]/(SpinalPheno$library_size * GLRALen))*10^9
# Spinal_Class_RPKM[9,] = (Spinal_Class_Raw[9,]/(SpinalPheno$library_size * TAGLen))*10^9
# Spinal_Class_RPKM[10,] = (Spinal_Class_Raw[10,]/(SpinalPheno$library_size * PCSKLen))*10^9

######################################################################################################