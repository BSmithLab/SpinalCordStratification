
#setwd("G:/SpinalCord/Publication/MetaData")

# NovaSeqMeta = read.csv("SpinalCord_NovaSeq_MetaData_ALS.csv")
# HiSeqMeta = read.csv("SpinalCord_HiSeq_MetaData_ALS.csv")
# 
# setwd("G:/SpinalCord/Publication/RawExpression/NovaSeq_NoGlia_JH")
# NovaExp = read.csv("NOVASEQ_SpinalCord_ALSCohort_ALS_MAD10k_NoGliaMarkers1k_GeneSYMBOL_TE_8-17-23.csv")
# setwd("G:/SpinalCord/Publication/RawExpression/HiSeq_NoGlia_JH")
# HiSeqExp = read.csv("HISEQ_SpinalCord_ALSCohort_NoGliaMarkers1k_MAD10k_GeneSYMBOL_TE_8-17-23.csv")
# 
# 

#Subtypes = read.csv("ALS_Cortex_SpinalCordNoOligo_SubtypeConcordance_8-14-23.csv")

setwd("G:/SpinalCord/Publication/SystemicAnalysis")
NovaMeta = read.csv("NovaSeq_PatientPheno_NoOligo_SpinalCord_SuppData_8-14-23.csv")
HiMeta = read.csv("HiSeq_PatientPheno_NoOligo_SpinalCord_SuppData_8-14-23.csv")

setwd("G:/SpinalCord/Publication")
Pheno = read.csv("CLINICAL_DATA_PRUDENCIO.csv")
Pheno$ExternalSampleId = gsub("-",".",Pheno$ExternalSampleId)

########################################## NovaSeq #############################################################

gind = NovaMeta$CervicalSample[which(NovaMeta$CervicalSubtype == "GLIA")]
oind = NovaMeta$CervicalSample[which(NovaMeta$CervicalSubtype == "OX")]
tind = NovaMeta$CervicalSample[which(NovaMeta$CervicalSubtype == "TD")]

hist(Pheno$RIN[Pheno$ExternalSampleId %in% gind],main="Glia Subtype RIN - Cervical Spinal Cord",col = "goldenrod1",xlim=c(0,10))
hist(Pheno$RIN[Pheno$ExternalSampleId %in% oind],main="OX Subtype RIN - Cervical Spinal Cord",col="navy",xlim=c(0,10))
hist(Pheno$RIN[Pheno$ExternalSampleId %in% tind],main="TD Subtype RIN - Cervical Spinal Cord",col="firebrick",xlim=c(0,10))



gind = NovaMeta$ThoracicSample[which(NovaMeta$ThoracicSubtype == "GLIA")]
oind = NovaMeta$ThoracicSample[which(NovaMeta$ThoracicSubtype == "OX")]
tind = NovaMeta$ThoracicSample[which(NovaMeta$ThoracicSubtype == "TD")]

gx = Pheno$RIN[Pheno$ExternalSampleId %in% gind]
ox = Pheno$RIN[Pheno$ExternalSampleId %in% oind]
tx = Pheno$RIN[Pheno$ExternalSampleId %in% tind]

hist(Pheno$RIN[Pheno$ExternalSampleId %in% gind],main="Glia Subtype RIN - Thoracic Spinal Cord",col = "goldenrod1",xlim=c(0,10))
hist(Pheno$RIN[Pheno$ExternalSampleId %in% oind],main="OX Subtype RIN - Thoracic Spinal Cord",col="navy",xlim=c(0,10))
hist(Pheno$RIN[Pheno$ExternalSampleId %in% tind],main="TD Subtype RIN - Thoracic Spinal Cord",col="firebrick",xlim=c(0,10))



gind = NovaMeta$LumbarSample[which(NovaMeta$LumbarSubtype == "GLIA")]
oind = NovaMeta$LumbarSample[which(NovaMeta$LumbarSubtype == "OX")]
tind = NovaMeta$LumbarSample[which(NovaMeta$LumbarSubtype == "TD")]

hist(Pheno$RIN[Pheno$ExternalSampleId %in% gind],main="Glia Subtype RIN - Lumbar Spinal Cord",col = "goldenrod1",xlim=c(0,10))
hist(Pheno$RIN[Pheno$ExternalSampleId %in% oind],main="OX Subtype RIN - Lumbar Spinal Cord",col="navy",xlim=c(0,10))
hist(Pheno$RIN[Pheno$ExternalSampleId %in% tind],main="TD Subtype RIN - Lumbar Spinal Cord",col="firebrick",xlim=c(0,10))



########################################## HiSeq #############################################################


gind = HiMeta$CervicalSample[which(HiMeta$CervicalSubtype == "GLIA")]
oind = HiMeta$CervicalSample[which(HiMeta$CervicalSubtype == "OX")]
tind = HiMeta$CervicalSample[which(HiMeta$CervicalSubtype == "TD")]

hist(Pheno$RIN[Pheno$ExternalSampleId %in% gind],main="Glia Subtype RIN - Cervical Spinal Cord - HiSeq",col = "goldenrod1",xlim=c(0,10))
hist(Pheno$RIN[Pheno$ExternalSampleId %in% oind],main="OX Subtype RIN - Cervical Spinal Cord - HiSeq",col="navy",xlim=c(0,10))
hist(Pheno$RIN[Pheno$ExternalSampleId %in% tind],main="TD Subtype RIN - Cervical Spinal Cord - HiSeq",col="firebrick",xlim=c(0,10))



gind = HiMeta$ThoracicSample[which(HiMeta$ThoracicSubtype == "GLIA")]
oind = HiMeta$ThoracicSample[which(HiMeta$ThoracicSubtype == "OX")]
tind = HiMeta$ThoracicSample[which(HiMeta$ThoracicSubtype == "TD")]

gx = Pheno$RIN[Pheno$ExternalSampleId %in% gind]
ox = Pheno$RIN[Pheno$ExternalSampleId %in% oind]
tx = Pheno$RIN[Pheno$ExternalSampleId %in% tind]

hist(Pheno$RIN[Pheno$ExternalSampleId %in% gind],main="Glia Subtype RIN - Thoracic Spinal Cord - HiSeq",col = "goldenrod1",xlim=c(0,10))
hist(Pheno$RIN[Pheno$ExternalSampleId %in% oind],main="OX Subtype RIN - Thoracic Spinal Cord - HiSeq",col="navy",xlim=c(0,10))
hist(Pheno$RIN[Pheno$ExternalSampleId %in% tind],main="TD Subtype RIN - Thoracic Spinal Cord - HiSeq",col="firebrick",xlim=c(0,10))



gind = HiMeta$LumbarSample[which(HiMeta$LumbarSubtype == "GLIA")]
oind = HiMeta$LumbarSample[which(HiMeta$LumbarSubtype == "OX")]
tind = HiMeta$LumbarSample[which(HiMeta$LumbarSubtype == "TD")]

hist(Pheno$RIN[Pheno$ExternalSampleId %in% gind],main="Glia Subtype RIN - Lumbar Spinal Cord - HiSeq",col = "goldenrod1",xlim=c(0,10))
hist(Pheno$RIN[Pheno$ExternalSampleId %in% oind],main="OX Subtype RIN - Lumbar Spinal Cord - HiSeq",col="navy",xlim=c(0,10))
hist(Pheno$RIN[Pheno$ExternalSampleId %in% tind],main="TD Subtype RIN - Lumbar Spinal Cord - HiSeq",col="firebrick",xlim=c(0,10))
