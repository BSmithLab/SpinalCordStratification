#############################  SQuIRE Transposable Element RNA-seq Post-Processing Part 2  ####################################################

#Written by: Jarrett Eshima
#Updated by: Raiyan Choudhury
#Date: June 4th, 2021
#For: Use by the Dr. Barbara Smith Lab at Arizona State University

## Description: This code is intended to help with the organization of SQuIRE post-processing count .RData files
# This code will provide a single output matrix of TPM counts as rows and RNA-seq samples as columns.
# Only the TEs found in all samples/subjects will be considered for downstream analysis.

#Important Note: This code is not fully automated. User is required to determine the appropriate number of "chunks" to load (dependent on system mem)

###############################################################################################################################################
############ USER PARAMETERS

# #Set working directory where .RData files are stored
# wd = "E:/NovaSeqTE/NovaSeq/ALS/Chunks"
# setwd(wd) 
# 
# #File Handle - the lettering scheme that precedes the integer and .RData file type
# handle = "SQuIRE_Post_Chunk"
# 
# #Number of "chunks" / .RData files
# nchunk = 10
# 
# #Chunk Length (number of subjects in each chunk - make sure to adjust for the last chunk)
# clength = 50
# 
# #If your last chunk is not the same length as other chunks - update arguments below
# allequal = F # change this to false if all chunks do not have the same file number
# totalfiles = 501
# 
# #Select directory for output file
# outputdir = "G:/MetaTE"

###############################################################################################################################################
### Read in .RDATA files (from Post-Processing Script) and identify overlapping features

#Set working directory where .RData files are stored
wd = "G:/SpinalCord/NovaSeqTE/NovaSeq/Controls/Chunks"
setwd(wd) 

#File Handle - the lettering scheme that precedes the integer and .RData file type
handle = "SQuIRE_Post_Chunk"

#Number of "chunks" / .RData files
nchunk = 2

#Select directory for output file
outputdir = "G:/SpinalCord/MetaTE"

#Name output Count Matrix (no file type handle needed)
savefile = "NovaSeq_SpinalCord_ControlCohort_TEExpression_8-15-23"

#Piecemeal Approach - .RData batches of 3
myrdata = paste(handle,seq(1,nchunk,1),".RData",sep = "")

#Load in all .RData files
for(y in 1:nchunk){
  load(myrdata[y])
  if((y %% 1) == 0) cat(".RData Chunk Read-In:",y,"\n")
}

#I have 9 chunks, so I am breaking it up into 1-3 and 4-6 and 7-9
#Comment the section corresponding to the chunks not being run
#tmp1 = get(paste("avalanche",clength*1,sep=""))
#tmp2 = get(paste("avalanche",clength*2,sep=""))
#tmp3 = get(paste("avalanche",clength*3,sep=""))
#TE_KEEP = tmp1[tmp1$TE_ID %in% tmp2$TE_ID,]
#TE_KEEP = TE_KEEP[TE_KEEP$TE_ID %in% tmp3$TE_ID,]
#write.csv(TE_KEEP,"TE_KEEP_Chunks1to3.csv")

######### NovaSeq ALS

# tmp1 = get("finalavalanche1")
# tmp2 = get("finalavalanche2")
# tmp3 = get("finalavalanche3")
# tmp4 = get("finalavalanche4")
# tmp5 = get("finalavalanche5")
# tmp6 = get("finalavalanche6")
# tmp7 = get("finalavalanche7")
# tmp8 = get("finalavalanche8")
# 
# 
# TE_KEEP_NS = tmp1[tmp1$TE_ID %in% tmp2$TE_ID,]
# TE_KEEP_NS = TE_KEEP_NS[TE_KEEP_NS$TE_ID %in% tmp3$TE_ID,]
# TE_KEEP_NS = TE_KEEP_NS[TE_KEEP_NS$TE_ID %in% tmp4$TE_ID,]
# TE_KEEP_NS = TE_KEEP_NS[TE_KEEP_NS$TE_ID %in% tmp5$TE_ID,]
# TE_KEEP_NS = TE_KEEP_NS[TE_KEEP_NS$TE_ID %in% tmp6$TE_ID,]
# TE_KEEP_NS = TE_KEEP_NS[TE_KEEP_NS$TE_ID %in% tmp7$TE_ID,]
# TE_KEEP_NS = TE_KEEP_NS[TE_KEEP_NS$TE_ID %in% tmp8$TE_ID,]

#Completed: NovaSeq
# setwd(outputdir)
# write.csv(TE_KEEP_NS,"NovaSeq_Final_TE_List.csv")

########### HiSeq ALS

# tmp1 = get("finalavalanche1")
# tmp2 = get("finalavalanche2")
# tmp3 = get("finalavalanche3")
# 
# 
# TE_KEEP_HS = tmp1[tmp1$TE_ID %in% tmp2$TE_ID,]
# TE_KEEP_HS = TE_KEEP_HS[TE_KEEP_HS$TE_ID %in% tmp3$TE_ID,]


#Completed: HiSeq
# setwd(outputdir)
# write.csv(TE_KEEP_HS,"HiSeq_Final_TE_List.csv")





#Read in file lists
# setwd(outputdir)
# TE_KEEP_NS = read.csv("NovaSeq_Final_TE_List.csv")
# TE_KEEP_HS = read.csv("HiSeq_Final_TE_List.csv")
# 
# FinalTE_List = TE_KEEP_NS[TE_KEEP_NS$TE_ID %in% TE_KEEP_HS$TE_ID,]
# colnames(FinalTE_List)[1] = "index"

#write.csv(FinalTE_List,"SpinalCord_FullCohort_Final_TE_List_6-2-23.csv")

# #Read in the Hi-Seq cohort and filter one last time - Doesnt work bc controls remove some must rerun 6-2-23
# HiSeq_Exp = read.csv("E:/Publication/RawExpression/HiSeq_SpinalCord_FullCohort_TECounts_HGND_4-8-23.csv")
# rownames(HiSeq_Exp) = HiSeq_Exp[,1]
# HiSeq_Exp = HiSeq_Exp[,-1]
# Clean_HiSeq_TE_Exp = HiSeq_Exp[rownames(HiSeq_Exp) %in% FinalTE_List$TE_ID,]
# write.csv(Clean_HiSeq_TE_Exp,"SpinalCord_HiSeqCohort_UpdatedTEs_6-2-23.csv")






###############################################################################################################################################

#Read in file lists
#outputdir = "E:/MetaTE"
setwd(outputdir)
TE_KEEP_NS = read.csv("NovaSeq_Final_TE_List.csv")
TE_KEEP_HS = read.csv("HiSeq_Final_TE_List.csv")

FinalTE_List = TE_KEEP_NS[TE_KEEP_NS$TE_ID %in% TE_KEEP_HS$TE_ID,]
IDs = FinalTE_List$TE_ID

#Completed: NovaSeq ALS, HiSeq ALS, HiSeq Control,

##############################################################################################################################################
#Generate a count matrix for the transposable elements
wd = "G:/SpinalCord/NovaSeqTE/NovaSeq/Controls"
setwd(wd)

FileList = list.files(path = wd, pattern = "_TEcounts.txt")
library(stringr)
accessions = str_sub(FileList,1,11)
TEFiles = names(table(accessions))
TEFinalList = TEFiles[!is.na(TEFiles)]


TECounts = matrix(NA,nrow=length(IDs),ncol=length(TEFinalList))
rownames(TECounts) = IDs
colnames(TECounts) = TEFinalList


for(i in 1:length(TEFinalList)){
  
  index = TEFinalList[i]
  SubjectRead = read.delim(FileList[i],fill = T,header = T,skipNul = T)
  FilteredSubjectRead = SubjectRead[, c("TE_ID", "score", "tot_counts")]
  SubjectName = paste("Subject",i, sep = "")
  
  assign(SubjectName,FilteredSubjectRead)
  tmp = get(SubjectName)
  #nam = substr(allfiles[index],1,11)
  #assign(nam,read.table(allfiles[index],header = T))
  
  #tmp = get(paste("Subject",i,sep=""))
  tmp = tmp[tmp$TE_ID %in% IDs,]
  
  
  for(j in 1:nrow(TECounts)){
    for(k in 1:nrow(tmp)){
      
      if(rownames(TECounts)[j] == tmp$TE_ID[k]){
        TECounts[j,i] = tmp$tot_counts[k]
      }
      
    }
  }
  
  if((i %% 1) == 0) cat("TE Counts Completed for Subject:",i,"of",length(TEFinalList),"\n")
}
setwd(outputdir)
save.image(paste(savefile,".RData"))
#load("NovaSeq_GSE153960_SpinalCord_TECounts_6-3-23.RData")
##############################################################################################################################################
#TE Count matrix combined with gene count matrix in the next script (ALSPatientStratification_DifferentialExpression_ALSPatients.R)

###############################################################################################################################################
###############NovaSeq

##Clean up accession names to alternative ID (CGND-HRA)

#Convert SRR to HGND
setwd("G:/SpinalCord/MetaTE")
lut = read.table("SpinalCord_MetaData_GSE153960_Full.csv",sep = ",",fill = T,header = T) #This file can be obtained from NCBI GEO repository (Accession: PRJNA644618)
cleanlut = lut[lut$Run %in% colnames(TECounts),]

#For novaseq cohort
dups = names(which(table(cleanlut$sample_id_alt)>1))
duplut = cleanlut[cleanlut$sample_id_alt %in% dups,]

nondups = names(which(table(cleanlut$sample_id_alt)==1))
nonduplut = cleanlut[cleanlut$sample_id_alt %in% nondups,]

keep = rep(NA,length(dups))
for(i in 1:length(dups)){
  currentdup = dups[i]
  sizes = duplut$Bases[which(duplut$sample_id_alt == currentdup)]
  tmpmax = max(sizes)
  
  for(j in 1:nrow(duplut)){
    
    
    
    if(duplut$sample_id_alt[j] == currentdup && duplut$Bases[j] == tmpmax){
      keep[i] = j
    }
    
  }
}

keep = keep[!is.na(keep)]
duplut2 = duplut[keep,]

table(colnames(duplut2) == colnames(nonduplut))

finallut = rbind(nonduplut,duplut2)

TECounts2 = TECounts[,colnames(TECounts) %in% finallut$Run]

CGND_IDs = finallut$sample_id_alt
SRR_IDs = finallut$Run

Convert = matrix(cbind(SRR_IDs,CGND_IDs),ncol=2)
TECountsF = TECounts2

for(i in 1:nrow(Convert)){
  for(j in 1:ncol(TECountsF)){
    if(Convert[i,1] == colnames(TECountsF)[j]){
      colnames(TECountsF)[j] = Convert[i,2]
    }
  }
}

setwd(outputdir)
#write.csv(TECountsF,"NovaSeq_SpinalCord_ControlCohort_TECounts_HGND_6-3-23_nodups.csv")
#write.csv(cleanlut,"NovaSeq_SpinalCord_ControlCohort_FinalMetaData_GSE153960_6-3-23.csv")


######### HiSeq & Controls

#Convert SRR to HGND
setwd("G:/SpinalCord/MetaTE")
lut = read.table("SpinalCord_MetaData_GSE153960_Full.csv",sep = ",",fill = T,header = T) #This file can be obtained from NCBI GEO repository (Accession: PRJNA644618)
cleanlut = lut[lut$Run %in% colnames(TECounts),]

CGND_IDs = cleanlut$sample_id_alt
SRR_IDs = cleanlut$Run

Convert = matrix(cbind(SRR_IDs,CGND_IDs),ncol=2)
TECountsF = TECounts

for(i in 1:nrow(Convert)){
  for(j in 1:ncol(TECounts)){
    if(Convert[i,1] == colnames(TECounts)[j]){
      colnames(TECountsF)[j] = Convert[i,2]
    }
  }
}

setwd(outputdir)
#write.csv(TECountsF,"HiSeq_SpinalCord_ControlCohort_TECounts_HGND_8-15-23.csv")
#write.csv(cleanlut,"HiSeq_SpinalCord_ControlCohort_FinalMetaData_GSE153960_8-15-23.csv")
write.csv(TECountsF,"NovaSeq_SpinalCord_ControlCohort_TECounts_HGND_8-15-23.csv")
write.csv(cleanlut,"NovaSeq_SpinalCord_ControlCohort_FinalMetaData_GSE153960_8-15-23.csv")
#load("GSE153960_TECounts_6-13-21.RData")
