##Use NMF 'feature scores' to determine the subset of MAD genes to use for supervised classification and GSEA

#This script generates the final feature set used in enrichment, network construction, and univariate analysis

#Written By: Jarrett Eshima
#For: Dr. Barbara Smith Lab

#Note 9/27/23: Covariate-dependent genes were not removed in the previous script. This ensures the top1000 features from both sequencing platform cohorts are considered. This fixes the issue of few overlapping genes between the sequencing platforms but leaves some dependency of the results on external covariates.

## Read in "totfeatures" files. These are obtained from SAKE by downloading all feature scores following NMF clustering.

NumReplicates = 11

#Read in the data
for(i in 1:NumReplicates){
  setwd(paste("G:/SpinalCord/Publication/Clustering/NovaSeq/4Covar_5k/RobustSubtypeAssignment/",seq(1,NumReplicates,1),sep="")[i])
  tmp = read.csv('tot_features.csv') #Generated from SAKE NMF output
  rownames(tmp) = make.unique(tmp$GeneCard)
  tmp = tmp[,-1]
  repname = paste("NovaRep",i,sep="")
  assign(repname,tmp)
}

NumReplicates = 11

for(i in 1:NumReplicates){
  setwd(paste("G:/SpinalCord/Publication/Clustering/HiSeq/4Covar_5k/RobustSubtypeAssignment/",seq(1,NumReplicates,1),sep="")[i])
  tmp = read.csv('tot_features.csv') #Generated from SAKE NMF output
  rownames(tmp) = make.unique(tmp$GeneCard)
  tmp = tmp[,-1]
  repname = paste("HiRep",i,sep="")
  assign(repname,tmp)
}


#subset out data based on averaged feature score across 11 runs
#both top x features for each subtype and top x features overall

#NovaSeq
AverageFeatureScoreNova = matrix(NA,nrow=nrow(NovaRep1),ncol = 2)
rownames(AverageFeatureScoreNova) = rownames(NovaRep1)
colnames(AverageFeatureScoreNova) = c("AverageScore","Rep1Group")
AverageFeatureScoreNova = data.frame(AverageFeatureScoreNova)

for(i in 1:nrow(AverageFeatureScoreNova)){
  r1 = which(rownames(NovaRep1) == rownames(AverageFeatureScoreNova)[i])
  r2 = which(rownames(NovaRep2) == rownames(AverageFeatureScoreNova)[i])
  r3 = which(rownames(NovaRep3) == rownames(AverageFeatureScoreNova)[i])
  r4 = which(rownames(NovaRep4) == rownames(AverageFeatureScoreNova)[i])
  r5 = which(rownames(NovaRep5) == rownames(AverageFeatureScoreNova)[i])
  r6 = which(rownames(NovaRep6) == rownames(AverageFeatureScoreNova)[i])
  r7 = which(rownames(NovaRep7) == rownames(AverageFeatureScoreNova)[i])
  r8 = which(rownames(NovaRep8) == rownames(AverageFeatureScoreNova)[i])
  r9 = which(rownames(NovaRep9) == rownames(AverageFeatureScoreNova)[i])
  r10 = which(rownames(NovaRep10) == rownames(AverageFeatureScoreNova)[i])
  r11 = which(rownames(NovaRep11) == rownames(AverageFeatureScoreNova)[i])
  
  tmp1 = mean(c(NovaRep1$featureScore[r1],NovaRep2$featureScore[r2],NovaRep3$featureScore[r3],NovaRep4$featureScore[r4],NovaRep5$featureScore[r5],NovaRep6$featureScore[r6],NovaRep7$featureScore[r7],NovaRep8$featureScore[r8],NovaRep9$featureScore[r9],NovaRep10$featureScore[r10]),NovaRep11$featureScore[r11])
  AverageFeatureScoreNova[i,1] = tmp1 #the average feature score across NMF replicates
  AverageFeatureScoreNova[i,2] = NovaRep1$Group[r1] #the NMF cluster of the feature from the first run only
  
  if((i %% 500) == 0) cat("% Done:",i/nrow(AverageFeatureScoreNova)*100,"\n")
}

orderNova = AverageFeatureScoreNova[order(AverageFeatureScoreNova$AverageScore,decreasing = T),]

head(orderNova)

#Does removing glia features have a strong impact on enrichment?
#No - less than 5% difference between the two feature sets. 
#test = read.csv("G:/SpinalCord/Publication/Clustering/NovaSeq/4Covar_wGlia_5k/RobustSubtypeAssignment/1/tot_features.csv")
# test = test[order(test$featureScore,decreasing = T),]
# head(test)
# 
# test1k = cbind(rownames(orderNova)[1:1000],test$GeneCard[1:1000])
# colnames(test1k) = c("Top1k_without","Top1k_WITH")
# tmp = data.frame(test1k)
# tmp$Top1k_without[tmp$Top1k_without %in% tmp$Top1k_WITH]
# #Percent of features that are the same between the two runs
# length(tmp$Top1k_without[tmp$Top1k_without %in% tmp$Top1k_WITH])/1000
# 
# test2k = cbind(rownames(orderNova)[1:2500],test$GeneCard[1:2500])
# colnames(test2k) = c("Top2k_without","Top2k_WITH")
# tmp = data.frame(test2k)
# tmp$Top2k_without[tmp$Top2k_without %in% tmp$Top2k_WITH]
# #Percent of features that are the same between the two runs
# length(tmp$Top2k_without[tmp$Top2k_without %in% tmp$Top2k_WITH])/2500
#
# test5k = cbind(rownames(orderNova)[1:5000],test$GeneCard[1:5000])
# colnames(test5k) = c("Top5k_without","Top5k_WITH")
# tmp = data.frame(test5k)
# tmp$Top5k_without[tmp$Top5k_without %in% tmp$Top5k_WITH]
# #Percent of features that are the same between the two runs
# length(tmp$Top5k_without[tmp$Top5k_without %in% tmp$Top5k_WITH])/5000

#HiSeq
AverageFeatureScoreHi = matrix(NA,nrow=nrow(HiRep1),ncol = 2)
rownames(AverageFeatureScoreHi) = rownames(HiRep1)
colnames(AverageFeatureScoreHi) = c("AverageScore","Rep1Group")
AverageFeatureScoreHi = data.frame(AverageFeatureScoreHi)

for(i in 1:nrow(AverageFeatureScoreHi)){
  r1 = which(rownames(HiRep1) == rownames(AverageFeatureScoreHi)[i])
  r2 = which(rownames(HiRep2) == rownames(AverageFeatureScoreHi)[i])
  r3 = which(rownames(HiRep3) == rownames(AverageFeatureScoreHi)[i])
  r4 = which(rownames(HiRep4) == rownames(AverageFeatureScoreHi)[i])
  r5 = which(rownames(HiRep5) == rownames(AverageFeatureScoreHi)[i])
  r6 = which(rownames(HiRep6) == rownames(AverageFeatureScoreHi)[i])
  r7 = which(rownames(HiRep7) == rownames(AverageFeatureScoreHi)[i])
  r8 = which(rownames(HiRep8) == rownames(AverageFeatureScoreHi)[i])
  r9 = which(rownames(HiRep9) == rownames(AverageFeatureScoreHi)[i])
  r10 = which(rownames(HiRep10) == rownames(AverageFeatureScoreHi)[i])
  r11 = which(rownames(HiRep11) == rownames(AverageFeatureScoreHi)[i])
  
  tmp1 = mean(c(HiRep1$featureScore[r1],HiRep2$featureScore[r2],HiRep3$featureScore[r3],HiRep4$featureScore[r4],HiRep5$featureScore[r5],HiRep6$featureScore[r6],HiRep7$featureScore[r7],HiRep8$featureScore[r8],HiRep9$featureScore[r9],HiRep10$featureScore[r10],HiRep11$featureScore[r11]))
  AverageFeatureScoreHi[i,1] = tmp1
  AverageFeatureScoreHi[i,2] = HiRep1$Group[r1]
  
  if((i %% 500) == 0) cat("% Done:",i/nrow(AverageFeatureScoreHi)*100,"\n")
}

orderHi = AverageFeatureScoreHi[order(AverageFeatureScoreHi$AverageScore,decreasing = T),]

head(orderHi)



#Top features overall - NovaSeq
Ntop1000ALL = rownames(orderNova)[1:1000] #Top 20% of MAD features
Ntop2500ALL = rownames(orderNova)[1:2500] #Top 50% of MAD features
NtopALL = rownames(orderNova)[1:5000] #ALL FEATURES

#Top features overall - HiSeq
Htop1000ALL = rownames(orderHi)[1:1000] #Top 20% of MAD features
Htop2500ALL = rownames(orderHi)[1:2500] #Top 50% of MAD features
HtopALL = rownames(orderHi)[1:5000] #ALL FEATURES

### Combine NovaSeq and HiSeq feature sets

#Combined All
Ctop1000ALL = names(table(c(Ntop1000ALL,Htop1000ALL)))
Ctop2500ALL = names(table(c(Ntop2500ALL,Htop2500ALL)))
Ctop5000ALL = names(table(c(NtopALL,HtopALL)))

#Write feature selection results
# wd = "G:/SpinalCord/Publication/FeatureSelection"
# setwd(wd)
# write.csv(orderNova,"NovaSeq_SpinalCord_nsNMF_AvgFeatureScores_9-27-23.csv")
# write.csv(orderHi,"HiSeq_SpinalCord_nsNMF_AvgFeatureScores_9-27-23.csv")
# write.table(Ctop1000ALL,"Enrichment_DE_SpinalCord_Top1k_nsNMFFeatureSet_9-27-23.txt",quote = F,row.names = F)
# write.table(Ctop2500ALL,"Enrichment_DE_SpinalCord_Top2k_nsNMFFeatureSet_9-27-23.txt",quote = F,row.names = F)

#####################################################################################################################################################################################
#ENSMBL LUT for VST Files (both VST files required because MAD10k and other processed files contain platform-specific genes)

#Get ENSEMBL ID and Symbol LUT for novaseq and hiseq platforms (and combine)
setwd("G:/SpinalCord/Publication/RawExpression/WithControls/NovaSeq_4Covar")

mad10000.ens = read.csv("NOVASEQ_SpinalCord_FULLCohort_VST_Gene_NoFiltering_TE_9-27-23.csv",header = T) #File generated in previous script
rownames(mad10000.ens) = mad10000.ens[,1]
mad10000.ens = mad10000.ens[,-1]
egenes = rownames(mad10000.ens)
newEG = sub("\\..*","",egenes) #Keep text to the left of the dot
ensembl_version = "https://apr2020.archive.ensembl.org" #Updated: https://jan2020.archive.ensembl.org
species="human"
ensembl = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host=ensembl_version)
gene_positions_nova <- biomaRt::getBM(filters = 'ensembl_gene_id',attributes=c('ensembl_gene_id','hgnc_symbol'), values = newEG, mart = ensembl)

# Commented: 9/27/23
# setwd("G:/SpinalCord/Publication/RawExpression/WithControls/HiSeq_4Covar")
# 
# mad10000.ens = read.csv("HISEQ_SpinalCord_FULLCohort_VST_Gene_NoFiltering_TE_9-27-23.csv",header=T) #File generated in previous script
# rownames(mad10000.ens) = mad10000.ens[,1] 
# mad10000.ens = mad10000.ens[,-1]
# egenes = rownames(mad10000.ens)
# newEG = sub("\\..*","",egenes) #Keep text to the left of the dot
# ensembl_version = "https://apr2020.archive.ensembl.org" #Updated: https://jan2020.archive.ensembl.org
# species="human"
# ensembl = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host=ensembl_version)
# gene_positions_hi <- biomaRt::getBM(filters = 'ensembl_gene_id',attributes=c('ensembl_gene_id','hgnc_symbol'), values = newEG, mart = ensembl)
# 
# ke = which(! gene_positions_nova$ensembl_gene_id %in% gene_positions_hi$ensembl_gene_id) #Only finds the genes found in the HiSeq Gene Symbol LUT but not the NovaSeq Gene Symbol LUT
# ke2 = gene_positions_nova[ke,]
# Combined_gene_positions = rbind(gene_positions_hi,ke2)

Combined_gene_positions = gene_positions_nova

#####################################################################################################################################################################################

##Parse VST Datasets for Classifier Development (run this section for both sequencing platform cohorts)
##NovaSeq
setwd("G:/SpinalCord/Publication/RawExpression/WithControls/NovaSeq_4Covar")
filen = "NOVASEQ_SpinalCord_FULLCohort_VST_Gene_NoFiltering_TE_9-27-23.csv" #File generated in previous script
vstcounts = read.csv(filen)
rownames(vstcounts) = vstcounts[,1]
vstcounts = vstcounts[,-1]


#Covert Ensembl rownames to gene symbol rownames
blank = rep(NA,nrow(vstcounts))
vstsymbol = cbind(blank,vstcounts)
vstsymbol[,1] = rownames(vstcounts)

#Vectorized for speed - 9/11/23
for(i in 1:nrow(vstsymbol)){
  
  matchind = which(gene_positions$ensembl_gene_id == sub("\\..*","",vstsymbol[i,1]))
  if(length(matchind)==1){
    vstsymbol[i,1] = gene_positions$hgnc_symbol[matchind]
  }else if(length(matchind) > 1){
    ORF1 = grep("^C.orf",gene_positions$hgnc_symbol[matchind])
    ORF2 = grep("^C..orf",gene_positions$hgnc_symbol[matchind])
    ORFind = as.numeric(c(ORF1,ORF2))
    if(length(ORFind)==1){
      matchind = matchind[-ORFind]
      vstsymbol[i,1] = gene_positions$hgnc_symbol[matchind]
    }else if(length(ORFind)>1){
      matchind = matchind[length(matchind)]
      vstsymbol[i,1] = gene_positions$hgnc_symbol[matchind]
    }
    
  }
  
  
  if((i %% 1000) == 0) cat("% Done:",i/nrow(vstsymbol)*100,"\n")
}

#Replace missing gene symbol with original ENSEMBL ID
for(i in 1:nrow(vstsymbol)){
  if(vstsymbol[i,1] == ""){
    vstsymbol[i,1] = rownames(vstsymbol)[i]
  }
}

novacheckpoint = vstsymbol
finalvstsymbol = vstsymbol

rownames(finalvstsymbol) = make.unique(finalvstsymbol[,1])
finalvstsymbol = finalvstsymbol[,-1]
rownames(finalvstsymbol) = make.unique(gsub('\\..*','',rownames(finalvstsymbol)))
finalvstsymbol =data.frame(finalvstsymbol)
FullNovaSymbol = finalvstsymbol


##HISEQ
setwd("G:/SpinalCord/Publication/RawExpression/WithControls/HiSeq_4Covar")
filen = "HISEQ_SpinalCord_FULLCohort_VST_Gene_NoFiltering_TE_9-27-23.csv"
vstcounts = read.csv(filen) #File generated in previous script
rownames(vstcounts) = vstcounts[,1]
vstcounts = vstcounts[,-1]


#Covert Ensembl rownames to gene symbol rownames
blank = rep(NA,nrow(vstcounts))
vstsymbol = cbind(blank,vstcounts)
vstsymbol[,1] = rownames(vstcounts)

#Vectorized for speed - 9/11/23
for(i in 1:nrow(vstsymbol)){
  
  matchind = which(Combined_gene_positions$ensembl_gene_id == sub("\\..*","",vstsymbol[i,1]))
  if(length(matchind)==1){
    vstsymbol[i,1] = Combined_gene_positions$hgnc_symbol[matchind]
  }else if(length(matchind) > 1){
    ORF1 = grep("^C.orf",Combined_gene_positions$hgnc_symbol[matchind])
    ORF2 = grep("^C..orf",Combined_gene_positions$hgnc_symbol[matchind])
    ORFind = as.numeric(c(ORF1,ORF2))
    if(length(ORFind)>0){
      matchind = matchind[-ORFind]
      vstsymbol[i,1] = Combined_gene_positions$hgnc_symbol[matchind]
    }else{
      matchind = matchind[length(matchind)]
      vstsymbol[i,1] = Combined_gene_positions$hgnc_symbol[matchind]
    }
    
  }
  
  
  if((i %% 1000) == 0) cat("% Done:",i/nrow(vstsymbol)*100,"\n")
}

#Replace missing gene symbol with original ENSEMBL ID
for(i in 1:nrow(vstsymbol)){
  if(vstsymbol[i,1] == ""){
    vstsymbol[i,1] = rownames(vstsymbol)[i]
  }
}

hicheckpoint = vstsymbol
finalvstsymbol = vstsymbol

rownames(finalvstsymbol) = make.unique(finalvstsymbol[,1])
finalvstsymbol = finalvstsymbol[,-1]
rownames(finalvstsymbol) = make.unique(gsub('\\..*','',rownames(finalvstsymbol)))
finalvstsymbol =data.frame(finalvstsymbol)
FullHiSymbol = finalvstsymbol


wd = 'G:/SpinalCord/Publication/Enrichment/FullCohort/4Covar'
setwd(wd)
save.image("SAKEfeatures4Enrichment_SpinalCord_withControls_4Covar_9-30-23.RData")
#load("G:/SpinalCord/Publication/Enrichment/FullCohort/SAKEfeatures4Enrichment_SpinalCord_withControls_NoOligo_8-15-23.RData")



##Correct all feature set lists for excel autoformat driven errors (MARCH3 --> 3-Mar)
myfeaturesets = c("Ctop1000ALL","Ctop2500ALL","Ctop5000ALL")
excelgenelist = c("1-Sep","2-Sep","3-Sep","4-Sep","5-Sep","6-Sep","7-Sep","8-Sep","9-Sep","10-Sep","11-Sep","12-Sep","1-Mar","2-Mar","3-Mar","4-Mar","5-Mar","6-Mar","7-Mar","8-Mar","9-Mar","10-Mar","11-Mar","12-Mar")
truegenelist = c("SEPT1","SEPT2","SEPT3","SEPT4","SEPT5","SEPT6","SEPT7","SEPT8","SEPT9","SEPT10","SEPT11","SEPT12","MARCH1","MARCH2","MARCH3","MARCH4","MARCH5","MARCH6","MARCH7","MARCH8","MARCH9","MARCH10","MARCH11","MARCH12")

ExcelerrorLUT = data.frame(cbind(excelgenelist,truegenelist))

for(i in 1:length(myfeaturesets)){
  
  tmp = get(myfeaturesets[i])
  
  for(j in 1:length(tmp)){
    for(k in 1:nrow(ExcelerrorLUT)){
      
      if(tmp[j] == ExcelerrorLUT$excelgenelist[k]){
        
        tmp[j] = ExcelerrorLUT$truegenelist[k]
        
      }
      
    }
  }
  nam = myfeaturesets[i]
  assign(paste(nam),tmp)
  if((i %% 1) == 0) cat("Feature Set Completed:",myfeaturesets[i],"\n")
} 


############################### Fix TE names (all symbols --> .) & most transcripts with a dash
#Top feature list - 1k
Ctop1000ALL[which(nchar(Ctop1000ALL)>20)] = gsub("\\|",".",Ctop1000ALL[which(nchar(Ctop1000ALL)>20)])
Ctop1000ALL[which(nchar(Ctop1000ALL)>20)] = gsub("\\:",".",Ctop1000ALL[which(nchar(Ctop1000ALL)>20)])
Ctop1000ALL[which(nchar(Ctop1000ALL)>20)] = gsub("-",".",Ctop1000ALL[which(nchar(Ctop1000ALL)>20)])
Ctop1000ALL[which(nchar(Ctop1000ALL)>20)] = gsub("\\+",".",Ctop1000ALL[which(nchar(Ctop1000ALL)>20)])
Ctop1000ALL[which(nchar(Ctop1000ALL)>20)] = gsub("\\_",".",Ctop1000ALL[which(nchar(Ctop1000ALL)>20)])
#Top feature list - 2.5k
Ctop2500ALL[which(nchar(Ctop2500ALL)>20)] = gsub("\\|",".",Ctop2500ALL[which(nchar(Ctop2500ALL)>20)])
Ctop2500ALL[which(nchar(Ctop2500ALL)>20)] = gsub("\\:",".",Ctop2500ALL[which(nchar(Ctop2500ALL)>20)])
Ctop2500ALL[which(nchar(Ctop2500ALL)>20)] = gsub("-",".",Ctop2500ALL[which(nchar(Ctop2500ALL)>20)])
Ctop2500ALL[which(nchar(Ctop2500ALL)>20)] = gsub("\\+",".",Ctop2500ALL[which(nchar(Ctop2500ALL)>20)])
Ctop2500ALL[which(nchar(Ctop2500ALL)>20)] = gsub("\\_",".",Ctop2500ALL[which(nchar(Ctop2500ALL)>20)])
#Top feature list - ALL
Ctop5000ALL[which(nchar(Ctop5000ALL)>20)] = gsub("\\|",".",Ctop5000ALL[which(nchar(Ctop5000ALL)>20)])
Ctop5000ALL[which(nchar(Ctop5000ALL)>20)] = gsub("\\:",".",Ctop5000ALL[which(nchar(Ctop5000ALL)>20)])
Ctop5000ALL[which(nchar(Ctop5000ALL)>20)] = gsub("-",".",Ctop5000ALL[which(nchar(Ctop5000ALL)>20)])
Ctop5000ALL[which(nchar(Ctop5000ALL)>20)] = gsub("\\+",".",Ctop5000ALL[which(nchar(Ctop5000ALL)>20)])
Ctop5000ALL[which(nchar(Ctop5000ALL)>20)] = gsub("\\_",".",Ctop5000ALL[which(nchar(Ctop5000ALL)>20)])
#NovaSeq
rownames(FullNovaSymbol)[which(nchar(rownames(FullNovaSymbol))>20)] = gsub("\\|",".",rownames(FullNovaSymbol)[which(nchar(rownames(FullNovaSymbol))>20)])
rownames(FullNovaSymbol)[which(nchar(rownames(FullNovaSymbol))>20)] = gsub("\\:",".",rownames(FullNovaSymbol)[which(nchar(rownames(FullNovaSymbol))>20)])
rownames(FullNovaSymbol)[which(nchar(rownames(FullNovaSymbol))>20)] = gsub("-",".",rownames(FullNovaSymbol)[which(nchar(rownames(FullNovaSymbol))>20)])
rownames(FullNovaSymbol)[which(nchar(rownames(FullNovaSymbol))>20)] = gsub("\\+",".",rownames(FullNovaSymbol)[which(nchar(rownames(FullNovaSymbol))>20)])
rownames(FullNovaSymbol)[which(nchar(rownames(FullNovaSymbol))>20)] = gsub("\\_",".",rownames(FullNovaSymbol)[which(nchar(rownames(FullNovaSymbol))>20)])
rownames(FullNovaSymbol) = gsub("-",".",rownames(FullNovaSymbol))
#HiSeq
rownames(FullHiSymbol)[which(nchar(rownames(FullHiSymbol))>20)] = gsub("\\|",".",rownames(FullHiSymbol)[which(nchar(rownames(FullHiSymbol))>20)])
rownames(FullHiSymbol)[which(nchar(rownames(FullHiSymbol))>20)] = gsub("\\:",".",rownames(FullHiSymbol)[which(nchar(rownames(FullHiSymbol))>20)])
rownames(FullHiSymbol)[which(nchar(rownames(FullHiSymbol))>20)] = gsub("-",".",rownames(FullHiSymbol)[which(nchar(rownames(FullHiSymbol))>20)])
rownames(FullHiSymbol)[which(nchar(rownames(FullHiSymbol))>20)] = gsub("\\+",".",rownames(FullHiSymbol)[which(nchar(rownames(FullHiSymbol))>20)])
rownames(FullHiSymbol)[which(nchar(rownames(FullHiSymbol))>20)] = gsub("\\_",".",rownames(FullHiSymbol)[which(nchar(rownames(FullHiSymbol))>20)])
rownames(FullHiSymbol) = gsub("-",".",rownames(FullHiSymbol))
#Check
rownames(FullNovaSymbol)[which(nchar(rownames(FullNovaSymbol))>20)]
rownames(FullHiSymbol)[which(nchar(rownames(FullHiSymbol))>20)]


######################Remove platform-specific covar-dependent genes - Top 1k
#NOVASEQ
# MAD10k = FullNovaSymbol
# platdif = which(! Ctop1000ALL %in% rownames(MAD10k))
# Ctop1000ALL2 = Ctop1000ALL[-platdif]
# #HISEQ
# MAD10k = FullHiSymbol
# h1 = which(! Ctop1000ALL2 %in% rownames(MAD10k))
# Ctop1000ALL3 = Ctop1000ALL2[-h1]

#Confirmed: Genes specific to a single sequencing platform were removed as covar-dependent, as determined by DESeq2
############################## 
#Ctop1000ALL[platdif] #missing in NovaSeq (covar-depedent)
#Ctop1000ALL2[h1] #missing in HiSeq (covar-depedent)

#library(stringr)
#which(str_sub(rownames(FullNovaSymbol),1,8) == "chr11.23")

# ind = which(Combined_gene_positions$hgnc_symbol == "TMEM114")
# Combined_gene_positions$ensembl_gene_id[ind]
# which(rownames(novacheckpoint) == Combined_gene_positions$ensembl_gene_id[ind])
# which(rownames(FullNovaSymbol) == "TMEM114")
############################## 

##NovaSeq
MAD10k = FullNovaSymbol
CA1000N = MAD10k[rownames(MAD10k) %in% Ctop1000ALL,]

##Hiseq
MAD10k = FullHiSymbol
CA1000H = MAD10k[rownames(MAD10k) %in% Ctop1000ALL,]

#Check that you can column bind the two platform VST count matrices
table(rownames(CA1000N) == rownames(CA1000H))


#Final Enrichment Matrix
CA1000 = cbind(CA1000N,CA1000H)


#Manually save feature set 
wd = 'G:/SpinalCord/Publication/Enrichment/FullCohort/4Covar'
setwd(wd)

filen2 = "SpinalCord_CombinedPlatform_NovaSeqSamples_Top1000.csv"
write.csv(CA1000N,filen2) 

filen2 = "SpinalCord_CombinedPlatform_HiSeqSamples_Top1000.csv"
write.csv(CA1000H,filen2) 

filen2 = "SpinalCord_CombinedPlatform_AllSamples_Top1000.csv"
write.csv(CA1000,filen2) 

#Full Feature Sets

filen2 = "SpinalCord_CombinedPlatform_NovaSeqSamples_VST_SYMBOL_Top2500.csv"
write.csv(FullNovaSymbol,filen2) 

filen2 = "SpinalCord_CombinedPlatform_HiSeqSamples_VST_SYMBOL_Top2500.csv"
write.csv(FullHiSymbol,filen2)

############Remove platform-specific sex-dependent genes - Top 2.5k
#NOVASEQ
#MAD10k = FullNovaSymbol
#platdif = which(! Ctop2500ALL %in% rownames(MAD10k))
#Ctop2500ALL2 = Ctop2500ALL[-platdif]
#HISEQ
#MAD10k = FullHiSymbol
#h1 = which(! Ctop2500ALL %in% rownames(MAD10k))
#Ctop2500ALL3 = Ctop2500ALL2[-h1]

#Ctop2500ALL[platdif] #missing in NovaSeq (covar-depedent)
#Ctop2500ALL2[h1] #missing in HiSeq (covar-depedent)

##NovaSeq
MAD10k = FullNovaSymbol
CA2500N = MAD10k[rownames(MAD10k) %in% Ctop2500ALL,]

##Hiseq
MAD10k = FullHiSymbol
CA2500H = MAD10k[rownames(MAD10k) %in% Ctop2500ALL,]

#Check that you can column bind the two platform VST count matrices
table(rownames(CA2500N) == rownames(CA2500H))


#Final Enrichment Matrix
CA2500 = cbind(CA2500N,CA2500H)


#Manually save feature set 
wd = 'G:/SpinalCord/Publication/Enrichment/FullCohort/4Covar'
setwd(wd)

filen2 = "SpinalCord_CombinedPlatform_NovaSeqSamples_Top2500_v2.csv"
write.csv(CA2500N,filen2) 

filen2 = "SpinalCord_CombinedPlatform_HiSeqSamples_Top2500_v2.csv"
write.csv(CA2500H,filen2) 

filen2 = "SpinalCord_CombinedPlatform_AllSamples_Top2500_v2.csv"
write.csv(CA2500,filen2) 

#############################################################################################

############Remove platform-specific sex-dependent genes - ALL
#NOVASEQ
#MAD10k = FullNovaSymbol
#platdif = which(! Ctop5000ALL %in% rownames(MAD10k))
#Ctop2500ALL2 = Ctop5000ALL[-platdif]
#HISEQ
#MAD10k = FullHiSymbol
#h1 = which(! Ctop5000ALL %in% rownames(MAD10k))
#Ctop2500ALL3 = Ctop2500ALL2[-h1]

#Ctop5000ALL[platdif] #missing in NovaSeq (covar-depedent)
#Ctop2500ALL2[h1] #missing in HiSeq (covar-depedent)

##NovaSeq
MAD10k = FullNovaSymbol
CA5000N = MAD10k[rownames(MAD10k) %in% Ctop5000ALL,]

##Hiseq
MAD10k = FullHiSymbol
CA5000H = MAD10k[rownames(MAD10k) %in% Ctop5000ALL,]

#Check that you can column bind the two platform VST count matrices
table(rownames(CA5000N) == rownames(CA5000H))


#Final Enrichment Matrix
CA5000 = cbind(CA5000N,CA5000H)


#Manually save feature set 
wd = 'G:/SpinalCord/Publication/Enrichment/FullCohort/4Covar'
setwd(wd)

filen2 = "SpinalCord_CombinedPlatform_NovaSeqSamples_Top5000.csv"
write.csv(CA5000N,filen2) 

filen2 = "SpinalCord_CombinedPlatform_HiSeqSamples_Top5000.csv"
write.csv(CA5000H,filen2) 

filen2 = "SpinalCord_CombinedPlatform_AllSamples_Top5000.csv"
write.csv(CA5000,filen2) 

#############################################################################################

#Phenotype list

PhenoLUT = data.frame(matrix(NA,nrow=ncol(CA1000),ncol=3))
colnames(PhenoLUT) = c("Sample","Subtype","Group")
PhenoLUT$Sample = colnames(CA1000)


NovaLabels = read.csv("G:/SpinalCord/Publication/Clustering/NovaSeq/4Covar_5k/RobustSubtypeAssignment/SpinalCord_NovaSeq_RobustSubtypeAssignment_4Covariate_9-18-23_majority.csv")
rownames(NovaLabels) = NovaLabels[,1];NovaLabels = NovaLabels[,-1]
colnames(NovaLabels) = sub("^[^_]*_", "", colnames(NovaLabels))#All text to the right of the underscore

HiLabels = read.csv("G:/SpinalCord/Publication/Clustering/HiSeq/4Covar_5k/RobustSubtypeAssignment/SpinalCord_HiSeq_RobustSubtypeAssignment_4Covariate_9-18-23_majority.csv")
rownames(HiLabels) = HiLabels[,1];HiLabels = HiLabels[,-1]
colnames(HiLabels) = sub("^[^_]*_", "", colnames(HiLabels))#All text to the right of the underscore

#Fill NovaSeq subtypes
for(i in 1:nrow(PhenoLUT)){
  for(j in 1:ncol(NovaLabels)){
    if(PhenoLUT$Sample[i] == colnames(NovaLabels)[j]){
      PhenoLUT$Subtype[i] = NovaLabels[12,j]
    }
  }
}

#Fill HiSeq subtypes
for(i in 1:nrow(PhenoLUT)){
  for(j in 1:ncol(HiLabels)){
    if(PhenoLUT$Sample[i] == colnames(HiLabels)[j]){
      PhenoLUT$Subtype[i] = HiLabels[12,j]
    }
  }
}

Clinical = read.csv("G:/SpinalCord/Publication/CLINICAL_DATA_PRUDENCIO.csv")
Clinical$ExternalSampleId = gsub("-",".",Clinical$ExternalSampleId)

for(i in 1:nrow(PhenoLUT)){
  for(j in 1:nrow(Clinical)){
    if(PhenoLUT$Sample[i] == Clinical$ExternalSampleId[j]){
      PhenoLUT$Group[i] = Clinical$Subject.Group[j]
    }
  }
}

table(colnames(CA1000) == PhenoLUT$Sample)
table(PhenoLUT$Group[which(is.na(PhenoLUT$Subtype))]) #All missing labels are controls 

PhenoLUT$Subtype[which(is.na(PhenoLUT$Subtype))] = "CTR"

PhenoLUT$Factor = NA

for(i in 1:nrow(PhenoLUT)){
  if(PhenoLUT$Subtype[i] == "TD"){
    PhenoLUT$Factor[i] = 1
  }else if(PhenoLUT$Subtype[i] =="GLIA"){
    PhenoLUT$Factor[i] = 2
  }else if(PhenoLUT$Subtype[i] == "OX"){
    PhenoLUT$Factor[i] = 3
  }else if(PhenoLUT$Subtype[i] == "CTR"){
    PhenoLUT$Factor[i] = 4
  }
}

wd = 'G:/SpinalCord/Publication/Enrichment/FullCohort/4Covar'
setwd(wd)
write.csv(PhenoLUT,"Samplewise_Subtypes_SuppData_9-27-23.csv")

#For .cls file; Manually edit in notepad++
#paste(PhenoLUT$Factor,sep="",collapse="")

writeLines(paste(PhenoLUT$Factor,sep="",collapse="\t"),"Pheno_4Covar.txt")

#Remove ENSG
tmp = CA1000
removeind1 = which(substr(rownames(tmp),1,4) == "ENSG")
hist(nchar(rownames(tmp))); rownames(tmp)[which(nchar(rownames(tmp))>25)]
removeind2 = which(nchar(rownames(tmp))>25)
removeind = as.numeric(names(table(c(removeind1,removeind2))))
tmp2 = tmp[-removeind,]
write.csv(tmp2,"SpinalCord_CombinedPlatform_AllSamples_Top1000_NoTE_NoENSG.csv")

tmp = CA2500
removeind1 = which(substr(rownames(tmp),1,4) == "ENSG")
hist(nchar(rownames(tmp))); rownames(tmp)[which(nchar(rownames(tmp))>25)]
removeind2 = which(nchar(rownames(tmp))>25)
removeind = as.numeric(names(table(c(removeind1,removeind2))))
tmp2 = tmp[-removeind,]
write.csv(tmp2,"SpinalCord_CombinedPlatform_AllSamples_Top2500_NoTE_NoENSG_9-30-23.csv")

tmp = CA5000
removeind1 = which(substr(rownames(tmp),1,4) == "ENSG")
hist(nchar(rownames(tmp))); rownames(tmp)[which(nchar(rownames(tmp))>25)]
removeind2 = which(nchar(rownames(tmp))>25)
removeind = as.numeric(names(table(c(removeind1,removeind2))))
tmp2 = tmp[-removeind,]
write.csv(tmp2,"SpinalCord_CombinedPlatform_AllSamples_Top5000_NoTE_NoENSG_9-30-23.csv")
