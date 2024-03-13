#######################################  IDENTIFY SEX-DEPENDENT GENES using DESeq2  #################################################

#Written by: Jarrett Eshima
#Date: April, 2021
#For: Use by the Dr. Barbara Smith Lab at Arizona State University

#Note: This code is not fully automated, as it has been designed to run in individual but related parts.
#Some parts can be skipped, depending on the downstream analysis you are interested in.

## Study Meta Information --
# n = 439 donors/subjects
# nsamples = 1659 RNA-seq files
# ntissues = 11
# GEO Series: GSE153960
# GEO SuperSeries: GSE137810

# n = 473 tissue-matched ALS patient samples
# n = 451 fully processed samples
# n = 140 from GSE124439
# n = 311 from GSE153960 
#######################################################################################################################################

#Load Libraries
library(dplyr)
library(Biobase)
library(limma)
library(DESeq2)
library(biomaRt)

########################################################################################################################################
#############################################  PART 1: LOAD COUNT DATA   #####################################################################
########################################################################################################################################
#SPECIFIC Gene Expression Omnibus (GEO) Study - GSE153960
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153960
#https://www.jci.org/articles/view/139741


setwd("G:/SpinalCord/Publication/RawExpression")
first = read.table(gzfile("GSE153960_Gene_counts_matrix_RSEM_Prudencio_et_al_2020.txt.gz"),sep="\t") #Available online at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153960
dim(first)


GSE15_RawCounts = first
backupcounts = GSE15_RawCounts
rownames(GSE15_RawCounts) = first[,1]
GSE15_RawCounts = GSE15_RawCounts[,-1]

colnames(GSE15_RawCounts) = substr(colnames(GSE15_RawCounts),1,14)

#which(substr(colnames(GSE15_RawCounts),1,14) == "CGND.HRA.00027")
#which(substr(colnames(GSE15_RawCounts),1,14) == "CGND.HRA.00530")

########################################################################################################################################
#############################################  PART 2: LOAD TE DATA   #####################################################################
########################################################################################################################################


TECounts = read.csv("NovaSeq_SpinalCord_ALSCohort_TECounts_HGND_6-3-23_nodups.csv") #Provided in Supplemental Dataset 1
rownames(TECounts) = TECounts[,1]
TECounts = TECounts[,-1]

ctrTECounts = read.csv("NovaSeq_SpinalCord_ControlCohort_TECounts_HGND_8-15-23.csv") #Provided in Supplemental Dataset 1
rownames(ctrTECounts) = ctrTECounts[,1]
ctrTECounts = ctrTECounts[,-1]

table(rownames(TECounts) == rownames(ctrTECounts))
TECounts = cbind(TECounts,ctrTECounts)


#Missing Sample checks
colnames(TECounts)[! colnames(TECounts) %in% colnames(GSE15_RawCounts)]# Which samples do we have TE counts for but not gene counts
#colnames(GSE15_RawCounts)[which(substr(colnames(GSE15_RawCounts),1,13) == "CGND.HRA.0053")]


HS_SpinalCord = GSE15_RawCounts[,colnames(GSE15_RawCounts) %in% colnames(TECounts)]
HS_TECounts = TECounts[,colnames(TECounts) %in% colnames(HS_SpinalCord)]

table(colnames(HS_SpinalCord) == colnames(HS_TECounts)) #Are our samples in the same order


#Reorder to allow row binding
tmp = data.frame(matrix(NA,nrow(HS_TECounts),ncol(HS_TECounts)))
rownames(tmp) = rownames(HS_TECounts)
colnames(tmp) = colnames(HS_SpinalCord)

for(i in 1:ncol(HS_TECounts)){
  
  for(j in 1:ncol(tmp)){
    
    if(colnames(HS_TECounts)[i] == colnames(tmp)[j]){
      
      tmp[,j] = HS_TECounts[,i]
      
    }
    
  }
  
}

HS_TECounts_org = tmp
table(colnames(HS_SpinalCord) == colnames(HS_TECounts_org))

SpinalCounts = rbind(HS_SpinalCord,HS_TECounts_org)
rownames(SpinalCounts) = make.unique(sub("\\..*","",rownames(SpinalCounts)))

########################################################################################################################################
#############################################  PART 3: LOAD PHENOTYPE DATA   #####################################################################
########################################################################################################################################
setwd("G:/SpinalCord/Publication")

#Coldata matrix for DESeq2
clinicaldata = read.csv("CLINICAL_DATA_PRUDENCIO.csv") 

##IMPORTANT: THERE ARE SAMPLES NOT INCLUDED IN THE GENE EXPRESSION MATRIX THAT WERE CONSIDERED BY PRUDENCIO AND HUMPHREYS
# THESE SAMPLES WILL NEED TO GO THROUGH STAR RSEM PIPELINE TOO
list = gsub("-",".",clinicaldata$ExternalSampleId)
table(list %in% colnames(TECounts))

coldata = clinicaldata
coldata$ExternalSampleId = gsub("-",".",coldata$ExternalSampleId)
coldata = coldata[coldata$ExternalSampleId %in% colnames(SpinalCounts),]

#Reorder coldata
tmp = data.frame(matrix(NA,nrow(coldata),ncol(coldata)))
rownames(tmp) = colnames(SpinalCounts)

for(i in 1:nrow(coldata)){
  
  for(j in 1:nrow(tmp)){
    
    if(coldata$ExternalSampleId[i] == rownames(tmp)[j]){
      
      tmp[j,] = coldata[i,]
      
    }
    
  }
  
}

colnames(tmp) = colnames(coldata)
coldata2 = tmp

setwd("G:/SpinalCord/Publication/MetaData")
Meta = read.csv("GSE153960Meta.txt") #Publicly available at: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA644618&o=acc_s%3Aa
setwd("G:/SpinalCord/Publication")
Clinical = read.csv("CLINICAL_DATA_PRUDENCIO.csv") #Provided by NYGC, by request
convertnames = gsub("-","\\.",Meta$sample_id_alt)
Meta$sample_id_alt = convertnames


FiltMeta = Meta[Meta$sample_id_alt %in% rownames(coldata2),]


reftable = table(FiltMeta$sample_id_alt)

removeind = rep(NA,nrow(FiltMeta)-length(reftable))
count = 1
for(i in 1:length(reftable)) {
  
  if(reftable[[i]] > 1){
    
    tmp = names(reftable[i])
    inds = which(FiltMeta$sample_id_alt == tmp)
    
    if(length(inds) == 2){
      removeind[count] = inds[1]
      count = count+1
    }else if(length(inds) == 3){
      removeind[seq(count,count+1,1)] = inds[1:2]
      count = count+2
    }else if(length(inds) > 3){
      cat("Problem")
    }
    
    
  }
  
}

SiteMeta = FiltMeta[-removeind,]



coldata2$Site = NA

for(i in 1:nrow(coldata2)){
  
  for(j in 1:nrow(SiteMeta)){
    
    if(rownames(coldata2)[i] == SiteMeta$sample_id_alt[j]){
      
      coldata2$Site[i] = SiteMeta$Project[j]
      
    }
    
  }
  
}

#Clean up site names
for(i in 1:nrow(coldata2)){
  
  if(coldata2$Site[i] == "NYGC ALS Consortium"){
    coldata2$Site[i] = "NYGC"
  }else{
    coldata2$Site[i] = "TargetALS"
  }
  
}

table(coldata2$Site)

colnames(coldata2)[11] = "Tissue"


#Check DESeq matrix requirements (all lines should be TRUE)
all(rownames(coldata2) == colnames(SpinalCounts))
ncol(SpinalCounts)==nrow(coldata2)
SpinalCounts = data.matrix(SpinalCounts)
all(is.numeric(SpinalCounts))

########################################################################################################################################
#############################################  PART 4: REMOVE GLIA MARKER GENES   #####################################################################
########################################################################################################################################

#Reference marker genes: https://www.nature.com/articles/s41593-022-01205-3

#Table S6
setwd("C:/Users/jeshima/Documents/Smith Lab/Fall 2023/Spinal Cord Project/MarkerGenes")
scSpinalMarkers = read.csv("JH_CNSMarkerGenes.csv")
colnames(scSpinalMarkers) = scSpinalMarkers[1,];scSpinalMarkers = scSpinalMarkers[-1,]
#table(scSpinalMarkers$term_id)

scolig = scSpinalMarkers$gene[which(scSpinalMarkers$term_id == "Oligo")]
scolig2 = scSpinalMarkers$gene[which(scSpinalMarkers$term_id == "Oligos")]
scolig3 = scSpinalMarkers$gene[which(scSpinalMarkers$term_id == "Oligodendrocytes")]
scastro = scSpinalMarkers$gene[which(scSpinalMarkers$term_id == "Astrocyte")]
scastro2 = scSpinalMarkers$gene[which(scSpinalMarkers$term_id == "Astrocytes")]
scmg = scSpinalMarkers$gene[which(scSpinalMarkers$term_id == "Microglia")]
scendo = scSpinalMarkers$gene[which(scSpinalMarkers$term_id == "Endothelial")]
scendo2 = scSpinalMarkers$gene[which(scSpinalMarkers$term_id == "Endothelial cells")]


gmarkers = names(table(c(scolig,scolig2,scolig3,scastro,scastro2,scmg,scendo,scendo2)))

newEG = gmarkers
ensembl_version = "https://apr2020.archive.ensembl.org"
species="human"
ensembl = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host=ensembl_version)
gene_positions <- biomaRt::getBM(filters = 'hgnc_symbol',attributes=c('ensembl_gene_id','hgnc_symbol'), values = newEG, mart = ensembl)

ensoligmarkers = c(gene_positions$ensembl_gene_id)

SpinalCounts_noolig = SpinalCounts[! rownames(SpinalCounts) %in% ensoligmarkers,]

dim(SpinalCounts);dim(SpinalCounts_noolig)
########################################################################################################################################
#############################################  PART 5: SEX-DEPENDENT GENE REMOVAL (DESeq2)   #####################################################################
########################################################################################################################################

setwd("G:/SpinalCord/Publication/RawExpression/WithControls/NovaSeq_4Covar")

###HiSeq Cohort

#Rounding non-integer counts recommended by RSEM authors for DESeq differential expression
rALSCountData = round(SpinalCounts_noolig,0)

dds = DESeqDataSetFromMatrix(countData = rALSCountData, colData = coldata2, design= ~ Sex, tidy=F)
dds$Sex = relevel(dds$Sex,ref = "Male")
dseq = DESeq(dds,betaPrior=T)
res = results(dseq)
sig1 = res[! is.na(res$padj) & res$padj<0.05,]
#write.csv(sig1,"NOVASEQ_SpinalCord_ALSCohort_NoGliaMarkers1k_NoRIN05_SexSigGenes_9-11-23.csv")


#Single sample with missing RIN must be removed (incomplete design equation)
FullPheno = coldata2[-which(is.na(coldata2$RIN)),]
rALSCountData_RIN = round(SpinalCounts_noolig[,-which(is.na(coldata2$RIN))],0)

FullPheno$RIN = scale(FullPheno$RIN,center = T)

#Covarate: RIN
dds = DESeqDataSetFromMatrix(countData = rALSCountData_RIN, colData = FullPheno, design= ~ RIN, tidy=F)
dseq = DESeq(dds,betaPrior=T)
res = results(dseq)
sig2 = res[! is.na(res$padj) & res$padj<0.05,]
#write.csv(sig2,"NOVASEQ_SpinalCord_ALSCohort_4Covar_RINSigGenes_9-15-23.csv")

#Covariate: Site of sample collection
dds = DESeqDataSetFromMatrix(countData = rALSCountData, colData = coldata2, design= ~ Site, tidy=F)
dds$Site = relevel(dds$Site,ref = "NYGC")
dseq = DESeq(dds,betaPrior=T)
res = results(dseq)
sig3 = res[! is.na(res$padj) & res$padj<0.05,]
#write.csv(sig3,"NOVASEQ_SpinalCord_ALSCohort_4Covar_SiteSigGenes_9-15-23.csv")

#Covariate: Tissue
dds = DESeqDataSetFromMatrix(countData = rALSCountData, colData = coldata2, design= ~ Tissue, tidy=F)
dseq = DESeq(dds,betaPrior=T)
res = results(dseq)
sig4 = res[! is.na(res$padj) & res$padj<0.05,]
#write.csv(sig4,"NOVASEQ_SpinalCord_ALSCohort_4Covar_TissueSigGenes_9-15-23.csv")


#sig = names(table(c(rownames(sig1),rownames(sig2),rownames(sig3),rownames(sig4))))


vsd = varianceStabilizingTransformation(dseq)
filtdata = assay(vsd)
#filtdata = vstcounts[! (rownames(vstcounts) %in% sig),]


file = "NOVASEQ_SpinalCord_FULLCohort_VST_Gene_NoFiltering_TE_9-27-23.csv"
write.csv(filtdata,file)
file = "NOVASEQ_SpinalCord_FULLCohort_VST_Gene_NoFiltering_TE_9-27-23.txt"
write.table(filtdata,file,sep="\t",quote=F,row.names=T,col.names=T)

madval = apply(as.matrix(filtdata),1,mad)
madval = madval[order(madval,decreasing=T)]

#Important Note: There are three "duplicated" gene names in the HiSEQ MAD10k gene set
#Two are genuine duplicates and are removed in ALSPatientStratification_ClusterPrep.R
# mad10000 = filtdata[match(names(madval)[1:10000],rownames(filtdata)),] #This gives slightly less than 10k genes
# file = "NOVASEQ_SpinalCord_FULLCohort_ALS_MAD10k_GeneENSEMBL_NoFiltering_TE_9-11-23.txt"
# write.table(mad10000,file,sep="\t",quote=F,row.names=T,col.names=T)
# file = "NOVASEQ_SpinalCord_FULLCohort_ALS_MAD10k_GeneENSEMBL_NoFiltering_TE_9-11-23.csv"
# write.csv(mad10000,file)


########################################################################################################################################
#############################################  PART 6: CONVERT ENSEMBLs --> HGNC Symbol   #####################################################################
########################################################################################################################################

#Convert ENSEMBL IDs to Gene Symbols
#Important Recommendation: Only the UTY and XIST ENSEMBL IDs are needed to complete Part 2
#You can wait until after filtering by median absolute deviation to convert Ensembl IDs to Gene Symbols (much quicker)
#Some small gene naming issues can occur if you use the Gene Symbol matrix generated in this Part of the code

MAD10k = filtdata
blank = rep(NA,nrow(MAD10k))
MAD10k = cbind(blank,MAD10k)
MAD10k = data.frame(MAD10k)

ensembl.genes = rownames(MAD10k) 
ensembl_version = "https://apr2020.archive.ensembl.org"
species="human"
ensembl = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host=ensembl_version)
gene_positions <- biomaRt::getBM(filters = 'ensembl_gene_id',attributes=c('ensembl_gene_id','hgnc_symbol'), values = ensembl.genes, mart = ensembl)

DONE = F
for(i in 1:nrow(MAD10k)){
  for(j in 1:nrow(gene_positions)){
    if(DONE == F){
      if(sub("\\..*","",rownames(MAD10k)[i]) == gene_positions$ensembl_gene_id[j]){
        MAD10k$blank[i] = gene_positions$hgnc_symbol[j]
        DONE = T
      }
    }
  }
  DONE = F
  if((i %% 100) == 0) cat("% Done:",i/nrow(MAD10k)*100,"\n")
}

#Replace missing gene symbol with original ENSEMBL ID
for(i in 1:nrow(MAD10k)){
  
  if(is.na(MAD10k$blank[i])){
    MAD10k$blank[i] = rownames(MAD10k)[i]
  }
  
  if(MAD10k$blank[i] == ""){
    MAD10k$blank[i] = rownames(MAD10k)[i]
  }
}


uniquern = make.names(MAD10k[,1],unique = T)
rownames(MAD10k) = uniquern
MAD10k_SAKE = MAD10k[,-1]


#Order rownames alphabetically
MAD10k_SAKE = MAD10k_SAKE[order(row.names(MAD10k_SAKE)),]
file = "NOVASEQ_SpinalCord_FULLCohort_GeneSYMBOL_NoFiltering_TE_9-27-23.csv"
write.csv(MAD10k_SAKE,file)

save.image("NovaSeq_MAD10k_Ready4Clustering_withControls_NoFiltering_9-27-23.RData")
#load("G:/SpinalCord/Publication/RawExpression/WithControls/NovaSeq_NoOligo/NovaSeq_MAD10k_Ready4Clustering_withControls_8-15-23.RData")


library(DESeq2)
tmp = estimateSizeFactors(dds)
DESeq_NormalizedCounts_60k = counts(tmp,normalized = T)
DESeq_NormalizedCounts_60k = DESeq_NormalizedCounts_60k[! (rownames(DESeq_NormalizedCounts_60k) %in% rownames(sig)),]
NormCounts = DESeq_NormalizedCounts_60k
write.csv(NormCounts,"NovaSeq_SpinalCord_withControls_MedianofRatios_ENSG_NoFiltering_9-27-23.csv")
