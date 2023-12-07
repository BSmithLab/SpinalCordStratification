##Univariate Analysis - Spinal Cord

#Written By: Jarrett Eshima
#For: Dr. Barbara Smith Lab
#Date: September 12th 2023

#######################################################################################################################################
#Load Libraries
library(dplyr)
library(Biobase)
library(limma)
library(DESeq2)
library(biomaRt)
library(ggplot2)
library(biomaRt)

RecoverDESeqSigGenes = function(DESeqres,nfeatures){
  if(typeof(DESeqres) == "S4"){
    resorder = DESeqres[order(DESeqres$padj,na.last = T,decreasing = F),]
    return(rownames(resorder)[1:nfeatures])
  }else{
    cat("The provided results are not in the expected format. /n")
    cat("Use this function to recover gene/transcript names associated with the top 'nfeatures' number of most differentially expressed features. /n")
  }
  
}

#######################################################################################################################################
#Read in filtered feature set
wd = 'G:/SpinalCord/Publication/Enrichment/FullCohort/4Covar'
setwd(wd)

PrognosticBM = read.csv("DE_FeatureSet_CombinedPlatform_10-11-23.csv",header = F) #Supplemental Table

##Get feature names from classifier matrix
Transcripts = as.character(unlist(PrognosticBM))
Transcripts = c(Transcripts,"TARDBP","TXN","OXR1","UBQLN1","UBQLN2","SOD1","BECN1","BECN2","UCP2") #Add a few genes previously associated with ALS

#Add in Cortex features from previous study. 
#Ref: https://www.nature.com/articles/s41467-022-35494-w
GliaTEs= c("chr1|98004628|98004776|AluJr:Alu:SINE|189|-","chr3|35763881|35764012|AluJb:Alu:SINE|114|+","chr4|21597081|21597444|THE1B:ERVL-MaLR:LTR|145|+","chr10|14299444|14299567|MIR:MIR:SINE|320|-","chr11|133150062|133150220|MIRb:MIR:SINE|304|-","chr12|79068732|79068841|AluJo:Alu:SINE|174|+","chr12|120642840|120642930|L2b:L2:LINE|367|-","chr21|40634798|40635090|L2a:L2:LINE|289|+")
OxTEs= c("chr2|130338399|130338546|L1ME4b:L1:LINE|212|+","chr6|49430916|49431136|LTR86A1:ERVL:LTR|291|-","chr6|116277660|116277934|AluSg:Alu:SINE|44|+","chr8|56958199|56958343|L2b:L2:LINE|303|-","chr14|62107151|62107446|AluJb:Alu:SINE|169|+","chr15|65891440|65891604|MIR3:MIR:SINE|247|+","chr19|46427065|46427223|L2c:L2:LINE|284|+","chr20|36652130|36652423|AluSx1:Alu:SINE|106|+")
TDTEs= c("chr1|35564978|35565131|MIR:MIR:SINE|248|-","chr1|244842569|244842704|L1ME2:L1:LINE|203|+","chr2|7016797|7016993|LTR78:ERV1:LTR|297|-","chr3|179902066|179902288|MIRb:MIR:SINE|355|-","chr6|68832632|68832783|MamSINE1:tRNA-RTE:SINE|351|+","chr13|66953336|66953415|L2a:L2:LINE|279|+","chr17|9935956|9936183|L1M4:L1:LINE|302|+","chrX|54815877|54816014|MER117:hAT-Charlie:DNA|248|-")
my36 = c("AIF1","APOC2","CD44","CHI3L2","CX3CR1","FOLH1","HLA-DRA","TLR7","TMEM125","TNC","TREM2","TYROBP","COL18A1","GABRA1","GAD2","GLRA3","HTR2A","OXR1","SERPINI1","SLC6A13","SLC17A6","TCIRG1","UBQLN2","UCP2","AGPAT4-IT1","CHKB-CPT1B","COL3A1","ENSG00000205041","ENSG00000258674","ENSG00000273151","GATA2-AS1","HSP90AB4P","LINC01347","MIR24-2","MIRLET7BHG","NANOGP4")
S8 = c("ALOX5AP","APOBR","APOC1","CCR5","CD68","CLEC7A","CR1","FPR3","MSR1","NCF2","NINJ2","ST6GALNAC2","TLR8","TNFRSF25","TREM1","VRK2")
S9 = c("B4GALT6","BECN1","COL4A6","COX4I2","CP","GABRA6","GPR22","MYH11","MYL9","NDUFA4L2","NOS3","NOTCH3","PCSK1","SOD1","TAGLN","UBQLN1")
S10 = c("ADAT3","COL6A3","EGLN1P1","ENSG00000263278","ENSG00000268670","ENSG00000279233","ITGBL1","KRT8P13","LINC00176","LINC00638","MIR219A2","NKX6-2","RPS20P22","SLX1B-SULT1A4","TP63","TUB-AS1")
S11 = c(GliaTEs,OxTEs,TDTEs)
S12 = c("AGER","AQP1","BECN2","C1D","CCDC154","ENSG00000260198","ENSG00000278434","ENSG00000279759","ENSG00000281969","HIST1H1T","HLA-DRB1","IFI30","IRF7","SELL","SERPINA1","SNX18P3","SOCS3","STH","TNRC6C-AS1","TUNAR")

CorFeat_Glia = c(my36[1:12],S8)
CorFeat_Ox = c(my36[13:24],S9)
CorFeat_TD = c(my36[25:36],S10)

gmarkers = read.csv("G:/SpinalCord/Publication/DifferentialExpression/4Covar/Overlapping_Independent_10-2-23.csv")[,2]

#allextra = c(CorFeat_Glia,CorFeat_Ox,CorFeat_TD)
allextra = names(table(c(CorFeat_Glia,CorFeat_Ox,CorFeat_TD,S12,gmarkers)))

TranscriptsPlus = c(Transcripts,allextra)
###############################################################################################################################

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

########################################################################################################################################
#############################################  PART 2: LOAD TE DATA   #####################################################################
########################################################################################################################################

### HiSeq

TECounts = read.csv("HiSeq_SpinalCord_ALSCohort_TECounts_HGND_6-3-23.csv")
rownames(TECounts) = TECounts[,1]
TECounts = TECounts[,-1]

ctrTECounts = read.csv("HiSeq_SpinalCord_ControlCohort_TECounts_HGND_8-15-23.csv")
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

SpinalCounts_HS = rbind(HS_SpinalCord,HS_TECounts_org)
rownames(SpinalCounts_HS) = make.unique(sub("\\..*","",rownames(SpinalCounts_HS)))


setwd("G:/SpinalCord/Publication")

#Coldata matrix for DESeq2
clinicaldata = read.csv("CLINICAL_DATA_PRUDENCIO.csv") 

##IMPORTANT: THERE ARE SAMPLES NOT INCLUDED IN THE GENE EXPRESSION MATRIX THAT WERE CONSIDERED BY PRUDENCIO AND HUMPHREYS
# THESE SAMPLES WILL NEED TO GO THROUGH STAR RSEM PIPELINE TOO
list = gsub("-",".",clinicaldata$ExternalSampleId)
table(list %in% colnames(TECounts))

coldata = clinicaldata
coldata$ExternalSampleId = gsub("-",".",coldata$ExternalSampleId)
coldata = coldata[coldata$ExternalSampleId %in% colnames(SpinalCounts_HS),]

#Reorder coldata
tmp = data.frame(matrix(NA,nrow(coldata),ncol(coldata)))
rownames(tmp) = colnames(SpinalCounts_HS)

for(i in 1:nrow(coldata)){
  
  for(j in 1:nrow(tmp)){
    
    if(coldata$ExternalSampleId[i] == rownames(tmp)[j]){
      
      tmp[j,] = coldata[i,]
      
    }
    
  }
  
}

colnames(tmp) = colnames(coldata)
HiSeq_coldata = tmp



#Check DESeq matrix requirements (all lines should be TRUE)
all(rownames(HiSeq_coldata) == colnames(SpinalCounts_HS))
ncol(SpinalCounts_HS)==nrow(HiSeq_coldata)
SpinalCounts_HS = data.matrix(SpinalCounts_HS)
all(is.numeric(SpinalCounts_HS))

setwd("G:/SpinalCord/Publication/RawExpression")

### NovaSeq

TECounts = read.csv("NovaSeq_SpinalCord_ALSCohort_TECounts_HGND_6-3-23_nodups.csv")
rownames(TECounts) = TECounts[,1]
TECounts = TECounts[,-1]

ctrTECounts = read.csv("NovaSeq_SpinalCord_ControlCohort_TECounts_HGND_8-15-23.csv")
rownames(ctrTECounts) = ctrTECounts[,1]
ctrTECounts = ctrTECounts[,-1]

table(rownames(TECounts) == rownames(ctrTECounts))
TECounts = cbind(TECounts,ctrTECounts)


#Missing Sample checks
colnames(TECounts)[! colnames(TECounts) %in% colnames(GSE15_RawCounts)]# Which samples do we have TE counts for but not gene counts
#colnames(GSE15_RawCounts)[which(substr(colnames(GSE15_RawCounts),1,13) == "CGND.HRA.0053")]


NS_SpinalCord = GSE15_RawCounts[,colnames(GSE15_RawCounts) %in% colnames(TECounts)]
NS_TECounts = TECounts[,colnames(TECounts) %in% colnames(NS_SpinalCord)]

table(colnames(NS_SpinalCord) == colnames(NS_TECounts)) #Are our samples in the same order


#Reorder to allow row binding
tmp = data.frame(matrix(NA,nrow(NS_TECounts),ncol(NS_TECounts)))
rownames(tmp) = rownames(NS_TECounts)
colnames(tmp) = colnames(NS_SpinalCord)

for(i in 1:ncol(NS_TECounts)){
  
  for(j in 1:ncol(tmp)){
    
    if(colnames(NS_TECounts)[i] == colnames(tmp)[j]){
      
      tmp[,j] = NS_TECounts[,i]
      
    }
    
  }
  
}

NS_TECounts_org = tmp
table(colnames(NS_SpinalCord) == colnames(NS_TECounts_org))

SpinalCounts_NS = rbind(NS_SpinalCord,NS_TECounts_org)
rownames(SpinalCounts_NS) = make.unique(sub("\\..*","",rownames(SpinalCounts_NS)))

setwd("G:/SpinalCord/Publication")

#Coldata matrix for DESeq2
clinicaldata = read.csv("CLINICAL_DATA_PRUDENCIO.csv") 

##IMPORTANT: THERE ARE SAMPLES NOT INCLUDED IN THE GENE EXPRESSION MATRIX THAT WERE CONSIDERED BY PRUDENCIO AND HUMPHREYS
# THESE SAMPLES WILL NEED TO GO THROUGH STAR RSEM PIPELINE TOO
list = gsub("-",".",clinicaldata$ExternalSampleId)
table(list %in% colnames(TECounts))

coldata = clinicaldata
coldata$ExternalSampleId = gsub("-",".",coldata$ExternalSampleId)
coldata = coldata[coldata$ExternalSampleId %in% colnames(SpinalCounts_NS),]

#Reorder coldata
tmp = data.frame(matrix(NA,nrow(coldata),ncol(coldata)))
rownames(tmp) = colnames(SpinalCounts_NS)

for(i in 1:nrow(coldata)){
  
  for(j in 1:nrow(tmp)){
    
    if(coldata$ExternalSampleId[i] == rownames(tmp)[j]){
      
      tmp[j,] = coldata[i,]
      
    }
    
  }
  
}

colnames(tmp) = colnames(coldata)
NovaSeq_coldata = tmp



#Check DESeq matrix requirements (all lines should be TRUE)
all(rownames(NovaSeq_coldata) == colnames(SpinalCounts_NS))
ncol(SpinalCounts_NS)==nrow(NovaSeq_coldata)
SpinalCounts_NS = data.matrix(SpinalCounts_NS)
all(is.numeric(SpinalCounts_NS))


#####################################################################################################################
########### DESeq2 for identification of subtype-specific significant genes #########################################
#####################################################################################################################

table(rownames(SpinalCounts_NS) == rownames(SpinalCounts_HS))

FullPheno = rbind(NovaSeq_coldata,HiSeq_coldata)
FullCount = cbind(SpinalCounts_NS,SpinalCounts_HS)

table(rownames(FullPheno) == colnames(FullCount))

for(i in 1:nrow(FullPheno)){
  if(FullPheno$Platform[i] == "NovaSeq"){
    
  }else{
    FullPheno$Platform[i] = "HiSeq"
  }
}


#Add RIN and Site to Pheno
setwd("G:/SpinalCord/Publication/MetaData")
Meta = read.csv("GSE153960Meta.txt") #Publicly available at: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA644618&o=acc_s%3Aa
setwd("G:/SpinalCord/Publication")
Clinical = read.csv("CLINICAL_DATA_PRUDENCIO.csv") #Provided by NYGC, by request
convertnames = gsub("-","\\.",Meta$sample_id_alt)
Meta$sample_id_alt = convertnames


FiltMeta = Meta[Meta$sample_id_alt %in% rownames(FullPheno),]


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



FullPheno$Site = NA

for(i in 1:nrow(FullPheno)){
  
  for(j in 1:nrow(SiteMeta)){
    
    if(rownames(FullPheno)[i] == SiteMeta$sample_id_alt[j]){
      
      FullPheno$Site[i] = SiteMeta$Project[j]
      
    }
    
  }
  
}

#Clean up site names
for(i in 1:nrow(FullPheno)){
  
  if(FullPheno$Site[i] == "NYGC ALS Consortium"){
    FullPheno$Site[i] = "NYGC"
  }else{
    FullPheno$Site[i] = "TargetALS"
  }
  
}

table(FullPheno$Site)

#Single sample with missing RIN must be removed (incomplete design equation)
table(rownames(FullPheno) == colnames(FullCount))
FullPheno_sr = FullPheno[-which(is.na(FullPheno$RIN)),]
FullCount_sr = FullCount[,-which(is.na(FullPheno$RIN))]

FullPheno_sr$RIN = scale(FullPheno_sr$RIN,center = T)


#Convert Names
a = which(substr(Transcripts,1,4) == "ENSG")
b = which(nchar(Transcripts)>25)
c = c(a,b)
filtTranscripts = Transcripts[-c]

egenes = filtTranscripts
newEG = sub("\\..*","",egenes) #Keep text to the left of the dot
ensembl_version = "https://apr2020.archive.ensembl.org" 
species="human"
ensembl = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host=ensembl_version)
gene_positions <- biomaRt::getBM(filters = 'hgnc_symbol',attributes=c('ensembl_gene_id','hgnc_symbol'), values = newEG, mart = ensembl)

setwd("G:/SpinalCord/Publication/DifferentialExpression")
Dec16archive = read.csv("Dec16Archive_SD2.csv") #maintain ENSG --> SYMBOL from previous paper

gene_positionsclean = gene_positions[! gene_positions$ensembl_gene_id %in% Dec16archive$ensembl_gene_id,]

gene_positions = rbind(gene_positionsclean,Dec16archive[Dec16archive$hgnc_symbol %in% allextra,])

#No conflicts
#gene_positions[gene_positions$ensembl_gene_id %in% Transcripts[a]]

blank = rep(NA,nrow(FullCount_sr))
FullCount_sr_symbol = cbind(blank,FullCount_sr)
FullCount_sr_symbol[,1] = rownames(FullCount_sr)

#Vectorized for speed - 9/11/23
for(i in 1:nrow(FullCount_sr_symbol)){
  
  matchind = which(gene_positions$ensembl_gene_id == sub("\\..*","",FullCount_sr_symbol[i,1]))
  if(length(matchind)==1){
    FullCount_sr_symbol[i,1] = gene_positions$hgnc_symbol[matchind]
  }else if(length(matchind) > 1){
    ORF1 = grep("^C.orf",gene_positions$hgnc_symbol[matchind])
    ORF2 = grep("^C..orf",gene_positions$hgnc_symbol[matchind])
    ORFind = as.numeric(c(ORF1,ORF2))
    if(length(ORFind)==1){
      matchind = matchind[-ORFind]
      FullCount_sr_symbol[i,1] = gene_positions$hgnc_symbol[matchind]
    }else if(length(ORFind)>1){
      matchind = matchind[length(matchind)]
      FullCount_sr_symbol[i,1] = gene_positions$hgnc_symbol[matchind]
    }
    
  }
  
  
  if((i %% 1000) == 0) cat("% Done:",i/nrow(FullCount_sr_symbol)*100,"\n")
}

#Is missing?
ismissing(FullCount_sr_symbol,my36,usefirstcol = T)


finalvstsymbol = FullCount_sr_symbol

rownames(finalvstsymbol) = make.unique(finalvstsymbol[,1])
finalvstsymbol = finalvstsymbol[,-1]
rownames(finalvstsymbol) = make.unique(gsub('\\..*','',rownames(finalvstsymbol)))
finalvstsymbol =data.frame(finalvstsymbol)

colnames(FullPheno_sr)[11] = "Tissue"

## Add Subtype

Clinical$ExternalSampleId = gsub('-','.',Clinical$ExternalSampleId)

setwd("G:/SpinalCord/Publication/Clustering/NovaSeq/4Covar_5k/RobustSubtypeAssignment")
NovaSeqlabs = read.csv("SpinalCord_NovaSeq_RobustSubtypeAssignment_4Covariate_9-18-23_majority.csv")
rownames(NovaSeqlabs) = NovaSeqlabs[,1]; NovaSeqlabs = NovaSeqlabs[,-1]
setwd("G:/SpinalCord/Publication/Clustering/HiSeq/4Covar_5k/RobustSubtypeAssignment")
HiSeqlabs = read.csv("SpinalCord_HiSeq_RobustSubtypeAssignment_4Covariate_9-18-23_majority.csv")
rownames(HiSeqlabs) = HiSeqlabs[,1]; HiSeqlabs = HiSeqlabs[,-1]

FullPheno_sr$Subtype = NA

for(i in 1:nrow(FullPheno_sr)){
  
  a = which(colnames(NovaSeqlabs) == FullPheno_sr$ExternalSampleId[i])
  if(length(a)>0){
    FullPheno_sr$Subtype[i] = NovaSeqlabs[12,a]
  }
  
  b = which(colnames(HiSeqlabs) == FullPheno_sr$ExternalSampleId[i])
  if(length(b)>0){
    FullPheno_sr$Subtype[i] = HiSeqlabs[12,b]
  }
  
}

for(i in 1:nrow(FullPheno_sr)){
  if(FullPheno_sr$disease_group[i] == "Control" && is.na(FullPheno_sr$Subtype[i])){
    FullPheno_sr$Subtype[i] = "CTR"
  }
}

#No missing subtype labels
#missinglabels = rownames(FullPheno_sr)[which(is.na(FullPheno_sr$Subtype))]


numericsymbol = matrix(as.numeric(unlist(finalvstsymbol)),nrow=nrow(finalvstsymbol))
rownames(numericsymbol) = rownames(finalvstsymbol); colnames(numericsymbol) = colnames(finalvstsymbol)
rCountData_rinsite = round(numericsymbol,0)

#design= ~ Platform + Site + RIN + Tissue + Subtype

dds_rinsite = DESeqDataSetFromMatrix(countData = rCountData_rinsite, colData = FullPheno_sr, design= ~ Platform + Site + RIN + Tissue + Subtype, tidy=F) #Subtype must be last (DESeq2 vignette)
dds_rinsite$Subtype = relevel(dds_rinsite$Subtype,ref = "CTR")
dseq_rinsite = DESeq(dds_rinsite,betaPrior=T)

setwd("G:/SpinalCord/Publication/DifferentialExpression/4Covar")

#Pairwise "contrast()"
glia.res = results(dseq_rinsite,contrast = c("Subtype","GLIA","CTR"))
write.csv(glia.res,"BothPlatform_AllSubjects_AllGenes_SubtypeSigGenes_10-18-23_GliaControl.csv")
#glia.sig = glia.res[! is.na(glia.res$padj) & glia.res$padj<0.05,]
filt.glia.sig = glia.res[rownames(glia.res) %in% Transcripts,]
write.csv(filt.glia.sig,"BothPlatform_AllSubjects_Top5000_SubtypeSigGenes_10-18-23_GliaControl.csv")

ox.res = results(dseq_rinsite,contrast = c("Subtype","OX","CTR"))
write.csv(ox.res,"BothPlatform_AllSubjects_AllGenes_SubtypeSigGenes_10-18-23_OxControl.csv")
#ox.sig = ox.res[! is.na(ox.res$padj) & ox.res$padj<0.05,]
filt.ox.sig = ox.res[rownames(ox.res) %in% Transcripts,]
write.csv(filt.ox.sig,"BothPlatform_AllSubjects_Top5000_SubtypeSigGenes_10-18-23_OxControl.csv")

TE.res = results(dseq_rinsite,contrast = c("Subtype","TD","CTR"))
write.csv(TE.res,"BothPlatform_AllSubjects_AllGenes_SubtypeSigGenes_10-18-23_TEControl.csv")
#TE.sig = TE.res[! is.na(TE.res$padj) & TE.res$padj<0.05,]
filt.TE.sig = TE.res[rownames(TE.res) %in% Transcripts,]
write.csv(filt.TE.sig,"BothPlatform_AllSubjects_Top5000_SubtypeSigGenes_10-18-23_TEControl.csv")

GT.res = results(dseq_rinsite,contrast = c("Subtype","GLIA","TD"))
write.csv(GT.res,"BothPlatform_AllSubjects_AllGenes_SubtypeSigGenes_10-18-23_GliaTE.csv")
#GT.sig = GT.res[! is.na(GT.res$padj) & GT.res$padj<0.05,]
filt.GT.sig = GT.res[rownames(GT.res) %in% Transcripts,]
write.csv(filt.GT.sig,"BothPlatform_AllSubjects_Top5000_SubtypeSigGenes_10-11-23_GliaTE.csv")

GO.res = results(dseq_rinsite,contrast = c("Subtype","GLIA","OX"))
write.csv(GO.res,"BothPlatform_AllSubjects_AllGenes_SubtypeSigGenes_10-18-23_GliaOX.csv")
#GO.sig = GO.res[! is.na(GO.res$padj) & GO.res$padj<0.05,]
filt.GO.sig = GO.res[rownames(GO.res) %in% Transcripts,]
write.csv(filt.GO.sig,"BothPlatform_AllSubjects_Top5000_SubtypeSigGenes_10-18-23_GliaOX.csv")

TO.res = results(dseq_rinsite,contrast = c("Subtype","TD","OX"))
write.csv(TO.res,"BothPlatform_AllSubjects_AllGenes_SubtypeSigGenes_10-18-23_TEOX.csv")
#TO.sig = TO.res[! is.na(TO.res$padj) & TO.res$padj<0.05,]
filt.TO.sig = TO.res[rownames(TO.res) %in% Transcripts,]
write.csv(filt.TO.sig,"BothPlatform_AllSubjects_Top5000_SubtypeSigGenes_10-18-23_TEOX.csv")

#Control Check - all non-neurological
table(FullPheno_sr$Subject.Group[which(FullPheno_sr$Subtype == "CTR")])

#All transcripts recovered
#Transcripts[! Transcripts %in% rownames(filt.glia.sig)]

filt.glia.sig2 = glia.res[rownames(glia.res) %in% TranscriptsPlus,]
filt.ox.sig2 = ox.res[rownames(ox.res) %in% TranscriptsPlus,]
filt.TE.sig2 = TE.res[rownames(TE.res) %in% TranscriptsPlus,]
filt.GT.sig2 = GT.res[rownames(GT.res) %in% TranscriptsPlus,]
filt.GO.sig2 = GO.res[rownames(GO.res) %in% TranscriptsPlus,]
filt.TO.sig2 = TO.res[rownames(TO.res) %in% TranscriptsPlus,]

write.csv(filt.glia.sig2,"BothPlatform_AllSubjects_Top5000_wCortexFeatures_SubtypeSigGenes_10-18-23_GliaControl.csv")
write.csv(filt.ox.sig2,"BothPlatform_AllSubjects_Top5000_wCortexFeatures_SubtypeSigGenes_10-18-23_OxControl.csv")
write.csv(filt.TE.sig2,"BothPlatform_AllSubjects_Top5000_wCortexFeatures_SubtypeSigGenes_10-18-23_TDControl.csv")
write.csv(filt.GT.sig2,"BothPlatform_AllSubjects_Top5000_wCortexFeatures_SubtypeSigGenes_10-18-23_GliaTD.csv")
write.csv(filt.GO.sig2,"BothPlatform_AllSubjects_Top5000_wCortexFeatures_SubtypeSigGenes_10-18-23_GliaOx.csv")
write.csv(filt.TO.sig2,"BothPlatform_AllSubjects_Top5000_wCortexFeatures_SubtypeSigGenes_10-18-23_TDOx.csv")

################################  DESeq2 Normalization  #########################################################################

tmp = estimateSizeFactors(dds_rinsite)
DESeq_NormalizedCounts_60k = counts(tmp,normalized = T)
DESeq_NormalizedCounts_60k_filtered = DESeq_NormalizedCounts_60k[rownames(DESeq_NormalizedCounts_60k) %in% TranscriptsPlus,]
NormCounts = DESeq_NormalizedCounts_60k_filtered

GliaIDs = FullPheno_sr$ExternalSampleId[which(FullPheno_sr$Subtype == "GLIA")]
OxIDs = FullPheno_sr$ExternalSampleId[which(FullPheno_sr$Subtype == "OX")]
TDIDs = FullPheno_sr$ExternalSampleId[which(FullPheno_sr$Subtype == "TD")]
CTRIDs = FullPheno_sr$ExternalSampleId[which(FullPheno_sr$Subtype == "CTR")]
samps = c(GliaIDs,OxIDs,TDIDs, CTRIDs)

#save.image("SpinalCord_DifferentialExpression_FULLCohort_Top5000_10-18-23.RData")
load("G:/SpinalCord/Publication/DifferentialExpression/4Covar/SpinalCord_DifferentialExpression_FULLCohort_Top5000_10-18-23.RData")


############################### HEATMAPs ################################################################

# #With the original cortex features
# HeatExpr = NormCounts[rownames(NormCounts) %in% my36,]
# zHeatExpr = zscore(HeatExpr,featuresarerows = T,removeNA = T,addlabels = T)
# MatrixtoHeatmap3(zHeatExpr,samplesarecols = T,featuresasrows = T,limits=c(-4,4),customfeats = rev(my36),customsamps = samps,sepindex = c(length(GliaIDs),length(GliaIDs)+length(OxIDs),length(GliaIDs)+length(OxIDs)+length(TDIDs)))
# 
# #Glia Supplemental
# HeatExpr = NormCounts[rownames(NormCounts) %in% S8,]
# zHeatExpr = zscore(HeatExpr,featuresarerows = T,removeNA = T,addlabels = T)
# MatrixtoHeatmap3(zHeatExpr,samplesarecols = T,featuresasrows = T,limits=c(-4,4),customfeats = rev(S8),customsamps = samps,sepindex = c(length(GliaIDs),length(GliaIDs)+length(OxIDs),length(GliaIDs)+length(OxIDs)+length(TDIDs)))
# 
# #Ox Supplemental
# HeatExpr = NormCounts[rownames(NormCounts) %in% S9,]
# zHeatExpr = zscore(HeatExpr,featuresarerows = T,removeNA = T,addlabels = T)
# MatrixtoHeatmap3(zHeatExpr,samplesarecols = T,featuresasrows = T,limits=c(-4,4),customfeats = rev(S9),customsamps = samps,sepindex = c(length(GliaIDs),length(GliaIDs)+length(OxIDs),length(GliaIDs)+length(OxIDs)+length(TDIDs)))
# 
# #TD Supplemental
# HeatExpr = NormCounts[rownames(NormCounts) %in% S10,]
# zHeatExpr = zscore(HeatExpr,featuresarerows = T,removeNA = T,addlabels = T)
# MatrixtoHeatmap3(zHeatExpr,samplesarecols = T,featuresasrows = T,limits=c(-4,4),customfeats = rev(S10),customsamps = samps,sepindex = c(length(GliaIDs),length(GliaIDs)+length(OxIDs),length(GliaIDs)+length(OxIDs)+length(TDIDs)))


##########Figure 4A

#Remove tentative features with fewer than 10 average counts - QC
# GliaFeatures = c("AHNAK","FKBP7","IMPA2","TMEM173","SLC10A6","SLC22A3","CHEK2","KCTD14","DLG2","B4GALNT2","CAPN3","CCL26","CCR4","C7","CCR9","CD59","CD300E","CFH","CLEC6A","IQGAP1","PPIC","CYFIP2","GLT8D2","IL1R1","IL20RB","LAMC1","MAP7","MYL9","RND3","SAMHD1","SELP","ST6GALNAC2","TAGLN","TYR","USP53","VIM")
# OxFeatures = c("ADAMTS20","ATXN8OS","B4GALT6","CPNE4","DDC","EFNA5","ENSG00000270953","ENSG00000286214","FGF13","GABRA1","GABRG3","GAD2","GALNT14","GLRA3","GRIA1","GRIN2A","GRM1","GRM7","HTR2A","KCNAB1","KCNC2","KCNH6","KCNS2","NMNAT2","NMS","NTS","PCLO","PCSK1","PKLR","RIMS2","SCN3A","SLC17A6","SLC35F4","SYN2","SYT1","TMEM35A","UBQLN2","UNC13C","UNC79","UNC80","UNCX","VSTM5")
# TDFeatures = c("APOBR","APOC1","CYP4X1","ENSG00000229771","ENSG00000237153","ENSG00000259601","ENSG00000261606","ENSG00000279149","ENSG00000280087","GBX1","GJB6","GTP2IP7","HMGB3P7","LINC01058","LINC01956","LINC02026","NKX6-2","NCF2","PTK6","POM121L7P","PTBP1P","SHH","SIGLEC15","TCF23")
# feats = c(GliaFeatures,OxFeatures,TDFeatures)
# HeatExpr = rCountData_rinsite[rownames(rCountData_rinsite) %in% feats,]
# rownames(HeatExpr)[which(rowMeans(HeatExpr,na.rm = T) < 10)]
# rownames(HeatExpr)[which(rowMeans(HeatExpr,na.rm = T) > 10)]


#Pub - Top 20 
GliaFeatures = c("AHNAK","B4GALNT2","CAPN3","CD59","CD300E","CFH","CHEK2","FKBP7","GLT8D2","IL1R1","IQGAP1","KCTD14","LAMC1","MAP7","PPIC","SAMHD1","SELP","TMEM173","USP53","VIM")
OxFeatures = c("CPNE4","EFNA5","FGF13","GABRG3","GALNT14","GRIA1","GRIN2A","GRM1","KCNH6","KCNS2","NMNAT2","NTS","PCLO","PKLR","RIMS2","SCN3A","SLC35F4","SYN2","SYT1","UNC13C")
TDFeatures = c("APOBR","APOC1","ENSG00000259601","ENSG00000185332","ENSG00000250608","ENSG00000275620","ENSG00000280087","ENSG00000285492","FPR3","LINC01091","ENSG00000259953","ENSG00000232310","M1AP","ENSG00000216802","NCF2","NLRP12","OTOAP1","ENSG00000234232","SLC28A1","TCF23")
feats = c(GliaFeatures,OxFeatures,TDFeatures)
HeatExpr = NormCounts[rownames(NormCounts) %in% feats,]
zHeatExpr = zscore(HeatExpr,featuresarerows = T,removeNA = T,addlabels = T)
#For plotting purposes adjust z-score values outside the range
zHeatExpr[which(zHeatExpr > 4)] = 4
zHeatExpr[which(zHeatExpr < -4)] = -4
#rownames(zHeatExpr)[which(rownames(zHeatExpr) == "ENSG00000259601")] == "DDX18P2"
MatrixtoHeatmap3(zHeatExpr,samplesarecols = T,featuresasrows = T,limits=c(-4,4),customfeats = rev(feats),customsamps = samps,sepindex = c(length(GliaIDs),length(GliaIDs)+length(OxIDs),length(GliaIDs)+length(OxIDs)+length(TDIDs)))

#Adjusted p-value matrix 

Adjpmat = data.frame(matrix(NA,nrow=length(feats),ncol=6))
rownames(Adjpmat) = feats
colnames(Adjpmat) = c("Glia_v_Ox","Glia_v_TD","Ox_v_TD","Glia_v_HC","Ox_v_HC","TD_v_HC")

for(i in 1:nrow(Adjpmat)){
  
  #Glia vs Ox
  GOind = which(rownames(filt.GO.sig2) == rownames(Adjpmat)[i])
  if(length(GOind)>0){
    Adjpmat$Glia_v_Ox[i] = filt.GO.sig2$padj[GOind]
  }else{
    Adjpmat$Glia_v_Ox[i] = NA
  }
  
  #Glia vs TD
  GTind = which(rownames(filt.GT.sig2) == rownames(Adjpmat)[i])
  if(length(GTind)>0){
    Adjpmat$Glia_v_TD[i] = filt.GT.sig2$padj[GTind]
  }else{
    Adjpmat$Glia_v_TD[i] = NA
  }
  
  #Ox vs TD
  OTind = which(rownames(filt.TO.sig2) == rownames(Adjpmat)[i])
  if(length(OTind)>0){
    Adjpmat$Ox_v_TD[i] = filt.TO.sig2$padj[OTind]
  }else{
    Adjpmat$Ox_v_TD[i] = NA
  }
  
  #Glia vs Control
  GCind = which(rownames(filt.glia.sig2) == rownames(Adjpmat)[i])
  if(length(GCind)>0){
    Adjpmat$Glia_v_HC[i] = filt.glia.sig2$padj[GCind]
  }else{
    Adjpmat$Glia_v_HC[i] = NA
  }
  
  #Ox vs Control
  OCind = which(rownames(filt.ox.sig2) == rownames(Adjpmat)[i])
  if(length(OCind)>0){
    Adjpmat$Ox_v_HC[i] = filt.ox.sig2$padj[OCind]
  }else{
    Adjpmat$Ox_v_HC[i] = NA
  }
  
  #TD vs Control
  TCind = which(rownames(filt.TE.sig2) == rownames(Adjpmat)[i])
  if(length(TCind)>0){
    Adjpmat$TD_v_HC[i] = filt.TE.sig2$padj[TCind]
  }else{
    Adjpmat$TD_v_HC[i] = NA
  }
  
}

#Pairwise differential expression - subtypes
MatrixtoHeatmap3(-log10(Adjpmat[,1:3]),samplesarecols = T,customfeats = rev(rownames(Adjpmat)),customsamps = colnames(Adjpmat)[1:3],colors = c("#2967d9","white","#cc8f02"),limits=c(-log10(0.05),max(-log10(Adjpmat[,1:3]))))
MatrixtoHeatmap3(-log10(Adjpmat[,1:3]),samplesarecols = T,customfeats = rev(rownames(Adjpmat)),customsamps = colnames(Adjpmat)[1:3],colors = c("#2967d9","white","#6941ba"),limits=c(-log10(0.05),max(-log10(Adjpmat[,1:3]))))

#Pairwise differential expression - controls
MatrixtoHeatmap3(-log10(Adjpmat[,4:6]),samplesarecols = T,customfeats = rev(rownames(Adjpmat)),customsamps = colnames(Adjpmat)[4:6],colors = c("#2967d9","white","#6941ba"),limits=c(-log10(0.05),max(-log10(Adjpmat[,4:6]))))
MatrixtoHeatmap3(-log10(Adjpmat[,4:6]),samplesarecols = T,customfeats = rev(rownames(Adjpmat)),customsamps = colnames(Adjpmat)[4:6],colors = c("#2967d9","white","#cc8f02"),limits=c(-log10(0.05),max(-log10(Adjpmat[,4:6]))))

###########################################################################

#Check covariate dependency on assigned subtype - Fig. S3
table(FullPheno_sr$Tissue)

#Tissue
table(FullPheno_sr$Subtype[which(FullPheno_sr$Tissue == "Spinal_Cord_Cervical")])
table(FullPheno_sr$Subtype[which(FullPheno_sr$Tissue == "Spinal_Cord_Thoracic")])
table(FullPheno_sr$Subtype[which(FullPheno_sr$Tissue == "Spinal_Cord_Lumbar")])

#Site
table(FullPheno_sr$Subtype[which(FullPheno_sr$Site == "NYGC")])
table(FullPheno_sr$Subtype[which(FullPheno_sr$Site == "TargetALS")])
#colors: #3c9955,#d49531

#Sex
table(FullPheno_sr$Subtype[which(FullPheno_sr$Sex == "Male")])
table(FullPheno_sr$Subtype[which(FullPheno_sr$Sex == "Female")])
#colors: skyblue1, pink3

#Above built in excel 

#RIN
par(mfrow=c(1,3))
hist(FullPheno_sr$RIN[which(FullPheno_sr$Subtype == "GLIA")],main = "ALS-Glia RIN",col="goldenrod1",ylim = c(0,25))
hist(FullPheno_sr$RIN[which(FullPheno_sr$Subtype == "TD")],main = "ALS-Glia RIN",col="firebrick",ylim = c(0,60))
hist(FullPheno_sr$RIN[which(FullPheno_sr$Subtype == "OX")],main = "ALS-Glia RIN",col="navy",ylim = c(0,40))

############################### VIOLINs ################################################################

#For plotting purposes (not used during statistical analysis), adjust zero count genes to 1

for(i in 1:nrow(NormCounts)){
  for(j in 1:ncol(NormCounts)){
    if(NormCounts[i,j] == 0){
      NormCounts[i,j] = 1
    }
  }
}

##Parse for plotting

#Subtype Normalized Count Matrix
GliaI = which(FullPheno_sr$Subtype == "GLIA")
GliaSamp = colnames(NormCounts)[GliaI]
OxI = which(FullPheno_sr$Subtype == "OX")
OxSamp = colnames(NormCounts)[OxI]
TEI = which(FullPheno_sr$Subtype == "TE")
TESamp = colnames(NormCounts)[TEI]
ONDI = which(FullPheno_sr$Subtype == "OND")
ONDSamp = colnames(NormCounts)[ONDI]
HCI = which(FullPheno_sr$Subtype == "Control")
HCSamp = colnames(NormCounts)[HCI]

################################### CNS LEVEL PLOTS #####################################################################


#Clean up naming

Subtype = FullPheno_sr$Subtype

for(i in 1:length(Subtype)){
  if(Subtype[i] == "GLIA"){
    Subtype[i] = "ALS-Glia"
  }else if(Subtype[i] == "OX"){
    Subtype[i] = "ALS-Ox"
  }else if(Subtype[i] == "TD"){
    Subtype[i] = "ALS-TD"
  }else{
    Subtype[i] = "Control"
  }
}

Subtype = factor(Subtype,levels = c("Control","ALS-Glia","ALS-Ox","ALS-TD"))


### AutoPlot
library(ggplot2)

#Violin at 3 levels, Cortex, Spinal Cord, Spinal Cord with Cortex patient labels
TDfeats = c("NKX6-2","TARDBP")
OXfeats = c("B4GALT6","GABRA1","GAD2","GLRA3","HTR2A","PCSK1","SLC17A6","UBQLN2")
Gfeats = c("MYL9","ST6GALNAC2","TAGLN")

#Standard Violins
#CortexCheck = c(my36,"CD68","TARDBP","IFI30")
#CortexSupp = c(S8,S9,S10)

PlotGene = c("CCR4","CR1","TTR","FCGR3B","CD300E","AGTR1","ANG","ANGPTL1","LILRA5","BDKRB2","SPIB","SERPINA5","CCR2","")
#Pick the reference subtype to add DE p-values contrasted with other groups 
#Options: Glia, Ox, TD, HC, All (case sensitive; note that 'All' option automatically selects the reference group... all comparisons are NOT performed) - this code does not support FTLD as the reference level
Focus = "Glia"

#Grab the index
for(i in 1:length(PlotGene)){
  tmpindex = which(rownames(NormCounts) == PlotGene[i])
  PlotCounts = NormCounts[tmpindex,]
  LogPlotCounts = log2(PlotCounts)
  logdat = data.frame(Subtype,LogPlotCounts)
  p = ggplot(logdat,aes(x=Subtype,y=LogPlotCounts,fill=Subtype)) + geom_boxplot()
  ttl = paste(PlotGene[i])
  p = p +ggtitle(ttl) + xlab("") + ylab("log2 Median-of-Ratios Counts")
  p = p+scale_x_discrete(limits=c("Control","ALS-Glia","ALS-Ox","ALS-TD"))
  p = p+scale_fill_manual(values = c("gray50","goldenrod1","navy","firebrick"))
  p = p+theme(axis.text = element_text(size=20), axis.title = element_text(size=20),plot.title = element_text(size=24))
  p = p+theme(legend.position = "none")
  p = p+theme(axis.title.x=element_text(vjust=-2))
  p = p+theme(axis.title.y=element_text(angle=90, vjust=6))
  p = p+ theme(plot.margin = unit(c(1,1,1,1), "cm"))
  p = p+theme(plot.title = element_text(hjust = 0.5))
  upperlim = max(LogPlotCounts)+max(LogPlotCounts)*0.65
  tmplim = min(LogPlotCounts)-min(LogPlotCounts)*0.15
  if(round(tmplim,0)==0){
    lowerlim = 0
  }else{
    lowerlim = tmplim
  }
  #Add p-values
  if(Focus == "Glia"){
    Controlp = filt.glia.sig #vs Control
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=2,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=1,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=2,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=1.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.15)
    }
    
    OXp = filt.GO.sig #vs OX
    tmpindex4 = which(rownames(OXp) == PlotGene[i])
    OXp = OXp$padj[tmpindex4]
    OXp = formatC(OXp,format = "e",digits = 2)
    if(OXp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.25),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.2,xend=2,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.2,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
      pstat = paste("FDR p-value:",OXp)
      p = p+annotate("text",label=pstat,x=2.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.30)
    }
    
    TEp = filt.GT.sig #vs TE
    tmpindex5 = which(rownames(TEp) == PlotGene[i])
    TEp = TEp$padj[tmpindex5]
    TEp = formatC(TEp,format = "e",digits = 2)
    if(TEp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.4,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.4),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.35,xend=2,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.35,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
      pstat = paste("FDR p-value:",TEp)
      p = p+annotate("text",label=pstat,x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.45)
    }
    
  }else if(Focus == "Ox"){
    Controlp = filt.ox.sig #vs Control
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=1,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.15)
    }
    
    OXp = filt.GO.sig #vs OX
    tmpindex4 = which(rownames(OXp) == PlotGene[i])
    OXp = OXp$padj[tmpindex4]
    OXp = formatC(OXp,format = "e",digits = 2)
    if(OXp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.25),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.2,xend=2,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.2,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
      pstat = paste("FDR p-value:",OXp)
      p = p+annotate("text",label=pstat,x=2.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.30)
    }
    
    TEp = filt.TO.sig #vs TE
    tmpindex5 = which(rownames(TEp) == PlotGene[i])
    TEp = TEp$padj[tmpindex5]
    TEp = formatC(TEp,format = "e",digits = 2)
    if(TEp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.4,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.4),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.35,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.35,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
      pstat = paste("FDR p-value:",TEp)
      p = p+annotate("text",label=pstat,x=3.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.45)
    }
    
  }else if(Focus == "TD"){
    Controlp = filt.TE.sig#vs Control
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=1,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=2.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.15)
    }
    
    Gliap = filt.GT.sig #vs Glia
    tmpindex4 = which(rownames(Gliap) == PlotGene[i])
    Gliap = Gliap$padj[tmpindex4]
    Gliap = formatC(Gliap,format = "e",digits = 2)
    if(Gliap>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.25),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.20,xend=2,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.20,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
      pstat = paste("FDR p-value:",Gliap)
      p = p+annotate("text",label=pstat,x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.30)
    }
    
    TEp = filt.TO.sig #vs OX
    tmpindex5 = which(rownames(TEp) == PlotGene[i])
    TEp = TEp$padj[tmpindex5]
    TEp = formatC(TEp,format = "e",digits = 2)
    if(TEp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.4,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.4),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.35,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.35,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
      pstat = paste("FDR p-value:",TEp)
      p = p+annotate("text",label=pstat,x=3.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.5)
    }
    
  }else if(Focus == "HC"){
    Controlp = filt.glia.sig #vs OND
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=2,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=1,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=2,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=1.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.15)
    }
    
    OXp = filt.ox.sig #vs OX
    tmpindex4 = which(rownames(OXp) == PlotGene[i])
    OXp = OXp$padj[tmpindex4]
    OXp = formatC(OXp,format = "e",digits = 2)
    if(OXp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.25),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.2,xend=1,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.2,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
      pstat = paste("FDR p-value:",OXp)
      p = p+annotate("text",label=pstat,x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.3)
    }
    
    TEp = filt.TE.sig #vs TE
    tmpindex5 = which(rownames(TEp) == PlotGene[i])
    TEp = TEp$padj[tmpindex5]
    TEp = formatC(TEp,format = "e",digits = 2)
    if(TEp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.4,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.4),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.35,xend=1,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.35,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
      pstat = paste("FDR p-value:",TEp)
      p = p+annotate("text",label=pstat,x=2.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.45)
    }
    
  }else{
    cat("Focus not specified, p-values not added to the plot")
  }
  print(p)
}




#Autoplot - Boxplot and Violin on log2 and norm count scales
for(i in 1:length(PlotGene)){
  
  #Grab the index
  tmpindex = which(rownames(NormCounts) == PlotGene[i])
  
  ##Boxplots
  
  #On the Normalized Count Scale
  PlotCounts = NormCounts[tmpindex,]
  dat = data.frame(Subtype,PlotCounts)
  p = ggplot(dat,aes(x=Subtype,y=PlotCounts,fill=Subtype)) + geom_boxplot()
  ttl = paste(PlotGene[i])
  p = p +ggtitle(ttl) + ylab("Median-of-Ratios Counts") + xlab("")
  p = p+scale_x_discrete(limits=c("Control","ALS-Glia","ALS-Ox","ALS-TD"))
  p = p+scale_fill_manual(values = c("gray50","goldenrod1","navy","firebrick"))
  p = p+theme(axis.text = element_text(size=20), axis.title = element_text(size=20),plot.title = element_text(size=24))
  p = p+theme(legend.position = "none")
  p = p+theme(axis.title.x=element_text(vjust=-2))
  p = p+theme(axis.title.y=element_text(angle=90, vjust=6))
  p = p+ theme(plot.margin = unit(c(1,1,1,1), "cm"))
  p = p+theme(plot.title = element_text(hjust = 0.5,size = 40))
  upperlim = max(PlotCounts)+max(PlotCounts)*0.65
  tmplim = min(PlotCounts)-min(PlotCounts)*0.15
  if(round(tmplim,0)==0){
    lowerlim = 0
  }else{
    lowerlim = tmplim
  }
  #Add p-values
  if(Focus == "Glia"){
    Controlp = filt.glia.sig #vs Control
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=1,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=1.5,y=max(PlotCounts)+max(PlotCounts)*0.15)
    }
    
    OXp = filt.GO.sig #vs OX
    tmpindex4 = which(rownames(OXp) == PlotGene[i])
    OXp = OXp$padj[tmpindex4]
    OXp = formatC(OXp,format = "e",digits = 2)
    if(OXp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.25,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.25),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.2,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.2,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      pstat = paste("FDR p-value:",OXp)
      p = p+annotate("text",label=pstat,x=2.5,y=max(PlotCounts)+max(PlotCounts)*0.30)
    }
    
    TEp = filt.GT.sig #vs TE
    tmpindex5 = which(rownames(TEp) == PlotGene[i])
    TEp = TEp$padj[tmpindex5]
    TEp = formatC(TEp,format = "e",digits = 2)
    if(TEp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.4,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.4),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      pstat = paste("FDR p-value:",TEp)
      p = p+annotate("text",label=pstat,x=3,y=max(PlotCounts)+max(PlotCounts)*0.45)
    }
    
  }else if(Focus == "Ox"){
    Controlp = filt.ox.sig #vs Control
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=1,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=2,y=max(PlotCounts)+max(PlotCounts)*0.15)
    }
    
    OXp = filt.GO.sig #vs OX
    tmpindex4 = which(rownames(OXp) == PlotGene[i])
    OXp = OXp$padj[tmpindex4]
    OXp = formatC(OXp,format = "e",digits = 2)
    if(OXp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.25,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.25),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.2,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.2,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      pstat = paste("FDR p-value:",OXp)
      p = p+annotate("text",label=pstat,x=2.5,y=max(PlotCounts)+max(PlotCounts)*0.30)
    }
    
    TEp = filt.TO.sig #vs TE
    tmpindex5 = which(rownames(TEp) == PlotGene[i])
    TEp = TEp$padj[tmpindex5]
    TEp = formatC(TEp,format = "e",digits = 2)
    if(TEp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.4,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.4),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      pstat = paste("FDR p-value:",TEp)
      p = p+annotate("text",label=pstat,x=3.5,y=max(PlotCounts)+max(PlotCounts)*0.45)
    }
    
  }else if(Focus == "TD"){
    Controlp = filt.TE.sig#vs Control
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=1,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=2.5,y=max(PlotCounts)+max(PlotCounts)*0.15)
    }
    
    Gliap = filt.GT.sig #vs Glia
    tmpindex4 = which(rownames(Gliap) == PlotGene[i])
    Gliap = Gliap$padj[tmpindex4]
    Gliap = formatC(Gliap,format = "e",digits = 2)
    if(Gliap>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.25,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.25),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.20,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(PlotCounts)+max(PlotCounts)*0.20,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      pstat = paste("FDR p-value:",Gliap)
      p = p+annotate("text",label=pstat,x=3,y=max(PlotCounts)+max(PlotCounts)*0.30)
    }
    
    TEp = filt.TO.sig #vs OX
    tmpindex5 = which(rownames(TEp) == PlotGene[i])
    TEp = TEp$padj[tmpindex5]
    TEp = formatC(TEp,format = "e",digits = 2)
    if(TEp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.4,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.4),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      pstat = paste("FDR p-value:",TEp)
      p = p+annotate("text",label=pstat,x=3.5,y=max(PlotCounts)+max(PlotCounts)*0.5)
    }
    
    }else if(Focus == "HC"){
      Controlp = filt.glia.sig #vs OND
      tmpindex2 = which(rownames(Controlp) == PlotGene[i])
      Controlp = Controlp$padj[tmpindex2]
      Controlp = formatC(Controlp,format = "e",digits = 2)
      if(Controlp>""){
        p = p+ylim(lowerlim,upperlim)
        p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.1),size=0.8)
        p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=1,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
        p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
        pstat = paste("FDR p-value:",Controlp)
        p = p+annotate("text",label=pstat,x=1.5,y=max(PlotCounts)+max(PlotCounts)*0.15)
      }
      
      OXp = filt.ox.sig #vs OX
      tmpindex4 = which(rownames(OXp) == PlotGene[i])
      OXp = OXp$padj[tmpindex4]
      OXp = formatC(OXp,format = "e",digits = 2)
      if(OXp>""){
        p = p+ylim(lowerlim,upperlim)
        p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.25,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.25),size=0.8)
        p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.2,xend=1,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
        p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.2,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
        pstat = paste("FDR p-value:",OXp)
        p = p+annotate("text",label=pstat,x=2,y=max(PlotCounts)+max(PlotCounts)*0.3)
      }
      
      TEp = filt.TE.sig #vs TE
      tmpindex5 = which(rownames(TEp) == PlotGene[i])
      TEp = TEp$padj[tmpindex5]
      TEp = formatC(TEp,format = "e",digits = 2)
      if(TEp>""){
        p = p+ylim(lowerlim,upperlim)
        p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.4,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.4),size=0.8)
        p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=1,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
        p = p+geom_segment(aes(x=4,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
        pstat = paste("FDR p-value:",TEp)
        p = p+annotate("text",label=pstat,x=2.5,y=max(PlotCounts)+max(PlotCounts)*0.45)
      }
    
  }else{
    cat("Focus not specified, p-values not added to the plot")
  }
  print(p)
  
  #On the Log2 Normalized Count Scale
  LogPlotCounts = log2(PlotCounts)
  logdat = data.frame(Subtype,LogPlotCounts)
  p = ggplot(logdat,aes(x=Subtype,y=LogPlotCounts,fill=Subtype)) + geom_boxplot()
  ttl = paste(PlotGene[i])
  p = p +ggtitle(ttl) + xlab("") + ylab("log2 Median-of-Ratios Counts")
  p = p+scale_x_discrete(limits=c("Control","ALS-Glia","ALS-Ox","ALS-TD"))
  p = p+scale_fill_manual(values = c("gray50","goldenrod1","navy","firebrick"))
  p = p+theme(axis.text = element_text(size=20), axis.title = element_text(size=20),plot.title = element_text(size=24))
  p = p+theme(legend.position = "none")
  p = p+theme(axis.title.x=element_text(vjust=-2))
  p = p+theme(axis.title.y=element_text(angle=90, vjust=6))
  p = p+ theme(plot.margin = unit(c(1,1,1,1), "cm"))
  p = p+theme(plot.title = element_text(hjust = 0.5))
  upperlim = max(LogPlotCounts)+max(LogPlotCounts)*0.65
  tmplim = min(LogPlotCounts)-min(LogPlotCounts)*0.15
  if(round(tmplim,0)==0){
    lowerlim = 0
  }else{
    lowerlim = tmplim
  }
  #Add p-values
  if(Focus == "Glia"){
    Controlp = filt.glia.sig #vs Control
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=1,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=1.5,y=max(PlotCounts)+max(PlotCounts)*0.15)
    }
    
    OXp = filt.GO.sig #vs OX
    tmpindex4 = which(rownames(OXp) == PlotGene[i])
    OXp = OXp$padj[tmpindex4]
    OXp = formatC(OXp,format = "e",digits = 2)
    if(OXp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.25,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.25),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.2,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.2,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      pstat = paste("FDR p-value:",OXp)
      p = p+annotate("text",label=pstat,x=2.5,y=max(PlotCounts)+max(PlotCounts)*0.30)
    }
    
    TEp = filt.GT.sig #vs TE
    tmpindex5 = which(rownames(TEp) == PlotGene[i])
    TEp = TEp$padj[tmpindex5]
    TEp = formatC(TEp,format = "e",digits = 2)
    if(TEp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.4,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.4),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      pstat = paste("FDR p-value:",TEp)
      p = p+annotate("text",label=pstat,x=3,y=max(PlotCounts)+max(PlotCounts)*0.45)
    }
    
  }else if(Focus == "Ox"){
    Controlp = filt.ox.sig #vs Control
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=1,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=2,y=max(PlotCounts)+max(PlotCounts)*0.15)
    }
    
    OXp = filt.GO.sig #vs OX
    tmpindex4 = which(rownames(OXp) == PlotGene[i])
    OXp = OXp$padj[tmpindex4]
    OXp = formatC(OXp,format = "e",digits = 2)
    if(OXp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.25,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.25),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.2,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.2,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      pstat = paste("FDR p-value:",OXp)
      p = p+annotate("text",label=pstat,x=2.5,y=max(PlotCounts)+max(PlotCounts)*0.30)
    }
    
    TEp = filt.TO.sig #vs TE
    tmpindex5 = which(rownames(TEp) == PlotGene[i])
    TEp = TEp$padj[tmpindex5]
    TEp = formatC(TEp,format = "e",digits = 2)
    if(TEp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.4,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.4),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      pstat = paste("FDR p-value:",TEp)
      p = p+annotate("text",label=pstat,x=3.5,y=max(PlotCounts)+max(PlotCounts)*0.45)
    }
    
  }else if(Focus == "TD"){
    Controlp = filt.TE.sig#vs Control
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=1,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=2.5,y=max(PlotCounts)+max(PlotCounts)*0.15)
    }
    
    Gliap = filt.GT.sig #vs Glia
    tmpindex4 = which(rownames(Gliap) == PlotGene[i])
    Gliap = Gliap$padj[tmpindex4]
    Gliap = formatC(Gliap,format = "e",digits = 2)
    if(Gliap>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.25,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.25),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.20,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(PlotCounts)+max(PlotCounts)*0.20,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      pstat = paste("FDR p-value:",Gliap)
      p = p+annotate("text",label=pstat,x=3,y=max(PlotCounts)+max(PlotCounts)*0.30)
    }
    
    TEp = filt.TO.sig #vs OX
    tmpindex5 = which(rownames(TEp) == PlotGene[i])
    TEp = TEp$padj[tmpindex5]
    TEp = formatC(TEp,format = "e",digits = 2)
    if(TEp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.4,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.4),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      pstat = paste("FDR p-value:",TEp)
      p = p+annotate("text",label=pstat,x=3.5,y=max(PlotCounts)+max(PlotCounts)*0.5)
    }
    
  }else if(Focus == "HC"){
    Controlp = filt.glia.sig #vs OND
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=1,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=1.5,y=max(PlotCounts)+max(PlotCounts)*0.15)
    }
    
    OXp = filt.ox.sig #vs OX
    tmpindex4 = which(rownames(OXp) == PlotGene[i])
    OXp = OXp$padj[tmpindex4]
    OXp = formatC(OXp,format = "e",digits = 2)
    if(OXp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.25,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.25),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.2,xend=1,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.2,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      pstat = paste("FDR p-value:",OXp)
      p = p+annotate("text",label=pstat,x=2,y=max(PlotCounts)+max(PlotCounts)*0.3)
    }
    
    TEp = filt.TE.sig #vs TE
    tmpindex5 = which(rownames(TEp) == PlotGene[i])
    TEp = TEp$padj[tmpindex5]
    TEp = formatC(TEp,format = "e",digits = 2)
    if(TEp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.4,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.4),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=1,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      pstat = paste("FDR p-value:",TEp)
      p = p+annotate("text",label=pstat,x=2.5,y=max(PlotCounts)+max(PlotCounts)*0.45)
    }
      
   }else{
    cat("Focus not specified, p-values not added to the plot")
  }
  print(p)
  
  ##Violin
  
  #Norm Count Scale
  p = ggplot(dat,aes(x=Subtype,y=PlotCounts,fill=Subtype)) + geom_violin()
  ttl = paste(PlotGene[i])
  p = p +ggtitle(ttl) + xlab("") + ylab("Median-of-Ratios Counts")
  p = p+scale_x_discrete(limits=c("Control","ALS-Glia","ALS-Ox","ALS-TD"))
  p = p+scale_fill_manual(values = c("gray50","goldenrod1","navy","firebrick"))
  p = p+theme(axis.text = element_text(size=20), axis.title = element_text(size=20),plot.title = element_text(size=24))
  p = p+theme(legend.position = "none")
  p = p+theme(panel.background = element_rect(fill = "gray90",colour = "gray80",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "gray80"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray80"))
  p = p+geom_dotplot(binaxis = 'y',stackdir = 'center',dotsize = 0.25,fill="white",stackratio = 0.75)
  p = p+theme(axis.title.x=element_text(vjust=-2))
  p = p+theme(axis.title.y=element_text(angle=90, vjust=6))
  p = p+ theme(plot.margin = unit(c(1,1,1,1), "cm"))
  p = p+theme(plot.title = element_text(hjust = 0.5))
  upperlim = max(PlotCounts)+max(PlotCounts)*0.65
  tmplim = min(PlotCounts)-min(PlotCounts)*0.15
  if(round(tmplim,0)==0){
    lowerlim = 0
  }else{
    lowerlim = tmplim
  }
  #Add p-values
  if(Focus == "Glia"){
    Controlp = filt.glia.sig #vs Control
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=1,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=1.5,y=max(PlotCounts)+max(PlotCounts)*0.15)
    }
    
    OXp = filt.GO.sig #vs OX
    tmpindex4 = which(rownames(OXp) == PlotGene[i])
    OXp = OXp$padj[tmpindex4]
    OXp = formatC(OXp,format = "e",digits = 2)
    if(OXp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.25,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.25),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.2,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.2,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      pstat = paste("FDR p-value:",OXp)
      p = p+annotate("text",label=pstat,x=2.5,y=max(PlotCounts)+max(PlotCounts)*0.30)
    }
    
    TEp = filt.GT.sig #vs TE
    tmpindex5 = which(rownames(TEp) == PlotGene[i])
    TEp = TEp$padj[tmpindex5]
    TEp = formatC(TEp,format = "e",digits = 2)
    if(TEp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.4,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.4),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      pstat = paste("FDR p-value:",TEp)
      p = p+annotate("text",label=pstat,x=3,y=max(PlotCounts)+max(PlotCounts)*0.45)
    }
    
  }else if(Focus == "Ox"){
    Controlp = filt.ox.sig #vs Control
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=1,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=2,y=max(PlotCounts)+max(PlotCounts)*0.15)
    }
    
    OXp = filt.GO.sig #vs OX
    tmpindex4 = which(rownames(OXp) == PlotGene[i])
    OXp = OXp$padj[tmpindex4]
    OXp = formatC(OXp,format = "e",digits = 2)
    if(OXp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.25,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.25),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.2,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.2,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      pstat = paste("FDR p-value:",OXp)
      p = p+annotate("text",label=pstat,x=2.5,y=max(PlotCounts)+max(PlotCounts)*0.30)
    }
    
    TEp = filt.TO.sig #vs TE
    tmpindex5 = which(rownames(TEp) == PlotGene[i])
    TEp = TEp$padj[tmpindex5]
    TEp = formatC(TEp,format = "e",digits = 2)
    if(TEp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.4,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.4),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      pstat = paste("FDR p-value:",TEp)
      p = p+annotate("text",label=pstat,x=3.5,y=max(PlotCounts)+max(PlotCounts)*0.45)
    }
    
  }else if(Focus == "TD"){
    Controlp = filt.TE.sig#vs Control
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=1,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=2.5,y=max(PlotCounts)+max(PlotCounts)*0.15)
    }
    
    Gliap = filt.GT.sig #vs Glia
    tmpindex4 = which(rownames(Gliap) == PlotGene[i])
    Gliap = Gliap$padj[tmpindex4]
    Gliap = formatC(Gliap,format = "e",digits = 2)
    if(Gliap>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.25,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.25),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.20,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(PlotCounts)+max(PlotCounts)*0.20,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      pstat = paste("FDR p-value:",Gliap)
      p = p+annotate("text",label=pstat,x=3,y=max(PlotCounts)+max(PlotCounts)*0.30)
    }
    
    TEp = filt.TO.sig #vs OX
    tmpindex5 = which(rownames(TEp) == PlotGene[i])
    TEp = TEp$padj[tmpindex5]
    TEp = formatC(TEp,format = "e",digits = 2)
    if(TEp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.4,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.4),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      pstat = paste("FDR p-value:",TEp)
      p = p+annotate("text",label=pstat,x=3.5,y=max(PlotCounts)+max(PlotCounts)*0.5)
    }
    
  }else if(Focus == "HC"){
    Controlp = filt.glia.sig #vs OND
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=1,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=1.5,y=max(PlotCounts)+max(PlotCounts)*0.15)
    }
    
    OXp = filt.ox.sig #vs OX
    tmpindex4 = which(rownames(OXp) == PlotGene[i])
    OXp = OXp$padj[tmpindex4]
    OXp = formatC(OXp,format = "e",digits = 2)
    if(OXp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.25,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.25),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.2,xend=1,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.2,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      pstat = paste("FDR p-value:",OXp)
      p = p+annotate("text",label=pstat,x=2,y=max(PlotCounts)+max(PlotCounts)*0.3)
    }
    
    TEp = filt.TE.sig #vs TE
    tmpindex5 = which(rownames(TEp) == PlotGene[i])
    TEp = TEp$padj[tmpindex5]
    TEp = formatC(TEp,format = "e",digits = 2)
    if(TEp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.4,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.4),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=1,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      pstat = paste("FDR p-value:",TEp)
      p = p+annotate("text",label=pstat,x=2.5,y=max(PlotCounts)+max(PlotCounts)*0.45)
    }
    
  }else{
    cat("Focus not specified, p-values not added to the plot")
  }
  print(p)
  
  #log2 scale
  p = ggplot(logdat,aes(x=Subtype,y=LogPlotCounts,fill=Subtype)) + geom_violin()
  ttl = paste(PlotGene[i])
  p = p +ggtitle(ttl) + xlab("") + ylab("log2 Median-of-Ratios Counts")
  p = p+scale_x_discrete(limits=c("Control","ALS-Glia","ALS-Ox","ALS-TD"))
  p = p+scale_fill_manual(values = c("gray50","goldenrod1","navy","firebrick"))
  p = p+theme(axis.text = element_text(size=20), axis.title = element_text(size=20),plot.title = element_text(size=24))
  p = p+theme(legend.position = "none")
  p = p+theme(panel.background = element_rect(fill = "gray90",colour = "gray80",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "gray80"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray80"))
  p = p+geom_dotplot(binaxis = 'y',stackdir = 'center',dotsize = 0.25,fill="white",stackratio = 0.75)
  p = p+theme(axis.title.x=element_text(vjust=-2))
  p = p+theme(axis.title.y=element_text(angle=90, vjust=6))
  p = p+theme(plot.margin = unit(c(1,1,1,1), "cm"))
  p = p+theme(plot.title = element_text(hjust = 0.5))
  upperlim = max(LogPlotCounts)+max(LogPlotCounts)*0.65
  tmplim = min(LogPlotCounts)-min(LogPlotCounts)*0.15
  if(round(tmplim,0)==0){
    lowerlim = 0
  }else{
    lowerlim = tmplim
  }
  #Add p-values
  if(Focus == "Glia"){
    Controlp = filt.glia.sig #vs Control
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=1,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=1.5,y=max(PlotCounts)+max(PlotCounts)*0.15)
    }
    
    OXp = filt.GO.sig #vs OX
    tmpindex4 = which(rownames(OXp) == PlotGene[i])
    OXp = OXp$padj[tmpindex4]
    OXp = formatC(OXp,format = "e",digits = 2)
    if(OXp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.25,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.25),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.2,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.2,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      pstat = paste("FDR p-value:",OXp)
      p = p+annotate("text",label=pstat,x=2.5,y=max(PlotCounts)+max(PlotCounts)*0.30)
    }
    
    TEp = filt.GT.sig #vs TE
    tmpindex5 = which(rownames(TEp) == PlotGene[i])
    TEp = TEp$padj[tmpindex5]
    TEp = formatC(TEp,format = "e",digits = 2)
    if(TEp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.4,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.4),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      pstat = paste("FDR p-value:",TEp)
      p = p+annotate("text",label=pstat,x=3,y=max(PlotCounts)+max(PlotCounts)*0.45)
    }
    
  }else if(Focus == "Ox"){
    Controlp = filt.ox.sig #vs Control
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=1,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=2,y=max(PlotCounts)+max(PlotCounts)*0.15)
    }
    
    OXp = filt.GO.sig #vs OX
    tmpindex4 = which(rownames(OXp) == PlotGene[i])
    OXp = OXp$padj[tmpindex4]
    OXp = formatC(OXp,format = "e",digits = 2)
    if(OXp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.25,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.25),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.2,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.2,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      pstat = paste("FDR p-value:",OXp)
      p = p+annotate("text",label=pstat,x=2.5,y=max(PlotCounts)+max(PlotCounts)*0.30)
    }
    
    TEp = filt.TO.sig #vs TE
    tmpindex5 = which(rownames(TEp) == PlotGene[i])
    TEp = TEp$padj[tmpindex5]
    TEp = formatC(TEp,format = "e",digits = 2)
    if(TEp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.4,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.4),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      pstat = paste("FDR p-value:",TEp)
      p = p+annotate("text",label=pstat,x=3.5,y=max(PlotCounts)+max(PlotCounts)*0.45)
    }
    
  }else if(Focus == "TD"){
    Controlp = filt.TE.sig#vs Control
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=1,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=2.5,y=max(PlotCounts)+max(PlotCounts)*0.15)
    }
    
    Gliap = filt.GT.sig #vs Glia
    tmpindex4 = which(rownames(Gliap) == PlotGene[i])
    Gliap = Gliap$padj[tmpindex4]
    Gliap = formatC(Gliap,format = "e",digits = 2)
    if(Gliap>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.25,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.25),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.20,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(PlotCounts)+max(PlotCounts)*0.20,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      pstat = paste("FDR p-value:",Gliap)
      p = p+annotate("text",label=pstat,x=3,y=max(PlotCounts)+max(PlotCounts)*0.30)
    }
    
    TEp = filt.TO.sig #vs OX
    tmpindex5 = which(rownames(TEp) == PlotGene[i])
    TEp = TEp$padj[tmpindex5]
    TEp = formatC(TEp,format = "e",digits = 2)
    if(TEp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.4,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.4),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      pstat = paste("FDR p-value:",TEp)
      p = p+annotate("text",label=pstat,x=3.5,y=max(PlotCounts)+max(PlotCounts)*0.5)
    }
    
  }else if(Focus == "HC"){
    Controlp = filt.glia.sig #vs OND
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=1,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(PlotCounts)+max(PlotCounts)*0.1,xend=2,yend=max(PlotCounts)+max(PlotCounts)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=1.5,y=max(PlotCounts)+max(PlotCounts)*0.15)
    }
    
    OXp = filt.ox.sig #vs OX
    tmpindex4 = which(rownames(OXp) == PlotGene[i])
    OXp = OXp$padj[tmpindex4]
    OXp = formatC(OXp,format = "e",digits = 2)
    if(OXp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.25,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.25),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.2,xend=1,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(PlotCounts)+max(PlotCounts)*0.2,xend=3,yend=max(PlotCounts)+max(PlotCounts)*0.2),size=0.8)
      pstat = paste("FDR p-value:",OXp)
      p = p+annotate("text",label=pstat,x=2,y=max(PlotCounts)+max(PlotCounts)*0.3)
    }
    
    TEp = filt.TE.sig #vs TE
    tmpindex5 = which(rownames(TEp) == PlotGene[i])
    TEp = TEp$padj[tmpindex5]
    TEp = formatC(TEp,format = "e",digits = 2)
    if(TEp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.4,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.4),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=1,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(PlotCounts)+max(PlotCounts)*0.35,xend=4,yend=max(PlotCounts)+max(PlotCounts)*0.35),size=0.8)
      pstat = paste("FDR p-value:",TEp)
      p = p+annotate("text",label=pstat,x=2.5,y=max(PlotCounts)+max(PlotCounts)*0.45)
    }
    
    
  }else{
    cat("Focus not specified, p-values not added to the plot")
  }
  print(p)
}

#Warnings can generally be ignored for this loop - there are some duplicated lines of code


#Plots for Publication
for(i in 1:length(PlotGene)){
  
  #Grab the index
  tmpindex = which(rownames(NormCounts) == PlotGene[i])
  
  #On the Normalized Count Scale
  PlotCounts = NormCounts[tmpindex,]
  dat = data.frame(Subtype,PlotCounts)
  
  LogPlotCounts = log2(PlotCounts)
  logdat = data.frame(Subtype,LogPlotCounts)
  
  #log2 scale
  p = ggplot(logdat,aes(x=Subtype,y=LogPlotCounts,fill=Subtype)) + geom_violin(kernel = "gaussian",scale="width")
  #p = p + stat_ydensity(adjust = 6,na.rm = T,scale = "count",bw=1)
  ttl = paste(PlotGene[i])
  p = p +ggtitle(ttl) + xlab("") + ylab("log2 Median-of-Ratios Counts")
  p = p+scale_x_discrete(limits=c("Control","FTLD","ALS-Glia","ALS-Ox","ALS-TD"))
  p = p+scale_fill_manual(values = c("gray50","gray30","goldenrod1","navy","firebrick"))
  p = p+theme(axis.text = element_text(size=20), axis.title = element_text(size=20),plot.title = element_text(size=24))
  p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "gray80"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray80"))
  p = p+theme(legend.position = "none")
  p = p+geom_dotplot(binaxis = 'y',stackdir = 'center',dotsize = 0.25,fill="white",stackratio = 0.75)
  p = p+theme(axis.title.x=element_text(vjust=-2))
  p = p+theme(axis.title.y=element_text(angle=90, vjust=4,size=24))
  p = p+theme(plot.margin = unit(c(1,1,1,1), "cm"))
  p = p+theme(axis.text.x = element_text(size = 30))
  p = p+theme(axis.text.y = element_text(size = 30))
  p = p+theme(plot.title = element_text(hjust = 0.5,size = 40))
  upperlim = max(LogPlotCounts)+max(LogPlotCounts)*0.65
  tmplim = min(LogPlotCounts)-min(LogPlotCounts)*0.15
  if(round(tmplim,0)==0){
    lowerlim = 0
  }else{
    lowerlim = tmplim
  }
  #Add p-values
  if(Focus == "Glia"){
    Controlp = filt.glia.sig #vs Control
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=1,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.15,size = 6)
    }
    
    ONDp = filt.glia.sig.ond #vs OND
    tmpindex3 = which(rownames(ONDp) == PlotGene[i])
    ONDp = ONDp$padj[tmpindex3]
    ONDp = formatC(ONDp,format = "e",digits = 2)
    if(ONDp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.25),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=2,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
      pstat = paste("FDR p-value:",ONDp)
      p = p+annotate("text",label=pstat,x=2.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.3,size = 6)
    }
    
    OXp = filt.GO.sig #vs OX
    tmpindex4 = which(rownames(OXp) == PlotGene[i])
    OXp = OXp$padj[tmpindex4]
    OXp = formatC(OXp,format = "e",digits = 2)
    if(OXp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.40,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.40),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.40,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.40,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
      pstat = paste("FDR p-value:",OXp)
      p = p+annotate("text",label=pstat,x=3.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.45,size = 6)
    }
    
    TEp = filt.GT.sig #vs TE
    tmpindex5 = which(rownames(TEp) == PlotGene[i])
    TEp = TEp$padj[tmpindex5]
    TEp = formatC(TEp,format = "e",digits = 2)
    if(TEp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.55),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.50),size=0.8)
      p = p+geom_segment(aes(x=5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.50),size=0.8)
      pstat = paste("FDR p-value:",TEp)
      p = p+annotate("text",label=pstat,x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.6,size = 6)
    }
    
  }else if(Focus == "Ox"){
    Controlp = filt.ox.sig #vs Control
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=1,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=2.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.15,size = 6)
    }
    
    ONDp = filt.ox.sig.ond #vs OND
    tmpindex3 = which(rownames(ONDp) == PlotGene[i])
    ONDp = ONDp$padj[tmpindex3]
    ONDp = formatC(ONDp,format = "e",digits = 2)
    if(ONDp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.25),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=2,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
      pstat = paste("FDR p-value:",ONDp)
      p = p+annotate("text",label=pstat,x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.3,size = 6)
    }
    
    OXp = filt.GO.sig #vs OX
    tmpindex4 = which(rownames(OXp) == PlotGene[i])
    OXp = OXp$padj[tmpindex4]
    OXp = formatC(OXp,format = "e",digits = 2)
    if(OXp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.40,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.40),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.40,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.40,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
      pstat = paste("FDR p-value:",OXp)
      p = p+annotate("text",label=pstat,x=3.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.45,size = 6)
    }
    
    TEp = filt.TO.sig #vs TE
    tmpindex5 = which(rownames(TEp) == PlotGene[i])
    TEp = TEp$padj[tmpindex5]
    TEp = formatC(TEp,format = "e",digits = 2)
    if(TEp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.55),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.50),size=0.8)
      p = p+geom_segment(aes(x=5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.50),size=0.8)
      pstat = paste("FDR p-value:",TEp)
      p = p+annotate("text",label=pstat,x=4.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.6,size = 6)
    }
    
  }else if(Focus == "TD"){
    Controlp = filt.TE.sig#vs Control
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=1,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
      p = p+geom_segment(aes(x=5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.15,size = 6)
    }
    
    ONDp = filt.TE.sig.ond #vs OND
    tmpindex3 = which(rownames(ONDp) == PlotGene[i])
    ONDp = ONDp$padj[tmpindex3]
    ONDp = formatC(ONDp,format = "e",digits = 2)
    if(ONDp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.25),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=2,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
      p = p+geom_segment(aes(x=5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
      pstat = paste("FDR p-value:",ONDp)
      p = p+annotate("text",label=pstat,x=3.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.3,size = 6)
    }
    
    Gliap = filt.GT.sig #vs Glia
    tmpindex4 = which(rownames(Gliap) == PlotGene[i])
    Gliap = Gliap$padj[tmpindex4]
    Gliap = formatC(Gliap,format = "e",digits = 2)
    if(Gliap>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.40,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.40),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.40,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
      p = p+geom_segment(aes(x=5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.40,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
      pstat = paste("FDR p-value:",Gliap)
      p = p+annotate("text",label=pstat,x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.45,size = 6)
    }
    
    TEp = filt.TO.sig #vs OX
    tmpindex5 = which(rownames(TEp) == PlotGene[i])
    TEp = TEp$padj[tmpindex5]
    TEp = formatC(TEp,format = "e",digits = 2)
    if(TEp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.55),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.50),size=0.8)
      p = p+geom_segment(aes(x=5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.50),size=0.8)
      pstat = paste("FDR p-value:",TEp)
      p = p+annotate("text",label=pstat,x=4.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.6,size = 6)
    }
    
  }else if(Focus == "HC"){
    Controlp = filt.COND.sig #vs OND
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=2,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=1,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=2,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=1.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.15,size = 6)
    }
    
    ONDp = filt.glia.sig #vs Glia
    tmpindex3 = which(rownames(ONDp) == PlotGene[i])
    ONDp = ONDp$padj[tmpindex3]
    ONDp = formatC(ONDp,format = "e",digits = 2)
    if(ONDp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.25),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=1,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
      pstat = paste("FDR p-value:",ONDp)
      p = p+annotate("text",label=pstat,x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.3,size = 6)
    }
    
    OXp = filt.ox.sig #vs OX
    tmpindex4 = which(rownames(OXp) == PlotGene[i])
    OXp = OXp$padj[tmpindex4]
    OXp = formatC(OXp,format = "e",digits = 2)
    if(OXp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.40,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.40),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.40,xend=1,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.40,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
      pstat = paste("FDR p-value:",OXp)
      p = p+annotate("text",label=pstat,x=2.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.45,size = 6)
    }
    
    TEp = filt.TE.sig #vs TE
    tmpindex5 = which(rownames(TEp) == PlotGene[i])
    TEp = TEp$padj[tmpindex5]
    TEp = formatC(TEp,format = "e",digits = 2)
    if(TEp>""){
      p = p+ylim(lowerlim,upperlim)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.55),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=1,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.50),size=0.8)
      p = p+geom_segment(aes(x=5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.50),size=0.8)
      pstat = paste("FDR p-value:",TEp)
      p = p+annotate("text",label=pstat,x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.6,size = 6)
    }
    
  }else if(Focus == "All"){
    
    #Reset axis limits
    upperlim = max(LogPlotCounts)+max(LogPlotCounts)*0.65
    tmplim = min(LogPlotCounts)-min(LogPlotCounts)*0.15
    if(round(tmplim,0)==0){
      lowerlim = 0
    }else{
      lowerlim = tmplim
    }
    
    #Get Subtype with max counts as reference level - all comparisons is too busy
    GCounts = NormCounts[,GliaI]
    OCounts = NormCounts[,OxI]
    TCounts = NormCounts[,TEI]
    HCounts = NormCounts[,HCI]
    ONDCounts = NormCounts[,ONDI]
    
    Gavg = mean(GCounts[tmpindex,])
    Oavg = mean(OCounts[tmpindex,])
    Tavg = mean(TCounts[tmpindex,])
    Havg = mean(HCounts[tmpindex,])
    ONDavg = mean(ONDCounts[tmpindex,])
    
    reflev = max(c(Gavg,Oavg,Tavg,Havg,ONDavg))
    
    if(Gavg == reflev){
      
      Controlp = filt.glia.sig #vs Control
      tmpindex2 = which(rownames(Controlp) == PlotGene[i])
      Controlp = Controlp$padj[tmpindex2]
      Controlp = formatC(Controlp,format = "e",digits = 2)
      if(Controlp>""){
        p = p+ylim(lowerlim,upperlim)
        p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.1),size=0.8)
        p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=1,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
        p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
        pstat = paste("FDR p-value:",Controlp)
        p = p+annotate("text",label=pstat,x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.15)
      }
      
      ONDp = filt.glia.sig.ond #vs OND
      tmpindex3 = which(rownames(ONDp) == PlotGene[i])
      ONDp = ONDp$padj[tmpindex3]
      ONDp = formatC(ONDp,format = "e",digits = 2)
      if(ONDp>""){
        p = p+ylim(lowerlim,upperlim)
        p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.25),size=0.8)
        p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=2,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
        p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
        pstat = paste("FDR p-value:",ONDp)
        p = p+annotate("text",label=pstat,x=2.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.3)
      }
      
      OXp = filt.GO.sig #vs OX
      tmpindex4 = which(rownames(OXp) == PlotGene[i])
      OXp = OXp$padj[tmpindex4]
      OXp = formatC(OXp,format = "e",digits = 2)
      if(OXp>""){
        p = p+ylim(lowerlim,upperlim)
        p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.40,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.40),size=0.8)
        p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.40,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
        p = p+geom_segment(aes(x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.40,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
        pstat = paste("FDR p-value:",OXp)
        p = p+annotate("text",label=pstat,x=3.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.45)
      }
      
      TEp = filt.GT.sig #vs TE
      tmpindex5 = which(rownames(TEp) == PlotGene[i])
      TEp = TEp$padj[tmpindex5]
      TEp = formatC(TEp,format = "e",digits = 2)
      if(TEp>""){
        p = p+ylim(lowerlim,upperlim)
        p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.55),size=0.8)
        p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.50),size=0.8)
        p = p+geom_segment(aes(x=5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.50),size=0.8)
        pstat = paste("FDR p-value:",TEp)
        p = p+annotate("text",label=pstat,x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.6)
      }
      
      
    }else if(Oavg == reflev){
      
      Controlp = filt.ox.sig #vs Control
      tmpindex2 = which(rownames(Controlp) == PlotGene[i])
      Controlp = Controlp$padj[tmpindex2]
      Controlp = formatC(Controlp,format = "e",digits = 2)
      if(Controlp>""){
        p = p+ylim(lowerlim,upperlim)
        p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.1),size=0.8)
        p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=1,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
        p = p+geom_segment(aes(x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
        pstat = paste("FDR p-value:",Controlp)
        p = p+annotate("text",label=pstat,x=2.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.15)
      }
      
      ONDp = filt.ox.sig.ond #vs OND
      tmpindex3 = which(rownames(ONDp) == PlotGene[i])
      ONDp = ONDp$padj[tmpindex3]
      ONDp = formatC(ONDp,format = "e",digits = 2)
      if(ONDp>""){
        p = p+ylim(lowerlim,upperlim)
        p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.25),size=0.8)
        p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=2,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
        p = p+geom_segment(aes(x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
        pstat = paste("FDR p-value:",ONDp)
        p = p+annotate("text",label=pstat,x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.3)
      }
      
      Gliap = filt.GO.sig #vs Glia
      tmpindex4 = which(rownames(Gliap) == PlotGene[i])
      Gliap = Gliap$padj[tmpindex4]
      Gliap = formatC(Gliap,format = "e",digits = 2)
      if(Gliap>""){
        p = p+ylim(lowerlim,upperlim)
        p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.40,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.40),size=0.8)
        p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.40,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
        p = p+geom_segment(aes(x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.40,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
        pstat = paste("FDR p-value:",Gliap)
        p = p+annotate("text",label=pstat,x=3.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.45)
      }
      
      TEp = filt.TO.sig #vs TE
      tmpindex5 = which(rownames(TEp) == PlotGene[i])
      TEp = TEp$padj[tmpindex5]
      TEp = formatC(TEp,format = "e",digits = 2)
      if(TEp>""){
        p = p+ylim(lowerlim,upperlim)
        p = p+geom_segment(aes(x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.55),size=0.8)
        p = p+geom_segment(aes(x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.50),size=0.8)
        p = p+geom_segment(aes(x=5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.50),size=0.8)
        pstat = paste("FDR p-value:",TEp)
        p = p+annotate("text",label=pstat,x=4.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.6)
      }
      
    }else if(Tavg == reflev){
      
      Controlp = filt.TE.sig #vs Control
      tmpindex2 = which(rownames(Controlp) == PlotGene[i])
      Controlp = Controlp$padj[tmpindex2]
      Controlp = formatC(Controlp,format = "e",digits = 2)
      if(Controlp>""){
        p = p+ylim(lowerlim,upperlim)
        p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.1),size=0.8)
        p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=1,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
        p = p+geom_segment(aes(x=5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
        pstat = paste("FDR p-value:",Controlp)
        p = p+annotate("text",label=pstat,x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.15)
      }
      
      ONDp = filt.TE.sig.ond #vs OND
      tmpindex3 = which(rownames(ONDp) == PlotGene[i])
      ONDp = ONDp$padj[tmpindex3]
      ONDp = formatC(ONDp,format = "e",digits = 2)
      if(ONDp>""){
        p = p+ylim(lowerlim,upperlim)
        p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.25),size=0.8)
        p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=2,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
        p = p+geom_segment(aes(x=5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
        pstat = paste("FDR p-value:",ONDp)
        p = p+annotate("text",label=pstat,x=3.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.3)
      }
      
      Gliap = filt.GT.sig #vs Glia
      tmpindex4 = which(rownames(Gliap) == PlotGene[i])
      Gliap = Gliap$padj[tmpindex4]
      Gliap = formatC(Gliap,format = "e",digits = 2)
      if(Gliap>""){
        p = p+ylim(lowerlim,upperlim)
        p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.40,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.40),size=0.8)
        p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.40,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
        p = p+geom_segment(aes(x=5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.40,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
        pstat = paste("FDR p-value:",Gliap)
        p = p+annotate("text",label=pstat,x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.45)
      }
      
      TEp = filt.TO.sig #vs OX
      tmpindex5 = which(rownames(TEp) == PlotGene[i])
      TEp = TEp$padj[tmpindex5]
      TEp = formatC(TEp,format = "e",digits = 2)
      if(TEp>""){
        p = p+ylim(lowerlim,upperlim)
        p = p+geom_segment(aes(x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.55),size=0.8)
        p = p+geom_segment(aes(x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.50),size=0.8)
        p = p+geom_segment(aes(x=5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.50),size=0.8)
        pstat = paste("FDR p-value:",TEp)
        p = p+annotate("text",label=pstat,x=4.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.6)
      }
      
    }else if(Havg == reflev){
      
      ONDp = filt.COND.sig #vs OND
      tmpindex3 = which(rownames(ONDp) == PlotGene[i])
      ONDp = ONDp$padj[tmpindex3]
      ONDp = formatC(ONDp,format = "e",digits = 2)
      if(ONDp>""){
        p = p+ylim(lowerlim,upperlim)
        p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=2,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.1),size=0.8)
        p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=1,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
        p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=2,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
        pstat = paste("FDR p-value:",ONDp)
        p = p+annotate("text",label=pstat,x=1.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.15)
      }
      
      Gliap = filt.glia.sig #vs Control
      tmpindex2 = which(rownames(Gliap) == PlotGene[i])
      Gliap = Gliap$padj[tmpindex2]
      Gliap = formatC(Gliap,format = "e",digits = 2)
      if(Gliap>""){
        p = p+ylim(lowerlim,upperlim)
        p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.25),size=0.8)
        p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=1,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
        p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
        pstat = paste("FDR p-value:",Gliap)
        p = p+annotate("text",label=pstat,x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.3)
      }
      
      OXp = filt.ox.sig #vs Control
      tmpindex2 = which(rownames(OXp) == PlotGene[i])
      OXp = OXp$padj[tmpindex2]
      OXp = formatC(OXp,format = "e",digits = 2)
      if(OXp>""){
        p = p+ylim(lowerlim,upperlim)
        p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.4,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.4),size=0.8)
        p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.4,xend=1,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
        p = p+geom_segment(aes(x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.4,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
        pstat = paste("FDR p-value:",OXp)
        p = p+annotate("text",label=pstat,x=2.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.45)
      }
      
      TEp = filt.TE.sig#vs Control
      tmpindex2 = which(rownames(TEp) == PlotGene[i])
      TEp = TEp$padj[tmpindex2]
      TEp = formatC(TEp,format = "e",digits = 2)
      if(TEp>""){
        p = p+ylim(lowerlim,upperlim)
        p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.55),size=0.8)
        p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=1,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.5),size=0.8)
        p = p+geom_segment(aes(x=5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.5),size=0.8)
        pstat = paste("FDR p-value:",TEp)
        p = p+annotate("text",label=pstat,x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.6)
      }
      
    }else if(ONDavg == reflev){
      
      Controlp = filt.COND.sig #vs OND
      tmpindex3 = which(rownames(Controlp) == PlotGene[i])
      Controlp = Controlp$padj[tmpindex3]
      Controlp = formatC(Controlp,format = "e",digits = 2)
      if(Controlp>""){
        p = p+ylim(lowerlim,upperlim)
        p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=2,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.1),size=0.8)
        p = p+geom_segment(aes(x=1,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=1,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
        p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.1,xend=2,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.05),size=0.8)
        pstat = paste("FDR p-value:",Controlp)
        p = p+annotate("text",label=pstat,x=1.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.15)
      }
      
      Gliap = filt.glia.sig.ond #vs Control
      tmpindex2 = which(rownames(Gliap) == PlotGene[i])
      Gliap = Gliap$padj[tmpindex2]
      Gliap = formatC(Gliap,format = "e",digits = 2)
      if(Gliap>""){
        p = p+ylim(lowerlim,upperlim)
        p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.25),size=0.8)
        p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=2,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
        p = p+geom_segment(aes(x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.25,xend=3,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.2),size=0.8)
        pstat = paste("FDR p-value:",Gliap)
        p = p+annotate("text",label=pstat,x=2.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.30)
      }
      
      OXp = filt.ox.sig.ond #vs Control
      tmpindex2 = which(rownames(OXp) == PlotGene[i])
      OXp = OXp$padj[tmpindex2]
      OXp = formatC(OXp,format = "e",digits = 2)
      if(OXp>""){
        p = p+ylim(lowerlim,upperlim)
        p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.4,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.4),size=0.8)
        p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.4,xend=2,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
        p = p+geom_segment(aes(x=4,y=max(LogPlotCounts)+max(LogPlotCounts)*0.4,xend=4,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.35),size=0.8)
        pstat = paste("FDR p-value:",OXp)
        p = p+annotate("text",label=pstat,x=3,y=max(LogPlotCounts)+max(LogPlotCounts)*0.45)
      }
      
      TEp = filt.TE.sig.ond #vs Control
      tmpindex2 = which(rownames(TEp) == PlotGene[i])
      TEp = TEp$padj[tmpindex2]
      TEp = formatC(TEp,format = "e",digits = 2)
      if(TEp>""){
        p = p+ylim(lowerlim,upperlim)
        p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.55),size=0.8)
        p = p+geom_segment(aes(x=2,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=2,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.5),size=0.8)
        p = p+geom_segment(aes(x=5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.55,xend=5,yend=max(LogPlotCounts)+max(LogPlotCounts)*0.5),size=0.8)
        pstat = paste("FDR p-value:",TEp)
        p = p+annotate("text",label=pstat,x=3.5,y=max(LogPlotCounts)+max(LogPlotCounts)*0.6)
      }
      
    }else(
      cat("Unknown Error")
    )
    
    
  }else{
    cat("Focus not specified, p-values not added to the plot")
  }
  print(p)
}



#######################################################################################################################
################################## SUPPLEMENTAL TE PLOTS ##############################################################
#######################################################################################################################

################ Subtype-specific Upregulated TE Expression #####################################################################

MyTEs = rownames(NormCounts)[which(nchar(rownames(NormCounts)) >20)]
SubtypeTE.glia = SubtypeTE.ox = SubtypeTE.te = SubtypeTE.ond = SubtypeTE.hc = rep(NA,length(MyTEs))
count1 = count2 = count3 = count4 = count5 = 1
for(i in 1:length(MyTEs)){
  
  ind = which(rownames(NormCounts) == MyTEs[i])
  GliaAvg = mean(NormCounts[ind,GliaI])
  OxAvg = mean(NormCounts[ind,OxI])
  TEAvg = mean(NormCounts[ind,TEI])
  ONDAvg = mean(NormCounts[ind,ONDI])
  HCAvg = mean(NormCounts[ind,HCI])
  
  pick = max(c(GliaAvg,OxAvg,TEAvg,ONDAvg,HCAvg))
  
  if(GliaAvg == pick){
    SubtypeTE.glia[count1] = strsplit(MyTEs[i],"\\|")[[1]][4]
    count1 = count1+1
  }else if(OxAvg == pick){
    SubtypeTE.ox[count2] = strsplit(MyTEs[i],"\\|")[[1]][4]
    count2 = count2+1
  }else if(TEAvg == pick){
    SubtypeTE.te[count3] = strsplit(MyTEs[i],"\\|")[[1]][4]
    count3 = count3+1
  }else if(ONDAvg == pick){
    SubtypeTE.ond[count4] = strsplit(MyTEs[i],"\\|")[[1]][4]
    count4 = count4+1
  }else{
    SubtypeTE.hc[count5] = strsplit(MyTEs[i],"\\|")[[1]][4]
    count5 = count5+1
  }
  
}

SubtypeTE.glia = SubtypeTE.glia[! is.na(SubtypeTE.glia)]
SubtypeTE.ox = SubtypeTE.ox[! is.na(SubtypeTE.ox)]
SubtypeTE.te = SubtypeTE.te[! is.na(SubtypeTE.te)]
SubtypeTE.ond = SubtypeTE.ond[! is.na(SubtypeTE.ond)]
SubtypeTE.hc = SubtypeTE.hc[! is.na(SubtypeTE.hc)]

Subfamily = names(table(c(names(table(SubtypeTE.glia)),names(table(SubtypeTE.ox)),names(table(SubtypeTE.te)),names(table(SubtypeTE.ond)),names(table(SubtypeTE.hc)))))
SubfamilyDF = matrix(NA,nrow=length(Subfamily),ncol=5)
rownames(SubfamilyDF) = Subfamily
colnames(SubfamilyDF) = c("ALS-TD","ALS-Ox","ALS-Glia","FTLD",'HC')

for(i in 1:nrow(SubfamilyDF)){
  
  TEcount = length(which(SubtypeTE.te == rownames(SubfamilyDF)[i]))
  Oxcount = length(which(SubtypeTE.ox == rownames(SubfamilyDF)[i]))
  Gliacount = length(which(SubtypeTE.glia == rownames(SubfamilyDF)[i]))
  ONDcount = length(which(SubtypeTE.ond == rownames(SubfamilyDF)[i]))
  HCcount = length(which(SubtypeTE.hc == rownames(SubfamilyDF)[i]))
  
  SubfamilyDF[i,1] = TEcount
  SubfamilyDF[i,2] = Oxcount
  SubfamilyDF[i,3] = Gliacount
  SubfamilyDF[i,4] = ONDcount
  SubfamilyDF[i,5] = HCcount
  
}

#Quick Plot
par(mar = c(6.1, 5.1, 4.1, 2.1)) 
mycols = rainbow(nrow(SubfamilyDF))
barplot(SubfamilyDF,ylim = c(0,300),ylab = "Number of Discriminatory TEs",xlab="Subtype",col = mycols,cex.axis = 1.2,cex.names = 1.5,cex.lab=1.5)
#legend(5.5,300,legend = rownames(SubfamilyDF),col = mycols,lty = 1,lwd=3,cex = 0.2)

#Legend without plot
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft",legend = rownames(SubfamilyDF),col = mycols,lty = 1,lwd=3,cex = 0.8,ncol=4)

#Source
# setwd("D:/Jarrett/Research/Nat Comm Reviews/Source")
# write.csv(SubfamilyDF,"FigS4A.csv")


#Pie chart version
par(mar = c(0, 1, 2, 1))
tmp1 = which(SubfamilyDF[,1] == 0)
SubfamilyDF2 = SubfamilyDF[-tmp1,]
mycols2 = mycols[-tmp1]
pct = SubfamilyDF2[,1]/sum(SubfamilyDF2[,1])
pie(SubfamilyDF2[,1],labels =  paste(rownames(SubfamilyDF2)," (",round(pct*100,1),"%)",sep=""),col=mycols2,main="ALS-TD Transposable Element Feature Breakdown",cex=0.5)

tmp1 = which(SubfamilyDF[,2] == 0)
SubfamilyDF2 = SubfamilyDF[-tmp1,]
mycols2 = mycols[-tmp1]
pct = SubfamilyDF2[,2]/sum(SubfamilyDF2[,2])
pie(SubfamilyDF2[,2],labels =  paste(rownames(SubfamilyDF2)," (",round(pct*100,1),"%)",sep=""),col=mycols2,main="ALS-Ox Transposable Element Feature Breakdown",cex=0.5)

tmp1 = which(SubfamilyDF[,3] == 0)
SubfamilyDF2 = SubfamilyDF[-tmp1,]
mycols2 = mycols[-tmp1]
pct = SubfamilyDF2[,3]/sum(SubfamilyDF2[,3])
pie(SubfamilyDF2[,3],labels =  paste(rownames(SubfamilyDF2)," (",round(pct*100,1),"%)",sep=""),col=mycols2,main="ALS-Glia Transposable Element Feature Breakdown",cex=1)


###################### Chromosome TE Pie Chart #####################################################################

a = table(substr(MyTEs,1,5))
for(i in 1:length(a)){
  if(substr(names(a)[i],5,5) == "|"){
    names(a)[i] = substr(names(a)[i],1,4)
  }
}
plot(a)

chroms = data.frame(matrix(NA,nrow=length(a),ncol = 2))
colnames(chroms) = c("chr","TEfeatures")
chroms$chr = paste("chr",seq(1,23),sep="")
chroms$chr[23] = "chrX"

for(i in 1:nrow(chroms)){
  for(j in 1:length(a)){
    if(chroms$chr[i] == names(a)[j]){
      chroms$TEfeatures[i] = a[[j]]
    }
    
  }
}

library(RColorBrewer)
mycols = colorRampPalette(c("cyan","dodgerblue4"))
#mycols = colorRampPalette(c("white","gray30"))
mycols = mycols(23)
chroms$pct = chroms$TEfeatures / sum(chroms$TEfeatures)
chroms$pct = chroms$pct * 100 
chroms$chr = factor(chroms$chr,levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"))
pie(chroms$TEfeatures,labels = paste(chroms$chr," (",round(chroms$pct,1),"%)",sep=""),col = mycols,main="chromosome vs number of discriminatory TE features")
