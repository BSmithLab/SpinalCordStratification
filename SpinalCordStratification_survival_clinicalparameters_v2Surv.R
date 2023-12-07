##Survival analysis for ALS subtypes defined using spinal cord transcriptomes

#Written By: Jarrett Eshima
#For: Dr. Barbara Smith Lab
#Date: October 5th 2023

library(survival)
library(survminer)
library(ranger)
library(ggplot2)
library(dplyr)

##########################################################################################################################################
################################    GSE153960 Spinal Cord Survival (n = 429)     #####################################################################
##########################################################################################################################################

wd = "G:/SpinalCord/Publication/Survival"
setwd(wd)
Pheno = read.csv("Table_S1.csv")

colnames(Pheno)[14:16] = c("AgeofOnset","AgeofDeath","Duration")

#Pie Charts
pie(table(Pheno$MetaCortex),col=c("gray50","goldenrod1","navy","firebrick"),init.angle = -10,main="Majority Agreement Cortex Labels")
pie(table(Pheno$MetaSpinal),col=c("gray50","goldenrod1","navy","firebrick"),init.angle = -20,main="Majority Agreement Spinal Labels")

FinalSurvivalMat = Pheno

#Clean up missing entries
for(i in 1:nrow(FinalSurvivalMat)){
  if(FinalSurvivalMat$AgeofOnset[i] == "Unknown"){
    FinalSurvivalMat$AgeofOnset[i] = NA
  }
  
  if(FinalSurvivalMat$Duration[i] == "Unknown"){
    FinalSurvivalMat$Duration[i] = NA
  }
}


FinalSurvival = FinalSurvivalMat
colnames(FinalSurvival)[16] = "time"
FinalSurvival$time = as.numeric(FinalSurvival$time)


## Majority assignment approach (matches cortex analysis)
km_fit = survfit(Surv(time) ~ MetaSpinal, data=FinalSurvival)
surv_pvalue(km_fit)
print(km_fit)
summary(km_fit, times = c(1,30,60,90*(1:10)))

par(mar = c(5, 5, 3, 2))
plot(km_fit,col = c("gray50","goldenrod","navy","firebrick"),main="ALS Subtype Kaplan-Meier Survival",xlab="Disease Duration (months)",ylab="Survival Probability",lwd=2.5,cex.axis = 1.5,cex.lab=1.5,cex.main=1.5,xaxt="n")
timeseq = seq(0,156,12)
axis(side = 1, at = timeseq,labels = T,tick = T,cex.axis = 1.5)
legend(130,1.02,legend=c("ALS-Glia","ALS-Ox","ALS-TD","Discordant"),lty = 1,lwd=5,col=c("goldenrod","navy","firebrick","gray50"),cex=1.4)
text(10,0.2,"p = 0.048",cex=1.5)

# segments(0,0.5,42,0.5,lty="dashed",lwd=2.5)
# segments(28,0,28,0.5,lty="dashed",lwd=2.5)
# segments(36,0,36,0.5,lty="dashed",lwd=2.5)
# segments(42,0,42,0.5,lty="dashed",lwd=2.5)

############ Considering only subjects with agreement between motor and frontal cortex samples ("Discordant" removed)
FinalSurvival_nodisc = FinalSurvival[-which(FinalSurvival$MetaSpinal == "Discordant"),]

km_fit = survfit(Surv(time) ~ MetaSpinal, data=FinalSurvival_nodisc)
surv_pvalue(km_fit)
print(km_fit)
summary(km_fit, times = c(1,30,60,90*(1:10)))
plot(km_fit,col = c("goldenrod","navy","firebrick"),main="ALS Patients Survival - No Discordant",xlab="Disease Duration (months)",ylab="Survival Probability",lwd=2.5,cex.axis = 1.5,cex.lab=1.5,cex.main=1.5,xaxt="n")
timeseq = seq(0,156,12)
axis(side = 1, at = timeseq,labels = T,tick = T,cex.axis = 1.5)
legend(138,1.02,legend=c("ALS-Glia","ALS-Ox","ALS-TD"),lty = 1,lwd=3.5,col=c("goldenrod","navy","firebrick"))
text(10,0.2,"p = 0.058",cex=1.5)

segments(0,0.5,35,0.5,lty="dashed",lwd=2.5)
segments(24,0,24,0.5,lty="dashed",lwd=2.5)
segments(35,0,35,0.5,lty="dashed",lwd=2.5)
segments(33,0,33,0.5,lty="dashed",lwd=2.5)

#############################################################################################################################

#Pairwise survival comparisons - WITH Discordant

tmp = which(FinalSurvival$MetaSpinal == "GLIA")
CleanSurvivalMat_OT = FinalSurvival[-tmp,]
tmp = which(CleanSurvivalMat_OT$MetaSpinal == "Discordant")
CleanSurvivalMat_OT = CleanSurvivalMat_OT[-tmp,]

tmp = which(FinalSurvival$MetaSpinal == "OX")
CleanSurvivalMat_TD = FinalSurvival[-tmp,]
tmp = which(CleanSurvivalMat_TD$MetaSpinal == "GLIA")
CleanSurvivalMat_TD = CleanSurvivalMat_TD[-tmp,]

tmp = which(FinalSurvival$MetaSpinal == "TD")
CleanSurvivalMat_GO = FinalSurvival[-tmp,]
tmp = which(CleanSurvivalMat_GO$MetaSpinal == "Discordant")
CleanSurvivalMat_GO = CleanSurvivalMat_GO[-tmp,]

tmp = which(FinalSurvival$MetaSpinal == "TD")
CleanSurvivalMat_GD = FinalSurvival[-tmp,]
tmp = which(CleanSurvivalMat_GD$MetaSpinal == "OX")
CleanSurvivalMat_GD = CleanSurvivalMat_GD[-tmp,]

tmp = which(FinalSurvival$MetaSpinal == "TD")
CleanSurvivalMat_OD = FinalSurvival[-tmp,]
tmp = which(CleanSurvivalMat_OD$MetaSpinal == "GLIA")
CleanSurvivalMat_OD = CleanSurvivalMat_OD[-tmp,]

tmp = which(FinalSurvival$MetaSpinal == "Discordant")
CleanSurvivalMat_GT = FinalSurvival[-tmp,]
tmp = which(CleanSurvivalMat_GT$MetaSpinal == "OX")
CleanSurvivalMat_GT = CleanSurvivalMat_GT[-tmp,]

#Glia vs Discordant - p = 0.032
km_fit = survfit(Surv(time) ~ MetaSpinal, data=CleanSurvivalMat_GD)
surv_pvalue(km_fit)
print(km_fit)

#Glia vs Ox - p = 0.023
km_fit = survfit(Surv(time) ~ MetaSpinal, data=CleanSurvivalMat_GO)
surv_pvalue(km_fit)
print(km_fit)

#Glia vs TD - p = 0.27
km_fit = survfit(Surv(time) ~ MetaSpinal, data=CleanSurvivalMat_GT)
surv_pvalue(km_fit)
print(km_fit)

#Ox vs TD - p = 0.12
km_fit = survfit(Surv(time) ~ MetaSpinal, data=CleanSurvivalMat_OT)
surv_pvalue(km_fit)
print(km_fit)

#Ox vs Discordant - p = 0.77
km_fit = survfit(Surv(time) ~ MetaSpinal, data=CleanSurvivalMat_OD)
surv_pvalue(km_fit)
print(km_fit)

#TD vs Discordant - p = 0.13
km_fit = survfit(Surv(time) ~ MetaSpinal, data=CleanSurvivalMat_TD)
surv_pvalue(km_fit)
print(km_fit)

#############################################################################################################################

#Pairwise survival comparisons - No Discordant

tmp = which(FinalSurvival_nodisc$MetaSpinal == "GLIA")
CleanSurvivalMat_NoGlia = FinalSurvival_nodisc[-tmp,]

tmp = which(FinalSurvival_nodisc$MetaSpinal == "OX")
CleanSurvivalMat_NoOx = FinalSurvival_nodisc[-tmp,]

tmp = which(FinalSurvival_nodisc$MetaSpinal == "TD")
CleanSurvivalMat_NoTD = FinalSurvival_nodisc[-tmp,]

#Ox vs TD without Discordant - p = 0.12
km_fit = survfit(Surv(time) ~ MetaSpinal, data=CleanSurvivalMat_NoGlia)
surv_pvalue(km_fit)
print(km_fit)


#Glia vs TD without Discordant - p = 0.27
km_fit = survfit(Surv(time) ~ MetaSpinal, data=CleanSurvivalMat_NoOx)
surv_pvalue(km_fit)
print(km_fit)


#Glia vs Ox without Discordant - p = 0.023
km_fit = survfit(Surv(time) ~ MetaSpinal, data=CleanSurvivalMat_NoTD)
surv_pvalue(km_fit)
print(km_fit)

############################################################################################################################################
################################ Define Summary Metrics for ALS Subtypes ###################################################################
############################################################################################################################################
#Fix missing data - confirmed with NYGC Metadata file

FinalSurvivalMat$Disease.Group[which(FinalSurvivalMat$Disease.Group == "#N/A")] = c("ALS-SOD1","ALS-SOD1","ALS-TDP")

#Clinical Parameters

G = which(FinalSurvivalMat$MetaSpinal == "GLIA")
O = which(FinalSurvivalMat$MetaSpinal == "OX")
TD = which(FinalSurvivalMat$MetaSpinal == "TD")
D = which(FinalSurvivalMat$MetaSpinal == "Discordant")
G = FinalSurvivalMat[G,]
O = FinalSurvivalMat[O,]
TD = FinalSurvivalMat[TD,]
D = FinalSurvivalMat[D,]

######################################### Subtype specific
#Age of onset - summary metrics
tmp = D$AgeofOnset
tmp = as.numeric(tmp)
tmp = tmp[!is.na(tmp)]
mean(tmp)
sd(tmp)/sqrt(length(tmp)) #standard error

#Age of death - summary metrics
tmp = D$AgeofDeath
tmp = as.numeric(tmp)
tmp = tmp[!is.na(tmp)]
mean(tmp)
sd(tmp)/sqrt(length(tmp)) #standard error

#Disease Duration - summary metrics
tmp = D$Duration
tmp = as.numeric(tmp)
tmp = tmp[!is.na(tmp)]
mean(tmp)
sd(tmp)/sqrt(length(tmp)) #standard error

#Site of onset - summary metrics
tmp = G$Site.of.Onset
table(tmp)
tmp = O$Site.of.Onset
table(tmp)
tmp = TD$Site.of.Onset
table(tmp)
tmp = D$Site.of.Onset
table(tmp)

#FTLD Comorbidity - summary metrics
#Site of onset
tmp = G$Disease.Group
table(tmp)
tmp = O$Disease.Group
table(tmp)
tmp = TD$Disease.Group
table(tmp)
tmp = D$Disease.Group
table(tmp)

######################################### Across all ALS Patients
#Age of onset
tmp = FinalSurvivalMat$AgeofOnset
tmp = as.numeric(tmp)
tmp = tmp[!is.na(tmp)]
mean(tmp)
sd(tmp)/sqrt(length(tmp)) #standard error 

#Age of death
tmp = FinalSurvivalMat$AgeofDeath
tmp = as.numeric(tmp)
tmp = tmp[!is.na(tmp)]
mean(tmp)
sd(tmp)/sqrt(length(tmp)) #standard error

#Disease Duration
tmp = FinalSurvivalMat$Duration
tmp = as.numeric(tmp)
tmp = tmp[!is.na(tmp)]
mean(tmp)
sd(tmp)/sqrt(length(tmp)) #standard error


#Healthy Controls

Expr = read.csv("G:/SpinalCord/Publication/FinalEnrichment/SpinalCord_CombinedPlatform_AllSamples_Top1000_MoR_NoTE_NoENSG_9-30-23.csv")
rownames(Expr) = Expr[,1];Expr = Expr[,-1]
allsamples = colnames(Expr)
clinical = read.csv("G:/SpinalCord/Publication/CLINICAL_DATA_PRUDENCIO.csv")
clinical$ExternalSampleId = gsub("-",".",clinical$ExternalSampleId)
Pheno = clinical[clinical$ExternalSampleId %in% allsamples,]
HC_Pheno = Pheno[which(Pheno$disease_group == "Control"),]
#write.csv(HC_Pheno,"SpinalCord_ControlPheno_SampleLevel_10-10-23.csv")

table(HC_Pheno$Sample.Source)

#Patient level
HCpatients = names(table(HC_Pheno$ExternalSubjectId))

HCPat_Pheno = data.frame(matrix(NA,length(HCpatients),ncol = ncol(HC_Pheno)))
colnames(HCPat_Pheno) = colnames(HC_Pheno)
HCPat_Pheno$ExternalSubjectId = HCpatients

for(i in 1:nrow(HCPat_Pheno)){
  ind = which(HC_Pheno$ExternalSubjectId == HCPat_Pheno$ExternalSubjectId[i])[1]
  HCPat_Pheno[i,] = HC_Pheno[ind,]
}

Clean_HCPat = HCPat_Pheno[,c(3,4,5,7,8,9,10,11,12,13,14,15,16,18)]
#write.csv(Clean_HCPat,"SpinalCord_ControlPheno_PatientLevel_10-10-23.csv")

table(Clean_HCPat$Sex)

Clean_HCPat$Age.at.Death[which(Clean_HCPat$Age.at.Death == "90 or Older")] = 90
table(Clean_HCPat$Age.at.Death)
AoD_HC = as.numeric(Clean_HCPat$Age.at.Death)
mean(AoD_HC)
sd(AoD_HC)/sqrt(length(AoD_HC)) #standard error


############################################################################################################################################
################################ Plot clinical parameters ##################################################################################
############################################################################################################################################

#This section considers clinical parameters without discordant subjects
#Discordant subjects are considered in the next section of this script

################################################# Age of onset

TDaoo = as.numeric(TD$AgeofOnset)
TDaoo = TDaoo[!is.na(TDaoo)]
OXaoo = as.numeric(O$AgeofOnset)
OXaoo = OXaoo[!is.na(OXaoo)]
GLaoo = as.numeric(G$AgeofOnset)
GLaoo = GLaoo[!is.na(GLaoo)]
Daoo = as.numeric(D$AgeofOnset)
Daoo = Daoo[!is.na(Daoo)]

par(mar = c(5, 5, 3, 2))
boxplot(GLaoo,OXaoo,TDaoo,Daoo,xaxt="n",main=c("Subtype Age of Onset"),cex.axis = 1.5,cex.main=1.5,col=c("goldenrod1","navy","firebrick","gray50"),xlab="Subtype",ylab="Age of Onset (years)",cex.lab=1.5,pch=20,cex=1.4,ylim=c(25,90))
axis(at=1:4,side=1,labels=c("ALS-Glia","ALS-Ox","ALS-TD","Discordant"),cex.axis=1.5)

#Quick ANOVA to check for post-hoc tests
Age = FinalSurvivalMat$AgeofOnset
ST = FinalSurvivalMat$MetaSpinal
anovadat = data.frame(cbind(Age,ST))
anovadat$ST = factor(anovadat$ST,ordered = T)
levels(anovadat$ST)

tmpindex = rep(NA,nrow(anovadat))
count=1
for(i in 1:nrow(anovadat)){
  if(anovadat$Age[i] == "Unknown"){
    tmpindex[count] = i
    count = count+1
  }else if(anovadat$Age[i] == "Not Applicable"){
    tmpindex[count] = i
    count = count+1
  }else if(anovadat$ST[i] == "Discordant"){
    tmpindex[count] = i
    count = count+1
  }
}
tmpindex = tmpindex[!is.na(tmpindex)]
NoDisc = anovadat[-tmpindex,]

oneway = aov(Age~ST,data=NoDisc)
summary(oneway) #suggests significant differences in age of onset dependent on subtype

oneway = aov(Age~ST,data=anovadat)
summary(oneway) #suggests significant differences in age of onset dependent on subtype

#Post hoc
PH1 = t.test(OXaoo,GLaoo) #p = 0.3457
PH2 = t.test(TDaoo,GLaoo) #p = 0.1563
PH3 = t.test(TDaoo,OXaoo) #p = 0.6554
PH4 = t.test(Daoo,GLaoo) #p = 0.0030
PH5 = t.test(Daoo,OXaoo) #p = 0.0541
PH6 = t.test(Daoo,TDaoo) #p = 0.1345

#MHTC 
pvals = c(PH1$p.value,PH2$p.value,PH3$p.value,PH4$p.value,PH5$p.value,PH6$p.value)
p.adjust(pvals,method = "fdr")

################################################# Site of onset (categories: Bulbar, Limb, Other)

#There is extra code at the end of the script to consider all categories: Axial, Axial-Limb, Axial-Bulbar, Bulbar, Limb, Bulbar-Limb, Generalized, Unknown
SiteofOnset = FinalSurvivalMat
for(i in 1:nrow(SiteofOnset)){
  if(SiteofOnset$Site.of.Onset[i] != "Limb" && SiteofOnset$Site.of.Onset[i] != "Bulbar"){
    SiteofOnset$Site.of.Onset[i] = "Other"
  }
}


tmpdata = SiteofOnset

ite = which(tmpdata$MetaSpinal == "TD")
iox = which(tmpdata$MetaSpinal == "OX")
igl = which(tmpdata$MetaSpinal == "GLIA")
idi = which(tmpdata$MetaSpinal == "Discordant")

tmpTE = tmpdata[ite,]
tmpOX = tmpdata[iox,]
tmpGL = tmpdata[igl,]
tmpDI = tmpdata[idi,]

TEsoo = tmpTE$Site.of.Onset
OXsoo = tmpOX$Site.of.Onset
GLsoo = tmpGL$Site.of.Onset
DIsoo = tmpDI$Site.of.Onset

a = table(TEsoo)
b = table(OXsoo)
c = table(GLsoo)
d = table(DIsoo)

barplot(a,main="ALS-TD Site of Onset",names.arg = c("Bulbar","Limb","Other"),cex.axis = 1.25,cex.names = 1.25,cex.lab=1.5,xlab="Subtype",ylab="Number of Patients",col="firebrick",ylim=c(0,60))
barplot(b,main="ALS-Ox Site of Onset",names.arg = c("Bulbar","Limb","Other"),cex.axis = 1.25,cex.names = 1.25,cex.lab=1.5,xlab="Subtype",ylab="Number of Patients",col="navy",ylim=c(0,60))
barplot(c,main="ALS-Glia Site of Onset",names.arg = c("Bulbar","Limb","Other"),cex.axis = 1.25,cex.names = 1.25,cex.lab=1.5,xlab="Subtype",ylab="Number of Patients",col="goldenrod1",ylim=c(0,60))
barplot(d,main="ALS-Discordant Site of Onset",names.arg = c("Bulbar","Limb","Other"),cex.axis = 1.25,cex.names = 1.25,cex.lab=1.5,xlab="Subtype",ylab="Number of Patients",col="gray50",ylim=c(0,60))

#Chi-squared test of independence (categorical/freq data)
CT = matrix(c(a,b,c),nrow = 3,ncol = 3)
rownames(CT)= c("Bulbar","Limb","Other")
colnames(CT)= c("ALS-TE","ALS-OX","ALS-GLIA")
CT = t(CT)

chisq.test(CT) # p = 0.75 - site of onset not associated with subtype

################################################# Age of death 

TDaod = as.numeric(TD$AgeofDeath)
TDaod = TEaod[!is.na(TEaod)]
OXaod = as.numeric(O$AgeofDeath)
OXaod = OXaod[!is.na(OXaod)]
GLaod = as.numeric(G$AgeofDeath)
GLaod = GLaod[!is.na(GLaod)]
Daod = as.numeric(D$AgeofDeath)
Daod = Daod[!is.na(Daod)]

par(mar = c(5, 5, 3, 2))
boxplot(GLaod,OXaod,TDaod,Daod,xaxt="n",main=c("Subtype Age of Death"),cex.axis = 1.5,cex.main=1.5,col=c("goldenrod1","navy","firebrick","gray50"),xlab="Subtype",ylab="Age of Death (years)",cex.lab=1.5,pch=20,cex=1.4,ylim=c(30,90))
axis(at=1:4,side=1,labels=c("ALS-Glia","ALS-Ox","ALS-TD","Discordant"),cex.axis=1.5)

#Quick ANOVA to check for post-hoc tests
Age = FinalSurvivalMat$AgeofDeath
ST = FinalSurvivalMat$MetaSpinal
anovadat = data.frame(cbind(Age,ST))
anovadat$ST = factor(anovadat$ST,ordered = T)
levels(anovadat$ST)

tmpindex = rep(NA,nrow(anovadat))
count=1
for(i in 1:nrow(anovadat)){
  if(anovadat$Age[i] == "Unknown"){
    tmpindex[count] = i
    count = count+1
  }else if(anovadat$Age[i] == "Not Applicable"){
    tmpindex[count] = i
    count = count+1
  }else if(anovadat$ST[i] == "Discordant"){
    tmpindex[count] = i
    count = count+1
  }
}
tmpindex = tmpindex[!is.na(tmpindex)]
NoDisc2 = anovadat[-tmpindex,]

oneway = aov(Age~ST,data=NoDisc2)
summary(oneway) #suggests no significant differences in age of death by subtype

oneway = aov(Age~ST,data=anovadat)
summary(oneway) #suggests no significant differences in age of death dependent on subtype

#Post hoc (not necessary)
PH1 = t.test(OXaod,GLaod) #p = 0.4802
PH2 = t.test(TDaod,GLaod) #p = 0.3412
PH3 = t.test(TDaod,OXaod) #p = 0.7385
PH4 = t.test(Daod,GLaod) #p = 0.0190
PH5 = t.test(Daod,OXaod) #p = 0.0532
PH6 = t.test(Daod,TDaod) #p = 0.1406

#MHTC 
pvals = c(PH1$p.value,PH2$p.value,PH3$p.value,PH4$p.value,PH5$p.value,PH6$p.value)
p.adjust(pvals,method = "fdr")

################################################# FTD Comorbidity

#Add in FTLD and Alzheimer's Comorbidity
table(G$Disease.Group)
table(O$Disease.Group)
table(TD$Disease.Group)
table(D$Disease.Group)

#Frontotemporal dementia comorbidity
gcomo = table(G$Disease.Group)[[3]]/length(G$Disease.Group) #FTLD
ocomo = table(O$Disease.Group)[[3]]/length(O$Disease.Group) #FTLD
tcomo = table(TD$Disease.Group)[[4]]/length(TD$Disease.Group) #FTLD
dcomo = table(D$Disease.Group)[[4]]/length(D$Disease.Group) #FTLD
vd = c(gcomo,ocomo,tcomo,dcomo)
barplot(vd*100,main="ALS Subtype FTLD Comorbidity",names.arg = c("ALS-Glia","ALS-Ox","ALS-TD","Discordant"),cex.axis = 1.5,cex.names = 1.5,cex.lab=1.5,xlab="Subtype",ylab="FTLD Comorbidity (%)",col=c("goldenrod1","navy","firebrick","gray50"),ylim = c(0,20))

#Alzheimers comorbidity
gcomo2 = table(G$Disease.Group)[[2]]/length(G$Disease.Group) #AD
ocomo2 = table(O$Disease.Group)[[2]]/length(O$Disease.Group) #AD
tcomo2 = table(TD$Disease.Group)[[3]]/length(TD$Disease.Group) #AD
dcomo2 = table(D$Disease.Group)[[3]]/length(D$Disease.Group) #AD
vd2 = c(gcomo2,ocomo2,tcomo2,dcomo2)
barplot(vd2*100,main="ALS Subtype Alzheimer's Comorbidity",names.arg = c("ALS-Glia","ALS-Ox","ALS-TD","Discordant"),cex.axis = 1.5,cex.names = 1.5,cex.lab=1.5,xlab="Subtype",ylab="Alzheimers Comorbidity (%)",col=c("goldenrod1","navy","firebrick","gray50"),ylim = c(0,15))

#Chi-Squared Test of Independence - FTLD

CTg = c(table(G$Disease.Group)[[3]],length(G$Disease.Group) - table(G$Disease.Group)[[3]])
CTo = c(table(O$Disease.Group)[[3]],length(O$Disease.Group) - table(O$Disease.Group)[[3]])
CTt = c(table(TD$Disease.Group)[[4]],length(TD$Disease.Group) - table(TD$Disease.Group)[[4]])
CTd = c(table(D$Disease.Group)[[4]],length(D$Disease.Group) - table(D$Disease.Group)[[4]])

CT = rbind(CTg,CTo,CTt,CTd)
colnames(CT) = c("FTLD","NotFTLD")
chisq.test(CT)

#Chi-Squared Test of Independence - AD

CTg = c(table(G$Disease.Group)[[2]],length(G$Disease.Group) - table(G$Disease.Group)[[2]])
CTo = c(table(O$Disease.Group)[[2]],length(O$Disease.Group) - table(O$Disease.Group)[[2]])
CTt = c(table(TD$Disease.Group)[[3]],length(TD$Disease.Group) - table(TD$Disease.Group)[[3]])
CTd = c(table(D$Disease.Group)[[3]],length(D$Disease.Group) - table(D$Disease.Group)[[3]])

CT = rbind(CTg,CTo,CTt,CTd)
colnames(CT) = c("AD","NotAD")
chisq.test(CT)