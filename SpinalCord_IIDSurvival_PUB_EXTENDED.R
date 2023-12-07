#Cox Proportional Hazard Regression Model for Survival (Frailty model)
#Written by: Jarrett Eshima
#Date: October 12th, 2023

####################################################################################################################################
#Load Packages
library(survival)
library(survminer)
library(beeswarm)
library(ggplot2)

#Load Data
setwd("G:/SpinalCord/Publication/Survival/IID/V2")
CNS_MEM = read.csv("Supplemental_Dataset3.csv")
#Prep Data
CNS_MEM$SubjectSex = as.factor(CNS_MEM$SubjectSex)
CNS_MEM$SubjectSex = relevel(CNS_MEM$SubjectSex,"Male")
CNS_MEM$disease_group = as.factor(CNS_MEM$disease_group)
CNS_MEM$disease_group = relevel(CNS_MEM$disease_group,"ALS-TDP")
CNS_MEM$Subtype = as.factor(CNS_MEM$Subtype)
CNS_MEM$Subtype = relevel(CNS_MEM$Subtype,"TD")
CNS_MEM$Site = as.factor(CNS_MEM$Site)
CNS_MEM$AgeOnset = as.numeric(CNS_MEM$AgeOnset)

########################################### COX PH MODEL ##################################################################

#EXTENDED MODEL - USING STEP FUNCTIONS
#Ref: https://stats.stackexchange.com/questions/144923/extended-cox-model-and-cox-zph/238964#238964

stratadata = survSplit(Surv(time) ~ SubjectSex + AgeOnset + Site + disease_group + Subtype, cut = 20, id = "PatientID",episode = "timegroup", data=CNS_MEM)

#Requires newdata argument - uses mean of covariate values - not biologically relevant / valid / useful
# coxfitadj = survfit(coxph(Surv(time) ~ SubjectSex:strata(tgroup) + AgeOnset + Site + disease_group:strata(tgroup) + strata(Subtype), cluster = PatientID, data=stratadata,model = T))
# ggsurvplot(coxfitadj,data = CNS_MEM,fun = "cumhaz",palette = c("firebrick","goldenrod1","navy"))
# ggsurvplot(coxfitadj,data = CNS_MEM,fun = "pct",palette = c("firebrick","goldenrod1","navy"))

# coxfit = coxph(Surv(time) ~ SubjectSex:strata(timegroup) + AgeOnset + Site + disease_group:strata(timegroup) + Subtype, cluster = PatientID, data=stratadata,model = T)
# summary(coxfit)
# cox.zph(coxfit,transform = "rank",terms=F) 

coxfit = coxph(Surv(time) ~ SubjectSex:strata(timegroup) + AgeOnset + Site + disease_group:strata(timegroup) + Subtype + cluster(PatientID), data=stratadata,model = T)
summary(coxfit)
cox.zph(coxfit,transform = "rank",terms=F) 

##############################################################################

#Manually construct a ggforest-like plot
#Hazard Ratio is the exp of the beta coefficient 

HRstats = summary(coxfit)$conf.int
HRstats = data.frame(HRstats)
colnames(HRstats) = c("HR","iHR","CIlow","CIhigh")

HRstats = HRstats[-which(is.na(rowSums(HRstats))),]

#Add in the reference level covariates
mylevels = c("AgeOnset","SiteTarget ALS PM core","SiteNYGC ALS Consortium","SubtypeGLIA","SubtypeOX","SubtypeTD","SubjectSexMale:strata(timegroup)timegroup=1","SubjectSexMale:strata(timegroup)timegroup=2","SubjectSexFemale:strata(timegroup)timegroup=1","SubjectSexFemale:strata(timegroup)timegroup=2","strata(timegroup)timegroup=1:disease_groupALS-SOD1","strata(timegroup)timegroup=2:disease_groupALS-SOD1","strata(timegroup)timegroup=1:disease_groupALS/AD","strata(timegroup)timegroup=2:disease_groupALS/AD","strata(timegroup)timegroup=1:disease_groupALS/FTLD","strata(timegroup)timegroup=2:disease_groupALS/FTLD","strata(timegroup)timegroup=1:disease_groupALS-TDP","strata(timegroup)timegroup=2:disease_groupALS-TDP")
reflevs = mylevels[! mylevels %in% rownames(HRstats)]
refmat = data.frame(matrix(NA,nrow = length(reflevs),ncol=4))
rownames(refmat) = reflevs; colnames(refmat) = colnames(HRstats)
refmat$HR = 1; refmat$iHR = exp(-1); refmat$CIlow = NA; refmat$CIhigh = NA
HRdat = rbind(HRstats,refmat)

HRdat$covars = rownames(HRdat)

#Reorder for plotting purposes
covarorder = c("SubjectSexMale:strata(timegroup)timegroup=1","SubjectSexMale:strata(timegroup)timegroup=2","SubjectSexFemale:strata(timegroup)timegroup=1","SubjectSexFemale:strata(timegroup)timegroup=2","strata(timegroup)timegroup=1:disease_groupALS-SOD1","strata(timegroup)timegroup=2:disease_groupALS-SOD1","strata(timegroup)timegroup=1:disease_groupALS/AD","strata(timegroup)timegroup=2:disease_groupALS/AD","strata(timegroup)timegroup=1:disease_groupALS/FTLD","strata(timegroup)timegroup=2:disease_groupALS/FTLD","strata(timegroup)timegroup=1:disease_groupALS-TDP","strata(timegroup)timegroup=2:disease_groupALS-TDP","SiteTarget ALS PM core","SiteNYGC ALS Consortium","SubtypeGLIA","SubtypeOX","SubtypeTD","AgeOnset")
HRdat2 = HRdat
rownames(HRdat2) = covarorder
for(i in 1:nrow(HRdat2)){
  ind = which(rownames(HRdat) == rownames(HRdat2)[i])
  HRdat2[i,] = HRdat[ind,]
  
}

#Get ggforest count numbers from coxfit$model
table(coxfit$model$SubjectSex)
table(coxfit$model$Subtype)
table(coxfit$model$Site)
table(coxfit$model$disease_group)

HRdat2$covars = factor(HRdat2$covars,levels = HRdat2$covars)

p = ggplot(HRdat2,aes(x=HR, y=covars)) + geom_point(size=5,shape=15)
p = p+geom_errorbar(aes(xmin = CIlow,xmax = CIhigh),size=1,width=0.1)
p = p+theme_bw()
p = p+xlim(-1,4)
p = p+geom_vline(xintercept = 1,linetype = "dashed")
p


############################################ P-value heatmap #############################################################################################################
#Build heatmap of mixed effect model results
pvals = summary(coxfit)$coefficients
pvals2 = unname(pvals[,6])

pmap = data.frame(pvals2)
rownames(pmap) = rownames(pvals); colnames(pmap) = "pval"

mylevels = c("AgeOnset","SiteTarget ALS PM core","SiteNYGC ALS Consortium","SubtypeGLIA","SubtypeOX","SubtypeTD","SubjectSexMale:strata(tgroup)tgroup=1","SubjectSexMale:strata(tgroup)tgroup=2","SubjectSexFemale:strata(tgroup)tgroup=1","SubjectSexFemale:strata(tgroup)tgroup=2","strata(tgroup)tgroup=1:disease_groupALS-SOD1","strata(tgroup)tgroup=2:disease_groupALS-SOD1","strata(tgroup)tgroup=1:disease_groupALS/AD","strata(tgroup)tgroup=2:disease_groupALS/AD","strata(tgroup)tgroup=1:disease_groupALS/FTLD","strata(tgroup)tgroup=2:disease_groupALS/FTLD","strata(tgroup)tgroup=1:disease_groupALS-TDP","strata(tgroup)tgroup=2:disease_groupALS-TDP")
reflevs = mylevels[! mylevels %in% rownames(pmap)]
refmat = data.frame(matrix(NA,nrow = length(reflevs),ncol=1))
rownames(refmat) = reflevs; colnames(refmat) = "pval"
refmat$pval = NA
pmap2 = rbind(pmap,refmat)

#Reorder for plotting purposes
covarorder = c("SubjectSexMale:strata(tgroup)tgroup=1","SubjectSexMale:strata(tgroup)tgroup=2","SubjectSexFemale:strata(tgroup)tgroup=1","SubjectSexFemale:strata(tgroup)tgroup=2","strata(tgroup)tgroup=1:disease_groupALS-SOD1","strata(tgroup)tgroup=2:disease_groupALS-SOD1","strata(tgroup)tgroup=1:disease_groupALS/AD","strata(tgroup)tgroup=2:disease_groupALS/AD","strata(tgroup)tgroup=1:disease_groupALS/FTLD","strata(tgroup)tgroup=2:disease_groupALS/FTLD","strata(tgroup)tgroup=1:disease_groupALS-TDP","strata(tgroup)tgroup=2:disease_groupALS-TDP","SiteTarget ALS PM core","SiteNYGC ALS Consortium","SubtypeGLIA","SubtypeOX","SubtypeTD","AgeOnset")
pmap3 = pmap2
rownames(pmap3) = covarorder
for(i in 1:nrow(pmap3)){
  ind = which(rownames(pmap2) == rownames(pmap3)[i])
  pmap3[i,] = pmap2[ind,]
}
MatrixtoHeatmap3(pmap3,limits=c(0,1),customfeats = rownames(pmap3),revcolors = T,pmap = T,customsamps = colnames(pmap3))

################################################# Diagnostics #############################################################################

#par(mfrow=c(2,2))
#Sex Residuals
boxplot(coxfit$residuals[which(coxfit$model$SubjectSex == "Male")],coxfit$residuals[which(coxfit$model$SubjectSex == "Female")],xaxt="n",ylab="Residuals",main="Sex Residuals",col=c("skyblue1","pink3"),outline = F)
axis(1,at=c(1,2),labels = c("Male","Female"))
beeswarm(at=1,coxfit$residuals[which(coxfit$model$SubjectSex == "Male")],pch=19,cex=0.75,add=T)
beeswarm(at=2,coxfit$residuals[which(coxfit$model$SubjectSex == "Female")],pch=19,cex=0.75,add=T)

#Site Residuals
boxplot(coxfit$residuals[which(coxfit$model$Site == "NYGC ALS Consortium")],coxfit$residuals[which(coxfit$model$Site == "Target ALS PM core")],xaxt="n",ylab="Residuals",main="Site Residuals",col=c("#3c9955","#d49531"),outline = F)
axis(1,at=c(1,2),labels = c("NYGC ALS Consortium","Target ALS PM core"))
beeswarm(at=1,coxfit$residuals[which(coxfit$model$Site == "NYGC ALS Consortium")],pch=19,cex=0.75,add=T)
beeswarm(at=2,coxfit$residuals[which(coxfit$model$Site == "Target ALS PM core")],pch=19,cex=0.75,add=T)

#Subtype Residuals
boxplot(coxfit$residuals[which(coxfit$model$Subtype == "GLIA")],coxfit$residuals[which(coxfit$model$Subtype == "OX")],coxfit$residuals[which(coxfit$model$Subtype == "TD")],xaxt="n",ylab="Residuals",main="Subtype Residuals",col=c("goldenrod1","navy","firebrick"),outline = F)
axis(1,at=c(1,2,3),labels = c("ALS-Glia","ALS-Ox","ALS-TD"))
beeswarm(at=1,coxfit$residuals[which(coxfit$model$Subtype == "GLIA")],pch=19,cex=0.75,add=T)
beeswarm(at=2,coxfit$residuals[which(coxfit$model$Subtype == "OX")],pch=19,cex=0.75,add=T)
beeswarm(at=3,coxfit$residuals[which(coxfit$model$Subtype == "TD")],pch=19,cex=0.75,add=T)

#Disease Group Residuals
boxplot(coxfit$residuals[which(coxfit$model$disease_group == "ALS-TDP")],coxfit$residuals[which(coxfit$model$disease_group == "ALS-SOD1")],coxfit$residuals[which(coxfit$model$disease_group == "ALS/AD")],coxfit$residuals[which(coxfit$model$disease_group == "ALS/FTLD")],xaxt="n",ylab="Residuals",main="Disease Group Residuals",col=c("#4db388","#4bb3af","#4d7db8","#4648bd"),outline = F)
axis(1,at=c(1,2,3,4),labels = c("ALS-TDP","ALS-SOD1","ALS/AD","ALS/FTLD"))
beeswarm(at=1,coxfit$residuals[which(coxfit$model$disease_group == "ALS-TDP")],pch=19,cex=0.5,add=T)
beeswarm(at=2,coxfit$residuals[which(coxfit$model$disease_group == "ALS-SOD1")],pch=19,cex=0.5,add=T)
beeswarm(at=3,coxfit$residuals[which(coxfit$model$disease_group == "ALS/AD")],pch=19,cex=0.5,add=T)
beeswarm(at=4,coxfit$residuals[which(coxfit$model$disease_group == "ALS/FTLD")],pch=19,cex=0.5,add=T)
#par(mfrow=c(1,1))

#Add point colors
customcols = rep(NA,nrow(coxfit$model))
for(i in 1:nrow(coxfit$model)){
  if(coxfit$model$Subtype[i] == "GLIA"){
    customcols[i] = "goldenrod1"
  }else if(coxfit$model$Subtype[i] == "OX"){
    customcols[i] = "navy"
  }else if(coxfit$model$Subtype[i] == "TD"){
    customcols[i] = "firebrick"
  }
}

#Check the proportional hazard assumption - for the most part, these meet the assumption (LOESS curve roughly 0 with slope 0)
ggcoxdiagnostics(coxfit,type="schoenfeld",linear.predictions = F,point.size = 2)
#Scaled Schoenfeld residuals - better diagnostic
ggcoxdiagnostics(coxfit,type="scaledsch",linear.predictions = F,point.size = 2) #At the covariate subgroup level

#At the covariate sublevel - ID troubling factor levels
par(mfrow=c(3,3))
plot(cox.zph(coxfit,transform = "rank",terms=F),col="#e63552") #this uses scaled schoenfeld residuals
cox.zph(coxfit,transform = "rank",terms=F)
par(mfrow=c(1,1))

#At the covariate level
par(mfrow=c(3,2))
plot(cox.zph(coxfit,transform = "rank",terms=T),col="#e63552")
cox.zph(coxfit,transform = "rank",terms=T)
par(mfrow=c(1,1))

#Different time scales for testing proportional hazard assumption
#Good reference: https://stats.stackexchange.com/questions/411075/what-does-km-transform-in-cox-zph-function-mean
# plot(cox.zph(coxfit,transform = "rank",terms=F))
# cox.zph(coxfit,transform = "rank",terms=F)
#
# plot(cox.zph(coxfit,transform = "km",terms=F))
# cox.zph(coxfit,transform = "km",terms=F)
# 
# plot(cox.zph(coxfit,transform = "identity",terms=F))
# cox.zph(coxfit,transform = "identity",terms=F)
# 
# plot(cox.zph(coxfit,transform = log,terms=F))
# cox.zph(coxfit,transform = log,terms=F)


#Check for outliers and influential points
ggcoxdiagnostics(coxfit,type="deviance",linear.predictions = F,point.col=customcols,point.size = 2)
ggcoxdiagnostics(coxfit,type="deviance",linear.predictions = T,point.col=customcols,point.size = 2)

#Score
ggcoxdiagnostics(coxfit,type="score",linear.predictions = F,point.col = rep(customcols,8),point.size = 2) #in the same order as the original matrix, according to R documentation

#dfbetas
ggcoxdiagnostics(coxfit,type="dfbetas",linear.predictions = F,point.size = 2)


#Residual vs Linear Predictor
plot(predict(coxfit,type="lp"),coxfit$residuals,main="Cortex Survival, Cox Proportional Hazard Regression - Residual vs Fitted",ylab="Residual",pch=19,col=customcols)
legend(1.15,-3.5,legend = c("ALS-Glia","ALS-Ox","ALS-TD"),col=c("goldenrod1","navy","firebrick"),pch = 20,pt.cex=2)

# #Survival vs LP
# plot(predict(coxfit,type="survival"),predict(coxfit,type="lp"),main="Cortex Survival, Cox Proportional Hazard Regression - Residual vs Fitted",pch=19,col=customcols)
# legend(1.15,-3.5,legend = c("ALS-Glia","ALS-Ox","ALS-TD"),col=c("goldenrod1","navy","firebrick"),pch = 20,pt.cex=2)
# plot(predict(coxfit,type="lp"),predict(coxfit,type="survival"),main="Cortex Survival, Cox Proportional Hazard Regression - Residual vs Fitted",pch=19,col=customcols)
# legend(1.15,-3.5,legend = c("ALS-Glia","ALS-Ox","ALS-TD"),col=c("goldenrod1","navy","firebrick"),pch = 20,pt.cex=2)

############################################# IID SURVIVAL - BY TISSUE REGION ###################################################################

CortexTissue = c("Frontal Cortex","Lateral Motor Cortex","Medial Motor Cortex","Other Motor Cortex")
Cortex_Pheno = CNS_MEM[which(CNS_MEM$tissue %in% CortexTissue),]

Survival = Cortex_Pheno

##################################### Frontal Cortex Survival (IID) ###########################################
FinalSurvival = Survival[which(Survival$tissue == "Frontal Cortex"),]

#Pairwise stats
GO = FinalSurvival[-which(FinalSurvival$Subtype == "TD"),]
km_fit = survfit(Surv(time) ~ Subtype, data=GO)
surv_pvalue(km_fit)

GT = FinalSurvival[-which(FinalSurvival$Subtype == "OX"),]
km_fit = survfit(Surv(time) ~ Subtype, data=GT)
surv_pvalue(km_fit)

TO = FinalSurvival[-which(FinalSurvival$Subtype == "GLIA"),]
km_fit = survfit(Surv(time) ~ Subtype, data=TO)
surv_pvalue(km_fit)


#Full Plot
km_fit = survfit(Surv(time) ~ Subtype, data=FinalSurvival)
surv_pvalue(km_fit)
print(km_fit)
summary(km_fit, times = c(1,30,60,90*(1:10)))

par(mar = c(5, 5, 3, 2))
plot(km_fit,col = c("firebrick","goldenrod1","navy"),main="IID ALS Survival - Frontal Cortex",xlab="Disease Duration (months)",ylab="Survival Probability",lwd=2.5,cex.axis = 1.5,cex.lab=1.5,cex.main=1.5,xaxt="n")
timeseq = seq(0,156,12)
axis(side = 1, at = timeseq,labels = T,tick = T,cex.axis = 1.5)
legend(120,1.02,legend=c("ALS-TD","ALS-Glia","ALS-Ox"),lty = 1,lwd=5,col=c("goldenrod","navy","firebrick"),cex=1.4)
text(135,0.75,"Glia vs Ox p = 0.238",cex=1.25)
text(135,0.7,"Glia vs TD p = 0.520",cex=1.25)
text(136,0.65,"Ox vs TD p = 0.619",cex=1.25)

##################################### Lateral Motor Cortex Survival (IID) #######################################
FinalSurvival = Survival[which(Survival$tissue == "Lateral Motor Cortex"),]

#Pairwise stats
GO = FinalSurvival[-which(FinalSurvival$Subtype == "TD"),]
km_fit = survfit(Surv(time) ~ Subtype, data=GO)
surv_pvalue(km_fit)

GT = FinalSurvival[-which(FinalSurvival$Subtype == "OX"),]
km_fit = survfit(Surv(time) ~ Subtype, data=GT)
surv_pvalue(km_fit)

TO = FinalSurvival[-which(FinalSurvival$Subtype == "GLIA"),]
km_fit = survfit(Surv(time) ~ Subtype, data=TO)
surv_pvalue(km_fit)


#Full Plot
km_fit = survfit(Surv(time) ~ Subtype, data=FinalSurvival)
surv_pvalue(km_fit)
print(km_fit)
summary(km_fit, times = c(1,30,60,90*(1:10)))

par(mar = c(5, 5, 3, 2))
plot(km_fit,col = c("firebrick","goldenrod1","navy"),main="IID ALS Survival - Lateral Motor Cortex",xlab="Disease Duration (months)",ylab="Survival Probability",lwd=2.5,cex.axis = 1.5,cex.lab=1.5,cex.main=1.5,xaxt="n")
timeseq = seq(0,156,12)
axis(side = 1, at = timeseq,labels = T,tick = T,cex.axis = 1.5)
legend(92,1.02,legend=c("ALS-TD","ALS-Glia","ALS-Ox"),lty = 1,lwd=5,col=c("goldenrod","navy","firebrick"),cex=1.4)
text(106,0.75,"Glia vs Ox p = 0.119",cex=1.25)
text(106,0.7,"Glia vs TD p = 0.294",cex=1.25)
text(107,0.65,"Ox vs TD p = 0.393",cex=1.25)

################################## Medial Motor Cortex Survival (IID) ########################################
FinalSurvival = Survival[which(Survival$tissue == "Medial Motor Cortex"),]

#Pairwise stats
GO = FinalSurvival[-which(FinalSurvival$Subtype == "TD"),]
km_fit = survfit(Surv(time) ~ Subtype, data=GO)
surv_pvalue(km_fit)

GT = FinalSurvival[-which(FinalSurvival$Subtype == "OX"),]
km_fit = survfit(Surv(time) ~ Subtype, data=GT)
surv_pvalue(km_fit)

TO = FinalSurvival[-which(FinalSurvival$Subtype == "GLIA"),]
km_fit = survfit(Surv(time) ~ Subtype, data=TO)
surv_pvalue(km_fit)


#Full Plot
km_fit = survfit(Surv(time) ~ Subtype, data=FinalSurvival)
surv_pvalue(km_fit)
print(km_fit)
summary(km_fit, times = c(1,30,60,90*(1:10)))

par(mar = c(5, 5, 3, 2))
plot(km_fit,col = c("firebrick","goldenrod1","navy"),main="IID ALS Survival - Medial Motor Cortex",xlab="Disease Duration (months)",ylab="Survival Probability",lwd=2.5,cex.axis = 1.5,cex.lab=1.5,cex.main=1.5,xaxt="n")
timeseq = seq(0,156,12)
axis(side = 1, at = timeseq,labels = T,tick = T,cex.axis = 1.5)
legend(120,1.02,legend=c("ALS-TD","ALS-Glia","ALS-Ox"),lty = 1,lwd=5,col=c("goldenrod","navy","firebrick"),cex=1.4)
text(135,0.75,"Glia vs Ox p = 0.081",cex=1.25)
text(135,0.7,"Glia vs TD p = 0.082",cex=1.25)
text(136,0.65,"Ox vs TD p = 0.726",cex=1.25)

##################################### Unspecified Motor Cortex Survival (IID) #####################################
FinalSurvival = Survival[which(Survival$tissue == "Other Motor Cortex"),]

#Pairwise stats
GO = FinalSurvival[-which(FinalSurvival$Subtype == "TD"),]
km_fit = survfit(Surv(time) ~ Subtype, data=GO)
surv_pvalue(km_fit)

GT = FinalSurvival[-which(FinalSurvival$Subtype == "OX"),]
km_fit = survfit(Surv(time) ~ Subtype, data=GT)
surv_pvalue(km_fit)

TO = FinalSurvival[-which(FinalSurvival$Subtype == "GLIA"),]
km_fit = survfit(Surv(time) ~ Subtype, data=TO)
surv_pvalue(km_fit)


#Full Plot
km_fit = survfit(Surv(time) ~ Subtype, data=FinalSurvival)
surv_pvalue(km_fit)
print(km_fit)
summary(km_fit, times = c(1,30,60,90*(1:10)))

par(mar = c(5, 5, 3, 2))
plot(km_fit,col = c("firebrick","goldenrod1","navy"),main="IID ALS Survival - Unspecified Motor Cortex",xlab="Disease Duration (months)",ylab="Survival Probability",lwd=2.5,cex.axis = 1.5,cex.lab=1.5,cex.main=1.5,xaxt="n")
timeseq = seq(0,156,12)
axis(side = 1, at = timeseq,labels = T,tick = T,cex.axis = 1.5)
legend(120,1.02,legend=c("ALS-TD","ALS-Glia","ALS-Ox"),lty = 1,lwd=5,col=c("goldenrod","navy","firebrick"),cex=1.4)
text(135,0.75,"Glia vs Ox p = 0.843",cex=1.25)
text(135,0.7,"Glia vs TD p = 0.959",cex=1.25)
text(136,0.65,"Ox vs TD p = 0.858",cex=1.25)

###############################################################

SpinalTissue = c("Spinal_Cord_Cervical","Spinal_Cord_Lumbar","Spinal_Cord_Thoracic")
Spinal_Pheno = CNS_MEM[which(CNS_MEM$tissue %in% SpinalTissue),]

Survival = Spinal_Pheno

##################################### Cervical Survival (IID) #############################################
FinalSurvival = Survival[which(Survival$tissue == "Spinal_Cord_Cervical"),]

#Pairwise stats
GO = FinalSurvival[-which(FinalSurvival$Subtype == "TD"),]
km_fit = survfit(Surv(time) ~ Subtype, data=GO)
surv_pvalue(km_fit)

GT = FinalSurvival[-which(FinalSurvival$Subtype == "OX"),]
km_fit = survfit(Surv(time) ~ Subtype, data=GT)
surv_pvalue(km_fit)

TO = FinalSurvival[-which(FinalSurvival$Subtype == "GLIA"),]
km_fit = survfit(Surv(time) ~ Subtype, data=TO)
surv_pvalue(km_fit)


#Full Plot
km_fit = survfit(Surv(time) ~ Subtype, data=FinalSurvival)
surv_pvalue(km_fit)
print(km_fit)
summary(km_fit, times = c(1,30,60,90*(1:10)))

par(mar = c(5, 5, 3, 2))
plot(km_fit,col = c("firebrick","goldenrod1","navy"),main="IID ALS Survival - Cervical Spinal Cord",xlab="Disease Duration (months)",ylab="Survival Probability",lwd=2.5,cex.axis = 1.5,cex.lab=1.5,cex.main=1.5,xaxt="n")
timeseq = seq(0,156,12)
axis(side = 1, at = timeseq,labels = T,tick = T,cex.axis = 1.5)
legend(120,1.02,legend=c("ALS-TD","ALS-Glia","ALS-Ox"),lty = 1,lwd=5,col=c("goldenrod","navy","firebrick"),cex=1.4)
text(135,0.75,"Glia vs Ox p = 0.029",cex=1.25)
text(135,0.7,"Glia vs TD p = 0.252",cex=1.25)
text(136,0.65,"Ox vs TD p = 0.176",cex=1.25)

####################################### Lumbar Survival (IID) #################################################
FinalSurvival = Survival[which(Survival$tissue == "Spinal_Cord_Lumbar"),]

#Pairwise stats
GO = FinalSurvival[-which(FinalSurvival$Subtype == "TD"),]
km_fit = survfit(Surv(time) ~ Subtype, data=GO)
surv_pvalue(km_fit)

GT = FinalSurvival[-which(FinalSurvival$Subtype == "OX"),]
km_fit = survfit(Surv(time) ~ Subtype, data=GT)
surv_pvalue(km_fit)

TO = FinalSurvival[-which(FinalSurvival$Subtype == "GLIA"),]
km_fit = survfit(Surv(time) ~ Subtype, data=TO)
surv_pvalue(km_fit)


#Full Plot
km_fit = survfit(Surv(time) ~ Subtype, data=FinalSurvival)
surv_pvalue(km_fit)
print(km_fit)
summary(km_fit, times = c(1,30,60,90*(1:10)))

par(mar = c(5, 5, 3, 2))
plot(km_fit,col = c("firebrick","goldenrod1","navy"),main="IID ALS Survival - Lumbar Spinal Cord",xlab="Disease Duration (months)",ylab="Survival Probability",lwd=2.5,cex.axis = 1.5,cex.lab=1.5,cex.main=1.5,xaxt="n")
timeseq = seq(0,156,12)
axis(side = 1, at = timeseq,labels = T,tick = T,cex.axis = 1.5)
legend(120,1.02,legend=c("ALS-TD","ALS-Glia","ALS-Ox"),lty = 1,lwd=5,col=c("goldenrod","navy","firebrick"),cex=1.4)
text(135,0.75,"Glia vs Ox p = 0.434",cex=1.25)
text(135,0.7,"Glia vs TD p = 0.476",cex=1.25)
text(136,0.65,"Ox vs TD p = 0.154",cex=1.25)


################################## Thoracic Survival (IID) ###########################################
FinalSurvival = Survival[which(Survival$tissue == "Spinal_Cord_Thoracic"),]

#Pairwise stats
GO = FinalSurvival[-which(FinalSurvival$Subtype == "TD"),]
km_fit = survfit(Surv(time) ~ Subtype, data=GO)
surv_pvalue(km_fit)

GT = FinalSurvival[-which(FinalSurvival$Subtype == "OX"),]
km_fit = survfit(Surv(time) ~ Subtype, data=GT)
surv_pvalue(km_fit)

TO = FinalSurvival[-which(FinalSurvival$Subtype == "GLIA"),]
km_fit = survfit(Surv(time) ~ Subtype, data=TO)
surv_pvalue(km_fit)


#Full Plot
km_fit = survfit(Surv(time) ~ Subtype, data=FinalSurvival)
surv_pvalue(km_fit)
print(km_fit)
summary(km_fit, times = c(1,30,60,90*(1:10)))

par(mar = c(5, 5, 3, 2))
plot(km_fit,col = c("firebrick","goldenrod1","navy"),main="IID ALS Survival - Thoracic Spinal Cord",xlab="Disease Duration (months)",ylab="Survival Probability",lwd=2.5,cex.axis = 1.5,cex.lab=1.5,cex.main=1.5,xaxt="n")
timeseq = seq(0,108,12)
axis(side = 1, at = timeseq,labels = T,tick = T,cex.axis = 1.5)
legend(85,1.02,legend=c("ALS-TD","ALS-Glia","ALS-Ox"),lty = 1,lwd=5,col=c("goldenrod","navy","firebrick"),cex=1.4)
text(95,0.75,"Glia vs Ox p = 0.110",cex=1.25)
text(95,0.7,"Glia vs TD p = 0.150",cex=1.25)
text(95.5,0.65,"Ox vs TD p = 0.199",cex=1.25)
