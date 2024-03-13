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

#Run independently in each tissue
table(CNS_MEM$tissue)
CNS_MEM_FC = CNS_MEM[which(CNS_MEM$tissue == "Frontal Cortex"),]
CNS_MEM_LMC = CNS_MEM[which(CNS_MEM$tissue == "Lateral Motor Cortex"),]
CNS_MEM_MMC = CNS_MEM[which(CNS_MEM$tissue == "Medial Motor Cortex"),]
CNS_MEM_UMC = CNS_MEM[which(CNS_MEM$tissue == "Other Motor Cortex"),]
CNS_MEM_Cerv = CNS_MEM[which(CNS_MEM$tissue == "Spinal_Cord_Cervical"),]
CNS_MEM_Thor = CNS_MEM[which(CNS_MEM$tissue == "Spinal_Cord_Thoracic"),]
CNS_MEM_Lumb = CNS_MEM[which(CNS_MEM$tissue == "Spinal_Cord_Lumbar"),]

########################################### COX Proportional Hazard MODEL ##################################################################
length(which(table(CNS_MEM_FC$PatientID)>1))
length(which(table(CNS_MEM_LMC$PatientID)>1))
length(which(table(CNS_MEM_MMC$PatientID)>1))
length(which(table(CNS_MEM_UMC$PatientID)>1))
length(which(table(CNS_MEM_Cerv$PatientID)>1))
length(which(table(CNS_MEM_Thor$PatientID)>1))
length(which(table(CNS_MEM_Lumb$PatientID)>1))
#No repeated measures, so +(1|PatientID) argument (random effect intercept) is not needed
#library(coxme)

#Frontal Cortex
coxfit = coxph(Surv(time) ~ SubjectSex + AgeOnset + disease_group + Subtype, data=CNS_MEM_FC,model = T)
summary(coxfit)
ggforest(coxfit,data = CNS_MEM_FC) #plot hazard ratios
cox.zph(coxfit) #check proportional hazard assumption

#Lateral Motor Cortex
coxfit = coxph(Surv(time) ~ SubjectSex + AgeOnset + disease_group + Subtype, data=CNS_MEM_LMC,model = T)
summary(coxfit)
ggforest(coxfit,data = CNS_MEM_LMC)
cox.zph(coxfit) #check proportional hazard assumption

#Medial Motor Cortex
coxfit = coxph(Surv(time) ~ SubjectSex + AgeOnset + disease_group + Subtype, data=CNS_MEM_MMC,model = T)
summary(coxfit)
ggforest(coxfit,data = CNS_MEM_MMC)
cox.zph(coxfit) #check proportional hazard assumption

#Other Motor Cortex
coxfit = coxph(Surv(time) ~ SubjectSex + AgeOnset + disease_group + Subtype, data=CNS_MEM_UMC,model = T)
summary(coxfit)
ggforest(coxfit,data = CNS_MEM_UMC) 
cox.zph(coxfit) #check proportional hazard assumption

#Cervical Spinal Cord
coxfit = coxph(Surv(time) ~ SubjectSex + AgeOnset + disease_group + Subtype, data=CNS_MEM_Cerv,model = T)
summary(coxfit)
ggforest(coxfit,data = CNS_MEM_Cerv)
cox.zph(coxfit) #check proportional hazard assumption

#Thoracic Spinal Cord
CNS_MEM_Thor$disease_group = factor(CNS_MEM_Thor$disease_group,levels = c("ALS-TDP","ALS/AD","ALS/FTLD"))
coxfit = coxph(Surv(time) ~ SubjectSex + AgeOnset + disease_group + Subtype, data=CNS_MEM_Thor,model = T)
summary(coxfit)
ggforest(coxfit,data = CNS_MEM_Thor) 
cox.zph(coxfit) #check proportional hazard assumption

#Lumbar Spinal Cord
coxfit = coxph(Surv(time) ~ SubjectSex + AgeOnset + disease_group + Subtype, data=CNS_MEM_Lumb,model = T)
summary(coxfit)
ggforest(coxfit,data = CNS_MEM_Lumb) 
cox.zph(coxfit) #check proportional hazard assumption


#################### Without disease group covariate

#Frontal Cortex
coxfit = coxph(Surv(time) ~ SubjectSex + AgeOnset + Subtype, data=CNS_MEM_FC,model = T)
summary(coxfit)
ggforest(coxfit,data = CNS_MEM_FC) #plot hazard ratios
cox.zph(coxfit) #check proportional hazard assumption

#Lateral Motor Cortex
coxfit = coxph(Surv(time) ~ SubjectSex + AgeOnset + Subtype, data=CNS_MEM_LMC,model = T)
summary(coxfit)
ggforest(coxfit,data = CNS_MEM_LMC)
cox.zph(coxfit) #check proportional hazard assumption

#Medial Motor Cortex
coxfit = coxph(Surv(time) ~ SubjectSex + AgeOnset + Subtype, data=CNS_MEM_MMC,model = T)
summary(coxfit)
ggforest(coxfit,data = CNS_MEM_MMC)
cox.zph(coxfit) #check proportional hazard assumption

#Other Motor Cortex
coxfit = coxph(Surv(time) ~ SubjectSex + AgeOnset + Subtype, data=CNS_MEM_UMC,model = T)
summary(coxfit)
ggforest(coxfit,data = CNS_MEM_UMC) 
cox.zph(coxfit) #check proportional hazard assumption

#Cervical Spinal Cord
coxfit = coxph(Surv(time) ~ SubjectSex + AgeOnset + Subtype, data=CNS_MEM_Cerv,model = T)
summary(coxfit)
ggforest(coxfit,data = CNS_MEM_Cerv)
cox.zph(coxfit) #check proportional hazard assumption

#Thoracic Spinal Cord
coxfit = coxph(Surv(time) ~ SubjectSex + AgeOnset + Subtype, data=CNS_MEM_Thor,model = T)
summary(coxfit)
ggforest(coxfit,data = CNS_MEM_Thor) 
cox.zph(coxfit) #check proportional hazard assumption

#Lumbar Spinal Cord
coxfit = coxph(Surv(time) ~ SubjectSex + AgeOnset + Subtype, data=CNS_MEM_Lumb,model = T)
summary(coxfit)
ggforest(coxfit,data = CNS_MEM_Lumb) 
cox.zph(coxfit) #check proportional hazard assumption

################################################# Frontal Cortex Diagnostics #############################################################################
coxfit = coxph(Surv(time) ~ SubjectSex + AgeOnset + disease_group + Subtype, data=CNS_MEM_FC,model = T)
summary(coxfit)

par(mfrow=c(1,3))
#Sex Residuals
boxplot(coxfit$residuals[which(coxfit$model$SubjectSex == "Male")],coxfit$residuals[which(coxfit$model$SubjectSex == "Female")],xaxt="n",ylab="Residuals",main="Sex Residuals",col=c("skyblue1","pink3"),outline = F,ylim=c(-3.5,1.5))
axis(1,at=c(1,2),labels = c("Male","Female"),cex.axis=1.5)
beeswarm(at=1,coxfit$residuals[which(coxfit$model$SubjectSex == "Male")],pch=19,cex=1,add=T)
beeswarm(at=2,coxfit$residuals[which(coxfit$model$SubjectSex == "Female")],pch=19,cex=1,add=T)

#Site Residuals
# boxplot(coxfit$residuals[which(coxfit$model$Site == "NYGC ALS Consortium")],coxfit$residuals[which(coxfit$model$Site == "Target ALS PM core")],xaxt="n",ylab="Residuals",main="Site Residuals",col=c("#3c9955","#d49531"),outline = F)
# axis(1,at=c(1,2),labels = c("NYGC ALS Consortium","Target ALS PM core"))
# beeswarm(at=1,coxfit$residuals[which(coxfit$model$Site == "NYGC ALS Consortium")],pch=19,cex=1,add=T)
# beeswarm(at=2,coxfit$residuals[which(coxfit$model$Site == "Target ALS PM core")],pch=19,cex=1,add=T)

#Subtype Residuals
boxplot(coxfit$residuals[which(coxfit$model$Subtype == "GLIA")],coxfit$residuals[which(coxfit$model$Subtype == "OX")],coxfit$residuals[which(coxfit$model$Subtype == "TD")],xaxt="n",ylab="Residuals",main="Subtype Residuals",col=c("goldenrod1","navy","firebrick"),outline = F,ylim=c(-3.5,1.5))
axis(1,at=c(1,2,3),labels = c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis=1.5)
beeswarm(at=1,coxfit$residuals[which(coxfit$model$Subtype == "GLIA")],pch=19,cex=1,add=T)
beeswarm(at=2,coxfit$residuals[which(coxfit$model$Subtype == "OX")],pch=19,cex=1,add=T)
beeswarm(at=3,coxfit$residuals[which(coxfit$model$Subtype == "TD")],pch=19,cex=1,add=T)

#Disease Group Residuals
boxplot(coxfit$residuals[which(coxfit$model$disease_group == "ALS-TDP")],coxfit$residuals[which(coxfit$model$disease_group == "ALS-SOD1")],coxfit$residuals[which(coxfit$model$disease_group == "ALS/AD")],coxfit$residuals[which(coxfit$model$disease_group == "ALS/FTLD")],xaxt="n",ylab="Residuals",main="Disease Group Residuals",col=c("#4db388","#4bb3af","#4d7db8","#4648bd"),outline = F,ylim=c(-3.5,1.5))
axis(1,at=c(1,2,3,4),labels = c("ALS-TDP","ALS-SOD1","ALS/AD","ALS/FTLD"),cex.axis=1.2)
beeswarm(at=1,coxfit$residuals[which(coxfit$model$disease_group == "ALS-TDP")],pch=19,cex=0.9,add=T)
beeswarm(at=2,coxfit$residuals[which(coxfit$model$disease_group == "ALS-SOD1")],pch=19,cex=0.9,add=T)
beeswarm(at=3,coxfit$residuals[which(coxfit$model$disease_group == "ALS/AD")],pch=19,cex=0.9,add=T)
beeswarm(at=4,coxfit$residuals[which(coxfit$model$disease_group == "ALS/FTLD")],pch=19,cex=0.9,add=T)
par(mfrow=c(1,1))


#Check the proportional hazard assumption - for the most part, these meet the assumption (LOESS curve roughly 0 with slope 0)
#ggcoxdiagnostics(coxfit,type="schoenfeld",linear.predictions = F,point.size = 2)
#Scaled Schoenfeld residuals - better diagnostic
ggcoxdiagnostics(coxfit,type="scaledsch",linear.predictions = F,point.size = 2) #At the covariate subgroup level



########## log(t) vs log(-log(S)) Diagnostics
#lines should be parallel to each other across time
#Citation: Kleinbaum DG, Klein M. Survival analysis a self-learning text. Springer; 1996.

#Stratify by subtype
coxfitadj = survfit(coxph(Surv(time) ~ SubjectSex + AgeOnset + disease_group + strata(Subtype), data=CNS_MEM_FC,model = T))
tmp = summary(coxfitadj)
plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Subtype")
linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
colnames(linedat) = c("time","S","strata")
linedat = linedat[order(linedat$time),]
lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
lines(linedat$time[which(linedat$strata == 3)],log(-log(linedat$S[which(linedat$strata == 3)])),col="#057a11")
legend(130,-3,legend = c("ALS-TD","ALS-Glia","ALS-Ox"),col = c("black","darkred","#057a11"),lty=1,lwd=4)

###Stratify by other covariates

#Sex
coxfitadj = survfit(coxph(Surv(time) ~ strata(SubjectSex) + AgeOnset + disease_group + Subtype, data=CNS_MEM_FC,model = T))
tmp = summary(coxfitadj)
plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Sex")
linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
colnames(linedat) = c("time","S","strata")
linedat = linedat[order(linedat$time),]
lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
legend(130,-3.5,legend = c("Male","Female"),col = c("black","darkred"),lty=1,lwd=4)

#Site
# coxfitadj = survfit(coxph(Surv(time) ~ SubjectSex + AgeOnset + strata(Site) + disease_group + Subtype, cluster = PatientID, data=CNS_MEM,model = T))
# tmp = summary(coxfitadj)
# plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Site")
# linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
# colnames(linedat) = c("time","S","strata")
# linedat = linedat[order(linedat$time),]
# lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
# lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
# legend(130,-4,legend = c("NYGC","Target ALS"),col = c("black","darkred"),lty=1,lwd=4)

#Disease
coxfitadj = survfit(coxph(Surv(time) ~ SubjectSex + AgeOnset + strata(disease_group) + Subtype, data=CNS_MEM_FC,model = T))
tmp = summary(coxfitadj)
plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Disease")
linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
colnames(linedat) = c("time","S","strata")
linedat = linedat[order(linedat$time),]
lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
lines(linedat$time[which(linedat$strata == 3)],log(-log(linedat$S[which(linedat$strata == 3)])),col="#057a11")
lines(linedat$time[which(linedat$strata == 4)],log(-log(linedat$S[which(linedat$strata == 4)])),col="#05557a")
legend(130,-2.5,legend = c("ALS-TDP","ALS-SOD1","ALS/AD","ALS/FTLD"),col = c("black","darkred","#057a11","#05557a"),lty=1,lwd=4)


################################################# Medial Motor Cortex Diagnostics #############################################################################
coxfit = coxph(Surv(time) ~ SubjectSex + AgeOnset + disease_group + Subtype, data=CNS_MEM_MMC,model = T)
summary(coxfit)

par(mfrow=c(1,3))
#Sex Residuals
boxplot(coxfit$residuals[which(coxfit$model$SubjectSex == "Male")],coxfit$residuals[which(coxfit$model$SubjectSex == "Female")],xaxt="n",ylab="Residuals",main="Sex Residuals",col=c("skyblue1","pink3"),outline = F,ylim=c(-3.5,1.5))
axis(1,at=c(1,2),labels = c("Male","Female"),cex.axis=1.5)
beeswarm(at=1,coxfit$residuals[which(coxfit$model$SubjectSex == "Male")],pch=19,cex=1,add=T)
beeswarm(at=2,coxfit$residuals[which(coxfit$model$SubjectSex == "Female")],pch=19,cex=1,add=T)

#Site Residuals
# boxplot(coxfit$residuals[which(coxfit$model$Site == "NYGC ALS Consortium")],coxfit$residuals[which(coxfit$model$Site == "Target ALS PM core")],xaxt="n",ylab="Residuals",main="Site Residuals",col=c("#3c9955","#d49531"),outline = F)
# axis(1,at=c(1,2),labels = c("NYGC ALS Consortium","Target ALS PM core"))
# beeswarm(at=1,coxfit$residuals[which(coxfit$model$Site == "NYGC ALS Consortium")],pch=19,cex=1,add=T)
# beeswarm(at=2,coxfit$residuals[which(coxfit$model$Site == "Target ALS PM core")],pch=19,cex=1,add=T)

#Subtype Residuals
boxplot(coxfit$residuals[which(coxfit$model$Subtype == "GLIA")],coxfit$residuals[which(coxfit$model$Subtype == "OX")],coxfit$residuals[which(coxfit$model$Subtype == "TD")],xaxt="n",ylab="Residuals",main="Subtype Residuals",col=c("goldenrod1","navy","firebrick"),outline = F,ylim=c(-3.5,1.5))
axis(1,at=c(1,2,3),labels = c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis=1.5)
beeswarm(at=1,coxfit$residuals[which(coxfit$model$Subtype == "GLIA")],pch=19,cex=1,add=T)
beeswarm(at=2,coxfit$residuals[which(coxfit$model$Subtype == "OX")],pch=19,cex=1,add=T)
beeswarm(at=3,coxfit$residuals[which(coxfit$model$Subtype == "TD")],pch=19,cex=1,add=T)

#Disease Group Residuals
boxplot(coxfit$residuals[which(coxfit$model$disease_group == "ALS-TDP")],coxfit$residuals[which(coxfit$model$disease_group == "ALS-SOD1")],coxfit$residuals[which(coxfit$model$disease_group == "ALS/AD")],coxfit$residuals[which(coxfit$model$disease_group == "ALS/FTLD")],xaxt="n",ylab="Residuals",main="Disease Group Residuals",col=c("#4db388","#4bb3af","#4d7db8","#4648bd"),outline = F,ylim=c(-3.5,1.5))
axis(1,at=c(1,2,3,4),labels = c("ALS-TDP","ALS-SOD1","ALS/AD","ALS/FTLD"),cex.axis=1.2)
beeswarm(at=1,coxfit$residuals[which(coxfit$model$disease_group == "ALS-TDP")],pch=19,cex=0.9,add=T)
beeswarm(at=2,coxfit$residuals[which(coxfit$model$disease_group == "ALS-SOD1")],pch=19,cex=0.9,add=T)
beeswarm(at=3,coxfit$residuals[which(coxfit$model$disease_group == "ALS/AD")],pch=19,cex=0.9,add=T)
beeswarm(at=4,coxfit$residuals[which(coxfit$model$disease_group == "ALS/FTLD")],pch=19,cex=0.9,add=T)
par(mfrow=c(1,1))

#Check the proportional hazard assumption - for the most part, these meet the assumption (LOESS curve roughly 0 with slope 0)
#ggcoxdiagnostics(coxfit,type="schoenfeld",linear.predictions = F,point.size = 2)
#Scaled Schoenfeld residuals - better diagnostic
ggcoxdiagnostics(coxfit,type="scaledsch",linear.predictions = F,point.size = 2) #At the covariate subgroup level



########## log(t) vs log(-log(S)) Diagnostics
#lines should be parallel to each other across time
#Citation: Kleinbaum DG, Klein M. Survival analysis a self-learning text. Springer; 1996.

#Stratify by subtype
coxfitadj = survfit(coxph(Surv(time) ~ SubjectSex + AgeOnset + disease_group + strata(Subtype), data=CNS_MEM_MMC,model = T))
tmp = summary(coxfitadj)
plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Subtype")
linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
colnames(linedat) = c("time","S","strata")
linedat = linedat[order(linedat$time),]
lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
lines(linedat$time[which(linedat$strata == 3)],log(-log(linedat$S[which(linedat$strata == 3)])),col="#057a11")
legend(130,-3,legend = c("ALS-TD","ALS-Glia","ALS-Ox"),col = c("black","darkred","#057a11"),lty=1,lwd=4)

###Stratify by other covariates

#Sex
coxfitadj = survfit(coxph(Surv(time) ~ strata(SubjectSex) + AgeOnset + disease_group + Subtype, data=CNS_MEM_MMC,model = T))
tmp = summary(coxfitadj)
plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Sex")
linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
colnames(linedat) = c("time","S","strata")
linedat = linedat[order(linedat$time),]
lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
legend(130,-3.5,legend = c("Male","Female"),col = c("black","darkred"),lty=1,lwd=4)

#Site
# coxfitadj = survfit(coxph(Surv(time) ~ SubjectSex + AgeOnset + strata(Site) + disease_group + Subtype, cluster = PatientID, data=CNS_MEM,model = T))
# tmp = summary(coxfitadj)
# plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Site")
# linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
# colnames(linedat) = c("time","S","strata")
# linedat = linedat[order(linedat$time),]
# lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
# lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
# legend(130,-4,legend = c("NYGC","Target ALS"),col = c("black","darkred"),lty=1,lwd=4)

#Disease
coxfitadj = survfit(coxph(Surv(time) ~ SubjectSex + AgeOnset + strata(disease_group) + Subtype, data=CNS_MEM_MMC,model = T))
tmp = summary(coxfitadj)
plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Disease")
linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
colnames(linedat) = c("time","S","strata")
linedat = linedat[order(linedat$time),]
lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
lines(linedat$time[which(linedat$strata == 3)],log(-log(linedat$S[which(linedat$strata == 3)])),col="#057a11")
lines(linedat$time[which(linedat$strata == 4)],log(-log(linedat$S[which(linedat$strata == 4)])),col="#05557a")
legend(130,-2.5,legend = c("ALS-TDP","ALS-SOD1","ALS/AD","ALS/FTLD"),col = c("black","darkred","#057a11","#05557a"),lty=1,lwd=4)

################################################# Lateral Motor Cortex Diagnostics #############################################################################
coxfit = coxph(Surv(time) ~ SubjectSex + AgeOnset + disease_group + Subtype, data=CNS_MEM_LMC,model = T)
summary(coxfit)

par(mfrow=c(1,3))
#Sex Residuals
boxplot(coxfit$residuals[which(coxfit$model$SubjectSex == "Male")],coxfit$residuals[which(coxfit$model$SubjectSex == "Female")],xaxt="n",ylab="Residuals",main="Sex Residuals",col=c("skyblue1","pink3"),outline = F,ylim=c(-3.5,1.5))
axis(1,at=c(1,2),labels = c("Male","Female"),cex.axis=1.5)
beeswarm(at=1,coxfit$residuals[which(coxfit$model$SubjectSex == "Male")],pch=19,cex=1,add=T)
beeswarm(at=2,coxfit$residuals[which(coxfit$model$SubjectSex == "Female")],pch=19,cex=1,add=T)

#Site Residuals
# boxplot(coxfit$residuals[which(coxfit$model$Site == "NYGC ALS Consortium")],coxfit$residuals[which(coxfit$model$Site == "Target ALS PM core")],xaxt="n",ylab="Residuals",main="Site Residuals",col=c("#3c9955","#d49531"),outline = F)
# axis(1,at=c(1,2),labels = c("NYGC ALS Consortium","Target ALS PM core"))
# beeswarm(at=1,coxfit$residuals[which(coxfit$model$Site == "NYGC ALS Consortium")],pch=19,cex=1,add=T)
# beeswarm(at=2,coxfit$residuals[which(coxfit$model$Site == "Target ALS PM core")],pch=19,cex=1,add=T)

#Subtype Residuals
boxplot(coxfit$residuals[which(coxfit$model$Subtype == "GLIA")],coxfit$residuals[which(coxfit$model$Subtype == "OX")],coxfit$residuals[which(coxfit$model$Subtype == "TD")],xaxt="n",ylab="Residuals",main="Subtype Residuals",col=c("goldenrod1","navy","firebrick"),outline = F,ylim=c(-3.5,1.5))
axis(1,at=c(1,2,3),labels = c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis=1.5)
beeswarm(at=1,coxfit$residuals[which(coxfit$model$Subtype == "GLIA")],pch=19,cex=1,add=T)
beeswarm(at=2,coxfit$residuals[which(coxfit$model$Subtype == "OX")],pch=19,cex=1,add=T)
beeswarm(at=3,coxfit$residuals[which(coxfit$model$Subtype == "TD")],pch=19,cex=1,add=T)

#Disease Group Residuals
boxplot(coxfit$residuals[which(coxfit$model$disease_group == "ALS-TDP")],coxfit$residuals[which(coxfit$model$disease_group == "ALS-SOD1")],coxfit$residuals[which(coxfit$model$disease_group == "ALS/AD")],coxfit$residuals[which(coxfit$model$disease_group == "ALS/FTLD")],xaxt="n",ylab="Residuals",main="Disease Group Residuals",col=c("#4db388","#4bb3af","#4d7db8","#4648bd"),outline = F,ylim=c(-3.5,1.5))
axis(1,at=c(1,2,3,4),labels = c("ALS-TDP","ALS-SOD1","ALS/AD","ALS/FTLD"),cex.axis=1.2)
beeswarm(at=1,coxfit$residuals[which(coxfit$model$disease_group == "ALS-TDP")],pch=19,cex=0.9,add=T)
beeswarm(at=2,coxfit$residuals[which(coxfit$model$disease_group == "ALS-SOD1")],pch=19,cex=0.9,add=T)
beeswarm(at=3,coxfit$residuals[which(coxfit$model$disease_group == "ALS/AD")],pch=19,cex=0.9,add=T)
beeswarm(at=4,coxfit$residuals[which(coxfit$model$disease_group == "ALS/FTLD")],pch=19,cex=0.9,add=T)
par(mfrow=c(1,1))


#Check the proportional hazard assumption - for the most part, these meet the assumption (LOESS curve roughly 0 with slope 0)
#ggcoxdiagnostics(coxfit,type="schoenfeld",linear.predictions = F,point.size = 2)
#Scaled Schoenfeld residuals - better diagnostic
ggcoxdiagnostics(coxfit,type="scaledsch",linear.predictions = F,point.size = 2) #At the covariate subgroup level



########## log(t) vs log(-log(S)) Diagnostics
#lines should be parallel to each other across time
#Citation: Kleinbaum DG, Klein M. Survival analysis a self-learning text. Springer; 1996.

#Stratify by subtype
coxfitadj = survfit(coxph(Surv(time) ~ SubjectSex + AgeOnset + disease_group + strata(Subtype), data=CNS_MEM_LMC,model = T))
tmp = summary(coxfitadj)
plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Subtype")
linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
colnames(linedat) = c("time","S","strata")
linedat = linedat[order(linedat$time),]
lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
lines(linedat$time[which(linedat$strata == 3)],log(-log(linedat$S[which(linedat$strata == 3)])),col="#057a11")
legend(130,-3,legend = c("ALS-TD","ALS-Glia","ALS-Ox"),col = c("black","darkred","#057a11"),lty=1,lwd=4)

###Stratify by other covariates

#Sex
coxfitadj = survfit(coxph(Surv(time) ~ strata(SubjectSex) + AgeOnset + disease_group + Subtype, data=CNS_MEM_LMC,model = T))
tmp = summary(coxfitadj)
plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Sex")
linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
colnames(linedat) = c("time","S","strata")
linedat = linedat[order(linedat$time),]
lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
legend(130,-3.5,legend = c("Male","Female"),col = c("black","darkred"),lty=1,lwd=4)

#Site
# coxfitadj = survfit(coxph(Surv(time) ~ SubjectSex + AgeOnset + strata(Site) + disease_group + Subtype, cluster = PatientID, data=CNS_MEM,model = T))
# tmp = summary(coxfitadj)
# plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Site")
# linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
# colnames(linedat) = c("time","S","strata")
# linedat = linedat[order(linedat$time),]
# lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
# lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
# legend(130,-4,legend = c("NYGC","Target ALS"),col = c("black","darkred"),lty=1,lwd=4)

#Disease
coxfitadj = survfit(coxph(Surv(time) ~ SubjectSex + AgeOnset + strata(disease_group) + Subtype, data=CNS_MEM_LMC,model = T))
tmp = summary(coxfitadj)
plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Disease")
linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
colnames(linedat) = c("time","S","strata")
linedat = linedat[order(linedat$time),]
lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
lines(linedat$time[which(linedat$strata == 3)],log(-log(linedat$S[which(linedat$strata == 3)])),col="#057a11")
lines(linedat$time[which(linedat$strata == 4)],log(-log(linedat$S[which(linedat$strata == 4)])),col="#05557a")
legend(130,-2.5,legend = c("ALS-TDP","ALS-SOD1","ALS/AD","ALS/FTLD"),col = c("black","darkred","#057a11","#05557a"),lty=1,lwd=4)

################################################# Unspec. Motor Cortex Diagnostics #############################################################################
coxfit = coxph(Surv(time) ~ SubjectSex + AgeOnset + disease_group + Subtype, data=CNS_MEM_UMC,model = T)
summary(coxfit)

par(mfrow=c(1,3))
#Sex Residuals
boxplot(coxfit$residuals[which(coxfit$model$SubjectSex == "Male")],coxfit$residuals[which(coxfit$model$SubjectSex == "Female")],xaxt="n",ylab="Residuals",main="Sex Residuals",col=c("skyblue1","pink3"),outline = F,ylim=c(-3.5,1.5))
axis(1,at=c(1,2),labels = c("Male","Female"),cex.axis=1.5)
beeswarm(at=1,coxfit$residuals[which(coxfit$model$SubjectSex == "Male")],pch=19,cex=1,add=T)
beeswarm(at=2,coxfit$residuals[which(coxfit$model$SubjectSex == "Female")],pch=19,cex=1,add=T)

#Site Residuals
# boxplot(coxfit$residuals[which(coxfit$model$Site == "NYGC ALS Consortium")],coxfit$residuals[which(coxfit$model$Site == "Target ALS PM core")],xaxt="n",ylab="Residuals",main="Site Residuals",col=c("#3c9955","#d49531"),outline = F)
# axis(1,at=c(1,2),labels = c("NYGC ALS Consortium","Target ALS PM core"))
# beeswarm(at=1,coxfit$residuals[which(coxfit$model$Site == "NYGC ALS Consortium")],pch=19,cex=1,add=T)
# beeswarm(at=2,coxfit$residuals[which(coxfit$model$Site == "Target ALS PM core")],pch=19,cex=1,add=T)

#Subtype Residuals
boxplot(coxfit$residuals[which(coxfit$model$Subtype == "GLIA")],coxfit$residuals[which(coxfit$model$Subtype == "OX")],coxfit$residuals[which(coxfit$model$Subtype == "TD")],xaxt="n",ylab="Residuals",main="Subtype Residuals",col=c("goldenrod1","navy","firebrick"),outline = F,ylim=c(-3.5,1.5))
axis(1,at=c(1,2,3),labels = c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis=1.5)
beeswarm(at=1,coxfit$residuals[which(coxfit$model$Subtype == "GLIA")],pch=19,cex=1,add=T)
beeswarm(at=2,coxfit$residuals[which(coxfit$model$Subtype == "OX")],pch=19,cex=1,add=T)
beeswarm(at=3,coxfit$residuals[which(coxfit$model$Subtype == "TD")],pch=19,cex=1,add=T)

#Disease Group Residuals
boxplot(coxfit$residuals[which(coxfit$model$disease_group == "ALS-TDP")],coxfit$residuals[which(coxfit$model$disease_group == "ALS-SOD1")],coxfit$residuals[which(coxfit$model$disease_group == "ALS/AD")],coxfit$residuals[which(coxfit$model$disease_group == "ALS/FTLD")],xaxt="n",ylab="Residuals",main="Disease Group Residuals",col=c("#4db388","#4bb3af","#4d7db8","#4648bd"),outline = F,ylim=c(-3.5,1.5))
axis(1,at=c(1,2,3,4),labels = c("ALS-TDP","ALS-SOD1","ALS/AD","ALS/FTLD"),cex.axis=1.2)
beeswarm(at=1,coxfit$residuals[which(coxfit$model$disease_group == "ALS-TDP")],pch=19,cex=0.9,add=T)
#beeswarm(at=2,coxfit$residuals[which(coxfit$model$disease_group == "ALS-SOD1")],pch=19,cex=0.9,add=T)
#beeswarm(at=3,coxfit$residuals[which(coxfit$model$disease_group == "ALS/AD")],pch=19,cex=0.9,add=T)
beeswarm(at=4,coxfit$residuals[which(coxfit$model$disease_group == "ALS/FTLD")],pch=19,cex=0.9,add=T)
par(mfrow=c(1,1))

#Check the proportional hazard assumption - for the most part, these meet the assumption (LOESS curve roughly 0 with slope 0)
#ggcoxdiagnostics(coxfit,type="schoenfeld",linear.predictions = F,point.size = 2)
#Scaled Schoenfeld residuals - better diagnostic
ggcoxdiagnostics(coxfit,type="scaledsch",linear.predictions = F,point.size = 2) #At the covariate subgroup level



########## log(t) vs log(-log(S)) Diagnostics
#lines should be parallel to each other across time
#Citation: Kleinbaum DG, Klein M. Survival analysis a self-learning text. Springer; 1996.

#Stratify by subtype
coxfitadj = survfit(coxph(Surv(time) ~ SubjectSex + AgeOnset + disease_group + strata(Subtype), data=CNS_MEM_UMC,model = T))
tmp = summary(coxfitadj)
plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Subtype")
linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
colnames(linedat) = c("time","S","strata")
linedat = linedat[order(linedat$time),]
lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
lines(linedat$time[which(linedat$strata == 3)],log(-log(linedat$S[which(linedat$strata == 3)])),col="#057a11")
legend(130,-2.5,legend = c("ALS-TD","ALS-Glia","ALS-Ox"),col = c("black","darkred","#057a11"),lty=1,lwd=4)

###Stratify by other covariates

#Sex
coxfitadj = survfit(coxph(Surv(time) ~ strata(SubjectSex) + AgeOnset + disease_group + Subtype, data=CNS_MEM_UMC,model = T))
tmp = summary(coxfitadj)
plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Sex")
linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
colnames(linedat) = c("time","S","strata")
linedat = linedat[order(linedat$time),]
lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
legend(130,-2.75,legend = c("Male","Female"),col = c("black","darkred"),lty=1,lwd=4)

#Site
# coxfitadj = survfit(coxph(Surv(time) ~ SubjectSex + AgeOnset + strata(Site) + disease_group + Subtype, cluster = PatientID, data=CNS_MEM,model = T))
# tmp = summary(coxfitadj)
# plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Site")
# linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
# colnames(linedat) = c("time","S","strata")
# linedat = linedat[order(linedat$time),]
# lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
# lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
# legend(130,-4,legend = c("NYGC","Target ALS"),col = c("black","darkred"),lty=1,lwd=4)

#Disease
coxfitadj = survfit(coxph(Surv(time) ~ SubjectSex + AgeOnset + Site + strata(disease_group) + Subtype, data=CNS_MEM_UMC,model = T))
tmp = summary(coxfitadj)
plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Disease")
linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
colnames(linedat) = c("time","S","strata")
linedat = linedat[order(linedat$time),]
lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
lines(linedat$time[which(linedat$strata == 3)],log(-log(linedat$S[which(linedat$strata == 3)])),col="#057a11")
lines(linedat$time[which(linedat$strata == 4)],log(-log(linedat$S[which(linedat$strata == 4)])),col="#05557a")
legend(130,-2.5,legend = c("ALS-TDP","ALS/FTLD"),col = c("black","darkred"),lty=1,lwd=4)

################################################# Cervical Spinal Diagnostics #############################################################################
coxfit = coxph(Surv(time) ~ SubjectSex + AgeOnset + disease_group + Subtype, data=CNS_MEM_Cerv,model = T)
summary(coxfit)

par(mfrow=c(1,3))
#Sex Residuals
boxplot(coxfit$residuals[which(coxfit$model$SubjectSex == "Male")],coxfit$residuals[which(coxfit$model$SubjectSex == "Female")],xaxt="n",ylab="Residuals",main="Sex Residuals",col=c("skyblue1","pink3"),outline = F,ylim=c(-3.5,1.5))
axis(1,at=c(1,2),labels = c("Male","Female"),cex.axis=1.5)
beeswarm(at=1,coxfit$residuals[which(coxfit$model$SubjectSex == "Male")],pch=19,cex=1,add=T)
beeswarm(at=2,coxfit$residuals[which(coxfit$model$SubjectSex == "Female")],pch=19,cex=1,add=T)

#Site Residuals
# boxplot(coxfit$residuals[which(coxfit$model$Site == "NYGC ALS Consortium")],coxfit$residuals[which(coxfit$model$Site == "Target ALS PM core")],xaxt="n",ylab="Residuals",main="Site Residuals",col=c("#3c9955","#d49531"),outline = F)
# axis(1,at=c(1,2),labels = c("NYGC ALS Consortium","Target ALS PM core"))
# beeswarm(at=1,coxfit$residuals[which(coxfit$model$Site == "NYGC ALS Consortium")],pch=19,cex=1,add=T)
# beeswarm(at=2,coxfit$residuals[which(coxfit$model$Site == "Target ALS PM core")],pch=19,cex=1,add=T)

#Subtype Residuals
boxplot(coxfit$residuals[which(coxfit$model$Subtype == "GLIA")],coxfit$residuals[which(coxfit$model$Subtype == "OX")],coxfit$residuals[which(coxfit$model$Subtype == "TD")],xaxt="n",ylab="Residuals",main="Subtype Residuals",col=c("goldenrod1","navy","firebrick"),outline = F,ylim=c(-3.5,1.5))
axis(1,at=c(1,2,3),labels = c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis=1.5)
beeswarm(at=1,coxfit$residuals[which(coxfit$model$Subtype == "GLIA")],pch=19,cex=1,add=T)
beeswarm(at=2,coxfit$residuals[which(coxfit$model$Subtype == "OX")],pch=19,cex=1,add=T)
beeswarm(at=3,coxfit$residuals[which(coxfit$model$Subtype == "TD")],pch=19,cex=1,add=T)

#Disease Group Residuals
boxplot(coxfit$residuals[which(coxfit$model$disease_group == "ALS-TDP")],coxfit$residuals[which(coxfit$model$disease_group == "ALS-SOD1")],coxfit$residuals[which(coxfit$model$disease_group == "ALS/AD")],coxfit$residuals[which(coxfit$model$disease_group == "ALS/FTLD")],xaxt="n",ylab="Residuals",main="Disease Group Residuals",col=c("#4db388","#4bb3af","#4d7db8","#4648bd"),outline = F,ylim=c(-3.5,1.5))
axis(1,at=c(1,2,3,4),labels = c("ALS-TDP","ALS-SOD1","ALS/AD","ALS/FTLD"),cex.axis=1.2)
beeswarm(at=1,coxfit$residuals[which(coxfit$model$disease_group == "ALS-TDP")],pch=19,cex=0.9,add=T)
beeswarm(at=2,coxfit$residuals[which(coxfit$model$disease_group == "ALS-SOD1")],pch=19,cex=0.9,add=T)
beeswarm(at=3,coxfit$residuals[which(coxfit$model$disease_group == "ALS/AD")],pch=19,cex=0.9,add=T)
beeswarm(at=4,coxfit$residuals[which(coxfit$model$disease_group == "ALS/FTLD")],pch=19,cex=0.9,add=T)
par(mfrow=c(1,1))

#Check the proportional hazard assumption - for the most part, these meet the assumption (LOESS curve roughly 0 with slope 0)
#ggcoxdiagnostics(coxfit,type="schoenfeld",linear.predictions = F,point.size = 2)
#Scaled Schoenfeld residuals - better diagnostic
ggcoxdiagnostics(coxfit,type="scaledsch",linear.predictions = F,point.size = 2) #At the covariate subgroup level



########## log(t) vs log(-log(S)) Diagnostics
#lines should be parallel to each other across time
#Citation: Kleinbaum DG, Klein M. Survival analysis a self-learning text. Springer; 1996.

#Stratify by subtype
coxfitadj = survfit(coxph(Surv(time) ~ SubjectSex + AgeOnset + disease_group + strata(Subtype), data=CNS_MEM_Cerv,model = T))
tmp = summary(coxfitadj)
plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Subtype")
linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
colnames(linedat) = c("time","S","strata")
linedat = linedat[order(linedat$time),]
lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
lines(linedat$time[which(linedat$strata == 3)],log(-log(linedat$S[which(linedat$strata == 3)])),col="#057a11")
legend(130,-3,legend = c("ALS-TD","ALS-Glia","ALS-Ox"),col = c("black","darkred","#057a11"),lty=1,lwd=4)

###Stratify by other covariates

#Sex
coxfitadj = survfit(coxph(Surv(time) ~ strata(SubjectSex) + AgeOnset + disease_group + Subtype, data=CNS_MEM_Cerv,model = T))
tmp = summary(coxfitadj)
plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Sex")
linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
colnames(linedat) = c("time","S","strata")
linedat = linedat[order(linedat$time),]
lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
legend(130,-3.5,legend = c("Male","Female"),col = c("black","darkred"),lty=1,lwd=4)

#Site
# coxfitadj = survfit(coxph(Surv(time) ~ SubjectSex + AgeOnset + strata(Site) + disease_group + Subtype, cluster = PatientID, data=CNS_MEM,model = T))
# tmp = summary(coxfitadj)
# plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Site")
# linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
# colnames(linedat) = c("time","S","strata")
# linedat = linedat[order(linedat$time),]
# lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
# lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
# legend(130,-4,legend = c("NYGC","Target ALS"),col = c("black","darkred"),lty=1,lwd=4)

#Disease
coxfitadj = survfit(coxph(Surv(time) ~ SubjectSex + AgeOnset + strata(disease_group) + Subtype, data=CNS_MEM_Cerv,model = T))
tmp = summary(coxfitadj)
plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Disease")
linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
colnames(linedat) = c("time","S","strata")
linedat = linedat[order(linedat$time),]
lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
lines(linedat$time[which(linedat$strata == 3)],log(-log(linedat$S[which(linedat$strata == 3)])),col="#057a11")
lines(linedat$time[which(linedat$strata == 4)],log(-log(linedat$S[which(linedat$strata == 4)])),col="#05557a")
legend(130,-2.5,legend = c("ALS-TDP","ALS-SOD1","ALS/AD","ALS/FTLD"),col = c("black","darkred","#057a11","#05557a"),lty=1,lwd=4)

################################################# Thoracic Spinal Diagnostics #############################################################################
coxfit = coxph(Surv(time) ~ SubjectSex + AgeOnset + disease_group + Subtype, data=CNS_MEM_Thor,model = T)
summary(coxfit)

par(mfrow=c(1,3))
#Sex Residuals
boxplot(coxfit$residuals[which(coxfit$model$SubjectSex == "Male")],coxfit$residuals[which(coxfit$model$SubjectSex == "Female")],xaxt="n",ylab="Residuals",main="Sex Residuals",col=c("skyblue1","pink3"),outline = F,ylim=c(-3.5,1.5))
axis(1,at=c(1,2),labels = c("Male","Female"),cex.axis=1.5)
beeswarm(at=1,coxfit$residuals[which(coxfit$model$SubjectSex == "Male")],pch=19,cex=1,add=T)
beeswarm(at=2,coxfit$residuals[which(coxfit$model$SubjectSex == "Female")],pch=19,cex=1,add=T)

#Site Residuals
# boxplot(coxfit$residuals[which(coxfit$model$Site == "NYGC ALS Consortium")],coxfit$residuals[which(coxfit$model$Site == "Target ALS PM core")],xaxt="n",ylab="Residuals",main="Site Residuals",col=c("#3c9955","#d49531"),outline = F)
# axis(1,at=c(1,2),labels = c("NYGC ALS Consortium","Target ALS PM core"))
# beeswarm(at=1,coxfit$residuals[which(coxfit$model$Site == "NYGC ALS Consortium")],pch=19,cex=1,add=T)
# beeswarm(at=2,coxfit$residuals[which(coxfit$model$Site == "Target ALS PM core")],pch=19,cex=1,add=T)

#Subtype Residuals
boxplot(coxfit$residuals[which(coxfit$model$Subtype == "GLIA")],coxfit$residuals[which(coxfit$model$Subtype == "OX")],coxfit$residuals[which(coxfit$model$Subtype == "TD")],xaxt="n",ylab="Residuals",main="Subtype Residuals",col=c("goldenrod1","navy","firebrick"),outline = F,ylim=c(-3.5,1.5))
axis(1,at=c(1,2,3),labels = c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis=1.5)
beeswarm(at=1,coxfit$residuals[which(coxfit$model$Subtype == "GLIA")],pch=19,cex=1,add=T)
beeswarm(at=2,coxfit$residuals[which(coxfit$model$Subtype == "OX")],pch=19,cex=1,add=T)
beeswarm(at=3,coxfit$residuals[which(coxfit$model$Subtype == "TD")],pch=19,cex=1,add=T)

#Disease Group Residuals
boxplot(coxfit$residuals[which(coxfit$model$disease_group == "ALS-TDP")],coxfit$residuals[which(coxfit$model$disease_group == "ALS-SOD1")],coxfit$residuals[which(coxfit$model$disease_group == "ALS/AD")],coxfit$residuals[which(coxfit$model$disease_group == "ALS/FTLD")],xaxt="n",ylab="Residuals",main="Disease Group Residuals",col=c("#4db388","#4bb3af","#4d7db8","#4648bd"),outline = F,ylim=c(-3.5,1.5))
axis(1,at=c(1,2,3,4),labels = c("ALS-TDP","ALS-SOD1","ALS/AD","ALS/FTLD"),cex.axis=1.2)
beeswarm(at=1,coxfit$residuals[which(coxfit$model$disease_group == "ALS-TDP")],pch=19,cex=0.9,add=T)
#beeswarm(at=2,coxfit$residuals[which(coxfit$model$disease_group == "ALS-SOD1")],pch=19,cex=0.9,add=T)
beeswarm(at=3,coxfit$residuals[which(coxfit$model$disease_group == "ALS/AD")],pch=19,cex=0.9,add=T)
beeswarm(at=4,coxfit$residuals[which(coxfit$model$disease_group == "ALS/FTLD")],pch=19,cex=0.9,add=T)
par(mfrow=c(1,1))

#Check the proportional hazard assumption - for the most part, these meet the assumption (LOESS curve roughly 0 with slope 0)
#ggcoxdiagnostics(coxfit,type="schoenfeld",linear.predictions = F,point.size = 2)
#Scaled Schoenfeld residuals - better diagnostic
ggcoxdiagnostics(coxfit,type="scaledsch",linear.predictions = F,point.size = 2) #At the covariate subgroup level



########## log(t) vs log(-log(S)) Diagnostics
#lines should be parallel to each other across time
#Citation: Kleinbaum DG, Klein M. Survival analysis a self-learning text. Springer; 1996.

#Stratify by subtype
coxfitadj = survfit(coxph(Surv(time) ~ SubjectSex + AgeOnset + disease_group + strata(Subtype), data=CNS_MEM_Thor,model = T))
tmp = summary(coxfitadj)
plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Subtype")
linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
colnames(linedat) = c("time","S","strata")
linedat = linedat[order(linedat$time),]
lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
lines(linedat$time[which(linedat$strata == 3)],log(-log(linedat$S[which(linedat$strata == 3)])),col="#057a11")
legend(90,-3,legend = c("ALS-TD","ALS-Glia","ALS-Ox"),col = c("black","darkred","#057a11"),lty=1,lwd=4)

###Stratify by other covariates

#Sex
coxfitadj = survfit(coxph(Surv(time) ~ strata(SubjectSex) + AgeOnset + disease_group + Subtype, data=CNS_MEM_Thor,model = T))
tmp = summary(coxfitadj)
plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Sex")
linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
colnames(linedat) = c("time","S","strata")
linedat = linedat[order(linedat$time),]
lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
legend(90,-2.5,legend = c("Male","Female"),col = c("black","darkred"),lty=1,lwd=4)

#Site
# coxfitadj = survfit(coxph(Surv(time) ~ SubjectSex + AgeOnset + strata(Site) + disease_group + Subtype, cluster = PatientID, data=CNS_MEM,model = T))
# tmp = summary(coxfitadj)
# plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Site")
# linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
# colnames(linedat) = c("time","S","strata")
# linedat = linedat[order(linedat$time),]
# lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
# lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
# legend(130,-4,legend = c("NYGC","Target ALS"),col = c("black","darkred"),lty=1,lwd=4)

#Disease
coxfitadj = survfit(coxph(Surv(time) ~ SubjectSex + AgeOnset  + strata(disease_group) + Subtype, data=CNS_MEM_Thor,model = T))
tmp = summary(coxfitadj)
plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Disease")
linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
colnames(linedat) = c("time","S","strata")
linedat = linedat[order(linedat$time),]
lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
lines(linedat$time[which(linedat$strata == 3)],log(-log(linedat$S[which(linedat$strata == 3)])),col="#057a11")
lines(linedat$time[which(linedat$strata == 4)],log(-log(linedat$S[which(linedat$strata == 4)])),col="#05557a")
legend(90,-2.5,legend = c("ALS-TDP","ALS/AD","ALS/FTLD"),col = c("black","darkred","#057a11"),lty=1,lwd=4)

################################################# Lumbar Spinal Diagnostics #############################################################################
coxfit = coxph(Surv(time) ~ SubjectSex + AgeOnset + disease_group + Subtype, data=CNS_MEM_Lumb,model = T)
summary(coxfit)

par(mfrow=c(1,3))
#Sex Residuals
boxplot(coxfit$residuals[which(coxfit$model$SubjectSex == "Male")],coxfit$residuals[which(coxfit$model$SubjectSex == "Female")],xaxt="n",ylab="Residuals",main="Sex Residuals",col=c("skyblue1","pink3"),outline = F,ylim=c(-3.5,1.5))
axis(1,at=c(1,2),labels = c("Male","Female"),cex.axis=1.5)
beeswarm(at=1,coxfit$residuals[which(coxfit$model$SubjectSex == "Male")],pch=19,cex=1,add=T)
beeswarm(at=2,coxfit$residuals[which(coxfit$model$SubjectSex == "Female")],pch=19,cex=1,add=T)

#Site Residuals
# boxplot(coxfit$residuals[which(coxfit$model$Site == "NYGC ALS Consortium")],coxfit$residuals[which(coxfit$model$Site == "Target ALS PM core")],xaxt="n",ylab="Residuals",main="Site Residuals",col=c("#3c9955","#d49531"),outline = F)
# axis(1,at=c(1,2),labels = c("NYGC ALS Consortium","Target ALS PM core"))
# beeswarm(at=1,coxfit$residuals[which(coxfit$model$Site == "NYGC ALS Consortium")],pch=19,cex=1,add=T)
# beeswarm(at=2,coxfit$residuals[which(coxfit$model$Site == "Target ALS PM core")],pch=19,cex=1,add=T)

#Subtype Residuals
boxplot(coxfit$residuals[which(coxfit$model$Subtype == "GLIA")],coxfit$residuals[which(coxfit$model$Subtype == "OX")],coxfit$residuals[which(coxfit$model$Subtype == "TD")],xaxt="n",ylab="Residuals",main="Subtype Residuals",col=c("goldenrod1","navy","firebrick"),outline = F,ylim=c(-3.5,1.5))
axis(1,at=c(1,2,3),labels = c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis=1.5)
beeswarm(at=1,coxfit$residuals[which(coxfit$model$Subtype == "GLIA")],pch=19,cex=1,add=T)
beeswarm(at=2,coxfit$residuals[which(coxfit$model$Subtype == "OX")],pch=19,cex=1,add=T)
beeswarm(at=3,coxfit$residuals[which(coxfit$model$Subtype == "TD")],pch=19,cex=1,add=T)

#Disease Group Residuals
boxplot(coxfit$residuals[which(coxfit$model$disease_group == "ALS-TDP")],coxfit$residuals[which(coxfit$model$disease_group == "ALS-SOD1")],coxfit$residuals[which(coxfit$model$disease_group == "ALS/AD")],coxfit$residuals[which(coxfit$model$disease_group == "ALS/FTLD")],xaxt="n",ylab="Residuals",main="Disease Group Residuals",col=c("#4db388","#4bb3af","#4d7db8","#4648bd"),outline = F,ylim=c(-3.5,1.5))
axis(1,at=c(1,2,3,4),labels = c("ALS-TDP","ALS-SOD1","ALS/AD","ALS/FTLD"),cex.axis=1.2)
beeswarm(at=1,coxfit$residuals[which(coxfit$model$disease_group == "ALS-TDP")],pch=19,cex=0.9,add=T)
beeswarm(at=2,coxfit$residuals[which(coxfit$model$disease_group == "ALS-SOD1")],pch=19,cex=0.9,add=T)
beeswarm(at=3,coxfit$residuals[which(coxfit$model$disease_group == "ALS/AD")],pch=19,cex=0.9,add=T)
beeswarm(at=4,coxfit$residuals[which(coxfit$model$disease_group == "ALS/FTLD")],pch=19,cex=0.9,add=T)
par(mfrow=c(1,1))

#Check the proportional hazard assumption - for the most part, these meet the assumption (LOESS curve roughly 0 with slope 0)
#ggcoxdiagnostics(coxfit,type="schoenfeld",linear.predictions = F,point.size = 2)
#Scaled Schoenfeld residuals - better diagnostic
ggcoxdiagnostics(coxfit,type="scaledsch",linear.predictions = F,point.size = 2) #At the covariate subgroup level



########## log(t) vs log(-log(S)) Diagnostics
#lines should be parallel to each other across time
#Citation: Kleinbaum DG, Klein M. Survival analysis a self-learning text. Springer; 1996.

#Stratify by subtype
coxfitadj = survfit(coxph(Surv(time) ~ SubjectSex + AgeOnset + disease_group + strata(Subtype), data=CNS_MEM_Lumb,model = T))
tmp = summary(coxfitadj)
plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Subtype")
linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
colnames(linedat) = c("time","S","strata")
linedat = linedat[order(linedat$time),]
lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
lines(linedat$time[which(linedat$strata == 3)],log(-log(linedat$S[which(linedat$strata == 3)])),col="#057a11")
legend(130,-3,legend = c("ALS-TD","ALS-Glia","ALS-Ox"),col = c("black","darkred","#057a11"),lty=1,lwd=4)

###Stratify by other covariates

#Sex
coxfitadj = survfit(coxph(Surv(time) ~ strata(SubjectSex) + AgeOnset + disease_group + Subtype, data=CNS_MEM_Lumb,model = T))
tmp = summary(coxfitadj)
plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Sex")
linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
colnames(linedat) = c("time","S","strata")
linedat = linedat[order(linedat$time),]
lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
legend(130,-3.5,legend = c("Male","Female"),col = c("black","darkred"),lty=1,lwd=4)

#Site
# coxfitadj = survfit(coxph(Surv(time) ~ SubjectSex + AgeOnset + strata(Site) + disease_group + Subtype, cluster = PatientID, data=CNS_MEM,model = T))
# tmp = summary(coxfitadj)
# plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Site")
# linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
# colnames(linedat) = c("time","S","strata")
# linedat = linedat[order(linedat$time),]
# lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
# lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
# legend(130,-4,legend = c("NYGC","Target ALS"),col = c("black","darkred"),lty=1,lwd=4)

#Disease
coxfitadj = survfit(coxph(Surv(time) ~ SubjectSex + AgeOnset + strata(disease_group) + Subtype, data=CNS_MEM_Lumb,model = T))
tmp = summary(coxfitadj)
plot(tmp$time,log(-log(tmp$surv)),col = tmp$strata,pch=20,cex=1,xlab = "log(t)",ylab = "log(-log(S))",main="Covariate: Disease")
linedat = data.frame(cbind(tmp$time,tmp$surv,tmp$strata))
colnames(linedat) = c("time","S","strata")
linedat = linedat[order(linedat$time),]
lines(linedat$time[which(linedat$strata == 1)],log(-log(linedat$S[which(linedat$strata == 1)])),col="black")
lines(linedat$time[which(linedat$strata == 2)],log(-log(linedat$S[which(linedat$strata == 2)])),col="darkred")
lines(linedat$time[which(linedat$strata == 3)],log(-log(linedat$S[which(linedat$strata == 3)])),col="#057a11")
lines(linedat$time[which(linedat$strata == 4)],log(-log(linedat$S[which(linedat$strata == 4)])),col="#05557a")
legend(130,-2.5,legend = c("ALS-TDP","ALS-SOD1","ALS/AD","ALS/FTLD"),col = c("black","darkred","#057a11","#05557a"),lty=1,lwd=4)

