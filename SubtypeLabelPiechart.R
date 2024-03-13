wd = "G:/SpinalCord/Publication/BonusCortex"
setwd(wd)

##############################################################################################################################################################################################################
#PATIENT IDs

SC_Nova = read.csv("NovaSeq_RobustSubtypeAssignment_10-11-21_majority.csv") #File generated in previous script

#colnames(SC_Nova) = str_sub(colnames(SC_Nova),6,nchar(colnames(SC_Nova))) 
SC_Nova = SC_Nova[-13:-20,]
rownames(SC_Nova) = SC_Nova[,1]; SC_Nova = SC_Nova[,-1]

#CGND.HRA.00290 and CGND.HRA.00042 are genuine duplicates (same patient, same tissue, different sequencing platforms)
#CGND.HRA.00290 is dropped so that all patient samples are analyzed on the same sequencing platform
#which(colnames(SC_Nova) == "CGND.HRA.00290")
#SC_Nova = SC_Nova[,-49]

#Spinal Cord Pie Chart
tmp = t(SC_Nova)
pie(table(tmp[,12]),col=c("goldenrod1","navy","firebrick"),main="NovaSeq Spinal Cord Labels")

SC_Hi = read.csv("HiSeq_RobustSubtypeAssignment_10-11-21_majority.csv") #File generated in previous script

#colnames(SC_Hi) = str_sub(colnames(SC_Hi),6,nchar(colnames(SC_Hi)))
SC_Hi = SC_Hi[-13,]
rownames(SC_Hi) = SC_Hi[,1]; SC_Hi = SC_Hi[,-1]

#Spinal Cord Pie Chart
tmp2 = t(SC_Hi)
pie(table(tmp2[,12]),col=c("goldenrod1","navy","firebrick"),main="HiSeq Spinal Cord Labels")

#Spinal Cord Pie Chart
pie(table(c(tmp[,12],tmp2[,12])),col=c("goldenrod1","navy","firebrick"),main="Combined Platform Spinal Cord Labels")
