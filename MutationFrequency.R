#Spinal Cord Mutation Frequency

#######Both files are available in Supplemental Dataset 2

#Run once
# MetaData = read.csv("G:/SpinalCord/Publication/Manuscript/Tables/Publicationv2/PatientPhenotype.csv")
# MetaData$Mutation = NA
# 
# samplepheno = read.csv("G:/SpinalCord/Publication/Manuscript/Tables/Publicationv2/SampleLevelPhenotype.csv")
# 
# for(i in 1:nrow(MetaData)){
#   
#   ind = which(samplepheno$ExternalSubjectId == MetaData$Patient[i])
#   
#   if(length(ind)>1){
#     ind = ind[1]
#   }
#   
#   MetaData$Mutation[i] = samplepheno$Reported.Genomic.Mutations..from.sites...NOT.in.any.way.associated.with.WGS.data.from.NYGC.[ind]
# }
# 
# setwd("C:/Users/jeshima/Documents/Smith Lab/Spring 2024/Dissertation Figures/Mutations")
# write.csv(MetaData,"Mutation_Metadata.csv")

####################### Known genetic associations (C9orf72,SOD1) at the patient level
library(stringr)
library(ggplot2)

setwd("C:/Users/jeshima/Documents/Smith Lab/Spring 2024/Dissertation Figures/Mutations")
MetaData = read.csv("Mutation_Metadata.csv")

#Reference values
nsubtype = as.numeric(length(table(MetaData$MetaSpinal)))
nmutations = 2 #c9orf72 and sod1
C9_char = as.numeric(nchar("C9orf72"))
SOD_char = as.numeric(nchar("SOD1"))
MySubtypes = names(table(MetaData$MetaSpinal))
nGlia = as.numeric(length(which(MetaData$MetaSpinal == "GLIA")))
nOx = as.numeric(length(which(MetaData$MetaSpinal == "OX")))
nTD = as.numeric(length(which(MetaData$MetaSpinal == "TD")))
nDisc = as.numeric(length(which(MetaData$MetaSpinal == "Discordant")))
npatients = as.numeric(nrow(MetaData))

#use string searching method 

for(i in 1:nsubtype){
  
  currentsub = MySubtypes[i]
  
  #Build container
  tmpmat = data.frame(matrix(0,nmutations,nmutations+1))
  colnames(tmpmat) = c("C9orf72","SOD1","Unknown")
  rownames(tmpmat) = c("Positive","Negative")
  
  
  for(j in 1:npatients){
    
    if(MetaData$MetaSpinal[j] == MySubtypes[i]){
      
      muts = str_split(MetaData$Mutation[j],",") #Manually converted ; to , in excel
      nmuts = length(muts[[1]])
      
      for(k in 1:nmuts){
        
        charlen = nchar(muts[[1]][k])
        tmpstring = muts[[1]][k]
        
        if(muts[[1]][k] == "Unknown"){
          tmpmat[1,3] = tmpmat[1,3]+1
        }
        
        #C9orf72
        for(m in 1:charlen-C9_char+1){
          
          if(str_sub(tmpstring,m,m+C9_char-1) == "C9orf72"){
            
            for(p in 1:charlen){
              if(str_sub(tmpstring,p,p+2) == "Pos" || str_sub(tmpstring,p,p+2) =="pos"){ #Ignores case sensitivity
                tmpmat[1,1] = tmpmat[1,1]+1
              }else if(str_sub(tmpstring,p,p+2) == "Neg" || str_sub(tmpstring,p,p+2) =="neg"){
                tmpmat[2,1] = tmpmat[2,1]+1
              }
            }
            
          }
          
        }
        
        #SOD1
        for(n in 1:charlen-SOD_char+1){
          
          if(str_sub(tmpstring,n,n+SOD_char-1) == "SOD1"){
            
            for(q in 1:charlen){
              
              if(str_sub(tmpstring,q,q+2) == "Pos" || str_sub(tmpstring,q,q+2) == "pos"){
                tmpmat[1,2] = tmpmat[1,2]+1
              }else if(str_sub(tmpstring,q,q+2) == "Neg" || str_sub(tmpstring,q,q+2) =="neg"){
                tmpmat[2,2] = tmpmat[2,2]+1
              }
              
            }
            
          }
          
        }
        
      }
      
    }
    
    
  }
  #Return results
  tmpname = paste(MySubtypes[i],"Mutations",sep="")
  assign(tmpname,tmpmat)
}


#Plot

CleanLabels = c("Discordant","ALS-Glia","ALS-Ox","ALS-TD")

#Generate plotting data frames
plotdatpos = data.frame(matrix(NA,3,4))
colnames(plotdatpos) = MySubtypes
rownames(plotdatpos) = colnames(GLIAMutations)
plotdatpos[,1] = as.numeric(DiscordantMutations[1,])
plotdatpos[,2] = as.numeric(GLIAMutations[1,])
plotdatpos[,3] = as.numeric(OXMutations[1,])
plotdatpos[,4] = as.numeric(TDMutations[1,])

plotdatneg = data.frame(matrix(NA,2,4))
colnames(plotdatneg) = MySubtypes
rownames(plotdatneg) = colnames(GLIAMutations)[1:2]
plotdatneg[,1] = as.numeric(DiscordantMutations[2,1:2])
plotdatneg[,2] = as.numeric(GLIAMutations[2,1:2])
plotdatneg[,3] = as.numeric(OXMutations[2,1:2])
plotdatneg[,4] = as.numeric(TDMutations[2,1:2])

plotdatpos = as.matrix(plotdatpos)
colnames(plotdatpos) = CleanLabels
rownames(plotdatpos) = colnames(GLIAMutations)
plotdatneg = as.matrix(plotdatneg)
colnames(plotdatneg) = CleanLabels
rownames(plotdatneg) = colnames(GLIAMutations)[1:2]

#Positive Mutations
barplot(plotdatpos,beside=T,col = c("palegreen2","plum2","ivory3"),main = "Positive Mutations by Subtype",ylab = "Frequency",ylim=c(0,35))
legend(1,35,legend = c("C9orf72","SOD1","Unknown"),col = c("palegreen2","plum2","ivory3"),pch=15,pt.cex=1.5)
abline(h=0)

#Negative Mutations
barplot(plotdatneg,beside=T,col = c("palegreen2","plum2"),main = "Negative Mutations by Subtype",ylab = "Frequency",ylim=c(0,50))
legend(1,50,legend = c("C9orf72","SOD1"),col = c("palegreen2","plum2"),pch=15,pt.cex=1.5)
abline(h=0)

#Discordant Mutations
DiscordantMutations = as.matrix(DiscordantMutations)
customcol = c(rep(c("palegreen2","plum2"),nmutations),"gray50")
barplot(DiscordantMutations,beside=T,col = customcol,main = "Discordant Mutations",ylab = "Frequency",ylim=c(0,20),xaxt = "n")
axis(1,c(2,5,7.5),labels = c("C9orf72","SOD1","Unknown"))
legend(7,20,legend = c("Positive","Negative","Unknown"),col = c("palegreen2","plum2","gray50"),pch=15,pt.cex=1.5)
abline(h=0)


#GLIA Mutations
GLIAMutations = as.matrix(GLIAMutations)
customcol = c(rep(c("palegreen2","plum2"),nmutations),"gray50")
barplot(GLIAMutations,beside=T,col = customcol,main = "ALS-Glia Mutations",ylab = "Frequency",ylim=c(0,30),xaxt = "n")
axis(1,c(2,5,7.5),labels = c("C9orf72","SOD1","Unknown"))
legend(7,30,legend = c("Positive","Negative","Unknown"),col = c("palegreen2","plum2","gray50"),pch=15,pt.cex=1.5)
abline(h=0)


#OX Mutations
OXMutations = as.matrix(OXMutations)
customcol = c(rep(c("palegreen2","plum2"),nmutations),"gray50")
barplot(OXMutations,beside=T,col = customcol,main = "ALS-Ox Mutations",ylab = "Frequency",ylim=c(0,50),xaxt = "n")
axis(1,c(2,5,7.5),labels = c("C9orf72","SOD1","Unknown"))
legend(7,50,legend = c("Positive","Negative","Unknown"),col = c("palegreen2","plum2","gray50"),pch=15,pt.cex=1.5)
abline(h=0)


#TD Mutations
TDMutations = as.matrix(TDMutations)
customcol = c(rep(c("palegreen2","plum2"),nmutations),"gray50")
barplot(TDMutations,beside=T,col = customcol,main = "ALS-TD Mutations",ylab = "Frequency",ylim=c(0,30),xaxt = "n")
axis(1,c(2,5,7.5),labels = c("C9orf72","SOD1","Unknown"))
legend(1,30,legend = c("Positive","Negative","Unknown"),col = c("palegreen2","plum2","gray50"),pch=15,pt.cex=1.5)
abline(h=0)


#Manually clean the matrices
GLIAMutations2 = GLIAMutations
tmp = rep(NA,ncol(GLIAMutations))
GLIAMutations2 = rbind(GLIAMutations2,tmp)
rownames(GLIAMutations2) = c("Positive","Negative","Unknown")
GLIAMutations2[3,3] = GLIAMutations2[1,3]
GLIAMutations2[1,3] = NA
GLIAMutations2[2,3] = NA

DiscordantMutations2 = DiscordantMutations
tmp = rep(NA,ncol(DiscordantMutations))
DiscordantMutations2 = rbind(DiscordantMutations2,tmp)
rownames(DiscordantMutations2) = c("Positive","Negative","Unknown")
DiscordantMutations2[3,3] = DiscordantMutations2[1,3]
DiscordantMutations2[1,3] = NA
DiscordantMutations2[2,3] = NA

OXMutations2 = OXMutations
tmp = rep(NA,ncol(OXMutations))
OXMutations2 = rbind(OXMutations2,tmp)
rownames(OXMutations2) = c("Positive","Negative","Unknown")
OXMutations2[3,3] = OXMutations2[1,3]
OXMutations2[1,3] = NA
OXMutations2[2,3] = NA

TDMutations2 = TDMutations
tmp = rep(NA,ncol(TDMutations))
TDMutations2 = rbind(TDMutations2,tmp)
rownames(TDMutations2) = c("Positive","Negative","Unknown")
TDMutations2[3,3] = TDMutations2[1,3]
TDMutations2[1,3] = NA
TDMutations2[2,3] = NA



#Grouped and stacked bar chart

#Create data frame

FullBar = data.frame(matrix(NA,20,4))
colnames(FullBar) = c("Facet","Group","Mutation","Value")
FullBar$Facet = c(rep("Discordant",5),rep("GLIA",5),rep("OX",5),rep("TD",5))
FullBar$Group = rep(c(rep(c("Positive","Negative"),2),"Unknown"),nsubtype)
FullBar$Mutation = rep(c(rep("C9orf72",2),rep("SOD1",2),"Unknown"),nsubtype)

#Fill data frame
for(i in 1:nsubtype){
  
  tmpdat = get(paste(MySubtypes[i],"Mutations2",sep=""))
  
  for(j in 1:nrow(tmpdat)){
    for(k in 1:ncol(tmpdat)){
      for(l in 1:nrow(FullBar)){
        if(FullBar$Facet[l] == MySubtypes[i]){
          if(rownames(tmpdat)[j] == FullBar$Group[l] && colnames(tmpdat)[k] == FullBar$Mutation[l]){
            FullBar$Value[l] = tmpdat[j,k]
          }
          
        }
      }
      
    }
    
  }
  
}

FullBar$Facet = c(rep("Discordant",5),rep("ALS-Glia",5),rep("ALS-Ox",5),rep("ALS-TD",5))

#mycol = rep(c(rep("lightsalmon2",2),rep("palegreen3",2),"thistle3"),nsubtype)

p = ggplot(FullBar,aes(x=Group,y=Value,fill=Mutation)) + geom_bar(stat = "identity",position = "stack") + facet_grid(~Facet)
p = p+scale_fill_manual(values = c("lightsalmon2","palegreen3","thistle3"))
p = p +ggtitle("Spinal Cord C9orf72 and SOD1 Mutations by Subtype") + xlab("") + ylab("Frequency")
p = p+theme(axis.text = element_text(size=12), axis.title = element_text(size=14),plot.title = element_text(size=22))
p = p+theme(axis.title.y=element_text(angle=90, vjust=2,size = 18))
p = p+theme(plot.title = element_text(hjust = 0.5))
p = p+theme(text = element_text(size=14))
p = p+theme(plot.margin = margin(t=10,r=10,b=10,l=10))
p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "gray80"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray80"))
p = p+theme(axis.text.x = element_text(size = 13))
p = p+theme(axis.text.y = element_text(size = 16))
p = p+theme(legend.text = element_text(size = 16))
p = p+theme(legend.title = element_text(size = 16))
p


#################chi-squared stats

#Contingency Tables (mutation-wise)

CSmat = data.frame(matrix(NA,nsubtype,nmutations+2))
colnames(CSmat) = c("Subtype","Positive_C9orf72","Negative_C9orf72","Unknown")
CSmat$Subtype = c("Discordant","ALS-Glia","ALS-Ox","ALS-TD")

CSmatSOD = data.frame(matrix(NA,nsubtype,nmutations+2))
colnames(CSmatSOD) = c("Subtype","Positive_SOD1","Negative_SOD1","Unknown")
CSmatSOD$Subtype = c("Discordant","ALS-Glia","ALS-Ox","ALS-TD")

for(i in 1:nrow(CSmat)){
  
  for(j in 1:nrow(FullBar)){
    
    if(CSmat$Subtype[i] == FullBar$Facet[j] && FullBar$Group[j] == "Positive" && FullBar$Mutation[j] == "C9orf72"){
      CSmat$Positive_C9orf72[i] = FullBar$Value[j]
    }else if(CSmat$Subtype[i] == FullBar$Facet[j] && FullBar$Group[j] == "Negative" && FullBar$Mutation[j] == "C9orf72"){
      CSmat$Negative_C9orf72[i] = FullBar$Value[j]
    }else if(CSmat$Subtype[i] == FullBar$Facet[j] && FullBar$Group[j] == "Unknown" && FullBar$Mutation[j] == "Unknown"){
      CSmat$Unknown[i] = FullBar$Value[j]
    }
    
    
    if(CSmatSOD$Subtype[i] == FullBar$Facet[j] && FullBar$Group[j] == "Positive" && FullBar$Mutation[j] == "SOD1"){
      CSmatSOD$Positive_SOD1[i] = FullBar$Value[j]
    }else if(CSmatSOD$Subtype[i] == FullBar$Facet[j] && FullBar$Group[j] == "Negative" && FullBar$Mutation[j] == "SOD1"){
      CSmatSOD$Negative_SOD1[i] = FullBar$Value[j]
    }else if(CSmatSOD$Subtype[i] == FullBar$Facet[j] && FullBar$Group[j] == "Unknown" && FullBar$Mutation[j] == "Unknown"){
      CSmatSOD$Unknown[i] = FullBar$Value[j]
    }
    
  }
  
}

#Clean up
rownames(CSmat) = CSmat[,1]
CSmat = CSmat[,-1]
rownames(CSmatSOD) = CSmatSOD[,1]
CSmatSOD = CSmatSOD[,-1]

#With unknown category
chisq.test(CSmat) #p = 0.006119 (significance driven by unknown category, see results below)
chisq.test(CSmatSOD) #p = 0.483

#Without unknown category
CSmat2 = CSmat[-which(colnames(CSmat) == "Unknown")]
CSmatSOD2 = CSmatSOD[-which(colnames(CSmatSOD) == "Unknown")]

chisq.test(CSmat2) #p = 0.4726
chisq.test(CSmatSOD2) #p = 0.2148
