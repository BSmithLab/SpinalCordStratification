
## Enrichr Prep

wd = "D:/SpinalCord/Publication/Enrichr"
setwd(wd)


Expr = read.csv("SpinalCord_CombinedPlatform_AllSamples_Top5000_MoR_NoTE_NoENSG_9-30-23.csv") #File found at: https://figshare.com/authors/Jarrett_Eshima/13813720
rownames(Expr) = Expr$Gene; Expr = Expr[,-1]

Pheno = read.csv("Samplewise_Subtypes_SuppData_9-27-23.csv") #Files found at: https://figshare.com/authors/Jarrett_Eshima/13813720


Glia = as.character(Pheno$Sample[which(Pheno$Subtype == "GLIA")])
Ox = as.character(Pheno$Sample[which(Pheno$Subtype == "OX")])
TD = as.character(Pheno$Sample[which(Pheno$Subtype == "TD")])
CTR = as.character(Pheno$Sample[which(Pheno$Subtype == "CTR")])
#length(which(Pheno$Subtype != "CTR"))


Glia_Expr = Expr[,colnames(Expr) %in% Glia]
Ox_Expr = Expr[,colnames(Expr) %in% Ox]
TD_Expr = Expr[,colnames(Expr) %in% TD]

Glia_Expr = matrix(as.numeric(unlist(Glia_Expr)),nrow(Glia_Expr),ncol(Glia_Expr))
Ox_Expr = matrix(as.numeric(unlist(Ox_Expr)),nrow(Ox_Expr),ncol(Ox_Expr))
TD_Expr = matrix(as.numeric(unlist(TD_Expr)),nrow(TD_Expr),ncol(TD_Expr))

#table(rownames(Glia_Expr) == rownames(TD_Expr))

GliaUP = GliaDOWN = OxUP = OxDOWN = TDUP = TDDOWN = rep(NA,nrow(Expr))
for(i in 1:nrow(Expr)){
  
  tmpGlia = median(Glia_Expr[i,],na.rm = T)
  tmpOx = median(Ox_Expr[i,],na.rm = T)
  tmpTD = median(TD_Expr[i,],na.rm = T)
  
  medians = c(tmpGlia,tmpOx,tmpTD)
  
  if(length(which(medians == max(medians))) == 1){
    if(tmpGlia == max(medians)){
      GliaUP[i] = rownames(Expr)[i]
    }
    if(tmpOx == max(medians)){
      OxUP[i] = rownames(Expr)[i]
    }
    if(tmpTD == max(medians)){
      TDUP[i] = rownames(Expr)[i]
    }
  }
  
  if(length(which(medians == min(medians))) == 1){
    if(tmpGlia == min(medians)){
      GliaDOWN[i] = rownames(Expr)[i]
    }
    if(tmpOx == min(medians)){
      OxDOWN[i] = rownames(Expr)[i]
    }
    if(tmpTD == min(medians)){
      TDDOWN[i] = rownames(Expr)[i]
    }
  }
  
  
  
  
}


GliaUP = GliaUP[!is.na(GliaUP)]
GliaDOWN = GliaDOWN[!is.na(GliaDOWN)]
OxUP = OxUP[!is.na(OxUP)]
OxDOWN = OxDOWN[!is.na(OxDOWN)]
TDUP = TDUP[!is.na(TDUP)]
TDDOWN = TDDOWN[!is.na(TDDOWN)]


#More than twice (error)
tmp = c(GliaUP,GliaDOWN,OxUP,OxDOWN,TDUP,TDDOWN)
names(table(tmp))[which(table(tmp)>2)]

#Use at input to Enrichr web-app
write.csv(GliaUP,"Enrichr_4Covar_Glia_Upregulated.csv")
write.csv(GliaDOWN,"Enrichr_4Covar_Glia_Downregulated.csv")
write.csv(OxUP,"Enrichr_4Covar_OX_Upregulated.csv")
write.csv(OxDOWN,"Enrichr_4Covar_OX_Downregulated.csv")
write.csv(TDUP,"Enrichr_4Covar_TD_Upregulated.csv")
write.csv(TDDOWN,"Enrichr_4Covar_TD_Downregulated.csv")

