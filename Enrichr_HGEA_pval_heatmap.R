
wd = "G:/SpinalCord/Publication/Enrichr"
setwd(wd)

#Result files downloaded from Enrichr as a CSV table
Glia_Up = read.table("GliaUpreg_Reactome_2022_table.txt",header = T,sep = "\t")
Glia_Down = read.table("GliaDownreg_Reactome_2022_table.txt",header = T,sep = "\t")
Ox_Up = read.table("OxUpreg_Reactome_2022_table.txt",header = T,sep = "\t")
Ox_Down = read.table("OxDownreg_Reactome_2022_table.txt",header = T,sep = "\t")
TD_Up = read.table("TDUpreg_Reactome_2022_table.txt",header = T,sep = "\t")
TD_Down = read.table("TDDownreg_Reactome_2022_table.txt",header = T,sep = "\t")


myenrichmentlist = c("Transmission Across Chemical Synapses R-HSA-112315","Potassium Channels R-HSA-1296071","Sialic Acid Metabolism R-HSA-4085001","Olfactory Signaling Pathway R-HSA-381753",
                     "Synthesis Of Substrates In N-glycan Biosythesis R-HSA-446219","mRNA Splicing R-HSA-72172","Regulation Of Innate Immune Responses To Cytosolic DNA R-HSA-3134975",
                     "Biological Oxidations R-HSA-211859","Adaptive Immune System R-HSA-1280218","Chemokine Receptors Bind Chemokines R-HSA-380108",
                     "MTOR Signaling R-HSA-165159","Complement Cascade R-HSA-166658", "Metabolism Of RNA R-HSA-8953854","SUMOylation R-HSA-2990846",
                     "Peptide Ligand-Binding Receptors R-HSA-375276","SLC Transporter Disorders R-HSA-5619102","Metabolism Of Vitamins And Cofactors R-HSA-196854")


Glia_Up_Filt = Glia_Up[Glia_Up$Term %in% myenrichmentlist,]
Ox_Up_Filt = Ox_Up[Ox_Up$Term %in% myenrichmentlist,]
TD_Up_Filt = TD_Up[TD_Up$Term %in% myenrichmentlist,]
Glia_Down_Filt = Glia_Down[Glia_Down$Term %in% myenrichmentlist,]
Ox_Down_Filt = Ox_Down[Ox_Down$Term %in% myenrichmentlist,]
TD_Down_Filt = TD_Down[TD_Down$Term %in% myenrichmentlist,]


EnrichrHeatmap = data.frame(matrix(NA,length(myenrichmentlist),ncol=6))
rownames(EnrichrHeatmap) = myenrichmentlist
colnames(EnrichrHeatmap) = c("ALS-Glia Up","ALS-Glia Down","ALS-Ox Up","ALS-Ox Down","ALS-TD Up","ALS-TD Down")

for(i in 1:nrow(EnrichrHeatmap)){
  
  ind = which(Glia_Up_Filt$Term == rownames(EnrichrHeatmap)[i])
  if(length(ind)>0){
    EnrichrHeatmap$`ALS-Glia Up`[i] = -log10(Glia_Up_Filt$Adjusted.P.value[ind])
  }
  
  
  ind2 = which(Glia_Down_Filt$Term == rownames(EnrichrHeatmap)[i])
  if(length(ind2)>0){
    EnrichrHeatmap$`ALS-Glia Down`[i] = log10(Glia_Down_Filt$Adjusted.P.value[ind2])
  }
  
  
  ind3 = which(Ox_Up_Filt$Term == rownames(EnrichrHeatmap)[i])
  if(length(ind3)>0){
    EnrichrHeatmap$`ALS-Ox Up`[i] = -log10(Ox_Up_Filt$Adjusted.P.value[ind3])
  }
  
  
  ind4 = which(Ox_Down_Filt$Term == rownames(EnrichrHeatmap)[i])
  if(length(ind4)>0){
    EnrichrHeatmap$`ALS-Ox Down`[i] = log10(Ox_Down_Filt$Adjusted.P.value[ind4])
  }
  
  
  ind5 = which(TD_Up_Filt$Term == rownames(EnrichrHeatmap)[i])
  if(length(ind5)>0){
    EnrichrHeatmap$`ALS-TD Up`[i] = -log10(TD_Up_Filt$Adjusted.P.value[ind5])
  }
  
  
  ind6 = which(TD_Down_Filt$Term == rownames(EnrichrHeatmap)[i])
  if(length(ind6)>0){
    EnrichrHeatmap$`ALS-TD Down`[i] = log10(TD_Down_Filt$Adjusted.P.value[ind6])
  }
  
  
}


#Cleanup nonsignificant values, set to 0 for plotting purposes
fEnrichrHeatmap = EnrichrHeatmap

for(i in 1:nrow(fEnrichrHeatmap)){
  for(j in 1:ncol(fEnrichrHeatmap)){
    if(is.na(fEnrichrHeatmap[i,j])){
      fEnrichrHeatmap[i,j] = 0
    }else if(fEnrichrHeatmap[i,j] > 0 && fEnrichrHeatmap[i,j] < -log10(0.05)){
      fEnrichrHeatmap[i,j] = 0
    }else if(fEnrichrHeatmap[i,j] < 0 && fEnrichrHeatmap[i,j] > log10(0.05)){
      fEnrichrHeatmap[i,j] = 0
    }
  }
}

#Combine negative and positive enrichment
FinalHeatmap = data.frame(matrix(NA,nrow(fEnrichrHeatmap),ncol=3))
rownames(FinalHeatmap) = rownames(fEnrichrHeatmap); colnames(FinalHeatmap) = c("ALS-Glia","ALS-Ox","ALS-TD")

for(i in 1:nrow(fEnrichrHeatmap)){
  
  ind = which(fEnrichrHeatmap[i,1:2] != 0)
  if(length(ind)>0){
    FinalHeatmap$`ALS-Glia`[i] = fEnrichrHeatmap[i,ind]
  }else{
    FinalHeatmap$`ALS-Glia`[i] = 0
  }
  
  ind2 = which(fEnrichrHeatmap[i,3:4] != 0)
  if(length(ind2)>0){
    FinalHeatmap$`ALS-Ox`[i] = fEnrichrHeatmap[i,ind2+2]
  }else{
    FinalHeatmap$`ALS-Ox`[i] = 0
  }
  
  ind3 = which(fEnrichrHeatmap[i,5:6] != 0)
  if(length(ind3)>0){
    FinalHeatmap$`ALS-TD`[i] = fEnrichrHeatmap[i,ind3+4]
  }else{
    FinalHeatmap$`ALS-TD`[i] = 0
  }
  
}


typeof(FinalHeatmap)
FinalHeatmap2 = matrix(as.numeric(unlist(FinalHeatmap)),nrow(FinalHeatmap),ncol(FinalHeatmap))
rownames(FinalHeatmap2) = rownames(FinalHeatmap); colnames(FinalHeatmap2) = colnames(FinalHeatmap)

MatrixtoHeatmap3(FinalHeatmap2,samplesarecols = T,featuresasrows = T,limits = c(-10,5),customfeats = rev(myenrichmentlist),customsamps = colnames(FinalHeatmap2))
