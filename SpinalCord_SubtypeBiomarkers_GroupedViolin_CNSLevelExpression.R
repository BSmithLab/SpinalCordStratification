
## Grouped Violin by CNS and Subtype
#Written by: Jarrett Eshima
#Date: October 23, 2023
#Smith Research Lab

library(ggplot2)

###########################################  LOAD DATA  ########################################################################
#Cortex Expression
load("G:/SpinalCord/Publication/DifferentialExpression/CNSPlots/ALSPatientStratification_UnivariateDatasets_RINSite_PeerReview.RData") #Available at: https://figshare.com/authors/Jarrett_Eshima/13813720

#Save Cortex Level Expression
Cortex_GO = filt.GO.sig
Cortex_GT = filt.GT.sig
Cortex_GC = filt.glia.sig
Cortex_TO = filt.TO.sig
Cortex_OC = filt.ox.sig
Cortex_TC = filt.TE.sig
Cortex_Gond = filt.glia.sig.ond
Cortex_Oond = filt.ox.sig.ond
Cortex_Tond = filt.TE.sig.ond
Cortex_Cond = filt.COND.sig

CortexExpr = NormCounts
CortexRaw = rCountData_rinsite
CortexPheno = FullPheno_sr


#Spinal Cord Expression
load("G:/SpinalCord/Publication/DifferentialExpression/4Covar/SpinalCord_DifferentialExpression_FULLCohort_Top5000_10-18-23.RData") #Available at: https://figshare.com/authors/Jarrett_Eshima/13813720

Spinal_GO = filt.GO.sig2
Spinal_GT = filt.GT.sig2
Spinal_GC = filt.glia.sig2
Spinal_TO = filt.TO.sig2
Spinal_OC = filt.ox.sig2
Spinal_TC = filt.TE.sig2

SpinalExpr = NormCounts
SpinalRaw = rCountData_rinsite
SpinalPheno = FullPheno_sr

###########################################  PREP  ########################################################################

### For plotting purposes (not used during statistical analysis), adjust zero count genes to 1

for(i in 1:nrow(CortexExpr)){
  for(j in 1:ncol(CortexExpr)){
    if(CortexExpr[i,j] == 0){
      CortexExpr[i,j] = 1
    }
  }
}

for(i in 1:nrow(SpinalExpr)){
  for(j in 1:ncol(SpinalExpr)){
    if(SpinalExpr[i,j] == 0){
      SpinalExpr[i,j] = 1
    }
  }
}

### Clean up naming

CortexSubtype = CortexPheno$Subtype

for(i in 1:length(CortexSubtype)){
  if(CortexSubtype[i] == "GLIA"){
    CortexSubtype[i] = "Cortex: ALS-Glia"
  }else if(CortexSubtype[i] == "OX"){
    CortexSubtype[i] = "Cortex: ALS-Ox"
  }else if(CortexSubtype[i] == "TE"){
    CortexSubtype[i] = "Cortex: ALS-TD"
  }else if(CortexSubtype[i] =="Control"){
    CortexSubtype[i] = "Cortex: Control"
  }else{
    CortexSubtype[i] = "Cortex: FTLD"
  }
}

#CortexSubtype = factor(CortexSubtype,levels = c("Control","FTLD","ALS-Glia","ALS-Ox","ALS-TD"))

SpinalSubtype = SpinalPheno$Subtype

for(i in 1:length(SpinalSubtype)){
  if(SpinalSubtype[i] == "GLIA"){
    SpinalSubtype[i] = "Spinal: ALS-Glia"
  }else if(SpinalSubtype[i] == "OX"){
    SpinalSubtype[i] = "Spinal: ALS-Ox"
  }else if(SpinalSubtype[i] == "TD"){
    SpinalSubtype[i] = "Spinal: ALS-TD"
  }else{
    SpinalSubtype[i] = "Spinal: Control"
  }
}

#SpinalSubtype = factor(SpinalSubtype,levels = c("Control","ALS-Glia","ALS-Ox","ALS-TD"))

Subtype = c(CortexSubtype,SpinalSubtype)
Subtype = factor(Subtype,levels = c("Cortex: Control","Cortex: FTLD","Cortex: ALS-Glia","Cortex: ALS-Ox","Cortex: ALS-TD","Spinal: Control","Spinal: ALS-Glia","Spinal: ALS-Ox","Spinal: ALS-TD"))


################################### CNS LEVEL VIOLIN PLOTS #####################################################################

#Violin at 3 levels, Cortex, Spinal Cord, Spinal Cord with Cortex patient labels
TDfeats = c("NKX6-2","TARDBP")
OXfeats = c("B4GALT6","GABRA1","GAD2","GLRA3","HTR2A","PCSK1","SLC17A6","UBQLN2")
Gfeats = c("MYL9","ST6GALNAC2","TAGLN")


PlotGene = "TARDBP"
#Pick the reference subtype to add DE p-values contrasted with other groups 
#Options: Glia, Ox, TD (case sensitive; note that 'All' option automatically selects the reference group... all comparisons are NOT performed) - this code does not support FTLD as the reference level
Focus = ""

#Looped Autoplot
for(i in 1:length(PlotGene)){
  #Cortex
  cortexindex = which(rownames(CortexExpr) == PlotGene[i])
  PlotCounts_cortex = CortexExpr[cortexindex,]
  LogPlotCounts_cortex = log2(PlotCounts_cortex)
  
  #Spinal Cord
  spinalindex = which(rownames(SpinalExpr) == PlotGene[i])
  PlotCounts_spinal = SpinalExpr[spinalindex,]
  LogPlotCounts_spinal = log2(PlotCounts_spinal)
  
  LogPlotCounts = c(LogPlotCounts_cortex,LogPlotCounts_spinal)
  
  #Generate Plot
  logdat = data.frame(Subtype,LogPlotCounts)
  p = ggplot(logdat,aes(x=Subtype,y=LogPlotCounts,fill=Subtype)) + geom_violin()
  ttl = paste(PlotGene[i])
  p = p +ggtitle(ttl) + xlab("") + ylab("log2 Median-of-Ratios Counts")
  p = p+scale_x_discrete(limits=c("Cortex: Control","Cortex: FTLD","Cortex: ALS-Glia","Cortex: ALS-Ox","Cortex: ALS-TD","Spinal: Control","Spinal: ALS-Glia","Spinal: ALS-Ox","Spinal: ALS-TD"))
  p = p+scale_fill_manual(values = c("gray50","gray15","goldenrod1","navy","firebrick","gray50","goldenrod1","navy","firebrick"))
  p = p+geom_dotplot(binaxis = 'y',stackdir = 'center',dotsize = 0.25,fill="white",stackratio = 1,binwidth = 0.3)
  #p = p+geom_beeswarm(cex=0.2,color = "white")
  p = p+theme(axis.text = element_text(size=14), axis.title = element_text(size=14),plot.title = element_text(size=24))
  p = p+theme(legend.position = "none")
  p = p+theme(axis.title.x=element_text(vjust=-2))
  p = p+theme(axis.title.y=element_text(angle=90, vjust=6))
  p = p+theme(plot.margin = unit(c(1,1,1,1), "cm"))
  p = p+theme(plot.title = element_text(hjust = 0.5))
  p = p+theme(panel.background = element_rect(fill = 'white', color = 'white'),panel.grid.major = element_line(color = 'gray75'),panel.grid.minor = element_line(color = 'gray75'))
  #p = p+theme_bw()
  
  #Set Limits for p-values
  upperlim = max(LogPlotCounts)+max(LogPlotCounts)*0.65
  lowerlim = 0
  
  #For TARDBP
  #upperlim = max(LogPlotCounts)+max(LogPlotCounts)*0.1
  #lowerlim = 9
  p = p+ylim(lowerlim,upperlim)
  
  #Adaptive lower limit
  # tmplim = min(LogPlotCounts)-min(LogPlotCounts)*0.15
  # if(round(tmplim,0)==0){
  #   lowerlim = 0
  # }else{
  #   lowerlim = tmplim
  # }
  #Add p-values
  if(Focus == "Glia"){
    
    ########## CORTEX
    
    Controlp = Cortex_GC #vs Control
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.1,xend=3,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.1,xend=1,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.05),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.1,xend=3,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=2,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.15)
    }
    
    ONDp = Cortex_Gond #vs Control
    tmpindex3 = which(rownames(ONDp) == PlotGene[i])
    ONDp = ONDp$padj[tmpindex3]
    ONDp = formatC(ONDp,format = "e",digits = 2)
    if(ONDp>""){
      p = p+geom_segment(aes(x=2,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.25,xend=3,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.25),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.25,xend=2,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.2),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.25,xend=3,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.2),size=0.8)
      pstat = paste("FDR p-value:",ONDp)
      p = p+annotate("text",label=pstat,x=2.5,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.30)
    }
    
    OXp = Cortex_GO #vs OX
    tmpindex4 = which(rownames(OXp) == PlotGene[i])
    OXp = OXp$padj[tmpindex4]
    OXp = formatC(OXp,format = "e",digits = 2)
    if(OXp>""){
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.4,xend=4,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.4),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.4,xend=3,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.35),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.4,xend=4,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.35),size=0.8)
      pstat = paste("FDR p-value:",OXp)
      p = p+annotate("text",label=pstat,x=3.5,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.45)
    }
    
    TEp = Cortex_GT #vs TE
    tmpindex5 = which(rownames(TEp) == PlotGene[i])
    TEp = TEp$padj[tmpindex5]
    TEp = formatC(TEp,format = "e",digits = 2)
    if(TEp>""){
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.55,xend=5,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.55),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.55,xend=3,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.5),size=0.8)
      p = p+geom_segment(aes(x=5,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.55,xend=5,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.5),size=0.8)
      pstat = paste("FDR p-value:",TEp)
      p = p+annotate("text",label=pstat,x=4,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.60)
    }
    
    ########## SPINAL CORD
    
    Controlp2 = Spinal_GC #vs Control
    tmpindex2 = which(rownames(Controlp2) == PlotGene[i])
    Controlp2 = Controlp2$padj[tmpindex2]
    Controlp2 = formatC(Controlp2,format = "e",digits = 2)
    if(Controlp2>""){
      p = p+geom_segment(aes(x=6,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.1,xend=7,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.1),size=0.8)
      p = p+geom_segment(aes(x=6,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.1,xend=6,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.05),size=0.8)
      p = p+geom_segment(aes(x=7,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.1,xend=7,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp2)
      p = p+annotate("text",label=pstat,x=6.5,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.15)
    }
    
    OXp2 = Spinal_GO #vs OX
    tmpindex4 = which(rownames(OXp2) == PlotGene[i])
    OXp2 = OXp2$padj[tmpindex4]
    OXp2 = formatC(OXp2,format = "e",digits = 2)
    if(OXp2>""){
      p = p+geom_segment(aes(x=7,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.25,xend=8,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.25),size=0.8)
      p = p+geom_segment(aes(x=7,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.25,xend=7,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.2),size=0.8)
      p = p+geom_segment(aes(x=8,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.25,xend=8,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.2),size=0.8)
      pstat = paste("FDR p-value:",OXp2)
      p = p+annotate("text",label=pstat,x=7.5,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.3)
    }
    
    TEp2 = Spinal_GT #vs TE
    tmpindex5 = which(rownames(TEp2) == PlotGene[i])
    TEp2 = TEp2$padj[tmpindex5]
    TEp2 = formatC(TEp2,format = "e",digits = 2)
    if(TEp2>""){
      p = p+geom_segment(aes(x=7,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.4,xend=9,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.4),size=0.8)
      p = p+geom_segment(aes(x=7,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.4,xend=7,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.35),size=0.8)
      p = p+geom_segment(aes(x=9,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.4,xend=9,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.35),size=0.8)
      pstat = paste("FDR p-value:",TEp2)
      p = p+annotate("text",label=pstat,x=8,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.45)
    }
    
  }else if(Focus == "Ox"){
    
    ########## CORTEX
    
    Controlp = Cortex_OC #vs Control
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.1,xend=4,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.1,xend=1,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.05),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.1,xend=4,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=2.5,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.15)
    }
    
    ONDp = Cortex_Oond #vs Control
    tmpindex3 = which(rownames(ONDp) == PlotGene[i])
    ONDp = ONDp$padj[tmpindex3]
    ONDp = formatC(ONDp,format = "e",digits = 2)
    if(ONDp>""){
      p = p+geom_segment(aes(x=2,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.25,xend=4,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.25),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.25,xend=2,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.2),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.25,xend=4,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.2),size=0.8)
      pstat = paste("FDR p-value:",ONDp)
      p = p+annotate("text",label=pstat,x=3,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.30)
    }
    
    Gliap = Cortex_GO #vs OX
    tmpindex4 = which(rownames(Gliap) == PlotGene[i])
    Gliap = Gliap$padj[tmpindex4]
    Gliap = formatC(Gliap,format = "e",digits = 2)
    if(Gliap>""){
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.4,xend=4,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.4),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.4,xend=3,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.35),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.4,xend=4,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.35),size=0.8)
      pstat = paste("FDR p-value:",Gliap)
      p = p+annotate("text",label=pstat,x=3.5,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.45)
    }
    
    TEp = Cortex_TO #vs TE
    tmpindex5 = which(rownames(TEp) == PlotGene[i])
    TEp = TEp$padj[tmpindex5]
    TEp = formatC(TEp,format = "e",digits = 2)
    if(TEp>""){
      p = p+geom_segment(aes(x=4,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.55,xend=5,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.55),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.55,xend=4,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.5),size=0.8)
      p = p+geom_segment(aes(x=5,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.55,xend=5,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.5),size=0.8)
      pstat = paste("FDR p-value:",TEp)
      p = p+annotate("text",label=pstat,x=4.5,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.6)
    }
    
    ########## SPINAL CORD
    
    Controlp2 = Spinal_OC #vs Control
    tmpindex2 = which(rownames(Controlp2) == PlotGene[i])
    Controlp2 = Controlp2$padj[tmpindex2]
    Controlp2 = formatC(Controlp2,format = "e",digits = 2)
    if(Controlp2>""){
      p = p+geom_segment(aes(x=6,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.1,xend=8,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.1),size=0.8)
      p = p+geom_segment(aes(x=6,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.1,xend=6,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.05),size=0.8)
      p = p+geom_segment(aes(x=8,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.1,xend=8,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp2)
      p = p+annotate("text",label=pstat,x=7,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.15)
    }
    
    Gliap2 = Spinal_GO #vs OX
    tmpindex4 = which(rownames(Gliap2) == PlotGene[i])
    Gliap2 = Gliap2$padj[tmpindex4]
    Gliap2 = formatC(Gliap2,format = "e",digits = 2)
    if(Gliap2>""){
      p = p+geom_segment(aes(x=7,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.25,xend=8,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.25),size=0.8)
      p = p+geom_segment(aes(x=7,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.25,xend=7,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.2),size=0.8)
      p = p+geom_segment(aes(x=8,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.25,xend=8,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.2),size=0.8)
      pstat = paste("FDR p-value:",Gliap2)
      p = p+annotate("text",label=pstat,x=7.5,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.3)
    }
    
    TEp2 = Spinal_TO #vs TE
    tmpindex5 = which(rownames(TEp2) == PlotGene[i])
    TEp2 = TEp2$padj[tmpindex5]
    TEp2 = formatC(TEp2,format = "e",digits = 2)
    if(TEp2>""){
      p = p+geom_segment(aes(x=8,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.4,xend=9,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.4),size=0.8)
      p = p+geom_segment(aes(x=8,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.4,xend=8,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.35),size=0.8)
      p = p+geom_segment(aes(x=9,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.4,xend=9,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.35),size=0.8)
      pstat = paste("FDR p-value:",TEp2)
      p = p+annotate("text",label=pstat,x=8.5,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.45)
    }
    
  }else if(Focus == "TD"){
    
    ########## CORTEX
    
    Controlp = Cortex_TC #vs Control
    tmpindex2 = which(rownames(Controlp) == PlotGene[i])
    Controlp = Controlp$padj[tmpindex2]
    Controlp = formatC(Controlp,format = "e",digits = 2)
    if(Controlp>""){
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.1,xend=5,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.1),size=0.8)
      p = p+geom_segment(aes(x=1,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.1,xend=1,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.05),size=0.8)
      p = p+geom_segment(aes(x=5,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.1,xend=5,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp)
      p = p+annotate("text",label=pstat,x=3,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.15)
    }
    
    ONDp = Cortex_Tond #vs Control
    tmpindex3 = which(rownames(ONDp) == PlotGene[i])
    ONDp = ONDp$padj[tmpindex3]
    ONDp = formatC(ONDp,format = "e",digits = 2)
    if(ONDp>""){
      p = p+geom_segment(aes(x=2,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.25,xend=5,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.25),size=0.8)
      p = p+geom_segment(aes(x=2,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.25,xend=2,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.2),size=0.8)
      p = p+geom_segment(aes(x=5,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.25,xend=5,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.2),size=0.8)
      pstat = paste("FDR p-value:",ONDp)
      p = p+annotate("text",label=pstat,x=3.5,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.30)
    }
    
    Gliap = Cortex_GT #vs Glia
    tmpindex4 = which(rownames(Gliap) == PlotGene[i])
    Gliap = Gliap$padj[tmpindex4]
    Gliap = formatC(Gliap,format = "e",digits = 2)
    if(Gliap>""){
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.4,xend=5,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.4),size=0.8)
      p = p+geom_segment(aes(x=3,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.4,xend=3,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.35),size=0.8)
      p = p+geom_segment(aes(x=5,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.4,xend=5,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.35),size=0.8)
      pstat = paste("FDR p-value:",Gliap)
      p = p+annotate("text",label=pstat,x=4,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.45)
    }
    
    OXp = Cortex_TO #vs OX
    tmpindex5 = which(rownames(OXp) == PlotGene[i])
    OXp = OXp$padj[tmpindex5]
    OXp = formatC(OXp,format = "e",digits = 2)
    if(OXp>""){
      p = p+geom_segment(aes(x=4,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.55,xend=5,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.55),size=0.8)
      p = p+geom_segment(aes(x=4,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.55,xend=4,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.5),size=0.8)
      p = p+geom_segment(aes(x=5,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.55,xend=5,yend=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.5),size=0.8)
      pstat = paste("FDR p-value:",OXp)
      p = p+annotate("text",label=pstat,x=4.5,y=max(LogPlotCounts_cortex)+max(LogPlotCounts_cortex)*0.6)
    }
    
    ########## SPINAL CORD
    
    Controlp2 = Spinal_TC #vs Control
    tmpindex2 = which(rownames(Controlp2) == PlotGene[i])
    Controlp2 = Controlp2$padj[tmpindex2]
    Controlp2 = formatC(Controlp2,format = "e",digits = 2)
    if(Controlp2>""){
      p = p+geom_segment(aes(x=6,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.1,xend=9,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.1),size=0.8)
      p = p+geom_segment(aes(x=6,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.1,xend=6,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.05),size=0.8)
      p = p+geom_segment(aes(x=9,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.1,xend=9,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.05),size=0.8)
      pstat = paste("FDR p-value:",Controlp2)
      p = p+annotate("text",label=pstat,x=7.5,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.15)
    }
    
    Gliap2 = Spinal_GT #vs OX
    tmpindex4 = which(rownames(Gliap2) == PlotGene[i])
    Gliap2 = Gliap2$padj[tmpindex4]
    Gliap2 = formatC(Gliap2,format = "e",digits = 2)
    if(Gliap2>""){
      p = p+geom_segment(aes(x=7,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.25,xend=9,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.25),size=0.8)
      p = p+geom_segment(aes(x=7,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.25,xend=7,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.2),size=0.8)
      p = p+geom_segment(aes(x=9,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.25,xend=9,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.2),size=0.8)
      pstat = paste("FDR p-value:",Gliap2)
      p = p+annotate("text",label=pstat,x=8,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.3)
    }
    
    OXp2 = Spinal_TO #vs TE
    tmpindex5 = which(rownames(OXp2) == PlotGene[i])
    OXp2 = OXp2$padj[tmpindex5]
    OXp2 = formatC(OXp2,format = "e",digits = 2)
    if(OXp2>""){
      p = p+geom_segment(aes(x=8,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.4,xend=9,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.4),size=0.8)
      p = p+geom_segment(aes(x=8,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.4,xend=8,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.35),size=0.8)
      p = p+geom_segment(aes(x=9,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.4,xend=9,yend=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.35),size=0.8)
      pstat = paste("FDR p-value:",OXp2)
      p = p+annotate("text",label=pstat,x=8.5,y=max(LogPlotCounts_spinal)+max(LogPlotCounts_spinal)*0.45)
    }
    
  }else{
    cat("Focus not specified or control reference level provided, p-values not added to the plot")
  }
  print(p)
}
