############################################## NOVASEQ ########################################################################################
#At the tissue level (no meta assigned subtype)

#NovaSeq: Without glial markers (Jack Humphrey paper) and RIN-dependent genes
CNSSubtype = read.csv("G:/SpinalCord/Publication/SystemicAnalysis/NoGliaNoRIN/NovaSeq_ConcordancePatientPheno_NoGliaNoRIN_SpinalCord_ALLTISSUE_9-5-23.csv")
CNSSubtype = CNSSubtype[-which(CNSSubtype$CortexSubtype == "Discordant"),]

colnames(CNSSubtype)[c(22,23,24,25)] = c("FrontalCortex","MedialMotorCortex","LateralMotorCortex","UnspecMotorCortex")

########################################## FRONTAL CORTEX ##################################################################################################

concountcerv = concountlumb = concountthor = discountcerv = discountlumb = discountthor = 0
miscountcerv = miscountlumb = miscountthor = 0 #Add in a counter to keep track of unavailable tissue samples
GtoT_cer = GtoO_cer = TtoO_cer = TtoG_cer = OtoT_cer = OtoG_cer = 0
GtoT_tho = GtoO_tho = TtoO_tho = TtoG_tho = OtoT_tho = OtoG_tho = 0
GtoT_lum = GtoO_lum = TtoO_lum = TtoG_lum = OtoT_lum = OtoG_lum = 0
GoodG_cer = GoodO_cer = GoodT_cer = GoodG_tho = GoodO_tho = GoodT_tho = GoodG_lum = GoodO_lum = GoodT_lum = 0
for(i in 1:nrow(CNSSubtype)){
  if(CNSSubtype$CervicalSubtype[i] == "" || CNSSubtype$FrontalCortex[i] == ""){
    miscountcerv = miscountcerv+1
  }else if(CNSSubtype$FrontalCortex[i] == CNSSubtype$CervicalSubtype[i]){
    concountcerv = concountcerv+1
    if(CNSSubtype$CervicalSubtype[i] =="GLIA"){
      GoodG_cer = GoodG_cer+1
    }else if(CNSSubtype$CervicalSubtype[i] == "OX"){
      GoodO_cer = GoodO_cer+1
    }else if(CNSSubtype$CervicalSubtype[i] == "TD"){
      GoodT_cer = GoodT_cer+1
    }
  }else{
    discountcerv = discountcerv+1
    
    if(CNSSubtype$FrontalCortex[i] == "GLIA" && CNSSubtype$CervicalSubtype[i] == "TD"){
      GtoT_cer = GtoT_cer+1
    }else if(CNSSubtype$FrontalCortex[i] == "GLIA" && CNSSubtype$CervicalSubtype[i] == "OX"){
      GtoO_cer = GtoO_cer+1
    }
    
    if(CNSSubtype$FrontalCortex[i] == "TD" && CNSSubtype$CervicalSubtype[i] == "GLIA"){
      TtoG_cer = TtoG_cer+1
    }else if(CNSSubtype$FrontalCortex[i] == "TD" && CNSSubtype$CervicalSubtype[i] == "OX"){
      TtoO_cer = TtoO_cer+1
    }
    
    if(CNSSubtype$FrontalCortex[i] == "OX" && CNSSubtype$CervicalSubtype[i] == "TD"){
      OtoT_cer = OtoT_cer+1
    }else if(CNSSubtype$FrontalCortex[i] == "OX" && CNSSubtype$CervicalSubtype[i] == "GLIA"){
      OtoG_cer = OtoG_cer+1
    }
    
  }
  
  if(CNSSubtype$ThoracicSubtype[i] == "" || CNSSubtype$FrontalCortex[i] == ""){
    miscountthor = miscountthor+1
  }else if(CNSSubtype$FrontalCortex[i] == CNSSubtype$ThoracicSubtype[i]){
    concountthor = concountthor+1
    if(CNSSubtype$ThoracicSubtype[i] =="GLIA"){
      GoodG_tho = GoodG_tho+1
    }else if(CNSSubtype$ThoracicSubtype[i] == "OX"){
      GoodO_tho = GoodO_tho+1
    }else if(CNSSubtype$ThoracicSubtype[i] == "TD"){
      GoodT_tho = GoodT_tho+1
    }
  }else{
    discountthor = discountthor+1
    
    if(CNSSubtype$FrontalCortex[i] == "GLIA" && CNSSubtype$ThoracicSubtype[i] == "TD"){
      GtoT_tho = GtoT_tho+1
    }else if(CNSSubtype$FrontalCortex[i] == "GLIA" && CNSSubtype$ThoracicSubtype[i] == "OX"){
      GtoO_tho = GtoO_tho+1
    }
    
    if(CNSSubtype$FrontalCortex[i] == "TD" && CNSSubtype$ThoracicSubtype[i] == "GLIA"){
      TtoG_tho = TtoG_tho+1
    }else if(CNSSubtype$FrontalCortex[i] == "TD" && CNSSubtype$ThoracicSubtype[i] == "OX"){
      TtoO_tho = TtoO_tho+1
    }
    
    if(CNSSubtype$FrontalCortex[i] == "OX" && CNSSubtype$ThoracicSubtype[i] == "TD"){
      OtoT_tho = OtoT_tho+1
    }else if(CNSSubtype$FrontalCortex[i] == "OX" && CNSSubtype$ThoracicSubtype[i] == "GLIA"){
      OtoG_tho = OtoG_tho+1
    }
    
  }
  
  if(CNSSubtype$LumbarSubtype[i] == "" || CNSSubtype$FrontalCortex[i] == ""){
    miscountlumb = miscountlumb+1
  }else if(CNSSubtype$FrontalCortex[i] == CNSSubtype$LumbarSubtype[i]){
    concountlumb = concountlumb+1
    if(CNSSubtype$LumbarSubtype[i] =="GLIA"){
      GoodG_lum = GoodG_lum+1
    }else if(CNSSubtype$LumbarSubtype[i] == "OX"){
      GoodO_lum = GoodO_lum+1
    }else if(CNSSubtype$LumbarSubtype[i] == "TD"){
      GoodT_lum = GoodT_lum+1
    }
  }else{
    discountlumb = discountlumb+1
    
    if(CNSSubtype$FrontalCortex[i] == "GLIA" && CNSSubtype$LumbarSubtype[i] == "TD"){
      GtoT_lum = GtoT_lum+1
    }else if(CNSSubtype$FrontalCortex[i] == "GLIA" && CNSSubtype$LumbarSubtype[i] == "OX"){
      GtoO_lum = GtoO_lum+1
    }
    
    if(CNSSubtype$FrontalCortex[i] == "TD" && CNSSubtype$LumbarSubtype[i] == "GLIA"){
      TtoG_lum = TtoG_lum+1
    }else if(CNSSubtype$FrontalCortex[i] == "TD" && CNSSubtype$LumbarSubtype[i] == "OX"){
      TtoO_lum = TtoO_lum+1
    }
    
    if(CNSSubtype$FrontalCortex[i] == "OX" && CNSSubtype$LumbarSubtype[i] == "TD"){
      OtoT_lum = OtoT_lum+1
    }else if(CNSSubtype$FrontalCortex[i] == "OX" && CNSSubtype$LumbarSubtype[i] == "GLIA"){
      OtoG_lum = OtoG_lum+1
    }
  }
  
  
}


#NoGlia Angle: 140
#NoGliaNoRin Angle: 110
ConcordantPercent_cerv = concountcerv/sum(concountcerv,discountcerv)
DiscordantPercent_cerv = 1-ConcordantPercent_cerv
pie(c(ConcordantPercent_cerv*100,DiscordantPercent_cerv*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 104,main="Subtype Concordance - Cervical Spinal Cord - FC")


#NoGlia Angle: 129
#NoGliaNoRin Angle: 115
ConcordantPercent_thor = concountthor/sum(concountthor,discountthor)
DiscordantPercent_thor = 1-ConcordantPercent_thor
pie(c(ConcordantPercent_thor*100,DiscordantPercent_thor*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 118,main="Subtype Concordance - Thoracic Spinal Cord - FC")


#NoGlia Angle: 127
#NoGliaNoRin Angle: 98
ConcordantPercent_lumb = concountlumb/sum(concountlumb,discountlumb)
DiscordantPercent_lumb = 1-ConcordantPercent_lumb
pie(c(ConcordantPercent_lumb*100,DiscordantPercent_lumb*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 97,main="Subtype Concordance - Lumbar Spinal Cord - FC")


################################# Cortex vs Cervical Spinal Cord

#ALS-Glia
a = GoodG_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
b = GtoO_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
c = GtoT_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
pie(c(a,b,c),main="ALS-Glia Clustering: Frontal Cortex vs Cervical Spinal Cord",col = c("goldenrod1","chartreuse3","darkorange1"),labels=c("Cortex Concordant","Cortex Discordant: Ox","Cortex Discordant: TD"),init.angle = 90)

#ALS-Ox
a = GoodO_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
b = OtoG_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
c = OtoT_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
pie(c(a,b,c),main="ALS-Ox Clustering: Frontal Cortex vs Cervical Spinal Cord",col = c("navy","chartreuse3","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: TD"),init.angle = 145)

#ALS-TD
a = GoodT_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
b = TtoG_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
c = TtoO_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
pie(c(a,b,c),main="ALS-TD Clustering: Frontal Cortex vs Cervical Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: Ox"),init.angle = 55)

################################# Cortex vs Thoracic Spinal Cord

#ALS-Glia - no ALS-Glia subtypes in the NovaSeq thoracic spinal cord 
# a = GoodG_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
# b = GtoO_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
# c = GtoT_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
# pie(c(a,b,c),main="ALS-Glia Clustering: Frontal Cortex vs Thoracic Spinal Cord",col = c("goldenrod1","chartreuse3","darkorange1"),labels=c("Cortex Concordant","Cortex Discordant: Ox","Cortex Discordant: TD"),init.angle = 0)
# pie(a,main="ALS-Glia Clustering: Frontal Cortex vs Thoracic Spinal Cord",col = "goldenrod1",labels=c("Cortex Concordant"),init.angle = 0)

#ALS-Ox
a = GoodO_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
b = OtoG_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
c = OtoT_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
pie(c(a,b,c),main="ALS-Ox Clustering: Frontal Cortex vs Thoracic Spinal Cord",col = c("navy","chartreuse3","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: TD"),init.angle = 107)

#ALS-TD
a = GoodT_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
b = TtoG_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
c = TtoO_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
pie(c(a,b,c),main="ALS-TD Clustering: Frontal Cortex vs Thoracic Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: Ox"),init.angle = 0)


################################# Cortex vs Lumbar Spinal Cord

#ALS-Glia
a = GoodG_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
b = GtoO_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
c = GtoT_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
pie(c(a,b,c),main="ALS-Glia Clustering: Frontal Cortex vs Lumbar Spinal Cord",col = c("goldenrod1","chartreuse3","darkorange1"),labels=c("Cortex Concordant","Cortex Discordant: Ox","Cortex Discordant: TD"),init.angle = 70)

#ALS-Ox
a = GoodO_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
b = OtoG_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
c = OtoT_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
pie(c(a,b,c),main="ALS-Ox Clustering: Frontal Cortex vs Lumbar Spinal Cord",col = c("navy","chartreuse3","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: TD"),init.angle = 110)

#ALS-TD
a = GoodT_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
b = TtoG_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
c = TtoO_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
pie(c(a,b,c),main="ALS-TD Clustering: Frontal Cortex vs Lumbar Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: Ox"),init.angle = 102)

Frontal_Concord = matrix(c(concountcerv,concountthor,concountlumb,discountcerv,discountthor,discountlumb),ncol=2)
colnames(Frontal_Concord) = c("Concordant","Discordant")
rownames(Frontal_Concord) = c("CervicalSpinal","ThoracicSpinal","LumbarSpinal")

F_C_Re = matrix(c(GoodG_cer,OtoG_cer,TtoG_cer,GtoO_cer,GoodO_cer,TtoO_cer,GtoT_cer,OtoT_cer,GoodT_cer),ncol=3)
colnames(F_C_Re) = c("Glia","Ox","TD"); rownames(F_C_Re) = c("Glia","Ox","TD")
F_T_Re = matrix(c(GoodG_tho,OtoG_tho,TtoG_tho,GtoO_tho,GoodO_tho,TtoO_tho,GtoT_tho,OtoT_tho,GoodT_tho),ncol=3)
colnames(F_T_Re) = c("Glia","Ox","TD"); rownames(F_T_Re) = c("Glia","Ox","TD")
F_L_Re = matrix(c(GoodG_lum,OtoG_lum,TtoG_lum,GtoO_lum,GoodO_lum,TtoO_lum,GtoT_lum,OtoT_lum,GoodT_lum),ncol=3)
colnames(F_L_Re) = c("Glia","Ox","TD"); rownames(F_L_Re) = c("Glia","Ox","TD")


###################################### Lateral MOTOR CORTEX ##################################################################################################

concountcerv = concountlumb = concountthor = discountcerv = discountlumb = discountthor = 0
miscountcerv = miscountlumb = miscountthor = 0 #Add in a counter to keep track of unavailable tissue samples
GtoT_cer = GtoO_cer = TtoO_cer = TtoG_cer = OtoT_cer = OtoG_cer = 0
GtoT_tho = GtoO_tho = TtoO_tho = TtoG_tho = OtoT_tho = OtoG_tho = 0
GtoT_lum = GtoO_lum = TtoO_lum = TtoG_lum = OtoT_lum = OtoG_lum = 0
GoodG_cer = GoodO_cer = GoodT_cer = GoodG_tho = GoodO_tho = GoodT_tho = GoodG_lum = GoodO_lum = GoodT_lum = 0
for(i in 1:nrow(CNSSubtype)){
  if(CNSSubtype$CervicalSubtype[i] == "" || CNSSubtype$LateralMotorCortex[i] == ""){
    miscountcerv = miscountcerv+1
  }else if(CNSSubtype$LateralMotorCortex[i] == CNSSubtype$CervicalSubtype[i]){
    concountcerv = concountcerv+1
    if(CNSSubtype$CervicalSubtype[i] =="GLIA"){
      GoodG_cer = GoodG_cer+1
    }else if(CNSSubtype$CervicalSubtype[i] == "OX"){
      GoodO_cer = GoodO_cer+1
    }else if(CNSSubtype$CervicalSubtype[i] == "TD"){
      GoodT_cer = GoodT_cer+1
    }
  }else{
    discountcerv = discountcerv+1
    
    if(CNSSubtype$LateralMotorCortex[i] == "GLIA" && CNSSubtype$CervicalSubtype[i] == "TD"){
      GtoT_cer = GtoT_cer+1
    }else if(CNSSubtype$LateralMotorCortex[i] == "GLIA" && CNSSubtype$CervicalSubtype[i] == "OX"){
      GtoO_cer = GtoO_cer+1
    }
    
    if(CNSSubtype$LateralMotorCortex[i] == "TD" && CNSSubtype$CervicalSubtype[i] == "GLIA"){
      TtoG_cer = TtoG_cer+1
    }else if(CNSSubtype$LateralMotorCortex[i] == "TD" && CNSSubtype$CervicalSubtype[i] == "OX"){
      TtoO_cer = TtoO_cer+1
    }
    
    if(CNSSubtype$LateralMotorCortex[i] == "OX" && CNSSubtype$CervicalSubtype[i] == "TD"){
      OtoT_cer = OtoT_cer+1
    }else if(CNSSubtype$LateralMotorCortex[i] == "OX" && CNSSubtype$CervicalSubtype[i] == "GLIA"){
      OtoG_cer = OtoG_cer+1
    }
    
  }
  
  if(CNSSubtype$ThoracicSubtype[i] == "" || CNSSubtype$LateralMotorCortex[i] == ""){
    miscountthor = miscountthor+1
  }else if(CNSSubtype$LateralMotorCortex[i] == CNSSubtype$ThoracicSubtype[i]){
    concountthor = concountthor+1
    if(CNSSubtype$ThoracicSubtype[i] =="GLIA"){
      GoodG_tho = GoodG_tho+1
    }else if(CNSSubtype$ThoracicSubtype[i] == "OX"){
      GoodO_tho = GoodO_tho+1
    }else if(CNSSubtype$ThoracicSubtype[i] == "TD"){
      GoodT_tho = GoodT_tho+1
    }
  }else{
    discountthor = discountthor+1
    
    if(CNSSubtype$LateralMotorCortex[i] == "GLIA" && CNSSubtype$ThoracicSubtype[i] == "TD"){
      GtoT_tho = GtoT_tho+1
    }else if(CNSSubtype$LateralMotorCortex[i] == "GLIA" && CNSSubtype$ThoracicSubtype[i] == "OX"){
      GtoO_tho = GtoO_tho+1
    }
    
    if(CNSSubtype$LateralMotorCortex[i] == "TD" && CNSSubtype$ThoracicSubtype[i] == "GLIA"){
      TtoG_tho = TtoG_tho+1
    }else if(CNSSubtype$LateralMotorCortex[i] == "TD" && CNSSubtype$ThoracicSubtype[i] == "OX"){
      TtoO_tho = TtoO_tho+1
    }
    
    if(CNSSubtype$LateralMotorCortex[i] == "OX" && CNSSubtype$ThoracicSubtype[i] == "TD"){
      OtoT_tho = OtoT_tho+1
    }else if(CNSSubtype$LateralMotorCortex[i] == "OX" && CNSSubtype$ThoracicSubtype[i] == "GLIA"){
      OtoG_tho = OtoG_tho+1
    }
    
  }
  
  if(CNSSubtype$LumbarSubtype[i] == "" || CNSSubtype$LateralMotorCortex[i] == ""){
    miscountlumb = miscountlumb+1
  }else if(CNSSubtype$LateralMotorCortex[i] == CNSSubtype$LumbarSubtype[i]){
    concountlumb = concountlumb+1
    if(CNSSubtype$LumbarSubtype[i] =="GLIA"){
      GoodG_lum = GoodG_lum+1
    }else if(CNSSubtype$LumbarSubtype[i] == "OX"){
      GoodO_lum = GoodO_lum+1
    }else if(CNSSubtype$LumbarSubtype[i] == "TD"){
      GoodT_lum = GoodT_lum+1
    }
  }else{
    discountlumb = discountlumb+1
    
    if(CNSSubtype$LateralMotorCortex[i] == "GLIA" && CNSSubtype$LumbarSubtype[i] == "TD"){
      GtoT_lum = GtoT_lum+1
    }else if(CNSSubtype$LateralMotorCortex[i] == "GLIA" && CNSSubtype$LumbarSubtype[i] == "OX"){
      GtoO_lum = GtoO_lum+1
    }
    
    if(CNSSubtype$LateralMotorCortex[i] == "TD" && CNSSubtype$LumbarSubtype[i] == "GLIA"){
      TtoG_lum = TtoG_lum+1
    }else if(CNSSubtype$LateralMotorCortex[i] == "TD" && CNSSubtype$LumbarSubtype[i] == "OX"){
      TtoO_lum = TtoO_lum+1
    }
    
    if(CNSSubtype$LateralMotorCortex[i] == "OX" && CNSSubtype$LumbarSubtype[i] == "TD"){
      OtoT_lum = OtoT_lum+1
    }else if(CNSSubtype$LateralMotorCortex[i] == "OX" && CNSSubtype$LumbarSubtype[i] == "GLIA"){
      OtoG_lum = OtoG_lum+1
    }
  }
  
  
}

#NoGlia Angle: 140
#NoGliaNoRin Angle: 110
ConcordantPercent_cerv = concountcerv/sum(concountcerv,discountcerv)
DiscordantPercent_cerv = 1-ConcordantPercent_cerv
pie(c(ConcordantPercent_cerv*100,DiscordantPercent_cerv*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 113,main="Subtype Concordance - Cervical Spinal Cord - LMC")


#NoGlia Angle: 129
#NoGliaNoRin Angle: 115
ConcordantPercent_thor = concountthor/sum(concountthor,discountthor)
DiscordantPercent_thor = 1-ConcordantPercent_thor
pie(c(ConcordantPercent_thor*100,DiscordantPercent_thor*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 120,main="Subtype Concordance - Thoracic Spinal Cord - LMC")


#NoGlia Angle: 127
#NoGliaNoRin Angle: 98
ConcordantPercent_lumb = concountlumb/sum(concountlumb,discountlumb)
DiscordantPercent_lumb = 1-ConcordantPercent_lumb
pie(c(ConcordantPercent_lumb*100,DiscordantPercent_lumb*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 98,main="Subtype Concordance - Lumbar Spinal Cord - LMC")


################################# Cortex vs Cervical Spinal Cord

#ALS-Glia
a = GoodG_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
b = GtoO_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
c = GtoT_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
pie(c(a,b,c),main="ALS-Glia Clustering: Lateral Motor Cortex vs Cervical Spinal Cord",col = c("goldenrod1","chartreuse3","darkorange1"),labels=c("Cortex Concordant","Cortex Discordant: Ox","Cortex Discordant: TD"),init.angle = 0)

#ALS-Ox
a = GoodO_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
b = OtoG_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
c = OtoT_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
pie(c(a,b,c),main="ALS-Ox Clustering: Lateral Motor Cortex vs Cervical Spinal Cord",col = c("navy","chartreuse3","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: TD"),init.angle = 135)

#ALS-TD
a = GoodT_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
b = TtoG_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
c = TtoO_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
pie(c(a,b,c),main="ALS-TD Clustering: Lateral Motor Cortex vs Cervical Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: Ox"),init.angle = 21)

################################# Cortex vs Thoracic Spinal Cord

#ALS-Glia
a = GoodG_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
b = GtoO_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
c = GtoT_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
pie(c(a,b,c),main="ALS-Glia Clustering: Lateral Motor Cortex vs Thoracic Spinal Cord",col = c("goldenrod1","chartreuse3","darkorange1"),labels=c("Cortex Concordant","Cortex Discordant: Ox","Cortex Discordant: TD"),init.angle = 0)

#ALS-Ox
a = GoodO_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
b = OtoG_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
c = OtoT_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
pie(c(a,b,c),main="ALS-Ox Clustering: Lateral Motor Cortex vs Thoracic Spinal Cord",col = c("navy","chartreuse3","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: TD"),init.angle = 108)

#ALS-TD: No Thoracic TD samples with matching lateral motor cortex subtype
# a = GoodT_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
# b = TtoG_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
# c = TtoO_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
# pie(c(a,b,c),main="ALS-TD Clustering: Lateral Motor Cortex vs Thoracic Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: Ox"),init.angle = 165)


################################# Cortex vs Lumbar Spinal Cord

#ALS-Glia
a = GoodG_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
b = GtoO_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
c = GtoT_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
pie(c(a,b,c),main="ALS-Glia Clustering: Lateral Motor Cortex vs Lumbar Spinal Cord",col = c("goldenrod1","chartreuse3","darkorange1"),labels=c("Cortex Concordant","Cortex Discordant: Ox","Cortex Discordant: TD"),init.angle = 90)

#ALS-Ox
a = GoodO_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
b = OtoG_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
c = OtoT_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
pie(c(a,b,c),main="ALS-Ox Clustering: Lateral Motor Cortex vs Lumbar Spinal Cord",col = c("navy","chartreuse3","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: TD"),init.angle = 100)

#ALS-TD
a = GoodT_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
b = TtoG_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
c = TtoO_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
pie(c(a,b,c),main="ALS-TD Clustering: Lateral Motor Cortex vs Lumbar Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: Ox"),init.angle = 75)


Lateral_Concord = matrix(c(concountcerv,concountthor,concountlumb,discountcerv,discountthor,discountlumb),ncol=2)
colnames(Lateral_Concord) = c("Concordant","Discordant")
rownames(Lateral_Concord) = c("CervicalSpinal","ThoracicSpinal","LumbarSpinal")

LM_C_Re = matrix(c(GoodG_cer,OtoG_cer,TtoG_cer,GtoO_cer,GoodO_cer,TtoO_cer,GtoT_cer,OtoT_cer,GoodT_cer),ncol=3)
colnames(LM_C_Re) = c("Glia","Ox","TD"); rownames(LM_C_Re) = c("Glia","Ox","TD")
LM_T_Re = matrix(c(GoodG_tho,OtoG_tho,TtoG_tho,GtoO_tho,GoodO_tho,TtoO_tho,GtoT_tho,OtoT_tho,GoodT_tho),ncol=3)
colnames(LM_T_Re) = c("Glia","Ox","TD"); rownames(LM_T_Re) = c("Glia","Ox","TD")
LM_L_Re = matrix(c(GoodG_lum,OtoG_lum,TtoG_lum,GtoO_lum,GoodO_lum,TtoO_lum,GtoT_lum,OtoT_lum,GoodT_lum),ncol=3)
colnames(LM_L_Re) = c("Glia","Ox","TD"); rownames(LM_L_Re) = c("Glia","Ox","TD")

################################## MEDIAL MOTOR CORTEX ##################################################################################################

concountcerv = concountlumb = concountthor = discountcerv = discountlumb = discountthor = 0
miscountcerv = miscountlumb = miscountthor = 0 #Add in a counter to keep track of unavailable tissue samples
GtoT_cer = GtoO_cer = TtoO_cer = TtoG_cer = OtoT_cer = OtoG_cer = 0
GtoT_tho = GtoO_tho = TtoO_tho = TtoG_tho = OtoT_tho = OtoG_tho = 0
GtoT_lum = GtoO_lum = TtoO_lum = TtoG_lum = OtoT_lum = OtoG_lum = 0
GoodG_cer = GoodO_cer = GoodT_cer = GoodG_tho = GoodO_tho = GoodT_tho = GoodG_lum = GoodO_lum = GoodT_lum = 0
for(i in 1:nrow(CNSSubtype)){
  if(CNSSubtype$CervicalSubtype[i] == "" || CNSSubtype$MedialMotorCortex[i] == ""){
    miscountcerv = miscountcerv+1
  }else if(CNSSubtype$MedialMotorCortex[i] == CNSSubtype$CervicalSubtype[i]){
    concountcerv = concountcerv+1
    if(CNSSubtype$CervicalSubtype[i] =="GLIA"){
      GoodG_cer = GoodG_cer+1
    }else if(CNSSubtype$CervicalSubtype[i] == "OX"){
      GoodO_cer = GoodO_cer+1
    }else if(CNSSubtype$CervicalSubtype[i] == "TD"){
      GoodT_cer = GoodT_cer+1
    }
  }else{
    discountcerv = discountcerv+1
    
    if(CNSSubtype$MedialMotorCortex[i] == "GLIA" && CNSSubtype$CervicalSubtype[i] == "TD"){
      GtoT_cer = GtoT_cer+1
    }else if(CNSSubtype$MedialMotorCortex[i] == "GLIA" && CNSSubtype$CervicalSubtype[i] == "OX"){
      GtoO_cer = GtoO_cer+1
    }
    
    if(CNSSubtype$MedialMotorCortex[i] == "TD" && CNSSubtype$CervicalSubtype[i] == "GLIA"){
      TtoG_cer = TtoG_cer+1
    }else if(CNSSubtype$MedialMotorCortex[i] == "TD" && CNSSubtype$CervicalSubtype[i] == "OX"){
      TtoO_cer = TtoO_cer+1
    }
    
    if(CNSSubtype$MedialMotorCortex[i] == "OX" && CNSSubtype$CervicalSubtype[i] == "TD"){
      OtoT_cer = OtoT_cer+1
    }else if(CNSSubtype$MedialMotorCortex[i] == "OX" && CNSSubtype$CervicalSubtype[i] == "GLIA"){
      OtoG_cer = OtoG_cer+1
    }
    
  }
  
  if(CNSSubtype$ThoracicSubtype[i] == "" || CNSSubtype$MedialMotorCortex[i] == ""){
    miscountthor = miscountthor+1
  }else if(CNSSubtype$MedialMotorCortex[i] == CNSSubtype$ThoracicSubtype[i]){
    concountthor = concountthor+1
    if(CNSSubtype$ThoracicSubtype[i] =="GLIA"){
      GoodG_tho = GoodG_tho+1
    }else if(CNSSubtype$ThoracicSubtype[i] == "OX"){
      GoodO_tho = GoodO_tho+1
    }else if(CNSSubtype$ThoracicSubtype[i] == "TD"){
      GoodT_tho = GoodT_tho+1
    }
  }else{
    discountthor = discountthor+1
    
    if(CNSSubtype$MedialMotorCortex[i] == "GLIA" && CNSSubtype$ThoracicSubtype[i] == "TD"){
      GtoT_tho = GtoT_tho+1
    }else if(CNSSubtype$MedialMotorCortex[i] == "GLIA" && CNSSubtype$ThoracicSubtype[i] == "OX"){
      GtoO_tho = GtoO_tho+1
    }
    
    if(CNSSubtype$MedialMotorCortex[i] == "TD" && CNSSubtype$ThoracicSubtype[i] == "GLIA"){
      TtoG_tho = TtoG_tho+1
    }else if(CNSSubtype$MedialMotorCortex[i] == "TD" && CNSSubtype$ThoracicSubtype[i] == "OX"){
      TtoO_tho = TtoO_tho+1
    }
    
    if(CNSSubtype$MedialMotorCortex[i] == "OX" && CNSSubtype$ThoracicSubtype[i] == "TD"){
      OtoT_tho = OtoT_tho+1
    }else if(CNSSubtype$MedialMotorCortex[i] == "OX" && CNSSubtype$ThoracicSubtype[i] == "GLIA"){
      OtoG_tho = OtoG_tho+1
    }
    
  }
  
  if(CNSSubtype$LumbarSubtype[i] == "" || CNSSubtype$MedialMotorCortex[i] == ""){
    miscountlumb = miscountlumb+1
  }else if(CNSSubtype$MedialMotorCortex[i] == CNSSubtype$LumbarSubtype[i]){
    concountlumb = concountlumb+1
    if(CNSSubtype$LumbarSubtype[i] =="GLIA"){
      GoodG_lum = GoodG_lum+1
    }else if(CNSSubtype$LumbarSubtype[i] == "OX"){
      GoodO_lum = GoodO_lum+1
    }else if(CNSSubtype$LumbarSubtype[i] == "TD"){
      GoodT_lum = GoodT_lum+1
    }
  }else{
    discountlumb = discountlumb+1
    
    if(CNSSubtype$MedialMotorCortex[i] == "GLIA" && CNSSubtype$LumbarSubtype[i] == "TD"){
      GtoT_lum = GtoT_lum+1
    }else if(CNSSubtype$MedialMotorCortex[i] == "GLIA" && CNSSubtype$LumbarSubtype[i] == "OX"){
      GtoO_lum = GtoO_lum+1
    }
    
    if(CNSSubtype$MedialMotorCortex[i] == "TD" && CNSSubtype$LumbarSubtype[i] == "GLIA"){
      TtoG_lum = TtoG_lum+1
    }else if(CNSSubtype$MedialMotorCortex[i] == "TD" && CNSSubtype$LumbarSubtype[i] == "OX"){
      TtoO_lum = TtoO_lum+1
    }
    
    if(CNSSubtype$MedialMotorCortex[i] == "OX" && CNSSubtype$LumbarSubtype[i] == "TD"){
      OtoT_lum = OtoT_lum+1
    }else if(CNSSubtype$MedialMotorCortex[i] == "OX" && CNSSubtype$LumbarSubtype[i] == "GLIA"){
      OtoG_lum = OtoG_lum+1
    }
  }
  
  
}

#NoGlia Angle: 140
#NoGliaNoRin Angle: 110
ConcordantPercent_cerv = concountcerv/sum(concountcerv,discountcerv)
DiscordantPercent_cerv = 1-ConcordantPercent_cerv
pie(c(ConcordantPercent_cerv*100,DiscordantPercent_cerv*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 114,main="Subtype Concordance - Cervical Spinal Cord - MMC")


#NoGlia Angle: 129
#NoGliaNoRin Angle: 115
ConcordantPercent_thor = concountthor/sum(concountthor,discountthor)
DiscordantPercent_thor = 1-ConcordantPercent_thor
pie(c(ConcordantPercent_thor*100,DiscordantPercent_thor*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 120,main="Subtype Concordance - Thoracic Spinal Cord - MMC")


#NoGlia Angle: 127
#NoGliaNoRin Angle: 98
ConcordantPercent_lumb = concountlumb/sum(concountlumb,discountlumb)
DiscordantPercent_lumb = 1-ConcordantPercent_lumb
pie(c(ConcordantPercent_lumb*100,DiscordantPercent_lumb*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 90,main="Subtype Concordance - Lumbar Spinal Cord - MMC")



################################# Cortex vs Cervical Spinal Cord

#ALS-Glia
a = GoodG_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
b = GtoO_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
c = GtoT_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
pie(c(a,b,c),main="ALS-Glia Clustering: Medial Motor Cortex vs Cervical Spinal Cord",col = c("goldenrod1","chartreuse3","darkorange1"),labels=c("Cortex Concordant","Cortex Discordant: Ox","Cortex Discordant: TD"),init.angle = 0)

#ALS-Ox
a = GoodO_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
b = OtoG_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
c = OtoT_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
pie(c(a,b,c),main="ALS-Ox Clustering: Medial Motor Cortex vs Cervical Spinal Cord",col = c("navy","chartreuse3","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: TD"),init.angle = 140)

#ALS-TD
a = GoodT_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
b = TtoG_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
c = TtoO_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
pie(c(a,b,c),main="ALS-TD Clustering: Medial Motor Cortex vs Cervical Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: Ox"),init.angle = 20)

################################# Cortex vs Thoracic Spinal Cord

#ALS-Glia: NA
# a = GoodG_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
# b = GtoO_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
# c = GtoT_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
# pie(c(a,b,c),main="ALS-Glia Clustering: Medial Motor Cortex vs Thoracic Spinal Cord",col = c("goldenrod1","chartreuse3","darkorange1"),labels=c("Cortex Concordant","Cortex Discordant: Ox","Cortex Discordant: TD"),init.angle = 70)

#ALS-Ox
a = GoodO_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
b = OtoG_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
c = OtoT_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
pie(c(a,b,c),main="ALS-Ox Clustering: Medial Motor Cortex vs Thoracic Spinal Cord",col = c("navy","chartreuse3","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: TD"),init.angle = 108)

#ALS-TD
a = GoodT_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
b = TtoG_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
c = TtoO_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
pie(c(a,b,c),main="ALS-TD Clustering: Medial Motor Cortex vs Thoracic Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: Ox"),init.angle = 0)


################################# Cortex vs Lumbar Spinal Cord

#ALS-Glia
a = GoodG_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
b = GtoO_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
c = GtoT_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
pie(c(a,b,c),main="ALS-Glia Clustering: Medial Motor Cortex vs Lumbar Spinal Cord",col = c("goldenrod1","chartreuse3","darkorange1"),labels=c("Cortex Concordant","Cortex Discordant: Ox","Cortex Discordant: TD"),init.angle = 0)

#ALS-Ox
a = GoodO_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
b = OtoG_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
c = OtoT_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
pie(c(a,b,c),main="ALS-Ox Clustering: Medial Motor Cortex vs Lumbar Spinal Cord",col = c("navy","chartreuse3","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: TD"),init.angle = 95)

#ALS-TD
a = GoodT_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
b = TtoG_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
c = TtoO_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
pie(c(a,b,c),main="ALS-TD Clustering: Medial Motor Cortex vs Lumbar Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: Ox"),init.angle = 50)


Medial_Concord = matrix(c(concountcerv,concountthor,concountlumb,discountcerv,discountthor,discountlumb),ncol=2)
colnames(Medial_Concord) = c("Concordant","Discordant")
rownames(Medial_Concord) = c("CervicalSpinal","ThoracicSpinal","LumbarSpinal")

MM_C_Re = matrix(c(GoodG_cer,OtoG_cer,TtoG_cer,GtoO_cer,GoodO_cer,TtoO_cer,GtoT_cer,OtoT_cer,GoodT_cer),ncol=3)
colnames(MM_C_Re) = c("Glia","Ox","TD"); rownames(MM_C_Re) = c("Glia","Ox","TD")
MM_T_Re = matrix(c(GoodG_tho,OtoG_tho,TtoG_tho,GtoO_tho,GoodO_tho,TtoO_tho,GtoT_tho,OtoT_tho,GoodT_tho),ncol=3)
colnames(MM_T_Re) = c("Glia","Ox","TD"); rownames(MM_T_Re) = c("Glia","Ox","TD")
MM_L_Re = matrix(c(GoodG_lum,OtoG_lum,TtoG_lum,GtoO_lum,GoodO_lum,TtoO_lum,GtoT_lum,OtoT_lum,GoodT_lum),ncol=3)
colnames(MM_L_Re) = c("Glia","Ox","TD"); rownames(MM_L_Re) = c("Glia","Ox","TD")

################################## Unspecified MOTOR CORTEX ##################################################################################################

concountcerv = concountlumb = concountthor = discountcerv = discountlumb = discountthor = 0
miscountcerv = miscountlumb = miscountthor = 0 #Add in a counter to keep track of unavailable tissue samples
GtoT_cer = GtoO_cer = TtoO_cer = TtoG_cer = OtoT_cer = OtoG_cer = 0
GtoT_tho = GtoO_tho = TtoO_tho = TtoG_tho = OtoT_tho = OtoG_tho = 0
GtoT_lum = GtoO_lum = TtoO_lum = TtoG_lum = OtoT_lum = OtoG_lum = 0
GoodG_cer = GoodO_cer = GoodT_cer = GoodG_tho = GoodO_tho = GoodT_tho = GoodG_lum = GoodO_lum = GoodT_lum = 0
for(i in 1:nrow(CNSSubtype)){
  if(CNSSubtype$CervicalSubtype[i] == "" || CNSSubtype$UnspecMotorCortex[i] == ""){
    miscountcerv = miscountcerv+1
  }else if(CNSSubtype$UnspecMotorCortex[i] == CNSSubtype$CervicalSubtype[i]){
    concountcerv = concountcerv+1
    if(CNSSubtype$CervicalSubtype[i] =="GLIA"){
      GoodG_cer = GoodG_cer+1
    }else if(CNSSubtype$CervicalSubtype[i] == "OX"){
      GoodO_cer = GoodO_cer+1
    }else if(CNSSubtype$CervicalSubtype[i] == "TD"){
      GoodT_cer = GoodT_cer+1
    }
  }else{
    discountcerv = discountcerv+1
    
    if(CNSSubtype$UnspecMotorCortex[i] == "GLIA" && CNSSubtype$CervicalSubtype[i] == "TD"){
      GtoT_cer = GtoT_cer+1
    }else if(CNSSubtype$UnspecMotorCortex[i] == "GLIA" && CNSSubtype$CervicalSubtype[i] == "OX"){
      GtoO_cer = GtoO_cer+1
    }
    
    if(CNSSubtype$UnspecMotorCortex[i] == "TD" && CNSSubtype$CervicalSubtype[i] == "GLIA"){
      TtoG_cer = TtoG_cer+1
    }else if(CNSSubtype$UnspecMotorCortex[i] == "TD" && CNSSubtype$CervicalSubtype[i] == "OX"){
      TtoO_cer = TtoO_cer+1
    }
    
    if(CNSSubtype$UnspecMotorCortex[i] == "OX" && CNSSubtype$CervicalSubtype[i] == "TD"){
      OtoT_cer = OtoT_cer+1
    }else if(CNSSubtype$UnspecMotorCortex[i] == "OX" && CNSSubtype$CervicalSubtype[i] == "GLIA"){
      OtoG_cer = OtoG_cer+1
    }
    
  }
  
  if(CNSSubtype$ThoracicSubtype[i] == "" || CNSSubtype$UnspecMotorCortex[i] == ""){
    miscountthor = miscountthor+1
  }else if(CNSSubtype$UnspecMotorCortex[i] == CNSSubtype$ThoracicSubtype[i]){
    concountthor = concountthor+1
    if(CNSSubtype$ThoracicSubtype[i] =="GLIA"){
      GoodG_tho = GoodG_tho+1
    }else if(CNSSubtype$ThoracicSubtype[i] == "OX"){
      GoodO_tho = GoodO_tho+1
    }else if(CNSSubtype$ThoracicSubtype[i] == "TD"){
      GoodT_tho = GoodT_tho+1
    }
  }else{
    discountthor = discountthor+1
    
    if(CNSSubtype$UnspecMotorCortex[i] == "GLIA" && CNSSubtype$ThoracicSubtype[i] == "TD"){
      GtoT_tho = GtoT_tho+1
    }else if(CNSSubtype$UnspecMotorCortex[i] == "GLIA" && CNSSubtype$ThoracicSubtype[i] == "OX"){
      GtoO_tho = GtoO_tho+1
    }
    
    if(CNSSubtype$UnspecMotorCortex[i] == "TD" && CNSSubtype$ThoracicSubtype[i] == "GLIA"){
      TtoG_tho = TtoG_tho+1
    }else if(CNSSubtype$UnspecMotorCortex[i] == "TD" && CNSSubtype$ThoracicSubtype[i] == "OX"){
      TtoO_tho = TtoO_tho+1
    }
    
    if(CNSSubtype$UnspecMotorCortex[i] == "OX" && CNSSubtype$ThoracicSubtype[i] == "TD"){
      OtoT_tho = OtoT_tho+1
    }else if(CNSSubtype$UnspecMotorCortex[i] == "OX" && CNSSubtype$ThoracicSubtype[i] == "GLIA"){
      OtoG_tho = OtoG_tho+1
    }
    
  }
  
  if(CNSSubtype$LumbarSubtype[i] == "" || CNSSubtype$UnspecMotorCortex[i] == ""){
    miscountlumb = miscountlumb+1
  }else if(CNSSubtype$UnspecMotorCortex[i] == CNSSubtype$LumbarSubtype[i]){
    concountlumb = concountlumb+1
    if(CNSSubtype$LumbarSubtype[i] =="GLIA"){
      GoodG_lum = GoodG_lum+1
    }else if(CNSSubtype$LumbarSubtype[i] == "OX"){
      GoodO_lum = GoodO_lum+1
    }else if(CNSSubtype$LumbarSubtype[i] == "TD"){
      GoodT_lum = GoodT_lum+1
    }
  }else{
    discountlumb = discountlumb+1
    
    if(CNSSubtype$UnspecMotorCortex[i] == "GLIA" && CNSSubtype$LumbarSubtype[i] == "TD"){
      GtoT_lum = GtoT_lum+1
    }else if(CNSSubtype$UnspecMotorCortex[i] == "GLIA" && CNSSubtype$LumbarSubtype[i] == "OX"){
      GtoO_lum = GtoO_lum+1
    }
    
    if(CNSSubtype$UnspecMotorCortex[i] == "TD" && CNSSubtype$LumbarSubtype[i] == "GLIA"){
      TtoG_lum = TtoG_lum+1
    }else if(CNSSubtype$UnspecMotorCortex[i] == "TD" && CNSSubtype$LumbarSubtype[i] == "OX"){
      TtoO_lum = TtoO_lum+1
    }
    
    if(CNSSubtype$UnspecMotorCortex[i] == "OX" && CNSSubtype$LumbarSubtype[i] == "TD"){
      OtoT_lum = OtoT_lum+1
    }else if(CNSSubtype$UnspecMotorCortex[i] == "OX" && CNSSubtype$LumbarSubtype[i] == "GLIA"){
      OtoG_lum = OtoG_lum+1
    }
  }
  
  
}

#NoGlia Angle: 140
#NoGliaNoRin Angle: 110
ConcordantPercent_cerv = concountcerv/sum(concountcerv,discountcerv)
DiscordantPercent_cerv = 1-ConcordantPercent_cerv
pie(c(ConcordantPercent_cerv*100,DiscordantPercent_cerv*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 80,main="Subtype Concordance - Cervical Spinal Cord - Unspec. MC")


#No comparisons to be made in NovaSeq cohort
# ConcordantPercent_thor = concountthor/sum(concountthor,discountthor)
# DiscordantPercent_thor = 1-ConcordantPercent_thor
# pie(c(ConcordantPercent_thor*100,DiscordantPercent_thor*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 0,main="Subtype Concordance - Thoracic Spinal Cord - Unspec. MC")


#NoGlia Angle: 127
#NoGliaNoRin Angle: 98
ConcordantPercent_lumb = concountlumb/sum(concountlumb,discountlumb)
DiscordantPercent_lumb = 1-ConcordantPercent_lumb
pie(c(ConcordantPercent_lumb*100,DiscordantPercent_lumb*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 85,main="Subtype Concordance - Lumbar Spinal Cord - Unspec. MC")


################################# Cortex vs Cervical Spinal Cord

#ALS-Glia
a = GoodG_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
b = GtoO_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
c = GtoT_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
pie(c(a,b,c),main="ALS-Glia Clustering: Unspecified Motor Cortex vs Cervical Spinal Cord",col = c("goldenrod1","chartreuse3","darkorange1"),labels=c("Cortex Concordant","Cortex Discordant: Ox","Cortex Discordant: TD"),init.angle = 78)

#ALS-Ox
a = GoodO_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
b = OtoG_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
c = OtoT_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
pie(c(a,b,c),main="ALS-Ox Clustering: Unspecified Motor Cortex vs Cervical Spinal Cord",col = c("navy","chartreuse3","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: TD"),init.angle = 100)

#ALS-TD
a = GoodT_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
b = TtoG_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
c = TtoO_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
pie(c(a,b,c),main="ALS-TD Clustering: Unspecified Motor Cortex vs Cervical Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: Ox"),init.angle = 46)

################################# Cortex vs Thoracic Spinal Cord

#ALS-Glia
# a = GoodG_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
# b = GtoO_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
# c = GtoT_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
# pie(c(a,b,c),main="ALS-Glia Clustering: Unspecified Motor Cortex vs Thoracic Spinal Cord",col = c("goldenrod1","chartreuse3","darkorange1"),labels=c("Cortex Concordant","Cortex Discordant: Ox","Cortex Discordant: TD"),init.angle = 30)

#ALS-Ox
# a = GoodO_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
# b = OtoG_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
# c = OtoT_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
# pie(c(a,b,c),main="ALS-Ox Clustering: Unspecified Motor Cortex vs Thoracic Spinal Cord",col = c("navy","chartreuse3","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: TD"),init.angle = 120)
# pie(b,main="ALS-Ox Clustering: Unspecified Motor Cortex vs Thoracic Spinal Cord",col = c("chartreuse3"),labels="Cortex Discordant: Glia",init.angle = 0)

#ALS-TD
# a = GoodT_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
# b = TtoG_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
# c = TtoO_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
# pie(c(a,b,c),main="ALS-TD Clustering: Unspecified Motor Cortex vs Thoracic Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: Ox"),init.angle = 150)


################################# Cortex vs Lumbar Spinal Cord

#ALS-Glia
a = GoodG_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
b = GtoO_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
c = GtoT_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
pie(c(a,b,c),main="ALS-Glia Clustering: Unspecified Motor Cortex vs Lumbar Spinal Cord",col = c("goldenrod1","chartreuse3","darkorange1"),labels=c("Cortex Concordant","Cortex Discordant: Ox","Cortex Discordant: TD"),init.angle = 75)

#ALS-Ox
a = GoodO_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
b = OtoG_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
c = OtoT_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
pie(c(a,b,c),main="ALS-Ox Clustering: Unspecified Motor Cortex vs Lumbar Spinal Cord",col = c("navy","chartreuse3","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: TD"),init.angle = 103)

#ALS-TD
a = GoodT_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
b = TtoG_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
c = TtoO_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
pie(c(a,b,c),main="ALS-TD Clustering: Unspecified Motor Cortex vs Lumbar Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: Ox"),init.angle = 90)


Unspec_Concord = matrix(c(concountcerv,concountthor,concountlumb,discountcerv,discountthor,discountlumb),ncol=2)
colnames(Medial_Concord) = c("Concordant","Discordant")
rownames(Medial_Concord) = c("CervicalSpinal","ThoracicSpinal","LumbarSpinal")

Un_C_Re = matrix(c(GoodG_cer,OtoG_cer,TtoG_cer,GtoO_cer,GoodO_cer,TtoO_cer,GtoT_cer,OtoT_cer,GoodT_cer),ncol=3)
colnames(Un_C_Re) = c("Glia","Ox","TD"); rownames(Un_C_Re) = c("Glia","Ox","TD")
Un_T_Re = matrix(c(GoodG_tho,OtoG_tho,TtoG_tho,GtoO_tho,GoodO_tho,TtoO_tho,GtoT_tho,OtoT_tho,GoodT_tho),ncol=3)
colnames(Un_T_Re) = c("Glia","Ox","TD"); rownames(Un_T_Re) = c("Glia","Ox","TD")
Un_L_Re = matrix(c(GoodG_lum,OtoG_lum,TtoG_lum,GtoO_lum,GoodO_lum,TtoO_lum,GtoT_lum,OtoT_lum,GoodT_lum),ncol=3)
colnames(Un_L_Re) = c("Glia","Ox","TD"); rownames(Un_L_Re) = c("Glia","Ox","TD")

###########################################################################################################################################################################################################
###########################################################################################################################################################################################################
###########################################################################################################################################################################################################

############################################## HISEQ ########################################################################################
#At the tissue level (no meta assigned subtype)

#HiSeq: Without glial markers (Jack Humphrey paper) and RIN-dependent genes
CNSSubtype = read.csv("G:/SpinalCord/Publication/SystemicAnalysis/NoGliaNoRIN/HiSeq_ConcordancePatientPheno_NoGliaNoRIN_SpinalCord_ALLTISSUE_9-5-23.csv")
CNSSubtype = CNSSubtype[-which(CNSSubtype$CortexSubtype == "Discordant"),]

colnames(CNSSubtype)[c(22,23,24,25)] = c("FrontalCortex","MedialMotorCortex","LateralMotorCortex","UnspecMotorCortex")


########################################## FRONTAL CORTEX ##################################################################################################

concountcerv = concountlumb = concountthor = discountcerv = discountlumb = discountthor = 0
miscountcerv = miscountlumb = miscountthor = 0 #Add in a counter to keep track of unavailable tissue samples
GtoT_cer = GtoO_cer = TtoO_cer = TtoG_cer = OtoT_cer = OtoG_cer = 0
GtoT_tho = GtoO_tho = TtoO_tho = TtoG_tho = OtoT_tho = OtoG_tho = 0
GtoT_lum = GtoO_lum = TtoO_lum = TtoG_lum = OtoT_lum = OtoG_lum = 0
GoodG_cer = GoodO_cer = GoodT_cer = GoodG_tho = GoodO_tho = GoodT_tho = GoodG_lum = GoodO_lum = GoodT_lum = 0
for(i in 1:nrow(CNSSubtype)){
  if(CNSSubtype$CervicalSubtype[i] == "" || CNSSubtype$FrontalCortex[i] == ""){
    miscountcerv = miscountcerv+1
  }else if(CNSSubtype$FrontalCortex[i] == CNSSubtype$CervicalSubtype[i]){
    concountcerv = concountcerv+1
    if(CNSSubtype$CervicalSubtype[i] =="GLIA"){
      GoodG_cer = GoodG_cer+1
    }else if(CNSSubtype$CervicalSubtype[i] == "OX"){
      GoodO_cer = GoodO_cer+1
    }else if(CNSSubtype$CervicalSubtype[i] == "TD"){
      GoodT_cer = GoodT_cer+1
    }
  }else{
    discountcerv = discountcerv+1
    
    if(CNSSubtype$FrontalCortex[i] == "GLIA" && CNSSubtype$CervicalSubtype[i] == "TD"){
      GtoT_cer = GtoT_cer+1
    }else if(CNSSubtype$FrontalCortex[i] == "GLIA" && CNSSubtype$CervicalSubtype[i] == "OX"){
      GtoO_cer = GtoO_cer+1
    }
    
    if(CNSSubtype$FrontalCortex[i] == "TD" && CNSSubtype$CervicalSubtype[i] == "GLIA"){
      TtoG_cer = TtoG_cer+1
    }else if(CNSSubtype$FrontalCortex[i] == "TD" && CNSSubtype$CervicalSubtype[i] == "OX"){
      TtoO_cer = TtoO_cer+1
    }
    
    if(CNSSubtype$FrontalCortex[i] == "OX" && CNSSubtype$CervicalSubtype[i] == "TD"){
      OtoT_cer = OtoT_cer+1
    }else if(CNSSubtype$FrontalCortex[i] == "OX" && CNSSubtype$CervicalSubtype[i] == "GLIA"){
      OtoG_cer = OtoG_cer+1
    }
    
  }
  
  if(CNSSubtype$ThoracicSubtype[i] == "" || CNSSubtype$FrontalCortex[i] == ""){
    miscountthor = miscountthor+1
  }else if(CNSSubtype$FrontalCortex[i] == CNSSubtype$ThoracicSubtype[i]){
    concountthor = concountthor+1
    if(CNSSubtype$ThoracicSubtype[i] =="GLIA"){
      GoodG_tho = GoodG_tho+1
    }else if(CNSSubtype$ThoracicSubtype[i] == "OX"){
      GoodO_tho = GoodO_tho+1
    }else if(CNSSubtype$ThoracicSubtype[i] == "TD"){
      GoodT_tho = GoodT_tho+1
    }
  }else{
    discountthor = discountthor+1
    
    if(CNSSubtype$FrontalCortex[i] == "GLIA" && CNSSubtype$ThoracicSubtype[i] == "TD"){
      GtoT_tho = GtoT_tho+1
    }else if(CNSSubtype$FrontalCortex[i] == "GLIA" && CNSSubtype$ThoracicSubtype[i] == "OX"){
      GtoO_tho = GtoO_tho+1
    }
    
    if(CNSSubtype$FrontalCortex[i] == "TD" && CNSSubtype$ThoracicSubtype[i] == "GLIA"){
      TtoG_tho = TtoG_tho+1
    }else if(CNSSubtype$FrontalCortex[i] == "TD" && CNSSubtype$ThoracicSubtype[i] == "OX"){
      TtoO_tho = TtoO_tho+1
    }
    
    if(CNSSubtype$FrontalCortex[i] == "OX" && CNSSubtype$ThoracicSubtype[i] == "TD"){
      OtoT_tho = OtoT_tho+1
    }else if(CNSSubtype$FrontalCortex[i] == "OX" && CNSSubtype$ThoracicSubtype[i] == "GLIA"){
      OtoG_tho = OtoG_tho+1
    }
    
  }
  
  if(CNSSubtype$LumbarSubtype[i] == "" || CNSSubtype$FrontalCortex[i] == ""){
    miscountlumb = miscountlumb+1
  }else if(CNSSubtype$FrontalCortex[i] == CNSSubtype$LumbarSubtype[i]){
    concountlumb = concountlumb+1
    if(CNSSubtype$LumbarSubtype[i] =="GLIA"){
      GoodG_lum = GoodG_lum+1
    }else if(CNSSubtype$LumbarSubtype[i] == "OX"){
      GoodO_lum = GoodO_lum+1
    }else if(CNSSubtype$LumbarSubtype[i] == "TD"){
      GoodT_lum = GoodT_lum+1
    }
  }else{
    discountlumb = discountlumb+1
    
    if(CNSSubtype$FrontalCortex[i] == "GLIA" && CNSSubtype$LumbarSubtype[i] == "TD"){
      GtoT_lum = GtoT_lum+1
    }else if(CNSSubtype$FrontalCortex[i] == "GLIA" && CNSSubtype$LumbarSubtype[i] == "OX"){
      GtoO_lum = GtoO_lum+1
    }
    
    if(CNSSubtype$FrontalCortex[i] == "TD" && CNSSubtype$LumbarSubtype[i] == "GLIA"){
      TtoG_lum = TtoG_lum+1
    }else if(CNSSubtype$FrontalCortex[i] == "TD" && CNSSubtype$LumbarSubtype[i] == "OX"){
      TtoO_lum = TtoO_lum+1
    }
    
    if(CNSSubtype$FrontalCortex[i] == "OX" && CNSSubtype$LumbarSubtype[i] == "TD"){
      OtoT_lum = OtoT_lum+1
    }else if(CNSSubtype$FrontalCortex[i] == "OX" && CNSSubtype$LumbarSubtype[i] == "GLIA"){
      OtoG_lum = OtoG_lum+1
    }
  }
  
  
}


#NoGlia Angle: 140
#NoGliaNoRin Angle: 110
ConcordantPercent_cerv = concountcerv/sum(concountcerv,discountcerv)
DiscordantPercent_cerv = 1-ConcordantPercent_cerv
pie(c(ConcordantPercent_cerv*100,DiscordantPercent_cerv*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 134,main="Subtype Concordance - Cervical Spinal Cord - FC")


#NoGlia Angle: 129
#NoGliaNoRin Angle: 115
ConcordantPercent_thor = concountthor/sum(concountthor,discountthor)
DiscordantPercent_thor = 1-ConcordantPercent_thor
pie(c(ConcordantPercent_thor*100,DiscordantPercent_thor*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 118,main="Subtype Concordance - Thoracic Spinal Cord - FC")


#NoGlia Angle: 127
#NoGliaNoRin Angle: 98
ConcordantPercent_lumb = concountlumb/sum(concountlumb,discountlumb)
DiscordantPercent_lumb = 1-ConcordantPercent_lumb
pie(c(ConcordantPercent_lumb*100,DiscordantPercent_lumb*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 108,main="Subtype Concordance - Lumbar Spinal Cord - FC")


################################# Cortex vs Cervical Spinal Cord

#ALS-Glia
a = GoodG_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
b = GtoO_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
c = GtoT_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
pie(c(a,b,c),main="ALS-Glia Clustering: Frontal Cortex vs Cervical Spinal Cord",col = c("goldenrod1","chartreuse3","darkorange1"),labels=c("Cortex Concordant","Cortex Discordant: Ox","Cortex Discordant: TD"),init.angle = 135)

#ALS-Ox
a = GoodO_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
b = OtoG_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
c = OtoT_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
pie(c(a,b,c),main="ALS-Ox Clustering: Frontal Cortex vs Cervical Spinal Cord",col = c("navy","chartreuse3","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: TD"),init.angle = 135)

#ALS-TD
a = GoodT_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
b = TtoG_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
c = TtoO_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
pie(c(a,b,c),main="ALS-TD Clustering: Frontal Cortex vs Cervical Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: Ox"),init.angle = 130)

################################# Cortex vs Thoracic Spinal Cord

#ALS-Glia - no ALS-Glia subtypes in the NovaSeq thoracic spinal cord 
a = GoodG_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
b = GtoO_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
c = GtoT_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
pie(c(a,b,c),main="ALS-Glia Clustering: Frontal Cortex vs Thoracic Spinal Cord",col = c("goldenrod1","chartreuse3","darkorange1"),labels=c("Cortex Concordant","Cortex Discordant: Ox","Cortex Discordant: TD"),init.angle = 0)
# pie(a,main="ALS-Glia Clustering: Frontal Cortex vs Thoracic Spinal Cord",col = "goldenrod1",labels=c("Cortex Concordant"),init.angle = 0)

#ALS-Ox
a = GoodO_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
b = OtoG_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
c = OtoT_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
pie(c(a,b,c),main="ALS-Ox Clustering: Frontal Cortex vs Thoracic Spinal Cord",col = c("navy","chartreuse3","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: TD"),init.angle = 130)

#ALS-TD
a = GoodT_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
b = TtoG_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
c = TtoO_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
pie(c(a,b,c),main="ALS-TD Clustering: Frontal Cortex vs Thoracic Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: Ox"),init.angle = 147)


################################# Cortex vs Lumbar Spinal Cord

#ALS-Glia
a = GoodG_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
b = GtoO_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
c = GtoT_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
pie(c(a,b,c),main="ALS-Glia Clustering: Frontal Cortex vs Lumbar Spinal Cord",col = c("goldenrod1","chartreuse3","darkorange1"),labels=c("Cortex Concordant","Cortex Discordant: Ox","Cortex Discordant: TD"),init.angle = 0)

#ALS-Ox
a = GoodO_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
b = OtoG_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
c = OtoT_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
pie(c(a,b,c),main="ALS-Ox Clustering: Frontal Cortex vs Lumbar Spinal Cord",col = c("navy","chartreuse3","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: TD"),init.angle = 110)

#ALS-TD
a = GoodT_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
b = TtoG_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
c = TtoO_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
pie(c(a,b,c),main="ALS-TD Clustering: Frontal Cortex vs Lumbar Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: Ox"),init.angle = 102)

Frontal_Concord = matrix(c(concountcerv,concountthor,concountlumb,discountcerv,discountthor,discountlumb),ncol=2)
colnames(Frontal_Concord) = c("Concordant","Discordant")
rownames(Frontal_Concord) = c("CervicalSpinal","ThoracicSpinal","LumbarSpinal")

F_C_Re = matrix(c(GoodG_cer,OtoG_cer,TtoG_cer,GtoO_cer,GoodO_cer,TtoO_cer,GtoT_cer,OtoT_cer,GoodT_cer),ncol=3)
colnames(F_C_Re) = c("Glia","Ox","TD"); rownames(F_C_Re) = c("Glia","Ox","TD")
F_T_Re = matrix(c(GoodG_tho,OtoG_tho,TtoG_tho,GtoO_tho,GoodO_tho,TtoO_tho,GtoT_tho,OtoT_tho,GoodT_tho),ncol=3)
colnames(F_T_Re) = c("Glia","Ox","TD"); rownames(F_T_Re) = c("Glia","Ox","TD")
F_L_Re = matrix(c(GoodG_lum,OtoG_lum,TtoG_lum,GtoO_lum,GoodO_lum,TtoO_lum,GtoT_lum,OtoT_lum,GoodT_lum),ncol=3)
colnames(F_L_Re) = c("Glia","Ox","TD"); rownames(F_L_Re) = c("Glia","Ox","TD")


###################################### Lateral MOTOR CORTEX ##################################################################################################

concountcerv = concountlumb = concountthor = discountcerv = discountlumb = discountthor = 0
miscountcerv = miscountlumb = miscountthor = 0 #Add in a counter to keep track of unavailable tissue samples
GtoT_cer = GtoO_cer = TtoO_cer = TtoG_cer = OtoT_cer = OtoG_cer = 0
GtoT_tho = GtoO_tho = TtoO_tho = TtoG_tho = OtoT_tho = OtoG_tho = 0
GtoT_lum = GtoO_lum = TtoO_lum = TtoG_lum = OtoT_lum = OtoG_lum = 0
GoodG_cer = GoodO_cer = GoodT_cer = GoodG_tho = GoodO_tho = GoodT_tho = GoodG_lum = GoodO_lum = GoodT_lum = 0
for(i in 1:nrow(CNSSubtype)){
  if(CNSSubtype$CervicalSubtype[i] == "" || CNSSubtype$LateralMotorCortex[i] == ""){
    miscountcerv = miscountcerv+1
  }else if(CNSSubtype$LateralMotorCortex[i] == CNSSubtype$CervicalSubtype[i]){
    concountcerv = concountcerv+1
    if(CNSSubtype$CervicalSubtype[i] =="GLIA"){
      GoodG_cer = GoodG_cer+1
    }else if(CNSSubtype$CervicalSubtype[i] == "OX"){
      GoodO_cer = GoodO_cer+1
    }else if(CNSSubtype$CervicalSubtype[i] == "TD"){
      GoodT_cer = GoodT_cer+1
    }
  }else{
    discountcerv = discountcerv+1
    
    if(CNSSubtype$LateralMotorCortex[i] == "GLIA" && CNSSubtype$CervicalSubtype[i] == "TD"){
      GtoT_cer = GtoT_cer+1
    }else if(CNSSubtype$LateralMotorCortex[i] == "GLIA" && CNSSubtype$CervicalSubtype[i] == "OX"){
      GtoO_cer = GtoO_cer+1
    }
    
    if(CNSSubtype$LateralMotorCortex[i] == "TD" && CNSSubtype$CervicalSubtype[i] == "GLIA"){
      TtoG_cer = TtoG_cer+1
    }else if(CNSSubtype$LateralMotorCortex[i] == "TD" && CNSSubtype$CervicalSubtype[i] == "OX"){
      TtoO_cer = TtoO_cer+1
    }
    
    if(CNSSubtype$LateralMotorCortex[i] == "OX" && CNSSubtype$CervicalSubtype[i] == "TD"){
      OtoT_cer = OtoT_cer+1
    }else if(CNSSubtype$LateralMotorCortex[i] == "OX" && CNSSubtype$CervicalSubtype[i] == "GLIA"){
      OtoG_cer = OtoG_cer+1
    }
    
  }
  
  if(CNSSubtype$ThoracicSubtype[i] == "" || CNSSubtype$LateralMotorCortex[i] == ""){
    miscountthor = miscountthor+1
  }else if(CNSSubtype$LateralMotorCortex[i] == CNSSubtype$ThoracicSubtype[i]){
    concountthor = concountthor+1
    if(CNSSubtype$ThoracicSubtype[i] =="GLIA"){
      GoodG_tho = GoodG_tho+1
    }else if(CNSSubtype$ThoracicSubtype[i] == "OX"){
      GoodO_tho = GoodO_tho+1
    }else if(CNSSubtype$ThoracicSubtype[i] == "TD"){
      GoodT_tho = GoodT_tho+1
    }
  }else{
    discountthor = discountthor+1
    
    if(CNSSubtype$LateralMotorCortex[i] == "GLIA" && CNSSubtype$ThoracicSubtype[i] == "TD"){
      GtoT_tho = GtoT_tho+1
    }else if(CNSSubtype$LateralMotorCortex[i] == "GLIA" && CNSSubtype$ThoracicSubtype[i] == "OX"){
      GtoO_tho = GtoO_tho+1
    }
    
    if(CNSSubtype$LateralMotorCortex[i] == "TD" && CNSSubtype$ThoracicSubtype[i] == "GLIA"){
      TtoG_tho = TtoG_tho+1
    }else if(CNSSubtype$LateralMotorCortex[i] == "TD" && CNSSubtype$ThoracicSubtype[i] == "OX"){
      TtoO_tho = TtoO_tho+1
    }
    
    if(CNSSubtype$LateralMotorCortex[i] == "OX" && CNSSubtype$ThoracicSubtype[i] == "TD"){
      OtoT_tho = OtoT_tho+1
    }else if(CNSSubtype$LateralMotorCortex[i] == "OX" && CNSSubtype$ThoracicSubtype[i] == "GLIA"){
      OtoG_tho = OtoG_tho+1
    }
    
  }
  
  if(CNSSubtype$LumbarSubtype[i] == "" || CNSSubtype$LateralMotorCortex[i] == ""){
    miscountlumb = miscountlumb+1
  }else if(CNSSubtype$LateralMotorCortex[i] == CNSSubtype$LumbarSubtype[i]){
    concountlumb = concountlumb+1
    if(CNSSubtype$LumbarSubtype[i] =="GLIA"){
      GoodG_lum = GoodG_lum+1
    }else if(CNSSubtype$LumbarSubtype[i] == "OX"){
      GoodO_lum = GoodO_lum+1
    }else if(CNSSubtype$LumbarSubtype[i] == "TD"){
      GoodT_lum = GoodT_lum+1
    }
  }else{
    discountlumb = discountlumb+1
    
    if(CNSSubtype$LateralMotorCortex[i] == "GLIA" && CNSSubtype$LumbarSubtype[i] == "TD"){
      GtoT_lum = GtoT_lum+1
    }else if(CNSSubtype$LateralMotorCortex[i] == "GLIA" && CNSSubtype$LumbarSubtype[i] == "OX"){
      GtoO_lum = GtoO_lum+1
    }
    
    if(CNSSubtype$LateralMotorCortex[i] == "TD" && CNSSubtype$LumbarSubtype[i] == "GLIA"){
      TtoG_lum = TtoG_lum+1
    }else if(CNSSubtype$LateralMotorCortex[i] == "TD" && CNSSubtype$LumbarSubtype[i] == "OX"){
      TtoO_lum = TtoO_lum+1
    }
    
    if(CNSSubtype$LateralMotorCortex[i] == "OX" && CNSSubtype$LumbarSubtype[i] == "TD"){
      OtoT_lum = OtoT_lum+1
    }else if(CNSSubtype$LateralMotorCortex[i] == "OX" && CNSSubtype$LumbarSubtype[i] == "GLIA"){
      OtoG_lum = OtoG_lum+1
    }
  }
  
  
}

#NoGlia Angle: 140
#NoGliaNoRin Angle: 110
ConcordantPercent_cerv = concountcerv/sum(concountcerv,discountcerv)
DiscordantPercent_cerv = 1-ConcordantPercent_cerv
pie(c(ConcordantPercent_cerv*100,DiscordantPercent_cerv*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 118,main="Subtype Concordance - Cervical Spinal Cord - LMC")


#NoGlia Angle: 129
#NoGliaNoRin Angle: 115
ConcordantPercent_thor = concountthor/sum(concountthor,discountthor)
DiscordantPercent_thor = 1-ConcordantPercent_thor
pie(c(ConcordantPercent_thor*100,DiscordantPercent_thor*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 124,main="Subtype Concordance - Thoracic Spinal Cord - LMC")


#NoGlia Angle: 127
#NoGliaNoRin Angle: 98
ConcordantPercent_lumb = concountlumb/sum(concountlumb,discountlumb)
DiscordantPercent_lumb = 1-ConcordantPercent_lumb
pie(c(ConcordantPercent_lumb*100,DiscordantPercent_lumb*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 115,main="Subtype Concordance - Lumbar Spinal Cord - LMC")


################################# Cortex vs Cervical Spinal Cord

#ALS-Glia
a = GoodG_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
b = GtoO_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
c = GtoT_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
pie(c(a,b,c),main="ALS-Glia Clustering: Lateral Motor Cortex vs Cervical Spinal Cord",col = c("goldenrod1","chartreuse3","darkorange1"),labels=c("Cortex Concordant","Cortex Discordant: Ox","Cortex Discordant: TD"),init.angle = 80)

#ALS-Ox
a = GoodO_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
b = OtoG_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
c = OtoT_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
pie(c(a,b,c),main="ALS-Ox Clustering: Lateral Motor Cortex vs Cervical Spinal Cord",col = c("navy","chartreuse3","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: TD"),init.angle = 140)

#ALS-TD
a = GoodT_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
b = TtoG_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
c = TtoO_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
pie(c(a,b,c),main="ALS-TD Clustering: Lateral Motor Cortex vs Cervical Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: Ox"),init.angle = 100)

################################# Cortex vs Thoracic Spinal Cord

#ALS-Glia
a = GoodG_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
b = GtoO_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
c = GtoT_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
pie(c(a,b,c),main="ALS-Glia Clustering: Lateral Motor Cortex vs Thoracic Spinal Cord",col = c("goldenrod1","chartreuse3","darkorange1"),labels=c("Cortex Concordant","Cortex Discordant: Ox","Cortex Discordant: TD"),init.angle = 45)

#ALS-Ox
a = GoodO_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
b = OtoG_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
c = OtoT_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
pie(c(a,b,c),main="ALS-Ox Clustering: Lateral Motor Cortex vs Thoracic Spinal Cord",col = c("navy","chartreuse3","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: TD"),init.angle = 128)

#ALS-TD: No Thoracic TD samples with matching lateral motor cortex subtype
a = GoodT_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
b = TtoG_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
c = TtoO_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
pie(c(a,b,c),main="ALS-TD Clustering: Lateral Motor Cortex vs Thoracic Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: Ox"),init.angle = 20)


################################# Cortex vs Lumbar Spinal Cord

#ALS-Glia
a = GoodG_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
b = GtoO_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
c = GtoT_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
pie(c(a,b,c),main="ALS-Glia Clustering: Lateral Motor Cortex vs Lumbar Spinal Cord",col = c("goldenrod1","chartreuse3","darkorange1"),labels=c("Cortex Concordant","Cortex Discordant: Ox","Cortex Discordant: TD"),init.angle = 120)

#ALS-Ox
a = GoodO_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
b = OtoG_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
c = OtoT_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
pie(c(a,b,c),main="ALS-Ox Clustering: Lateral Motor Cortex vs Lumbar Spinal Cord",col = c("navy","chartreuse3","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: TD"),init.angle = 124)

#ALS-TD
a = GoodT_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
b = TtoG_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
c = TtoO_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
pie(c(a,b,c),main="ALS-TD Clustering: Lateral Motor Cortex vs Lumbar Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: Ox"),init.angle = 90)


Lateral_Concord = matrix(c(concountcerv,concountthor,concountlumb,discountcerv,discountthor,discountlumb),ncol=2)
colnames(Lateral_Concord) = c("Concordant","Discordant")
rownames(Lateral_Concord) = c("CervicalSpinal","ThoracicSpinal","LumbarSpinal")

LM_C_Re = matrix(c(GoodG_cer,OtoG_cer,TtoG_cer,GtoO_cer,GoodO_cer,TtoO_cer,GtoT_cer,OtoT_cer,GoodT_cer),ncol=3)
colnames(LM_C_Re) = c("Glia","Ox","TD"); rownames(LM_C_Re) = c("Glia","Ox","TD")
LM_T_Re = matrix(c(GoodG_tho,OtoG_tho,TtoG_tho,GtoO_tho,GoodO_tho,TtoO_tho,GtoT_tho,OtoT_tho,GoodT_tho),ncol=3)
colnames(LM_T_Re) = c("Glia","Ox","TD"); rownames(LM_T_Re) = c("Glia","Ox","TD")
LM_L_Re = matrix(c(GoodG_lum,OtoG_lum,TtoG_lum,GtoO_lum,GoodO_lum,TtoO_lum,GtoT_lum,OtoT_lum,GoodT_lum),ncol=3)
colnames(LM_L_Re) = c("Glia","Ox","TD"); rownames(LM_L_Re) = c("Glia","Ox","TD")

################################## MEDIAL MOTOR CORTEX ##################################################################################################

concountcerv = concountlumb = concountthor = discountcerv = discountlumb = discountthor = 0
miscountcerv = miscountlumb = miscountthor = 0 #Add in a counter to keep track of unavailable tissue samples
GtoT_cer = GtoO_cer = TtoO_cer = TtoG_cer = OtoT_cer = OtoG_cer = 0
GtoT_tho = GtoO_tho = TtoO_tho = TtoG_tho = OtoT_tho = OtoG_tho = 0
GtoT_lum = GtoO_lum = TtoO_lum = TtoG_lum = OtoT_lum = OtoG_lum = 0
GoodG_cer = GoodO_cer = GoodT_cer = GoodG_tho = GoodO_tho = GoodT_tho = GoodG_lum = GoodO_lum = GoodT_lum = 0
for(i in 1:nrow(CNSSubtype)){
  if(CNSSubtype$CervicalSubtype[i] == "" || CNSSubtype$MedialMotorCortex[i] == ""){
    miscountcerv = miscountcerv+1
  }else if(CNSSubtype$MedialMotorCortex[i] == CNSSubtype$CervicalSubtype[i]){
    concountcerv = concountcerv+1
    if(CNSSubtype$CervicalSubtype[i] =="GLIA"){
      GoodG_cer = GoodG_cer+1
    }else if(CNSSubtype$CervicalSubtype[i] == "OX"){
      GoodO_cer = GoodO_cer+1
    }else if(CNSSubtype$CervicalSubtype[i] == "TD"){
      GoodT_cer = GoodT_cer+1
    }
  }else{
    discountcerv = discountcerv+1
    
    if(CNSSubtype$MedialMotorCortex[i] == "GLIA" && CNSSubtype$CervicalSubtype[i] == "TD"){
      GtoT_cer = GtoT_cer+1
    }else if(CNSSubtype$MedialMotorCortex[i] == "GLIA" && CNSSubtype$CervicalSubtype[i] == "OX"){
      GtoO_cer = GtoO_cer+1
    }
    
    if(CNSSubtype$MedialMotorCortex[i] == "TD" && CNSSubtype$CervicalSubtype[i] == "GLIA"){
      TtoG_cer = TtoG_cer+1
    }else if(CNSSubtype$MedialMotorCortex[i] == "TD" && CNSSubtype$CervicalSubtype[i] == "OX"){
      TtoO_cer = TtoO_cer+1
    }
    
    if(CNSSubtype$MedialMotorCortex[i] == "OX" && CNSSubtype$CervicalSubtype[i] == "TD"){
      OtoT_cer = OtoT_cer+1
    }else if(CNSSubtype$MedialMotorCortex[i] == "OX" && CNSSubtype$CervicalSubtype[i] == "GLIA"){
      OtoG_cer = OtoG_cer+1
    }
    
  }
  
  if(CNSSubtype$ThoracicSubtype[i] == "" || CNSSubtype$MedialMotorCortex[i] == ""){
    miscountthor = miscountthor+1
  }else if(CNSSubtype$MedialMotorCortex[i] == CNSSubtype$ThoracicSubtype[i]){
    concountthor = concountthor+1
    if(CNSSubtype$ThoracicSubtype[i] =="GLIA"){
      GoodG_tho = GoodG_tho+1
    }else if(CNSSubtype$ThoracicSubtype[i] == "OX"){
      GoodO_tho = GoodO_tho+1
    }else if(CNSSubtype$ThoracicSubtype[i] == "TD"){
      GoodT_tho = GoodT_tho+1
    }
  }else{
    discountthor = discountthor+1
    
    if(CNSSubtype$MedialMotorCortex[i] == "GLIA" && CNSSubtype$ThoracicSubtype[i] == "TD"){
      GtoT_tho = GtoT_tho+1
    }else if(CNSSubtype$MedialMotorCortex[i] == "GLIA" && CNSSubtype$ThoracicSubtype[i] == "OX"){
      GtoO_tho = GtoO_tho+1
    }
    
    if(CNSSubtype$MedialMotorCortex[i] == "TD" && CNSSubtype$ThoracicSubtype[i] == "GLIA"){
      TtoG_tho = TtoG_tho+1
    }else if(CNSSubtype$MedialMotorCortex[i] == "TD" && CNSSubtype$ThoracicSubtype[i] == "OX"){
      TtoO_tho = TtoO_tho+1
    }
    
    if(CNSSubtype$MedialMotorCortex[i] == "OX" && CNSSubtype$ThoracicSubtype[i] == "TD"){
      OtoT_tho = OtoT_tho+1
    }else if(CNSSubtype$MedialMotorCortex[i] == "OX" && CNSSubtype$ThoracicSubtype[i] == "GLIA"){
      OtoG_tho = OtoG_tho+1
    }
    
  }
  
  if(CNSSubtype$LumbarSubtype[i] == "" || CNSSubtype$MedialMotorCortex[i] == ""){
    miscountlumb = miscountlumb+1
  }else if(CNSSubtype$MedialMotorCortex[i] == CNSSubtype$LumbarSubtype[i]){
    concountlumb = concountlumb+1
    if(CNSSubtype$LumbarSubtype[i] =="GLIA"){
      GoodG_lum = GoodG_lum+1
    }else if(CNSSubtype$LumbarSubtype[i] == "OX"){
      GoodO_lum = GoodO_lum+1
    }else if(CNSSubtype$LumbarSubtype[i] == "TD"){
      GoodT_lum = GoodT_lum+1
    }
  }else{
    discountlumb = discountlumb+1
    
    if(CNSSubtype$MedialMotorCortex[i] == "GLIA" && CNSSubtype$LumbarSubtype[i] == "TD"){
      GtoT_lum = GtoT_lum+1
    }else if(CNSSubtype$MedialMotorCortex[i] == "GLIA" && CNSSubtype$LumbarSubtype[i] == "OX"){
      GtoO_lum = GtoO_lum+1
    }
    
    if(CNSSubtype$MedialMotorCortex[i] == "TD" && CNSSubtype$LumbarSubtype[i] == "GLIA"){
      TtoG_lum = TtoG_lum+1
    }else if(CNSSubtype$MedialMotorCortex[i] == "TD" && CNSSubtype$LumbarSubtype[i] == "OX"){
      TtoO_lum = TtoO_lum+1
    }
    
    if(CNSSubtype$MedialMotorCortex[i] == "OX" && CNSSubtype$LumbarSubtype[i] == "TD"){
      OtoT_lum = OtoT_lum+1
    }else if(CNSSubtype$MedialMotorCortex[i] == "OX" && CNSSubtype$LumbarSubtype[i] == "GLIA"){
      OtoG_lum = OtoG_lum+1
    }
  }
  
  
}

#NoGlia Angle: 140
#NoGliaNoRin Angle: 110
ConcordantPercent_cerv = concountcerv/sum(concountcerv,discountcerv)
DiscordantPercent_cerv = 1-ConcordantPercent_cerv
pie(c(ConcordantPercent_cerv*100,DiscordantPercent_cerv*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 114,main="Subtype Concordance - Cervical Spinal Cord - MMC")


#NoGlia Angle: 129
#NoGliaNoRin Angle: 115
ConcordantPercent_thor = concountthor/sum(concountthor,discountthor)
DiscordantPercent_thor = 1-ConcordantPercent_thor
pie(c(ConcordantPercent_thor*100,DiscordantPercent_thor*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 120,main="Subtype Concordance - Thoracic Spinal Cord - MMC")


#NoGlia Angle: 127
#NoGliaNoRin Angle: 98
ConcordantPercent_lumb = concountlumb/sum(concountlumb,discountlumb)
DiscordantPercent_lumb = 1-ConcordantPercent_lumb
pie(c(ConcordantPercent_lumb*100,DiscordantPercent_lumb*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 90,main="Subtype Concordance - Lumbar Spinal Cord - MMC")



################################# Cortex vs Cervical Spinal Cord

#ALS-Glia
a = GoodG_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
b = GtoO_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
c = GtoT_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
pie(c(a,b,c),main="ALS-Glia Clustering: Medial Motor Cortex vs Cervical Spinal Cord",col = c("goldenrod1","chartreuse3","darkorange1"),labels=c("Cortex Concordant","Cortex Discordant: Ox","Cortex Discordant: TD"),init.angle = 100)

#ALS-Ox
a = GoodO_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
b = OtoG_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
c = OtoT_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
pie(c(a,b,c),main="ALS-Ox Clustering: Medial Motor Cortex vs Cervical Spinal Cord",col = c("navy","chartreuse3","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: TD"),init.angle = 123)

#ALS-TD
a = GoodT_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
b = TtoG_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
c = TtoO_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
pie(c(a,b,c),main="ALS-TD Clustering: Medial Motor Cortex vs Cervical Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: Ox"),init.angle = 108)

################################# Cortex vs Thoracic Spinal Cord

#ALS-Glia: NA
a = GoodG_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
b = GtoO_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
c = GtoT_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
pie(c(a,b,c),main="ALS-Glia Clustering: Medial Motor Cortex vs Thoracic Spinal Cord",col = c("goldenrod1","chartreuse3","darkorange1"),labels=c("Cortex Concordant","Cortex Discordant: Ox","Cortex Discordant: TD"),init.angle = 70)

#ALS-Ox
a = GoodO_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
b = OtoG_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
c = OtoT_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
pie(c(a,b,c),main="ALS-Ox Clustering: Medial Motor Cortex vs Thoracic Spinal Cord",col = c("navy","chartreuse3","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: TD"),init.angle = 130)

#ALS-TD
a = GoodT_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
b = TtoG_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
c = TtoO_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
pie(c(a,b,c),main="ALS-TD Clustering: Medial Motor Cortex vs Thoracic Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: Ox"),init.angle = 160)


################################# Cortex vs Lumbar Spinal Cord

#ALS-Glia
a = GoodG_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
b = GtoO_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
c = GtoT_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
pie(c(a,b,c),main="ALS-Glia Clustering: Medial Motor Cortex vs Lumbar Spinal Cord",col = c("goldenrod1","chartreuse3","darkorange1"),labels=c("Cortex Concordant","Cortex Discordant: Ox","Cortex Discordant: TD"),init.angle = 135)

#ALS-Ox
a = GoodO_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
b = OtoG_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
c = OtoT_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
pie(c(a,b,c),main="ALS-Ox Clustering: Medial Motor Cortex vs Lumbar Spinal Cord",col = c("navy","chartreuse3","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: TD"),init.angle = 90)

#ALS-TD
a = GoodT_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
b = TtoG_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
c = TtoO_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
pie(c(a,b,c),main="ALS-TD Clustering: Medial Motor Cortex vs Lumbar Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: Ox"),init.angle = 50)


Medial_Concord = matrix(c(concountcerv,concountthor,concountlumb,discountcerv,discountthor,discountlumb),ncol=2)
colnames(Medial_Concord) = c("Concordant","Discordant")
rownames(Medial_Concord) = c("CervicalSpinal","ThoracicSpinal","LumbarSpinal")

MM_C_Re = matrix(c(GoodG_cer,OtoG_cer,TtoG_cer,GtoO_cer,GoodO_cer,TtoO_cer,GtoT_cer,OtoT_cer,GoodT_cer),ncol=3)
colnames(MM_C_Re) = c("Glia","Ox","TD"); rownames(MM_C_Re) = c("Glia","Ox","TD")
MM_T_Re = matrix(c(GoodG_tho,OtoG_tho,TtoG_tho,GtoO_tho,GoodO_tho,TtoO_tho,GtoT_tho,OtoT_tho,GoodT_tho),ncol=3)
colnames(MM_T_Re) = c("Glia","Ox","TD"); rownames(MM_T_Re) = c("Glia","Ox","TD")
MM_L_Re = matrix(c(GoodG_lum,OtoG_lum,TtoG_lum,GtoO_lum,GoodO_lum,TtoO_lum,GtoT_lum,OtoT_lum,GoodT_lum),ncol=3)
colnames(MM_L_Re) = c("Glia","Ox","TD"); rownames(MM_L_Re) = c("Glia","Ox","TD")

################################## Unspecified MOTOR CORTEX ##################################################################################################

concountcerv = concountlumb = concountthor = discountcerv = discountlumb = discountthor = 0
miscountcerv = miscountlumb = miscountthor = 0 #Add in a counter to keep track of unavailable tissue samples
GtoT_cer = GtoO_cer = TtoO_cer = TtoG_cer = OtoT_cer = OtoG_cer = 0
GtoT_tho = GtoO_tho = TtoO_tho = TtoG_tho = OtoT_tho = OtoG_tho = 0
GtoT_lum = GtoO_lum = TtoO_lum = TtoG_lum = OtoT_lum = OtoG_lum = 0
GoodG_cer = GoodO_cer = GoodT_cer = GoodG_tho = GoodO_tho = GoodT_tho = GoodG_lum = GoodO_lum = GoodT_lum = 0
for(i in 1:nrow(CNSSubtype)){
  if(CNSSubtype$CervicalSubtype[i] == "" || CNSSubtype$UnspecMotorCortex[i] == ""){
    miscountcerv = miscountcerv+1
  }else if(CNSSubtype$UnspecMotorCortex[i] == CNSSubtype$CervicalSubtype[i]){
    concountcerv = concountcerv+1
    if(CNSSubtype$CervicalSubtype[i] =="GLIA"){
      GoodG_cer = GoodG_cer+1
    }else if(CNSSubtype$CervicalSubtype[i] == "OX"){
      GoodO_cer = GoodO_cer+1
    }else if(CNSSubtype$CervicalSubtype[i] == "TD"){
      GoodT_cer = GoodT_cer+1
    }
  }else{
    discountcerv = discountcerv+1
    
    if(CNSSubtype$UnspecMotorCortex[i] == "GLIA" && CNSSubtype$CervicalSubtype[i] == "TD"){
      GtoT_cer = GtoT_cer+1
    }else if(CNSSubtype$UnspecMotorCortex[i] == "GLIA" && CNSSubtype$CervicalSubtype[i] == "OX"){
      GtoO_cer = GtoO_cer+1
    }
    
    if(CNSSubtype$UnspecMotorCortex[i] == "TD" && CNSSubtype$CervicalSubtype[i] == "GLIA"){
      TtoG_cer = TtoG_cer+1
    }else if(CNSSubtype$UnspecMotorCortex[i] == "TD" && CNSSubtype$CervicalSubtype[i] == "OX"){
      TtoO_cer = TtoO_cer+1
    }
    
    if(CNSSubtype$UnspecMotorCortex[i] == "OX" && CNSSubtype$CervicalSubtype[i] == "TD"){
      OtoT_cer = OtoT_cer+1
    }else if(CNSSubtype$UnspecMotorCortex[i] == "OX" && CNSSubtype$CervicalSubtype[i] == "GLIA"){
      OtoG_cer = OtoG_cer+1
    }
    
  }
  
  if(CNSSubtype$ThoracicSubtype[i] == "" || CNSSubtype$UnspecMotorCortex[i] == ""){
    miscountthor = miscountthor+1
  }else if(CNSSubtype$UnspecMotorCortex[i] == CNSSubtype$ThoracicSubtype[i]){
    concountthor = concountthor+1
    if(CNSSubtype$ThoracicSubtype[i] =="GLIA"){
      GoodG_tho = GoodG_tho+1
    }else if(CNSSubtype$ThoracicSubtype[i] == "OX"){
      GoodO_tho = GoodO_tho+1
    }else if(CNSSubtype$ThoracicSubtype[i] == "TD"){
      GoodT_tho = GoodT_tho+1
    }
  }else{
    discountthor = discountthor+1
    
    if(CNSSubtype$UnspecMotorCortex[i] == "GLIA" && CNSSubtype$ThoracicSubtype[i] == "TD"){
      GtoT_tho = GtoT_tho+1
    }else if(CNSSubtype$UnspecMotorCortex[i] == "GLIA" && CNSSubtype$ThoracicSubtype[i] == "OX"){
      GtoO_tho = GtoO_tho+1
    }
    
    if(CNSSubtype$UnspecMotorCortex[i] == "TD" && CNSSubtype$ThoracicSubtype[i] == "GLIA"){
      TtoG_tho = TtoG_tho+1
    }else if(CNSSubtype$UnspecMotorCortex[i] == "TD" && CNSSubtype$ThoracicSubtype[i] == "OX"){
      TtoO_tho = TtoO_tho+1
    }
    
    if(CNSSubtype$UnspecMotorCortex[i] == "OX" && CNSSubtype$ThoracicSubtype[i] == "TD"){
      OtoT_tho = OtoT_tho+1
    }else if(CNSSubtype$UnspecMotorCortex[i] == "OX" && CNSSubtype$ThoracicSubtype[i] == "GLIA"){
      OtoG_tho = OtoG_tho+1
    }
    
  }
  
  if(CNSSubtype$LumbarSubtype[i] == "" || CNSSubtype$UnspecMotorCortex[i] == ""){
    miscountlumb = miscountlumb+1
  }else if(CNSSubtype$UnspecMotorCortex[i] == CNSSubtype$LumbarSubtype[i]){
    concountlumb = concountlumb+1
    if(CNSSubtype$LumbarSubtype[i] =="GLIA"){
      GoodG_lum = GoodG_lum+1
    }else if(CNSSubtype$LumbarSubtype[i] == "OX"){
      GoodO_lum = GoodO_lum+1
    }else if(CNSSubtype$LumbarSubtype[i] == "TD"){
      GoodT_lum = GoodT_lum+1
    }
  }else{
    discountlumb = discountlumb+1
    
    if(CNSSubtype$UnspecMotorCortex[i] == "GLIA" && CNSSubtype$LumbarSubtype[i] == "TD"){
      GtoT_lum = GtoT_lum+1
    }else if(CNSSubtype$UnspecMotorCortex[i] == "GLIA" && CNSSubtype$LumbarSubtype[i] == "OX"){
      GtoO_lum = GtoO_lum+1
    }
    
    if(CNSSubtype$UnspecMotorCortex[i] == "TD" && CNSSubtype$LumbarSubtype[i] == "GLIA"){
      TtoG_lum = TtoG_lum+1
    }else if(CNSSubtype$UnspecMotorCortex[i] == "TD" && CNSSubtype$LumbarSubtype[i] == "OX"){
      TtoO_lum = TtoO_lum+1
    }
    
    if(CNSSubtype$UnspecMotorCortex[i] == "OX" && CNSSubtype$LumbarSubtype[i] == "TD"){
      OtoT_lum = OtoT_lum+1
    }else if(CNSSubtype$UnspecMotorCortex[i] == "OX" && CNSSubtype$LumbarSubtype[i] == "GLIA"){
      OtoG_lum = OtoG_lum+1
    }
  }
  
  
}

#NoGlia Angle: 140
#NoGliaNoRin Angle: 110
ConcordantPercent_cerv = concountcerv/sum(concountcerv,discountcerv)
DiscordantPercent_cerv = 1-ConcordantPercent_cerv
pie(c(ConcordantPercent_cerv*100,DiscordantPercent_cerv*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 0,main="Subtype Concordance - Cervical Spinal Cord - Unspec. MC")


#No comparisons to be made in NovaSeq cohort
ConcordantPercent_thor = concountthor/sum(concountthor,discountthor)
DiscordantPercent_thor = 1-ConcordantPercent_thor
pie(c(ConcordantPercent_thor*100,DiscordantPercent_thor*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 0,main="Subtype Concordance - Thoracic Spinal Cord - Unspec. MC")


#NoGlia Angle: 127
#NoGliaNoRin Angle: 98
ConcordantPercent_lumb = concountlumb/sum(concountlumb,discountlumb)
DiscordantPercent_lumb = 1-ConcordantPercent_lumb
pie(c(ConcordantPercent_lumb*100,DiscordantPercent_lumb*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 0,main="Subtype Concordance - Lumbar Spinal Cord - Unspec. MC")


################################# Cortex vs Cervical Spinal Cord

#ALS-Glia
# a = GoodG_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
# b = GtoO_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
# c = GtoT_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
# pie(c(a,b,c),main="ALS-Glia Clustering: Unspecified Motor Cortex vs Cervical Spinal Cord",col = c("goldenrod1","chartreuse3","darkorange1"),labels=c("Cortex Concordant","Cortex Discordant: Ox","Cortex Discordant: TD"),init.angle = 78)

#ALS-Ox
a = GoodO_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
b = OtoG_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
c = OtoT_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
pie(c(a,b,c),main="ALS-Ox Clustering: Unspecified Motor Cortex vs Cervical Spinal Cord",col = c("navy","chartreuse3","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: TD"),init.angle = 90)

#ALS-TD
a = GoodT_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
b = TtoG_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
c = TtoO_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
pie(c(a,b,c),main="ALS-TD Clustering: Unspecified Motor Cortex vs Cervical Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: Ox"),init.angle = 90)

################################# Cortex vs Thoracic Spinal Cord

#ALS-Glia
# a = GoodG_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
# b = GtoO_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
# c = GtoT_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
# pie(c(a,b,c),main="ALS-Glia Clustering: Unspecified Motor Cortex vs Thoracic Spinal Cord",col = c("goldenrod1","chartreuse3","darkorange1"),labels=c("Cortex Concordant","Cortex Discordant: Ox","Cortex Discordant: TD"),init.angle = 30)

#ALS-Ox
a = GoodO_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
b = OtoG_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
c = OtoT_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
pie(c(a,b,c),main="ALS-Ox Clustering: Unspecified Motor Cortex vs Thoracic Spinal Cord",col = c("navy","chartreuse3","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: TD"),init.angle = 0)
# pie(b,main="ALS-Ox Clustering: Unspecified Motor Cortex vs Thoracic Spinal Cord",col = c("chartreuse3"),labels="Cortex Discordant: Glia",init.angle = 0)

#ALS-TD
# a = GoodT_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
# b = TtoG_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
# c = TtoO_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
# pie(c(a,b,c),main="ALS-TD Clustering: Unspecified Motor Cortex vs Thoracic Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: Ox"),init.angle = 150)


################################# Cortex vs Lumbar Spinal Cord

#ALS-Glia
# a = GoodG_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
# b = GtoO_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
# c = GtoT_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
# pie(c(a,b,c),main="ALS-Glia Clustering: Unspecified Motor Cortex vs Lumbar Spinal Cord",col = c("goldenrod1","chartreuse3","darkorange1"),labels=c("Cortex Concordant","Cortex Discordant: Ox","Cortex Discordant: TD"),init.angle = 75)

#ALS-Ox
a = GoodO_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
b = OtoG_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
c = OtoT_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
pie(c(a,b,c),main="ALS-Ox Clustering: Unspecified Motor Cortex vs Lumbar Spinal Cord",col = c("navy","chartreuse3","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: TD"),init.angle = 0)

#ALS-TD
a = GoodT_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
b = TtoG_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
c = TtoO_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
pie(c(a,b,c),main="ALS-TD Clustering: Unspecified Motor Cortex vs Lumbar Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Cortex Concordant","Cortex Discordant: Glia","Cortex Discordant: Ox"),init.angle = 0)


Unspec_Concord = matrix(c(concountcerv,concountthor,concountlumb,discountcerv,discountthor,discountlumb),ncol=2)
colnames(Medial_Concord) = c("Concordant","Discordant")
rownames(Medial_Concord) = c("CervicalSpinal","ThoracicSpinal","LumbarSpinal")

Un_C_Re = matrix(c(GoodG_cer,OtoG_cer,TtoG_cer,GtoO_cer,GoodO_cer,TtoO_cer,GtoT_cer,OtoT_cer,GoodT_cer),ncol=3)
colnames(Un_C_Re) = c("Glia","Ox","TD"); rownames(Un_C_Re) = c("Glia","Ox","TD")
Un_T_Re = matrix(c(GoodG_tho,OtoG_tho,TtoG_tho,GtoO_tho,GoodO_tho,TtoO_tho,GtoT_tho,OtoT_tho,GoodT_tho),ncol=3)
colnames(Un_T_Re) = c("Glia","Ox","TD"); rownames(Un_T_Re) = c("Glia","Ox","TD")
Un_L_Re = matrix(c(GoodG_lum,OtoG_lum,TtoG_lum,GtoO_lum,GoodO_lum,TtoO_lum,GtoT_lum,OtoT_lum,GoodT_lum),ncol=3)
colnames(Un_L_Re) = c("Glia","Ox","TD"); rownames(Un_L_Re) = c("Glia","Ox","TD")