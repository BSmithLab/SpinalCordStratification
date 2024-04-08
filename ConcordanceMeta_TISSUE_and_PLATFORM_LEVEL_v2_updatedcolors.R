########################################### TISSUE-SPECIFIC CONCORDANCE ####################################################################################
############################################################################################################################################################
#At the tissue level (no meta assigned subtype)

#NovaSeq: Without glial markers (Jack Humphrey paper) and RIN-dependent genes
#CNSSubtype = read.csv("G:/SpinalCord/Publication/SystemicAnalysis/NoGliaNoRIN/NovaSeq_ConcordancePatientPheno_NoGliaNoRIN_SpinalCord_ALLTISSUE_9-5-23.csv")

#Sex, Site, Tissue, and RIN dependent genes expression removal
CNSSubtype = read.csv("G:/SpinalCord/Publication/Manuscript/Tables/ALS_Cortex_SpinalCord_4Covar_SubtypeConcordance_AssignedSubtype_9-20-23.csv") #Files found at: https://figshare.com/authors/Jarrett_Eshima/13813720

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


ConcordantPercent_cerv = concountcerv/sum(concountcerv,discountcerv)
DiscordantPercent_cerv = 1-ConcordantPercent_cerv
pie(c(ConcordantPercent_cerv*100,DiscordantPercent_cerv*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 106,main="Subtype Concordance - Cervical Spinal Cord - FC")
text(-0.5,0,"70",cex=4)
text(0.5,0,"100",cex=4)

ConcordantPercent_thor = concountthor/sum(concountthor,discountthor)
DiscordantPercent_thor = 1-ConcordantPercent_thor
pie(c(ConcordantPercent_thor*100,DiscordantPercent_thor*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 120,main="Subtype Concordance - Thoracic Spinal Cord - FC")
text(-0.5,0,"18",cex=4)
text(0.5,0,"35",cex=4)

ConcordantPercent_lumb = concountlumb/sum(concountlumb,discountlumb)
DiscordantPercent_lumb = 1-ConcordantPercent_lumb
pie(c(ConcordantPercent_lumb*100,DiscordantPercent_lumb*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 104,main="Subtype Concordance - Lumbar Spinal Cord - FC")
text(-0.5,0,"68",cex=4)
text(0.5,0,"91",cex=4)

################################# Cortex vs Cervical Spinal Cord

#ALS-Glia
a = GoodG_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
b = GtoO_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
c = GtoT_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
pie(c(a,b,c),main="ALS-Glia Clustering: Frontal Cortex vs Cervical Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 104)
text(-0.5,0,"13",cex=6,col="white")
text(0.35,0.4,"9",cex=6,col="white")
text(0.35,-0.4,"9",cex=6,col="white")

#ALS-Ox
a = GoodO_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
b = OtoG_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
c = OtoT_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
pie(c(a,b,c),main="ALS-Ox Clustering: Frontal Cortex vs Cervical Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 128)
text(-0.5,0,"26",cex=6,col="white")
text(0.35,0.25,"48",cex=6,col="white")
text(-0.075,-0.5,"14",cex=6,col="white")

#ALS-TD
a = GoodT_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
b = TtoG_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
c = TtoO_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
pie(c(a,b,c),main="ALS-TD Clustering: Frontal Cortex vs Cervical Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 72)
text(-0.45,0,"31",cex=6,col="white")
text(0.4,0.2,"13",cex=6,col="white")
text(.4,-0.4,"7",cex=6,col="white")

################################# Cortex vs Thoracic Spinal Cord

#ALS-Glia - no ALS-Glia subtypes in the NovaSeq thoracic spinal cord 
a = GoodG_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
b = GtoO_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
c = GtoT_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
pie(c(a,b,c),main="ALS-Glia Clustering: Frontal Cortex vs Thoracic Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 110)
text(-0.5,0,"2",cex=6,col="white")
text(0.5,0,"3",cex=6,col="white")

#ALS-Ox
a = GoodO_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
b = OtoG_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
c = OtoT_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
pie(c(a,b,c),main="ALS-Ox Clustering: Frontal Cortex vs Thoracic Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 130)
text(-0.5,0,"9",cex=6,col="white")
text(0.4,0.25,"18",cex=6,col="white")
text(-0.075,-0.5,"6",cex=6,col="white")

#ALS-TD
a = GoodT_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
b = TtoG_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
c = TtoO_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
pie(c(a,b,c),main="ALS-TD Clustering: Frontal Cortex vs Thoracic Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 94)
text(-0.5,0,"7",cex=6,col="white")
text(0.275,0.425,"3",cex=6,col="white")
text(.4,-0.3,"5",cex=6,col="white")

################################# Cortex vs Lumbar Spinal Cord

#ALS-Glia
a = GoodG_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
b = GtoO_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
c = GtoT_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
pie(c(a,b,c),main="ALS-Glia Clustering: Frontal Cortex vs Lumbar Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 96)
text(-0.5,0,"12",cex=6,col="white")
text(0.075,0.6,"2",cex=6,col="white")
text(0.45,-0.1,"12",cex=6,col="white")

#ALS-Ox
a = GoodO_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
b = OtoG_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
c = OtoT_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
pie(c(a,b,c),main="ALS-Ox Clustering: Frontal Cortex vs Lumbar Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 103)
text(-0.45,0,"38",cex=6,col="white")
text(0.35,0.4,"25",cex=6,col="white")
text(.35,-0.35,"22",cex=6,col="white")

#ALS-TD
a = GoodT_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
b = TtoG_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
c = TtoO_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
pie(c(a,b,c),main="ALS-TD Clustering: Frontal Cortex vs Lumbar Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 110)
text(-0.5,0,"18",cex=6,col="white")
text(0.35,0.30,"18",cex=6,col="white")
text(.15,-0.45,"12",cex=6,col="white")

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


ConcordantPercent_cerv = concountcerv/sum(concountcerv,discountcerv)
DiscordantPercent_cerv = 1-ConcordantPercent_cerv
pie(c(ConcordantPercent_cerv*100,DiscordantPercent_cerv*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 115,main="Subtype Concordance - Cervical Spinal Cord - LMC")
text(-0.5,0,"33",cex=4)
text(0.5,0,"60",cex=4)

ConcordantPercent_thor = concountthor/sum(concountthor,discountthor)
DiscordantPercent_thor = 1-ConcordantPercent_thor
pie(c(ConcordantPercent_thor*100,DiscordantPercent_thor*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 108,main="Subtype Concordance - Thoracic Spinal Cord - LMC")
text(-0.5,0,"19",cex=4)
text(0.5,0,"29",cex=4)

ConcordantPercent_lumb = concountlumb/sum(concountlumb,discountlumb)
DiscordantPercent_lumb = 1-ConcordantPercent_lumb
pie(c(ConcordantPercent_lumb*100,DiscordantPercent_lumb*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 108,main="Subtype Concordance - Lumbar Spinal Cord - LMC")
text(-0.5,0,"33",cex=4)
text(0.5,0,"51",cex=4)

################################# Cortex vs Cervical Spinal Cord

#ALS-Glia
a = GoodG_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
b = GtoO_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
c = GtoT_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
pie(c(a,b,c),main="ALS-Glia Clustering: Lateral Motor Cortex vs Cervical Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 124)
text(-0.5,0,"4",cex=6,col="white")
text(0.3,0.35,"5",cex=6,col="white")
text(0.175,-0.45,"4",cex=6,col="white")

#ALS-Ox
a = GoodO_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
b = OtoG_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
c = OtoT_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
pie(c(a,b,c),main="ALS-Ox Clustering: Lateral Motor Cortex vs Cervical Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 120)
text(-0.5,0,"19",cex=6,col="white")
text(0.45,0.1,"34",cex=6,col="white")
text(-0.15,-0.55,"5",cex=6,col="white")

#ALS-TD
a = GoodT_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
b = TtoG_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
c = TtoO_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
pie(c(a,b,c),main="ALS-TD Clustering: Lateral Motor Cortex vs Cervical Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 98)
text(-0.5,0,"10",cex=6,col="white")
text(0.225,0.45,"4",cex=6,col="white")
text(.4,-0.25,"8",cex=6,col="white")

################################# Cortex vs Thoracic Spinal Cord

#ALS-Glia
a = GoodG_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
b = GtoO_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
c = GtoT_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
pie(c(a,b,c),main="ALS-Glia Clustering: Lateral Motor Cortex vs Thoracic Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 98)
text(-0.5,0,"5",cex=6,col="white")
text(0.4,0.25,"4",cex=6,col="white")
text(0.225,-0.45,"2",cex=6,col="white")

#ALS-Ox
a = GoodO_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
b = OtoG_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
c = OtoT_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
pie(c(a,b,c),main="ALS-Ox Clustering: Lateral Motor Cortex vs Thoracic Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 125)
text(-0.5,0,"8",cex=6,col="white")
text(0.4,0.1,"15",cex=6,col="white")
text(-0.15,-0.55,"2",cex=6,col="white")


a = GoodT_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
b = TtoG_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
c = TtoO_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
pie(c(a,b,c),main="ALS-TD Clustering: Lateral Motor Cortex vs Thoracic Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 90)
text(-0.50,0,"6",cex=6,col="white")
text(0.275,0.45,"2",cex=6,col="white")
text(.45,-0.25,"4",cex=6,col="white")

################################# Cortex vs Lumbar Spinal Cord

#ALS-Glia
a = GoodG_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
b = GtoO_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
c = GtoT_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
pie(c(a,b,c),main="ALS-Glia Clustering: Lateral Motor Cortex vs Lumbar Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 108)
text(-0.5,0,"4",cex=6,col="white")
text(0.175,0.5,"2",cex=6,col="white")
text(0.4,-0.25,"4",cex=6,col="white")

#ALS-Ox
a = GoodO_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
b = OtoG_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
c = OtoT_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
pie(c(a,b,c),main="ALS-Ox Clustering: Lateral Motor Cortex vs Lumbar Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 100)
text(-0.45,0,"24",cex=6,col="white")
text(0.3,0.35,"16",cex=6,col="white")
text(0.3,-0.4,"14",cex=6,col="white")

#ALS-TD
a = GoodT_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
b = TtoG_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
c = TtoO_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
pie(c(a,b,c),main="ALS-TD Clustering: Lateral Motor Cortex vs Lumbar Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 135)
text(-0.55,0,"5",cex=6,col="white")
text(0.1,0.5,"6",cex=6,col="white")
text(.3,-0.35,"9",cex=6,col="white")

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


ConcordantPercent_cerv = concountcerv/sum(concountcerv,discountcerv)
DiscordantPercent_cerv = 1-ConcordantPercent_cerv
pie(c(ConcordantPercent_cerv*100,DiscordantPercent_cerv*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 108,main="Subtype Concordance - Cervical Spinal Cord - MMC")
text(-0.5,0,"38",cex=4)
text(0.5,0,"58",cex=4)

ConcordantPercent_thor = concountthor/sum(concountthor,discountthor)
DiscordantPercent_thor = 1-ConcordantPercent_thor
pie(c(ConcordantPercent_thor*100,DiscordantPercent_thor*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 116,main="Subtype Concordance - Thoracic Spinal Cord - MMC")
text(-0.5,0,"18",cex=4)
text(0.5,0,"32",cex=4)

ConcordantPercent_lumb = concountlumb/sum(concountlumb,discountlumb)
DiscordantPercent_lumb = 1-ConcordantPercent_lumb
pie(c(ConcordantPercent_lumb*100,DiscordantPercent_lumb*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 100,main="Subtype Concordance - Lumbar Spinal Cord - MMC")
text(-0.5,0,"38",cex=4)
text(0.5,0,"47",cex=4)


################################# Cortex vs Cervical Spinal Cord

#ALS-Glia
a = GoodG_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
b = GtoO_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
c = GtoT_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
pie(c(a,b,c),main="ALS-Glia Clustering: Medial Motor Cortex vs Cervical Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 148)
text(-0.5,0,"2",cex=6,col="white")
text(0.4,0.2,"7",cex=6,col="white")
text(-0.2,-0.45,"2",cex=6,col="white")

#ALS-Ox
a = GoodO_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
b = OtoG_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
c = OtoT_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
pie(c(a,b,c),main="ALS-Ox Clustering: Medial Motor Cortex vs Cervical Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 120)
text(-0.5,0,"20",cex=6,col="white")
text(0.4,0.2,"32",cex=6,col="white")
text(-0.05,-0.525,"10",cex=6,col="white")

#ALS-TD
a = GoodT_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
b = TtoG_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
c = TtoO_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
pie(c(a,b,c),main="ALS-TD Clustering: Medial Motor Cortex vs Cervical Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 55)
text(-0.45,0,"16",cex=6,col="white")
text(0.5,0.3,"3",cex=6,col="white")
text(.525,-0.2,"4",cex=6,col="white")

################################# Cortex vs Thoracic Spinal Cord

#ALS-Glia
a = GoodG_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
b = GtoO_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
c = GtoT_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
pie(c(a,b,c),main="ALS-Glia Clustering: Medial Motor Cortex vs Thoracic Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 160)
text(-0.6,0,"1",cex=6,col="white")
text(0.4,0.25,"6",cex=6,col="white")
text(-0.25,-0.45,"2",cex=6,col="white")

#ALS-Ox
a = GoodO_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
b = OtoG_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
c = OtoT_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
pie(c(a,b,c),main="ALS-Ox Clustering: Medial Motor Cortex vs Thoracic Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 116)
text(-0.5,0,"9",cex=6,col="white")
text(0.4,0.25,"12",cex=6,col="white")
text(0.05,-0.5,"4",cex=6,col="white")

#ALS-TD
a = GoodT_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
b = TtoG_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
c = TtoO_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
pie(c(a,b,c),main="ALS-TD Clustering: Medial Motor Cortex vs Thoracic Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 90)
text(-0.5,0,"8",cex=6,col="white")
text(0.115,0.6,"1",cex=6,col="white")
text(.5,-0.1,"7",cex=6,col="white")

################################# Cortex vs Lumbar Spinal Cord

#ALS-Glia
a = GoodG_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
b = GtoO_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
c = GtoT_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
pie(c(a,b,c),main="ALS-Glia Clustering: Medial Motor Cortex vs Lumbar Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 115)
text(-0.5,0,"4",cex=6,col="white")
text(0.35,0.35,"4",cex=6,col="white")
text(0.225,-0.4,"3",cex=6,col="white")

#ALS-Ox
a = GoodO_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
b = OtoG_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
c = OtoT_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
pie(c(a,b,c),main="ALS-Ox Clustering: Medial Motor Cortex vs Lumbar Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 95)
text(-0.45,0,"26",cex=6,col="white")
text(0.325,0.35,"15",cex=6,col="white")
text(0.325,-0.35,"15",cex=6,col="white")

#ALS-TD
a = GoodT_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
b = TtoG_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
c = TtoO_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
pie(c(a,b,c),main="ALS-TD Clustering: Medial Motor Cortex vs Lumbar Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 100)
text(-0.5,0,"8",cex=6,col="white")
text(0.1,0.6,"2",cex=6,col="white")
text(.45,-0.15,"8",cex=6,col="white")

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


ConcordantPercent_cerv = concountcerv/sum(concountcerv,discountcerv)
DiscordantPercent_cerv = 1-ConcordantPercent_cerv
pie(c(ConcordantPercent_cerv*100,DiscordantPercent_cerv*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 100,main="Subtype Concordance - Cervical Spinal Cord - Unspec. MC")
text(-0.5,0,"21",cex=4)
text(0.5,0,"27",cex=4)

ConcordantPercent_thor = concountthor/sum(concountthor,discountthor)
DiscordantPercent_thor = 1-ConcordantPercent_thor
pie(c(ConcordantPercent_thor*100,DiscordantPercent_thor*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 0,main="Subtype Concordance - Thoracic Spinal Cord - Unspec. MC")
text(-0.5,0,"2",cex=4)


ConcordantPercent_lumb = concountlumb/sum(concountlumb,discountlumb)
DiscordantPercent_lumb = 1-ConcordantPercent_lumb
pie(c(ConcordantPercent_lumb*100,DiscordantPercent_lumb*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 87,main="Subtype Concordance - Lumbar Spinal Cord - Unspec. MC")
text(-0.5,0,"22",cex=4)
text(0.5,0,"21",cex=4)

################################# Cortex vs Cervical Spinal Cord

#ALS-Glia
a = GoodG_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
b = GtoO_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
c = GtoT_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
pie(c(a,b,c),main="ALS-Glia Clustering: Unspecified Motor Cortex vs Cervical Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 100)
text(-0.5,0,"8",cex=6,col="white")
text(0.275,0.45,"4",cex=6,col="white")
text(0.4,-0.3,"6",cex=6,col="white")

#ALS-Ox
a = GoodO_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
b = OtoG_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
c = OtoT_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
pie(c(a,b,c),main="ALS-Ox Clustering: Unspecified Motor Cortex vs Cervical Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 112)
text(-0.5,0,"7",cex=6,col="white")
text(0.35,0.325,"7",cex=6,col="white")
text(0.2,-0.45,"5",cex=6,col="white")

#ALS-TD
a = GoodT_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
b = TtoG_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
c = TtoO_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
pie(c(a,b,c),main="ALS-TD Clustering: Unspecified Motor Cortex vs Cervical Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 82)
text(-0.5,0,"6",cex=6,col="white")
text(0.35,0.4,"2",cex=6,col="white")
text(0.4,-0.25,"3",cex=6,col="white")

################################# Cortex vs Thoracic Spinal Cord

#ALS-Glia
a = GoodG_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
b = GtoO_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
c = GtoT_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
pie(c(a,b,c),main="ALS-Glia Clustering: Unspecified Motor Cortex vs Thoracic Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 0)
text(-0.5,0,"1",cex=6,col="white")

#ALS-Ox
a = GoodO_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
b = OtoG_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
c = OtoT_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
pie(c(a,b,c),main="ALS-Ox Clustering: Unspecified Motor Cortex vs Thoracic Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 0)
text(-0.5,0,"1",cex=6,col="white")

#ALS-TD - NA
# a = GoodT_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
# b = TtoG_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
# c = TtoO_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
# pie(c(a,b,c),main="ALS-TD Clustering: Unspecified Motor Cortex vs Thoracic Spinal Cord",col = c("firebrick","darkorange1","#8d2ca3"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 150)


################################# Cortex vs Lumbar Spinal Cord

#ALS-Glia
a = GoodG_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
b = GtoO_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
c = GtoT_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
pie(c(a,b,c),main="ALS-Glia Clustering: Unspecified Motor Cortex vs Lumbar Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 85)
text(-0.5,0,"9",cex=6,col="white")
text(0.165,0.6,"1",cex=6,col="white")
text(0.5,-0.1,"7",cex=6,col="white")

#ALS-Ox
a = GoodO_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
b = OtoG_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
c = OtoT_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
pie(c(a,b,c),main="ALS-Ox Clustering: Unspecified Motor Cortex vs Lumbar Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 65)
text(-0.45,0,"10",cex=6,col="white")
text(0.42,0.4,"2",cex=6,col="white")
text(0.5,-0.2,"4",cex=6,col="white")

#ALS-TD
a = GoodT_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
b = TtoG_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
c = TtoO_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
pie(c(a,b,c),main="ALS-TD Clustering: Unspecified Motor Cortex vs Lumbar Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 125)
text(-0.5,0,"3",cex=6,col="white")
text(0.2,0.5,"3",cex=6,col="white")
text(0.3,-0.4,"4",cex=6,col="white")

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
####################################################  PLATFORM-LEVEL #######################################################################################################################################################
###########################################################################################################################################################################################################


###### CHECK DEPENDENCY ON SEQUENCING PLATFORM

#NovaSeq Pheno
NovaSeqPheno = read.csv("G:/SpinalCord/Publication/SystemicAnalysis/4Covar/SpinalCord_Concordance_NovaSeq_9-18-23.csv") #File generated in previous script, refer to Github readme


#HiSeq Pheno
HiSeqPheno = read.csv("G:/SpinalCord/Publication/SystemicAnalysis/4Covar/SpinalCord_Concordance_HiSeq_9-18-23.csv") #File generated in previous script, refer to Github readme


############################################## NOVASEQ ########################################################################################
#BREAK

CNSSubtype = NovaSeqPheno
rownames(CNSSubtype) = CNSSubtype$Patient
CNSSubtype = CNSSubtype[,3:9]
#colnames(CNSSubtype)[c(30,31,32,33)] = c("FrontalCortex","MedialMotorCortex","LateralMotorCortex","UnspecMotorCortex")

for(i in 1:nrow(CNSSubtype)){
  for(j in 1:ncol(CNSSubtype)){
    if(is.na(CNSSubtype[i,j])){
      CNSSubtype[i,j] = ""
    }
  }
}

for(i in 1:nrow(CNSSubtype)){
  for(j in 1:ncol(CNSSubtype)){
    if(! is.na(CNSSubtype[i,j]) && CNSSubtype[i,j] == "TE"){
      CNSSubtype[i,j] = "TD"
    }
  }
}


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


ConcordantPercent_cerv = concountcerv/sum(concountcerv,discountcerv)
DiscordantPercent_cerv = 1-ConcordantPercent_cerv
pie(c(ConcordantPercent_cerv*100,DiscordantPercent_cerv*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 102,main="NovaSeq Subtype Concordance - Cervical Spinal Cord - FC")
text(-0.5,0,"50",cex=4)
text(0.5,0,"66",cex=4)

ConcordantPercent_thor = concountthor/sum(concountthor,discountthor)
DiscordantPercent_thor = 1-ConcordantPercent_thor
pie(c(ConcordantPercent_thor*100,DiscordantPercent_thor*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 130,main="NovaSeq Subtype Concordance - Thoracic Spinal Cord - FC")
text(-0.55,0,"5",cex=4)
text(0.5,0,"13",cex=4)

ConcordantPercent_lumb = concountlumb/sum(concountlumb,discountlumb)
DiscordantPercent_lumb = 1-ConcordantPercent_lumb
pie(c(ConcordantPercent_lumb*100,DiscordantPercent_lumb*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 100,main="NovaSeq Subtype Concordance - Lumbar Spinal Cord - FC")
text(-0.5,0,"50",cex=4)
text(0.5,0,"62",cex=4)

################################# Cortex vs Cervical Spinal Cord

#ALS-Glia
a = GoodG_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
b = GtoO_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
c = GtoT_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
pie(c(a,b,c),main="NovaSeq ALS-Glia Clustering: Frontal Cortex vs Cervical Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 105)
text(-0.5,0,"11",cex=6,col="white")
text(0.3,0.45,"7",cex=6,col="white")
text(0.3,-0.3,"9",cex=6,col="white")

#ALS-Ox
a = GoodO_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
b = OtoG_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
c = OtoT_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
pie(c(a,b,c),main="NovaSeq ALS-Ox Clustering: Frontal Cortex vs Cervical Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 138)
text(-0.5,0,"12",cex=6,col="white")
text(0.4,0.15,"32",cex=6,col="white")
text(-0.225,-0.5,"6",cex=6,col="white")

#ALS-TD
a = GoodT_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
b = TtoG_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
c = TtoO_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
pie(c(a,b,c),main="NovaSeq ALS-TD Clustering: Frontal Cortex vs Cervical Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 55)
text(-0.45,0,"27",cex=6,col="white")
text(0.5,0.1,"10",cex=6,col="white")
text(.42,-0.43,"2",cex=6,col="white")

################################# Cortex vs Thoracic Spinal Cord

#ALS-Glia
a = GoodG_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
b = GtoO_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
c = GtoT_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
pie(c(a,b,c),main="NovaSeq ALS-Glia Clustering: Frontal Cortex vs Thoracic Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 0)
text(-0.5,0,"1",cex=6,col="white")

#ALS-Ox
a = GoodO_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
b = OtoG_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
c = OtoT_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
pie(c(a,b,c),main="NovaSeq ALS-Ox Clustering: Frontal Cortex vs Thoracic Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 160)
text(-0.6,0,"1",cex=6,col="white")
text(0.425,0.15,"7",cex=6,col="white")
text(-0.45,-0.375,"1",cex=6,col="white")

#ALS-TD
a = GoodT_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
b = TtoG_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
c = TtoO_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
pie(c(a,b,c),main="NovaSeq ALS-TD Clustering: Frontal Cortex vs Thoracic Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 90)
text(-0.5,0,"4",cex=6,col="white")
text(0.35,.35,"2",cex=6,col="white")
text(0.35,-.35,"2",cex=6,col="white")

################################# Cortex vs Lumbar Spinal Cord

#ALS-Glia
a = GoodG_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
b = GtoO_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
c = GtoT_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
pie(c(a,b,c),main="NovaSeq ALS-Glia Clustering: Frontal Cortex vs Lumbar Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 95)
text(-0.45,0,"12",cex=6,col="white")
text(0.45,-0.05,"12",cex=6,col="white")
text(0.025,0.6,"1",cex=6,col="white")

#ALS-Ox
a = GoodO_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
b = OtoG_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
c = OtoT_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
pie(c(a,b,c),main="NovaSeq ALS-Ox Clustering: Frontal Cortex vs Lumbar Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 100)
text(-0.5,0,"22",cex=6,col="white")
text(0.3,0.4,"14",cex=6,col="white")
text(.3,-0.4,"14",cex=6,col="white")

#ALS-TD
a = GoodT_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
b = TtoG_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
c = TtoO_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
pie(c(a,b,c),main="NovaSeq ALS-TD Clustering: Frontal Cortex vs Lumbar Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 102)
text(-0.5,0,"16",cex=6,col="white")
text(0.4,0.2,"15",cex=6,col="white")
text(.175,-0.5,"6",cex=6,col="white")

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


ConcordantPercent_cerv = concountcerv/sum(concountcerv,discountcerv)
DiscordantPercent_cerv = 1-ConcordantPercent_cerv
pie(c(ConcordantPercent_cerv*100,DiscordantPercent_cerv*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 113,main="NovaSeq Subtype Concordance - Cervical Spinal Cord - LMC")
text(-0.5,0,"18",cex=6)
text(0.5,0,"30",cex=6)

ConcordantPercent_thor = concountthor/sum(concountthor,discountthor)
DiscordantPercent_thor = 1-ConcordantPercent_thor
pie(c(ConcordantPercent_thor*100,DiscordantPercent_thor*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 120,main="NovaSeq Subtype Concordance - Thoracic Spinal Cord - LMC")
text(-0.55,0,"5",cex=6)
text(0.5,0,"10",cex=6)

ConcordantPercent_lumb = concountlumb/sum(concountlumb,discountlumb)
DiscordantPercent_lumb = 1-ConcordantPercent_lumb
pie(c(ConcordantPercent_lumb*100,DiscordantPercent_lumb*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 92,main="NovaSeq Subtype Concordance - Lumbar Spinal Cord - LMC")
text(-0.5,0,"21",cex=6)
text(0.5,0,"22",cex=6)

################################# Cortex vs Cervical Spinal Cord

#ALS-Glia
a = GoodG_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
b = GtoO_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
c = GtoT_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
pie(c(a,b,c),main="ALS-Glia Clustering: Lateral Motor Cortex vs Cervical Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 0)
text(-0.5,0,"2",cex=6,col="white")


#ALS-Ox
a = GoodO_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
b = OtoG_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
c = OtoT_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
pie(c(a,b,c),main="ALS-Ox Clustering: Lateral Motor Cortex vs Cervical Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 130)
text(-0.5,0,"10",cex=6,col="white")
text(0.45,0.1,"23",cex=6,col="white")
text(-0.3,-0.55,"2",cex=6,col="white")

#ALS-TD
a = GoodT_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
b = TtoG_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
c = TtoO_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
pie(c(a,b,c),main="ALS-TD Clustering: Lateral Motor Cortex vs Cervical Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 48)
text(-0.5,0,"8",cex=6,col="white")
text(0.5,0,"3",cex=6,col="white")

################################# Cortex vs Thoracic Spinal Cord

#ALS-Glia
a = GoodG_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
b = GtoO_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
c = GtoT_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
pie(c(a,b,c),main="ALS-Glia Clustering: Lateral Motor Cortex vs Thoracic Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 0)
text(-0.5,0,"1",cex=6,col="white")

#ALS-Ox
a = GoodO_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
b = OtoG_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
c = OtoT_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
pie(c(a,b,c),main="ALS-Ox Clustering: Lateral Motor Cortex vs Thoracic Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 145)
text(-0.5,0,"2",cex=6,col="white")
text(0.45,0.15,"7",cex=6,col="white")
text(-0.325,-0.45,"1",cex=6,col="white")

#ALS-TD 
a = GoodT_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
b = TtoG_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
c = TtoO_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
pie(c(a,b,c),main="ALS-TD Clustering: Lateral Motor Cortex vs Thoracic Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 45)
text(-0.5,0,"3",cex=6,col="white")
text(0.5,0,"1",cex=6,col="white")

#View(cbind(CNSSubtype$LateralMotorCortex,CNSSubtype$ThoracicSubtype))

################################# Cortex vs Lumbar Spinal Cord

#ALS-Glia
a = GoodG_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
b = GtoO_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
c = GtoT_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
pie(c(a,b,c),main="ALS-Glia Clustering: Lateral Motor Cortex vs Lumbar Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 90)
text(-0.5,0,"1",cex=6,col="white")
text(0.5,0,"1",cex=6,col="white")

#ALS-Ox
a = GoodO_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
b = OtoG_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
c = OtoT_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
pie(c(a,b,c),main="ALS-Ox Clustering: Lateral Motor Cortex vs Lumbar Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 90)
text(-0.5,0,"15",cex=6,col="white")
text(0.35,0.35,"7",cex=6,col="white")
text(0.35,-0.35,"9",cex=6,col="white")

#ALS-TD
a = GoodT_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
b = TtoG_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
c = TtoO_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
pie(c(a,b,c),main="ALS-TD Clustering: Lateral Motor Cortex vs Lumbar Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 90)
text(-0.5,0,"5",cex=6,col="white")
text(0.3,0.45,"2",cex=6,col="white")
text(.4,-0.25,"3",cex=6,col="white")

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


ConcordantPercent_cerv = concountcerv/sum(concountcerv,discountcerv)
DiscordantPercent_cerv = 1-ConcordantPercent_cerv
pie(c(ConcordantPercent_cerv*100,DiscordantPercent_cerv*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 105,main="NovaSeq Subtype Concordance - Cervical Spinal Cord - MMC")
text(-0.5,0,"22",cex=4)
text(0.5,0,"32",cex=4)

ConcordantPercent_thor = concountthor/sum(concountthor,discountthor)
DiscordantPercent_thor = 1-ConcordantPercent_thor
pie(c(ConcordantPercent_thor*100,DiscordantPercent_thor*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 116,main="NovaSeq Subtype Concordance - Thoracic Spinal Cord - MMC")
text(-0.5,0,"6",cex=4)
text(0.5,0,"11",cex=4)

ConcordantPercent_lumb = concountlumb/sum(concountlumb,discountlumb)
DiscordantPercent_lumb = 1-ConcordantPercent_lumb
pie(c(ConcordantPercent_lumb*100,DiscordantPercent_lumb*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 94,main="NovaSeq Subtype Concordance - Lumbar Spinal Cord - MMC")
text(-0.5,0,"22",cex=4)
text(0.5,0,"24",cex=4)


################################# Cortex vs Cervical Spinal Cord

#ALS-Glia
a = GoodG_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
b = GtoO_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
c = GtoT_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
pie(c(a,b,c),main="ALS-Glia Clustering: Medial Motor Cortex vs Cervical Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 0)
text(-0.5,0,"1",cex=6,col="white")

#ALS-Ox
a = GoodO_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
b = OtoG_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
c = OtoT_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
pie(c(a,b,c),main="ALS-Ox Clustering: Medial Motor Cortex vs Cervical Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 136)
text(-0.5,0,"10",cex=6,col="white")
text(0.45,0.15,"26",cex=6,col="white")
text(-0.25,-0.5,"4",cex=6,col="white")

#ALS-TD
a = GoodT_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
b = TtoG_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
c = TtoO_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
pie(c(a,b,c),main="ALS-TD Clustering: Medial Motor Cortex vs Cervical Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 15)
text(-0.45,0,"12",cex=6,col="white")
text(0.65,0,"1",cex=6,col="white")

################################# Cortex vs Thoracic Spinal Cord

#ALS-Glia - No ALS-Glia samples in the Medial Motor Cortex Novaseq Cohort
# a = GoodG_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
# b = GtoO_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
# c = GtoT_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
# pie(c(a,b,c),main="ALS-Glia Clustering: Medial Motor Cortex vs Thoracic Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 160)
# text(-0.55,0,"1",cex=6)
# text(0.45,0.35,"6",cex=6)
# text(-0.25,-0.5,"2",cex=6)


#ALS-Ox
a = GoodO_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
b = OtoG_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
c = OtoT_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
pie(c(a,b,c),main="ALS-Ox Clustering: Medial Motor Cortex vs Thoracic Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 145)
text(-0.5,0,"2",cex=6,col="white")
text(0.45,0.15,"7",cex=6,col="white")
text(-0.35,-0.5,"1",cex=6,col="white")

#ALS-TD
a = GoodT_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
b = TtoG_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
c = TtoO_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
pie(c(a,b,c),main="ALS-TD Clustering: Medial Motor Cortex vs Thoracic Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 77)
text(-0.5,0,"4",cex=6,col="white")
text(0.45,-0.2,"2",cex=6,col="white")
text(0.35,0.425,"1",cex=6,col="white")

################################# Cortex vs Lumbar Spinal Cord

#ALS-Glia
a = GoodG_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
b = GtoO_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
c = GtoT_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
pie(c(a,b,c),main="ALS-Glia Clustering: Medial Motor Cortex vs Lumbar Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 0)
text(-0.5,0,"1",cex=6,col="white")

#ALS-Ox
a = GoodO_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
b = OtoG_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
c = OtoT_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
pie(c(a,b,c),main="ALS-Ox Clustering: Medial Motor Cortex vs Lumbar Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 101)
text(-0.5,0,"15",cex=6,col="white")
text(0.3,0.4,"9",cex=6,col="white")
text(0.3,-0.375,"10",cex=6,col="white")

#ALS-TD
a = GoodT_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
b = TtoG_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
c = TtoO_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
pie(c(a,b,c),main="ALS-TD Clustering: Medial Motor Cortex vs Lumbar Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 80)
text(-0.5,0,"6",cex=6,col="white")
text(0.35,0.4,"2",cex=6,col="white")
text(.4,-0.25,"3",cex=6,col="white")

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


ConcordantPercent_cerv = concountcerv/sum(concountcerv,discountcerv)
DiscordantPercent_cerv = 1-ConcordantPercent_cerv
pie(c(ConcordantPercent_cerv*100,DiscordantPercent_cerv*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 94,main="Subtype Concordance - Cervical Spinal Cord - Unspec. MC")
text(-0.5,0,"20",cex=6)
text(0.5,0,"22",cex=6)

#There are no patients with unspecified motor cortex and thoracic spinal samples
# ConcordantPercent_thor = concountthor/sum(concountthor,discountthor)
# DiscordantPercent_thor = 1-ConcordantPercent_thor
# pie(c(ConcordantPercent_thor*100,DiscordantPercent_thor*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 0,main="Subtype Concordance - Thoracic Spinal Cord - Unspec. MC")
# text(-0.5,0,"2",cex=6)

ConcordantPercent_lumb = concountlumb/sum(concountlumb,discountlumb)
DiscordantPercent_lumb = 1-ConcordantPercent_lumb
pie(c(ConcordantPercent_lumb*100,DiscordantPercent_lumb*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 81,main="Subtype Concordance - Lumbar Spinal Cord - Unspec. MC")
text(-0.5,0,"22",cex=6)
text(0.5,0,"18",cex=6)

################################# Cortex vs Cervical Spinal Cord

#ALS-Glia
a = GoodG_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
b = GtoO_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
c = GtoT_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
pie(c(a,b,c),main="ALS-Glia Clustering: Unspecified Motor Cortex vs Cervical Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 96)
text(-0.5,0,"8",cex=6,col="white")
text(0.25,0.475,"3",cex=6,col="white")
text(0.4,-0.25,"6",cex=6,col="white")

#ALS-Ox
a = GoodO_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
b = OtoG_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
c = OtoT_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
pie(c(a,b,c),main="ALS-Ox Clustering: Unspecified Motor Cortex vs Cervical Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 102)
text(-0.5,0,"7",cex=6,col="white")
text(0.35,0.35,"5",cex=6,col="white")
text(0.3,-0.4,"4",cex=6,col="white")

#ALS-TD
a = GoodT_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
b = TtoG_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
c = TtoO_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
pie(c(a,b,c),main="ALS-TD Clustering: Unspecified Motor Cortex vs Cervical Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 82)
text(-0.5,0,"5",cex=6,col="white")
text(0.4,0.35,"2",cex=6,col="white")
text(0.4,-0.3,"2",cex=6,col="white")

################################# Cortex vs Thoracic Spinal Cord - There are no samples pairings in the NovaSeq Cohort
#View(cbind(CNSSubtype$UnspecMotorCortex,CNSSubtype$ThoracicSubtype))
#ALS-Glia
# a = GoodG_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
# b = GtoO_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
# c = GtoT_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
# pie(c(a,b,c),main="ALS-Glia Clustering: Unspecified Motor Cortex vs Thoracic Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 0)
# text(-0.55,0,"1",cex=6)

#ALS-Ox
# a = GoodO_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
# b = OtoG_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
# c = OtoT_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
# pie(c(a,b,c),main="ALS-Ox Clustering: Unspecified Motor Cortex vs Thoracic Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 0)
# text(-0.55,0,"1",cex=6)

#ALS-TD - NA
#a = GoodT_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
#b = TtoG_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
#c = TtoO_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
# pie(c(a,b,c),main="ALS-TD Clustering: Unspecified Motor Cortex vs Thoracic Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 150)


################################# Cortex vs Lumbar Spinal Cord

#ALS-Glia
a = GoodG_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
b = GtoO_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
c = GtoT_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
pie(c(a,b,c),main="ALS-Glia Clustering: Unspecified Motor Cortex vs Lumbar Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 78)
text(-0.5,0,"9",cex=6,col="white")
text(0.5,0,"7",cex=6,col="white")

#ALS-Ox
a = GoodO_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
b = OtoG_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
c = OtoT_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
pie(c(a,b,c),main="ALS-Ox Clustering: Unspecified Motor Cortex vs Lumbar Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 58)
text(-0.5,0,"10",cex=6,col="white")
text(0.425,0.475,"1",cex=6,col="white")
text(0.5,-0.1,"4",cex=6,col="white")

#ALS-TD
a = GoodT_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
b = TtoG_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
c = TtoO_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
pie(c(a,b,c),main="ALS-TD Clustering: Unspecified Motor Cortex vs Lumbar Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 122)
text(-0.5,0,"3",cex=6,col="white")
text(0.25,0.45,"3",cex=6,col="white")
text(0.25,-0.4,"3",cex=6,col="white")


############################################## HISEQ ########################################################################################
#BREAK

CNSSubtype = HiSeqPheno
rownames(CNSSubtype) = CNSSubtype$Patient
CNSSubtype = CNSSubtype[,3:9]
#colnames(CNSSubtype)[c(30,31,32,33)] = c("FrontalCortex","MedialMotorCortex","LateralMotorCortex","UnspecMotorCortex")

for(i in 1:nrow(CNSSubtype)){
  for(j in 1:ncol(CNSSubtype)){
    if(is.na(CNSSubtype[i,j])){
      CNSSubtype[i,j] = ""
    }
  }
}

for(i in 1:nrow(CNSSubtype)){
  for(j in 1:ncol(CNSSubtype)){
    if(! is.na(CNSSubtype[i,j]) && CNSSubtype[i,j] == "TE"){
      CNSSubtype[i,j] = "TD"
    }
  }
}

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


ConcordantPercent_cerv = concountcerv/sum(concountcerv,discountcerv)
DiscordantPercent_cerv = 1-ConcordantPercent_cerv
pie(c(ConcordantPercent_cerv*100,DiscordantPercent_cerv*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 114,main="HiSeq Subtype Concordance - Cervical Spinal Cord - FC")
text(-0.5,0,"20",cex=6)
text(0.5,0,"34",cex=6)

ConcordantPercent_thor = concountthor/sum(concountthor,discountthor)
DiscordantPercent_thor = 1-ConcordantPercent_thor
pie(c(ConcordantPercent_thor*100,DiscordantPercent_thor*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 115,main="HiSeq Subtype Concordance - Thoracic Spinal Cord - FC")
text(-0.5,0,"13",cex=6)
text(0.5,0,"22",cex=6)

ConcordantPercent_lumb = concountlumb/sum(concountlumb,discountlumb)
DiscordantPercent_lumb = 1-ConcordantPercent_lumb
pie(c(ConcordantPercent_lumb*100,DiscordantPercent_lumb*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 110,main="HiSeq Subtype Concordance - Lumbar Spinal Cord - FC")
text(-0.5,0,"18",cex=6)
text(0.5,0,"29",cex=6)

################################# Cortex vs Cervical Spinal Cord

#ALS-Glia
a = GoodG_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
b = GtoO_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
c = GtoT_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
pie(c(a,b,c),main="ALS-Glia Clustering: Frontal Cortex vs Cervical Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 90)
text(-0.5,0,"2",cex=6,col="white")
text(0.5,0,"2",cex=6,col="white")

#ALS-Ox
a = GoodO_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
b = OtoG_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
c = OtoT_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
pie(c(a,b,c),main="ALS-Ox Clustering: Frontal Cortex vs Cervical Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 114)
text(-0.5,0,"14",cex=6,col="white")
text(0.35,0.25,"16",cex=6,col="white")
text(0.125,-0.45,"8",cex=6,col="white")

#ALS-TD
a = GoodT_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
b = TtoG_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
c = TtoO_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
pie(c(a,b,c),main="ALS-TD Clustering: Frontal Cortex vs Cervical Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 120)
text(-0.5,0,"4",cex=6,col="white")
text(0.15,0.5,"3",cex=6,col="white")
text(.35,-0.35,"5",cex=6,col="white")

################################# Cortex vs Thoracic Spinal Cord

#ALS-Glia 
a = GoodG_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
b = GtoO_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
c = GtoT_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
pie(c(a,b,c),main="ALS-Glia Clustering: Frontal Cortex vs Thoracic Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 90)
text(-0.5,0,"2",cex=6,col="white")
text(0.5,0,"2",cex=6,col="white")

#ALS-Ox
a = GoodO_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
b = OtoG_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
c = OtoT_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
pie(c(a,b,c),main="ALS-Ox Clustering: Frontal Cortex vs Thoracic Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 120)
text(-0.5,0,"8",cex=6,col="white")
text(0.35,0.25,"11",cex=6,col="white")
text(0.07,-0.475,"5",cex=6,col="white")

#ALS-TD
a = GoodT_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
b = TtoG_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
c = TtoO_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
pie(c(a,b,c),main="ALS-TD Clustering: Frontal Cortex vs Thoracic Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 104)
text(-0.5,0,"3",cex=6,col="white")
text(0.1,0.55,"1",cex=6,col="white")
text(.4,-0.2,"3",cex=6,col="white")

################################# Cortex vs Lumbar Spinal Cord

#ALS-Glia
a = GoodG_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
b = GtoO_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
c = GtoT_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
pie(c(a,b,c),main="ALS-Glia Clustering: Frontal Cortex vs Lumbar Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 0)
text(-0.5,0,"1",cex=6,col="white")

#ALS-Ox
a = GoodO_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
b = OtoG_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
c = OtoT_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
pie(c(a,b,c),main="ALS-Ox Clustering: Frontal Cortex vs Lumbar Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 98)
text(-0.5,0,"16",cex=6,col="white")
text(0.35,0.3,"11",cex=6,col="white")
text(.3,-0.4,"8",cex=6,col="white")

#ALS-TD
a = GoodT_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
b = TtoG_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
c = TtoO_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
pie(c(a,b,c),main="ALS-TD Clustering: Frontal Cortex vs Lumbar Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 148)
text(-0.5,0,"2",cex=6,col="white")
text(-0.05,0.5,"3",cex=6,col="white")
text(.3,-0.3,"6",cex=6,col="white")

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


ConcordantPercent_cerv = concountcerv/sum(concountcerv,discountcerv)
DiscordantPercent_cerv = 1-ConcordantPercent_cerv
pie(c(ConcordantPercent_cerv*100,DiscordantPercent_cerv*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 120,main="HiSeq Subtype Concordance - Cervical Spinal Cord - LMC")
text(-0.5,0,"15",cex=6)
text(0.5,0,"30",cex=6)

ConcordantPercent_thor = concountthor/sum(concountthor,discountthor)
DiscordantPercent_thor = 1-ConcordantPercent_thor
pie(c(ConcordantPercent_thor*100,DiscordantPercent_thor*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 103,main="HiSeq Subtype Concordance - Thoracic Spinal Cord - LMC")
text(-0.5,0,"14",cex=6)
text(0.5,0,"19",cex=6)

ConcordantPercent_lumb = concountlumb/sum(concountlumb,discountlumb)
DiscordantPercent_lumb = 1-ConcordantPercent_lumb
pie(c(ConcordantPercent_lumb*100,DiscordantPercent_lumb*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 125,main="HiSeq Subtype Concordance - Lumbar Spinal Cord - LMC")
text(-0.525,0,"12",cex=6)
text(0.5,0,"29",cex=6)

################################# Cortex vs Cervical Spinal Cord

#ALS-Glia
a = GoodG_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
b = GtoO_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
c = GtoT_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
pie(c(a,b,c),main="ALS-Glia Clustering: Lateral Motor Cortex vs Cervical Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 115)
text(-0.5,0,"4",cex=6,col="white")
text(0.4,0.25,"5",cex=6,col="white")
text(0.075,-0.45,"2",cex=6,col="white")

#ALS-Ox
a = GoodO_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
b = OtoG_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
c = OtoT_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
pie(c(a,b,c),main="ALS-Ox Clustering: Lateral Motor Cortex vs Cervical Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 110)
text(-0.5,0,"9",cex=6,col="white")
text(0.4,0.15,"11",cex=6,col="white")
text(0.035,-0.5,"3",cex=6,col="white")

#ALS-TD
a = GoodT_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
b = TtoG_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
c = TtoO_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
pie(c(a,b,c),main="ALS-TD Clustering: Lateral Motor Cortex vs Cervical Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 148)
text(-0.5,0,"2",cex=6,col="white")
text(0.05,0.45,"4",cex=6,col="white")
text(.2,-0.35,"5",cex=6,col="white")

################################# Cortex vs Thoracic Spinal Cord

#ALS-Glia
a = GoodG_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
b = GtoO_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
c = GtoT_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
pie(c(a,b,c),main="ALS-Glia Clustering: Lateral Motor Cortex vs Thoracic Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 90)
text(-0.5,0,"5",cex=6,col="white")
text(0.45,0.15,"4",cex=6,col="white")
text(0.18,-0.525,"1",cex=6,col="white")

#ALS-Ox
a = GoodO_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
b = OtoG_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
c = OtoT_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
pie(c(a,b,c),main="ALS-Ox Clustering: Lateral Motor Cortex vs Thoracic Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 108)
text(-0.5,0,"6",cex=6,col="white")
text(0.45,0.1,"8",cex=6,col="white")
text(-0.075,-0.6,"1",cex=6,col="white")


a = GoodT_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
b = TtoG_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
c = TtoO_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
pie(c(a,b,c),main="ALS-TD Clustering: Lateral Motor Cortex vs Thoracic Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 112)
text(-0.5,0,"3",cex=6,col="white")
text(0.2,0.5,"2",cex=6,col="white")
text(.35,-0.35,"3",cex=6,col="white")

################################# Cortex vs Lumbar Spinal Cord

#ALS-Glia
a = GoodG_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
b = GtoO_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
c = GtoT_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
pie(c(a,b,c),main="ALS-Glia Clustering: Lateral Motor Cortex vs Lumbar Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 112)
text(-0.5,0,"3",cex=6,col="white")
text(0.2,0.45,"2",cex=6,col="white")
text(0.35,-0.35,"3",cex=6,col="white")

#ALS-Ox
a = GoodO_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
b = OtoG_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
c = OtoT_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
pie(c(a,b,c),main="ALS-Ox Clustering: Lateral Motor Cortex vs Lumbar Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 110)
text(-0.5,0,"9",cex=6,col="white")
text(0.35,0.3,"9",cex=6,col="white")
text(0.185,-0.45,"5",cex=6,col="white")

#ALS-TD
a = GoodT_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
b = TtoG_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
c = TtoO_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
pie(c(a,b,c),main="ALS-TD Clustering: Lateral Motor Cortex vs Lumbar Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 73)
text(-0.5,0,"6",cex=6,col="white")
text(0.5,0,"4",cex=6,col="white")

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


ConcordantPercent_cerv = concountcerv/sum(concountcerv,discountcerv)
DiscordantPercent_cerv = 1-ConcordantPercent_cerv
pie(c(ConcordantPercent_cerv*100,DiscordantPercent_cerv*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 110,main="HiSeq Subtype Concordance - Cervical Spinal Cord - MMC")
text(-0.5,0,"16",cex=6)
text(0.5,0,"26",cex=6)

ConcordantPercent_thor = concountthor/sum(concountthor,discountthor)
DiscordantPercent_thor = 1-ConcordantPercent_thor
pie(c(ConcordantPercent_thor*100,DiscordantPercent_thor*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 113,main="HiSeq Subtype Concordance - Thoracic Spinal Cord - MMC")
text(-0.5,0,"12",cex=6)
text(0.5,0,"21",cex=6)

ConcordantPercent_lumb = concountlumb/sum(concountlumb,discountlumb)
DiscordantPercent_lumb = 1-ConcordantPercent_lumb
pie(c(ConcordantPercent_lumb*100,DiscordantPercent_lumb*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 105,main="HiSeq Subtype Concordance - Lumbar Spinal Cord - MMC")
text(-0.5,0,"16",cex=6)
text(0.5,0,"23",cex=6)


################################# Cortex vs Cervical Spinal Cord

#ALS-Glia
a = GoodG_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
b = GtoO_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
c = GtoT_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
pie(c(a,b,c),main="ALS-Glia Clustering: Medial Motor Cortex vs Cervical Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 145)
text(-0.5,0,"2",cex=6,col="white")
text(0.45,0.125,"7",cex=6,col="white")
text(-0.325,-0.45,"1",cex=6,col="white")

#ALS-Ox
a = GoodO_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
b = OtoG_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
c = OtoT_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
pie(c(a,b,c),main="ALS-Ox Clustering: Medial Motor Cortex vs Cervical Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 98)
text(-0.5,0,"10",cex=6,col="white")
text(0.35,0.35,"6",cex=6,col="white")
text(0.35,-0.35,"6",cex=6,col="white")

#ALS-TD
a = GoodT_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
b = TtoG_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
c = TtoO_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
pie(c(a,b,c),main="ALS-TD Clustering: Medial Motor Cortex vs Cervical Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 110)
text(-0.5,0,"4",cex=6,col="white")
text(0.15,0.5,"2",cex=6,col="white")
text(.4,-0.25,"4",cex=6,col="white")

################################# Cortex vs Thoracic Spinal Cord

#ALS-Glia
a = GoodG_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
b = GtoO_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
c = GtoT_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
pie(c(a,b,c),main="ALS-Glia Clustering: Medial Motor Cortex vs Thoracic Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 160)
text(-0.6,0,"1",cex=6,col="white")
text(0.35,0.25,"6",cex=6,col="white")
text(-0.225,-0.4,"2",cex=6,col="white")

#ALS-Ox
a = GoodO_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
b = OtoG_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
c = OtoT_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
pie(c(a,b,c),main="ALS-Ox Clustering: Medial Motor Cortex vs Thoracic Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 94)
text(-0.5,0,"7",cex=6,col="white")
text(0.4,0.25,"5",cex=6,col="white")
text(0.25,-0.425,"3",cex=6,col="white")

#ALS-TD
a = GoodT_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
b = TtoG_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
c = TtoO_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
pie(c(a,b,c),main="ALS-TD Clustering: Medial Motor Cortex vs Thoracic Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 100)
text(-0.5,0,"4",cex=6,col="white")
text(0.5,0,"5",cex=6,col="white")

################################# Cortex vs Lumbar Spinal Cord

#ALS-Glia
a = GoodG_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
b = GtoO_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
c = GtoT_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
pie(c(a,b,c),main="ALS-Glia Clustering: Medial Motor Cortex vs Lumbar Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 126)
text(-0.5,0,"3",cex=6,col="white")
text(0.3,0.4,"4",cex=6,col="white")
text(0.175,-0.45,"3",cex=6,col="white")

#ALS-Ox
a = GoodO_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
b = OtoG_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
c = OtoT_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
pie(c(a,b,c),main="ALS-Ox Clustering: Medial Motor Cortex vs Lumbar Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 90)
text(-0.45,0,"11",cex=6,col="white")
text(0.35,0.325,"6",cex=6,col="white")
text(0.35,-0.375,"5",cex=6,col="white")

#ALS-TD
a = GoodT_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
b = TtoG_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
c = TtoO_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
pie(c(a,b,c),main="ALS-TD Clustering: Medial Motor Cortex vs Lumbar Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 129)
text(-0.5,0,"2",cex=6,col="white")
text(0.5,0,"5",cex=6,col="white")

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


ConcordantPercent_cerv = concountcerv/sum(concountcerv,discountcerv)
DiscordantPercent_cerv = 1-ConcordantPercent_cerv
pie(c(ConcordantPercent_cerv*100,DiscordantPercent_cerv*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 150,main="HiSeq Subtype Concordance - Cervical Spinal Cord - Unspec. MC")
text(-0.5,0,"1",cex=6)
text(0.5,0,"5",cex=6)

ConcordantPercent_thor = concountthor/sum(concountthor,discountthor)
DiscordantPercent_thor = 1-ConcordantPercent_thor
pie(c(ConcordantPercent_thor*100,DiscordantPercent_thor*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 0,main="HiSeq Subtype Concordance - Thoracic Spinal Cord - Unspec. MC")
text(-0.5,0,"2",cex=6)

ConcordantPercent_lumb = concountlumb/sum(concountlumb,discountlumb)
DiscordantPercent_lumb = 1-ConcordantPercent_lumb
pie(c(ConcordantPercent_lumb*100,DiscordantPercent_lumb*100),col=c("#4fcdf0","#d65e78"),labels = c("Concordant","Discordant"),init.angle = 0,main="HiSeq Subtype Concordance - Lumbar Spinal Cord - Unspec. MC")
text(-0.5,0,"3",cex=6)

################################# Cortex vs Cervical Spinal Cord

#ALS-Glia
a = GoodG_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
b = GtoO_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
c = GtoT_cer/sum(GoodG_cer,GtoO_cer,GtoT_cer)
pie(c(a,b,c),main="ALS-Glia Clustering: Unspecified Motor Cortex vs Cervical Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 0)
text(-0.5,0,"1",cex=6,col="white")

#ALS-Ox
a = GoodO_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
b = OtoG_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
c = OtoT_cer/sum(GoodO_cer,OtoG_cer,OtoT_cer)
pie(c(a,b,c),main="ALS-Ox Clustering: Unspecified Motor Cortex vs Cervical Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 120)
text(-0.5,0,"1",cex=6,col="white")
text(0.5,0,"2",cex=6,col="white")

#ALS-TD
a = GoodT_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
b = TtoG_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
c = TtoO_cer/sum(GoodT_cer,TtoG_cer,TtoO_cer)
pie(c(a,b,c),main="ALS-TD Clustering: Unspecified Motor Cortex vs Cervical Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 90)
text(-0.5,0,"1",cex=6,col="white")
text(0.5,0,"1",cex=6,col="white")

################################# Cortex vs Thoracic Spinal Cord

#ALS-Glia
a = GoodG_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
b = GtoO_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
c = GtoT_tho/sum(GoodG_tho,GtoO_tho,GtoT_tho)
pie(c(a,b,c),main="ALS-Glia Clustering: Unspecified Motor Cortex vs Thoracic Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 0)
text(-0.5,0,"1",cex=6,col="white")

#ALS-Ox
a = GoodO_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
b = OtoG_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
c = OtoT_tho/sum(GoodO_tho,OtoG_tho,OtoT_tho)
pie(c(a,b,c),main="ALS-Ox Clustering: Unspecified Motor Cortex vs Thoracic Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 0)
text(-0.5,0,"1",cex=6,col = "white")

#ALS-TD - There are no ALS-TD Unspecified motor samples in the HiSeq Cohort
# a = GoodT_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
# b = TtoG_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
# c = TtoO_tho/sum(GoodT_tho,TtoG_tho,TtoO_tho)
# pie(c(a,b,c),main="ALS-TD Clustering: Unspecified Motor Cortex vs Thoracic Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 150)


################################# Cortex vs Lumbar Spinal Cord

#ALS-Glia
a = GoodG_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
b = GtoO_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
c = GtoT_lum/sum(GoodG_lum,GtoO_lum,GtoT_lum)
pie(c(a,b,c),main="ALS-Glia Clustering: Unspecified Motor Cortex vs Lumbar Spinal Cord",col = c("goldenrod1","navy","firebrick"),labels=c("Concordant","Glia (Cortex) to Ox (Spinal)","Glia (Cortex) to TD (Spinal)"),init.angle = 0)
text(-0.5,0,"1",cex=6,col="white")

#ALS-Ox
a = GoodO_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
b = OtoG_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
c = OtoT_lum/sum(GoodO_lum,OtoG_lum,OtoT_lum)
pie(c(a,b,c),main="ALS-Ox Clustering: Unspecified Motor Cortex vs Lumbar Spinal Cord",col = c("navy","goldenrod1","firebrick"),labels=c("Concordant","Ox (Cortex) to Glia (Spinal)","Ox (Cortex) to TD (Spinal)"),init.angle = 0)
text(-0.5,0,"1",cex=6,col="white")

#ALS-TD
a = GoodT_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
b = TtoG_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
c = TtoO_lum/sum(GoodT_lum,TtoG_lum,TtoO_lum)
pie(c(a,b,c),main="ALS-TD Clustering: Unspecified Motor Cortex vs Lumbar Spinal Cord",col = c("firebrick","goldenrod1","navy"),labels=c("Concordant","TD (Cortex) to Glia (Spinal)","TD (Cortex) to Ox (Spinal)"),init.angle = 0)
text(-0.5,0,"1",cex=6,col="white")
