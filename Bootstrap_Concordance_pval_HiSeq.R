
nboot = 10000
PatientPheno = read.csv("C:/Users/Jarrett Eshima/Documents/Research/Dissertation/Concordance/SpinalCord_Concordance_HiSeq_9-18-23.csv")
rownames(PatientPheno) = PatientPheno$Patient
PatientPheno = PatientPheno[,2:8]

for(i in 1:nrow(PatientPheno)){
  for(j in 1:ncol(PatientPheno)){
    if(!is.na(PatientPheno[i,j]) && PatientPheno[i,j] == "TE"){
      PatientPheno[i,j] = "TD"
    }
  }
}

subtypes = c("OX","GLIA","TD")
cortexprobs = c(3/5.5,1/5.5,1.5/5.5)
spinalprobs = c(2/4.5,1/4.5,1.5/4.5)

################################## FRONTAL CORTEX ########################################

FC_Cervical_nConcordant = FC_Cervical_nDiscordant = FC_Thoracic_nConcordant = FC_Thoracic_nDiscordant = FC_Lumbar_nConcordant = FC_Lumbar_nDiscordant = rep(NA,nboot)

for(b in 1:nboot){
  
  currentPheno = PatientPheno
  
  for(j in 1:nrow(currentPheno)){
    for(k in 1:ncol(currentPheno)){
      
      if(!is.na(currentPheno[j,k]) && k <= 4 ){
        currentPheno[j,k] = sample(subtypes,1,prob = cortexprobs)
      }else if(! is.na(currentPheno[j,k]) && k > 4){
        currentPheno[j,k] = sample(subtypes,1,prob = spinalprobs)
      }else if(is.na(currentPheno[j,k])){
        currentPheno[j,k] = ""
      }
      
    }
  }
  
  CNSSubtype = currentPheno
  
  #Now compare tissue level concordance, use previous code.
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
  
  FC_Cervical_nConcordant[b] = concountcerv
  FC_Cervical_nDiscordant[b] = discountcerv
  
  FC_Thoracic_nConcordant[b] = concountthor
  FC_Thoracic_nDiscordant[b] = discountthor
  
  FC_Lumbar_nConcordant[b] = concountlumb
  FC_Lumbar_nDiscordant[b] = discountlumb
  
  if((b %% 100) == 0) cat("% Done:",(b/nboot)*100,"\n")
  
}


par(mfrow=c(1,2))
##Cervical Spinal Cord
hist(FC_Cervical_nConcordant,main="Frontal Cortex vs Cervical Spinal",breaks = 20)
abline(v=20,col="red") #This is the actual number of concordant sample pairings
#get the pvalue associated with the zscore from the bootstrapped distribution. Distribution is roughly normal for this nboot. One tailed bc expected to be greater number of concordant obs. than by chance
pnorm(q=(20-mean(FC_Cervical_nConcordant))/sd(FC_Cervical_nConcordant) ,lower.tail = F)

hist(FC_Cervical_nDiscordant,main="Frontal Cortex vs Cervical Spinal",breaks = 20)
abline(v=34,col="red")
#get the pvalue associated with the zscore from the bootstrapped distribution. One tailed bc expected to be lower number of discordant obs. than by chance
pnorm(q=(34-mean(FC_Cervical_nDiscordant))/sd(FC_Cervical_nDiscordant) ,lower.tail = T)

### Thoracic Spinal Cord
hist(FC_Thoracic_nConcordant,main="Frontal Cortex vs Thoracic Spinal",breaks = 20)
abline(v=13,col="red")
pnorm(q=(13-mean(FC_Thoracic_nConcordant))/sd(FC_Thoracic_nConcordant) ,lower.tail = F)
hist(FC_Thoracic_nDiscordant,main="Frontal Cortex vs Thoracic Spinal",breaks = 20)
abline(v=22,col="red")
pnorm(q=(22-mean(FC_Thoracic_nDiscordant))/sd(FC_Thoracic_nDiscordant) ,lower.tail = T)

### Lumbar Spinal Cord
hist(FC_Lumbar_nConcordant,main="Frontal Cortex vs Lumbar Spinal",breaks = 20)
abline(v=18,col="red")
pnorm(q=(18-mean(FC_Lumbar_nConcordant))/sd(FC_Lumbar_nConcordant) ,lower.tail = F)
hist(FC_Lumbar_nDiscordant,main="Frontal Cortex vs Lumbar Spinal",breaks = 20)
abline(v=29,col="red")
pnorm(q=(29-mean(FC_Lumbar_nDiscordant))/sd(FC_Lumbar_nDiscordant) ,lower.tail = T)

par(mfrow=c(1,1))


################################## MEDIAL MOTOR CORTEX ########################################
nboot = 10000
PatientPheno = read.csv("C:/Users/Jarrett Eshima/Documents/Research/Dissertation/Concordance/SpinalCord_Concordance_HiSeq_9-18-23.csv")
rownames(PatientPheno) = PatientPheno$Patient
PatientPheno = PatientPheno[,2:8]

for(i in 1:nrow(PatientPheno)){
  for(j in 1:ncol(PatientPheno)){
    if(!is.na(PatientPheno[i,j]) && PatientPheno[i,j] == "TE"){
      PatientPheno[i,j] = "TD"
    }
  }
}

subtypes = c("OX","GLIA","TD")
cortexprobs = c(3/5.5,1/5.5,1.5/5.5)
spinalprobs = c(2/4.5,1/4.5,1.5/4.5)
MM_Cervical_nConcordant = MM_Cervical_nDiscordant = MM_Thoracic_nConcordant = MM_Thoracic_nDiscordant = MM_Lumbar_nConcordant = MM_Lumbar_nDiscordant = rep(NA,nboot)

for(b in 1:nboot){
  
  currentPheno = PatientPheno
  
  for(j in 1:nrow(currentPheno)){
    for(k in 1:ncol(currentPheno)){
      
      if(!is.na(currentPheno[j,k]) && k <= 4 ){
        currentPheno[j,k] = sample(subtypes,1,prob = cortexprobs)
      }else if(! is.na(currentPheno[j,k]) && k > 4){
        currentPheno[j,k] = sample(subtypes,1,prob = spinalprobs)
      }else if(is.na(currentPheno[j,k])){
        currentPheno[j,k] = ""
      }
      
    }
  }
  
  CNSSubtype = currentPheno
  
  #Now compare tissue level concordance, use previous code.
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
  
  MM_Cervical_nConcordant[b] = concountcerv
  MM_Cervical_nDiscordant[b] = discountcerv
  
  MM_Thoracic_nConcordant[b] = concountthor
  MM_Thoracic_nDiscordant[b] = discountthor
  
  MM_Lumbar_nConcordant[b] = concountlumb
  MM_Lumbar_nDiscordant[b] = discountlumb
  
  if((b %% 100) == 0) cat("% Done:",(b/nboot)*100,"\n")
  
}


par(mfrow=c(1,2))
hist(MM_Cervical_nConcordant,main="Medial Motor Cortex vs Cervical Spinal",breaks = 20)
abline(v=16,col="red")
pnorm(q=(16-mean(MM_Cervical_nConcordant))/sd(MM_Cervical_nConcordant) ,lower.tail = F)
hist(MM_Cervical_nDiscordant,main="Medial Motor Cortex vs Cervical Spinal",breaks = 20)
abline(v=26,col="red")
pnorm(q=(26-mean(MM_Cervical_nDiscordant))/sd(MM_Cervical_nDiscordant) ,lower.tail = T)

hist(MM_Thoracic_nConcordant,main="Medial Motor Cortex vs Thoracic Spinal",breaks = 20)
abline(v=12,col="red")
pnorm(q=(12-mean(MM_Thoracic_nConcordant))/sd(MM_Thoracic_nConcordant) ,lower.tail = F)
hist(MM_Thoracic_nDiscordant,main="Medial Motor Cortex vs Thoracic Spinal",breaks = 20)
abline(v=21,col="red")
pnorm(q=(21-mean(MM_Thoracic_nDiscordant))/sd(MM_Thoracic_nDiscordant) ,lower.tail = T)

hist(MM_Lumbar_nConcordant,main="Medial Motor Cortex vs Lumbar Spinal",breaks = 20)
abline(v=16,col="red")
pnorm(q=(16-mean(MM_Lumbar_nConcordant))/sd(MM_Lumbar_nConcordant) ,lower.tail = F)
hist(MM_Lumbar_nDiscordant,main="Medial Motor Cortex vs Lumbar Spinal",breaks = 20)
abline(v=23,col="red")
pnorm(q=(23-mean(MM_Lumbar_nDiscordant))/sd(MM_Lumbar_nDiscordant) ,lower.tail = T)

par(mfrow=c(1,1))

################################## LATERAL MOTOR CORTEX ########################################

nboot = 10000
PatientPheno = read.csv("C:/Users/Jarrett Eshima/Documents/Research/Dissertation/Concordance/SpinalCord_Concordance_HiSeq_9-18-23.csv")
rownames(PatientPheno) = PatientPheno$Patient
PatientPheno = PatientPheno[,2:8]

for(i in 1:nrow(PatientPheno)){
  for(j in 1:ncol(PatientPheno)){
    if(!is.na(PatientPheno[i,j]) && PatientPheno[i,j] == "TE"){
      PatientPheno[i,j] = "TD"
    }
  }
}

subtypes = c("OX","GLIA","TD")
cortexprobs = c(3/5.5,1/5.5,1.5/5.5)
spinalprobs = c(2/4.5,1/4.5,1.5/4.5)
LM_Cervical_nConcordant = LM_Cervical_nDiscordant = LM_Thoracic_nConcordant = LM_Thoracic_nDiscordant = LM_Lumbar_nConcordant = LM_Lumbar_nDiscordant = rep(NA,nboot)

for(b in 1:nboot){
  
  currentPheno = PatientPheno
  
  for(j in 1:nrow(currentPheno)){
    for(k in 1:ncol(currentPheno)){
      
      if(!is.na(currentPheno[j,k]) && k <= 4 ){
        currentPheno[j,k] = sample(subtypes,1,prob = cortexprobs)
      }else if(! is.na(currentPheno[j,k]) && k > 4){
        currentPheno[j,k] = sample(subtypes,1,prob = spinalprobs)
      }else if(is.na(currentPheno[j,k])){
        currentPheno[j,k] = ""
      }
      
    }
  }
  
  CNSSubtype = currentPheno
  
  #Now compare tissue level concordance, use previous code.
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
  
  LM_Cervical_nConcordant[b] = concountcerv
  LM_Cervical_nDiscordant[b] = discountcerv
  
  LM_Thoracic_nConcordant[b] = concountthor
  LM_Thoracic_nDiscordant[b] = discountthor
  
  LM_Lumbar_nConcordant[b] = concountlumb
  LM_Lumbar_nDiscordant[b] = discountlumb
  
  if((b %% 100) == 0) cat("% Done:",(b/nboot)*100,"\n")
  
}


par(mfrow=c(1,2))
hist(LM_Cervical_nConcordant,main="Lateral Motor Cortex vs Cervical Spinal",breaks = 20)
abline(v=15,col="red")
pnorm(q=(15-mean(LM_Cervical_nConcordant))/sd(LM_Cervical_nConcordant) ,lower.tail = F)
hist(LM_Cervical_nDiscordant,main="Lateral Motor Cortex vs Cervical Spinal",breaks = 20)
abline(v=30,col="red")
pnorm(q=(30-mean(LM_Cervical_nDiscordant))/sd(LM_Cervical_nDiscordant) ,lower.tail = T)

hist(LM_Thoracic_nConcordant,main="Lateral Motor Cortex vs Thoracic Spinal",breaks = 20)
abline(v=14,col="red")
pnorm(q=(14-mean(LM_Thoracic_nConcordant))/sd(LM_Thoracic_nConcordant) ,lower.tail = F)
hist(LM_Thoracic_nDiscordant,main="Lateral Motor Cortex vs Thoracic Spinal",breaks = 20)
abline(v=19,col="red")
pnorm(q=(19-mean(LM_Thoracic_nDiscordant))/sd(LM_Thoracic_nDiscordant) ,lower.tail = T)

hist(LM_Lumbar_nConcordant,main="Lateral Motor Cortex vs Lumbar Spinal",breaks = 20)
abline(v=12,col="red")
pnorm(q=(12-mean(LM_Lumbar_nConcordant))/sd(LM_Lumbar_nConcordant) ,lower.tail = F)
hist(LM_Lumbar_nDiscordant,main="Lateral Motor Cortex vs Lumbar Spinal",breaks = 20)
abline(v=29,col="red")
pnorm(q=(29-mean(LM_Lumbar_nDiscordant))/sd(LM_Lumbar_nDiscordant) ,lower.tail = T)

par(mfrow=c(1,1))

################################## Unspecified MOTOR CORTEX ########################################

nboot = 10000
PatientPheno = read.csv("C:/Users/Jarrett Eshima/Documents/Research/Dissertation/Concordance/SpinalCord_Concordance_HiSeq_9-18-23.csv")
rownames(PatientPheno) = PatientPheno$Patient
PatientPheno = PatientPheno[,2:8]

for(i in 1:nrow(PatientPheno)){
  for(j in 1:ncol(PatientPheno)){
    if(!is.na(PatientPheno[i,j]) && PatientPheno[i,j] == "TE"){
      PatientPheno[i,j] = "TD"
    }
  }
}

subtypes = c("OX","GLIA","TD")
cortexprobs = c(3/5.5,1/5.5,1.5/5.5)
spinalprobs = c(2/4.5,1/4.5,1.5/4.5)
UM_Cervical_nConcordant = UM_Cervical_nDiscordant = UM_Thoracic_nConcordant = UM_Thoracic_nDiscordant = UM_Lumbar_nConcordant = UM_Lumbar_nDiscordant = rep(NA,nboot)

for(b in 1:nboot){
  
  currentPheno = PatientPheno
  
  for(j in 1:nrow(currentPheno)){
    for(k in 1:ncol(currentPheno)){
      
      if(!is.na(currentPheno[j,k]) && k <= 4 ){
        currentPheno[j,k] = sample(subtypes,1,prob = cortexprobs)
      }else if(! is.na(currentPheno[j,k]) && k > 4){
        currentPheno[j,k] = sample(subtypes,1,prob = spinalprobs)
      }else if(is.na(currentPheno[j,k])){
        currentPheno[j,k] = ""
      }
      
    }
  }
  
  CNSSubtype = currentPheno
  
  #Now compare tissue level concordance, use previous code.
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
  
  UM_Cervical_nConcordant[b] = concountcerv
  UM_Cervical_nDiscordant[b] = discountcerv
  
  UM_Thoracic_nConcordant[b] = concountthor
  UM_Thoracic_nDiscordant[b] = discountthor
  
  UM_Lumbar_nConcordant[b] = concountlumb
  UM_Lumbar_nDiscordant[b] = discountlumb
  
  if((b %% 100) == 0) cat("% Done:",(b/nboot)*100,"\n")
  
}


par(mfrow=c(1,2))
hist(UM_Cervical_nConcordant,main="Unspecified Motor Cortex vs Cervical Spinal",breaks = 20)
abline(v=1,col="red")
pnorm(q=(1-mean(UM_Cervical_nConcordant))/sd(UM_Cervical_nConcordant) ,lower.tail = F)
hist(UM_Cervical_nDiscordant,main="Unspecified Motor Cortex vs Cervical Spinal",breaks = 20)
abline(v=5,col="red")
pnorm(q=(5-mean(UM_Cervical_nDiscordant))/sd(UM_Cervical_nDiscordant) ,lower.tail = T)

hist(UM_Thoracic_nConcordant,main="Unspecified Motor Cortex vs Thoracic Spinal",breaks = 20)
abline(v=0,col="red")
pnorm(q=(0-mean(UM_Thoracic_nConcordant))/sd(UM_Thoracic_nConcordant) ,lower.tail = F)
hist(UM_Thoracic_nDiscordant,main="Unspecified Motor Cortex vs Thoracic Spinal",breaks = 20)
abline(v=2,col="red")
pnorm(q=(2-mean(UM_Thoracic_nDiscordant))/sd(UM_Thoracic_nDiscordant) ,lower.tail = T)

hist(UM_Lumbar_nConcordant,main="Unspecified Motor Cortex vs Lumbar Spinal",breaks = 20)
abline(v=0,col="red")
pnorm(q=(0-mean(UM_Lumbar_nConcordant))/sd(UM_Lumbar_nConcordant) ,lower.tail = F)
hist(UM_Lumbar_nDiscordant,main="Unspecified Motor Cortex vs Lumbar Spinal",breaks = 20)
abline(v=3,col="red")
pnorm(q=(3-mean(UM_Lumbar_nDiscordant))/sd(UM_Lumbar_nDiscordant) ,lower.tail = T)

par(mfrow=c(1,1))

#setwd("C:/Users/Jarrett Eshima/Documents/Research/Dissertation")
#save.image("Bootstrap_Concordance_10kiter_1-19-24_HiSeq.RData")
load("C:/Users/Jarrett Eshima/Documents/Research/Dissertation/Concordance/Bootstrap_Concordance_10kiter_1-19-24_HiSeq.RData")
