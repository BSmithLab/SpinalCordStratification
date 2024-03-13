#User Input - SQuIRE Post-Processing

wd = "E:/NovaSeqTE/NovaSeq/Controls"
setwd(wd)

FileList = list.files(path = wd, pattern = "_TEcounts.txt") #Files obtained from SQuIRE output

chunknumber = 2
start=56
stop=109
outputdir = "E:/NovaSeqTE/NovaSeq/Controls/Chunks"

file = paste("SQuIRE_Post_Chunk",chunknumber,".RData",sep="")

#################################################################################

range = length(start:stop)

library(stringr)
accessions = str_sub(FileList,1,11)
TEFiles = names(table(accessions))
TEFinalList = TEFiles[!is.na(TEFiles)]

##Read in Subjects and Filter
lcounter = 1
for (i in start:stop){
  TECounts = read.delim(FileList[i],fill = T,header = T,skipNul = T)
  FilteredTECounts = TECounts[, c("TE_ID", "score", "tot_counts")]
  SubjectName = paste("Subject",i, sep = "")
  
  assign(SubjectName,FilteredTECounts)
  Subject = get(SubjectName)
  Remove = rep(NA,nrow(Subject))
  count = 1
  for (k in 1:nrow(Subject)) {
    if(Subject$score[k] < 99){
      Remove[count] = k
      count = count + 1
    }
  }
  Remove = Remove[!is.na(Remove)]
  assign(SubjectName,Subject[-Remove,])
  if((i %% 1) == 0) cat("Read-In % Completed:", (lcounter/range)*100,"\n",sep="")
  lcounter = lcounter+1
}


lcounter = 1
for(i in start:stop){
  #If this is the first time we are executing this loop
  if(i == start) {
    tmp = get(paste("Subject",i, sep = ""))
    tmp2 = get(paste("Subject",i+1,sep=""))
    snow = tmp[tmp$TE_ID %in% tmp2$TE_ID,]
  }else if(i > start && i < stop){
    tmp = snow 
    tmp2 = get(paste("Subject",i+1,sep=""))
    snow = tmp[tmp$TE_ID %in% tmp2$TE_ID,]
  }
  
  if(i == stop){
    assign(paste("finalavalanche",chunknumber,sep=""),snow)
  }
  
  if((i %% 10) == 0) assign(paste("avalanche",i,sep=""),snow)
  if((i %% 2) == 0) cat("Avalanche % Completed:", (lcounter/range)*100,"\n",sep="")
  lcounter = lcounter+1
}



setwd(outputdir)
save.image(file=file)