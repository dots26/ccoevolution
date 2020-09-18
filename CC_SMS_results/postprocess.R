rm(list=ls())
meanHV <- NULL
for(ID in c(2)){
  sumHV <- 0
  count <- 0
  instanceCount <-  length(list.files(path = paste0("global_stepsize/WFG2/test",ID)))
  for(instance in 1){
    print(paste(ID,instance))
    load(paste0("global_stepsize/WFG2/test",ID,"/instance",instance,".Rdata"))
    sumHV <- sumHV + abc$hv_log
    count <- count + 1
  }
  meanHV <- append(meanHV,sumHV/count)
}
