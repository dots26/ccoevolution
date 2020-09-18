meanHV <- NULL
for(ID in 7){
  sumHV <- 0
  for(instance in 1:10){
    print(paste(ID,instance))
    fna <- paste0("1iter_ID",ID,"_instance",instance,".Rdata")
    load(paste0("scaled_down/",fna))
    plot(0,0,xlim=c(0,6),ylim=c(0,10),type="l")
    title(paste(ID,instance))
    points(abc$y)
    sumHV <- sumHV + abc$hv_log
  }
  meanHV <- append(meanHV,sumHV/10)
}
