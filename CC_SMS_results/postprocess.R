rm(list=ls())
meanHV <- NULL
nIterA <- c(1,250,500,1000)

nPoint <- 1000
# paretoSetSample <- matrix((randtoolbox::sobol(n = nPoint,dim = k,scrambling = 3))*(1:k)*2,nrow=k)
# paretoSetSample <- t(matrix((randtoolbox::sobol(k*nPoint))*(1:k)*2,ncol=k))
# paretoSetSample <- matrix((runif(k*nPoint)^50)*(1:k)*2,nrow=k)
#paretoSetSample <- matrix((runif(k*nPoint))*(1:k)*2,nrow=k)
paretoSetSample <- matrix(1:nPoint * (1/nPoint),nrow=1)*(1:k)*2 # wfg2
aa <- pracma::repmat(rep(0,nPoint),k-1,1)
paretoSetSample <- rbind(paretoSetSample,aa)
bb <- pracma::repmat(matrix(0.7*((k+1):nVar),ncol=1),1,nPoint)
paretoSetSample <- rbind(paretoSetSample,bb)
paretoFrontSample <- func(individual = paretoSetSample,nObj = 2,k = 20)
xx <- (0:1000)/1000
yy <- (1-xx*xx)^0.5
xx <- paretoFrontSample[1,order(paretoFrontSample[1,])]
yy <- paretoFrontSample[2,order(paretoFrontSample[1,])]

for (nIter in 250)
for(ID in c(3)){
  sumHV <- 0
  count <- 0
  # instanceCount <-  length(list.files(path = paste0("global_stepsize/WFG2/test",ID)))
  for(instance in 1:10){
    # print(paste(ID,instance))
    load(paste0("WFG2/fixed_step/",nIter,"iter_ID",ID,"_instance",instance,".Rdata"))
    sumHV <- sumHV + abc$hv_log
    print(abc$hv_log)
    count <- count + 1
    plot(xx,yy,xlim=c(0,6),ylim=c(0,8),type="l")
    points(abc$y)
    title(paste(ID,instance))
  }
  meanHV <- append(meanHV,sumHV/count)
}
View(meanHV)
