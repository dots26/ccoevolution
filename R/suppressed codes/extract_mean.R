load('data_97000.Rdata')
dataLength <- length(res[1,])
sumVal <- 0
sumConv <- double(length(res[,1]$conv))
for(i in 1:dataLength){
  sumVal <- res[,i]$y + sumVal
  sumConv <- sumConv + res[,i]$conv
}

meanVal <- sumVal/dataLength
meanConv <- sumConv/dataLength

FE <- (1:10)*50000
FE <- c(1000,FE)

plot(FE,meanConv,type = 'line',col='red')


load('data_98000.Rdata')
dataLength <- length(res[1,])
sumVal <- 0
sumConv <- double(length(res[,1]$conv))
for(i in 1:dataLength){
  sumVal <- res[,i]$y + sumVal
  sumConv <- sumConv + res[,i]$conv
}

meanVal <- sumVal/dataLength
meanConv <- sumConv/dataLength

line(FE,meanConv,type = 'line',col='blue')

load('data_99000.Rdata')
dataLength <- length(res[1,])
sumVal <- 0
sumConv <- double(length(res[,1]$conv))
for(i in 1:dataLength){
  sumVal <- res[,i]$y + sumVal
  sumConv <- sumConv + res[,i]$conv
}

meanVal <- sumVal/dataLength
meanConv <- sumConv/dataLength

line(FE,meanConv,type = 'line',col='green')

load('data_100000.Rdata')
dataLength <- length(res[1,])
sumVal <- 0
sumConv <- double(length(res[,1]$conv))
for(i in 1:dataLength){
  sumVal <- res[,i]$y + sumVal
  sumConv <- sumConv + res[,i]$conv
}

meanVal <- sumVal/dataLength
meanConv <- sumConv/dataLength

line(FE,meanConv,type = 'line',col='black')
