# r<- 20
# a <- sensitivity::morris(model=tfun,
#                          factor=nVar,
#                          r = r,
#                          design = list(type='oat',levels=8,grid.jump=1),
#                          binf=lbound,
#                          bsup=ubound,
#                          scale=F,mainfun=fun,...)
# bestPopIndex <- which.min(a$y)
# bestPop <- a$X[bestPopIndex,]
# bestObj <- min(a$y)
#
# contextVector <- bestPop
# nEval <- nEval + r*(nVar+1)
# mu.star <- apply(a$ee, 2, function(a) mean(abs(a)))
# sigma <- apply(a$ee, 2, sd)
#
# print('Using k-means...') # clustering
# assignedGroup <- kmeans(cbind((mu.star),(sigma)),nLevel)$cluster
#
# EF <- NULL
# group <- NULL
# sigmaMuAll <- sum(mu.star)
# for(groupIndex in 1:nLevel){
#   sigmaMu <- sum(mu.star[which(assignedGroup==groupIndex)])
#   EF <- append(EF,sigmaMu/sigmaMuAll)
#   group <- append(group,list(which(assignedGroup==groupIndex)))
# }
# clusterOrder <- order(EF,decreasing = T)
# save(group,file='MOFBVE_group.Rdata')
# # }
# print('Groups assigned')
# print(group[clusterOrder])
# error checking on groups
