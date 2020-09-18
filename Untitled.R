myfun <- function(individual, nObj,k = nObj-1){
  M <- nObj
  if(is.vector(individual))
    individual <- matrix(individual,ncol=1)
  nIndividual <- ncol(individual)
  n <- nrow(individual) # number of variables
  l <- n-k
  
  individual1 <- individual
  individual2 <- individual
  individual3 <- individual
  x <- pracma::zeros(M,nIndividual)
  h <- x
  individual <- individual/seq(2,2*n,2)
  
  # first transformation
  for(i in 1:k){
    individual1[i, ] <- individual[i, ]
  }
  for(i in (k+1):n){
    individual1[i, ] <- MaOEA:::s_linear(individual[i, ],0.35)
  }
  print(individual1)
  # second transform
  for(i in 1:k){
    individual2[i, ] <- individual1[i, ]
  }
  for(i in (k+1):n){
    individual2[i, ] <- MaOEA:::b_flat(individual1[i, ],0.8,0.75,0.85)
  }
  
  # third transform
  for(i in 1:n){
    individual3[i, ] <- MaOEA:::b_poly(individual2[i, ],0.02)
  }
  # fourth transform
  for(i in 1:(M-1)){
    rsumMinIndex <- (i-1)*k/(M-1)+1
    rsumMaxIndex <- (i*k)/(M-1)
    weightMin <- 2*rsumMinIndex
    weightMax <- 2*rsumMaxIndex
    
    weightVector <- seq(weightMin,weightMax,2)
    
    x[i, ] <- MaOEA:::r_sum(individual3[rsumMinIndex:rsumMaxIndex,,drop=F],weightVector)
  }
  rsumMinIndex <- k+1
  rsumMaxIndex <- n
  weightMin <- 2*rsumMinIndex
  weightMax <- 2*rsumMaxIndex
  
  weightVector <- seq(weightMin,weightMax,2)
  
  x[M, ] <- MaOEA:::r_sum(individual3[rsumMinIndex:rsumMaxIndex,,drop=F],weightVector)
  
  # shape function
  for(i in 1:(M-1)){
    h[i, ] <- MaOEA:::shape_convex(M,i,x)
  }
  h[M, ] <- MaOEA:::shape_mixed(M,x,1,5)
  
  S <- seq(2,2*M,2)
  
  obj_val <- x[M,] + h*S
  return(obj_val)
}


WFG2 <- function(individual, nObj,k = nObj-1){
  M <- nObj
  if(is.vector(individual))
    individual <- matrix(individual,ncol=1)
  nIndividual <- ncol(individual)
  n <- nrow(individual) # number of variables
  l <- n-k
  
  
  individual1 <- individual
  individual2 <- individual
  individual3 <- individual
  x <- pracma::zeros(M,nIndividual)
  h <- x
  
  individual <- individual/seq(2,2*n,2)
  # first transformation
  for(i in 1:k){
    individual1[i, ] <- individual[i, ]
  }
  for(i in (k+1):n){
    individual1[i, ] <- MaOEA:::s_linear(individual[i, ],0.35)
  }
  
  # second transform
  for(i in 1:k){
    individual2[i, ] <- individual1[i, ]
  }
  for(i in (k+1):(k+ (l/2) ) ){
    individual2[i, ] <- MaOEA:::r_nonsep(c(individual1[k+2*(i-k)-1, ],individual1[k+2*(i-k), ]),2)
  }
  
  # third transform
  for(i in 1:(M-1)){
    rsumMinIndex <- (i-1)*k/(M-1)+1
    rsumMaxIndex <- (i*k)/(M-1)
    weightVector <- rep(1,rsumMaxIndex-rsumMinIndex+1)
    
    x[i, ] <- MaOEA:::r_sum(individual2[rsumMinIndex:rsumMaxIndex,,drop=F],weightVector)
  }
  rsumMinIndex <- k+1
  rsumMaxIndex <- k+l/2
  
  weightVector <- rep(1,rsumMaxIndex-rsumMinIndex+1)
  x[M, ] <- MaOEA:::r_sum(individual2[rsumMinIndex:rsumMaxIndex,,drop=F],weightVector)
  
  
  # shape function
  for(i in 1:(M-1)){
    h[i, ] <- MaOEA:::shape_convex(M,i,x)
  }
  h[M, ] <- MaOEA:::shape_disconnected(M,x,1,1,5)
  
  S <- seq(2,2*M,2)
  
  obj_val <- x[M,] + h*S
  return(obj_val)
}

b_flat <- function(y,A,B,C){
  fun1 <- function(y1,A,B,C){
    return( min(c(0,floor(y1-B)))*(A-(A*y1/B)))
  }
  fun2 <- function(y1,A,B,C){
    return( min(c(0,floor(C-y1)))*(1-A)*(y1-C)/(1-C))
  }
  
  xx <- sapply(X = y,FUN = fun1,A=A,B=B,C=C)
  yy <- sapply(X = y,FUN = fun2,A=A,B=B,C=C)
  
  x <- A + xx - yy
  return(x)
}
