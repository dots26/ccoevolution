f_elliptic <- function(x){
  if(is.vector(x))
    x <- matrix(x,ncol=length(x))

  D <- ncol(x)
  i <- 1:D
  conditioning <- matrix(10^(6*(i-1)/(D-1)),nrow=1)

  if(nrow(x)>1){
    res <- (conditioning%*%t(x*x))
  }else{
    res <- sum(conditioning*x*x)
  }
  return(res)
}

f_ackley <- function(z){
  if(is.vector(z))
    z <- matrix(z,ncol=length(z))

  D <- ncol(z)
  if(nrow(z)>1){
    ackley <- -20 * exp(-0.2 * (1/D*rowSums(z*z))^0.5 ) - exp(1/D*rowSums(cos(2*pi*z)) ) + 20 + exp(1)
  }else{
    ackley <- -20 * exp(-0.2 * (1/D*sum(z*z))^0.5 ) - exp(1/D*sum(cos(2*pi*z)) ) + 20 + exp(1)
  }
  return(ackley)
}

f_schewefel_1_2 <- function(x){
  if(is.vector(x))
    x <- matrix(x,ncol=length(x))

  D <- ncol(x)

  summand <- numeric(nrow(x))
  prevSum <- summand
  if(nrow(x)>1){
    for(i in 1:D){
      prevSum <- prevSum + x[,i]
      summand <- summand + ( prevSum* prevSum)
    }
  }else{
    for(i in 1:D){
      prevSum <- prevSum + x[,i]
      summand <- summand + (  prevSum* prevSum)
    }
  }
  return(summand)
}

t_osz <- function(x){
  xlog <- log(abs(x))* as.integer(x!=0)

  xlog[which(is.nan(xlog))] <- 0

  c1 <- 10*as.integer(x>0) +  5.5*as.integer(x<=0)
  c2 <- 7.9*as.integer(x>0) +  3.1*as.integer(x<=0)

  t <- sign(x) * exp( xlog + 0.0049*sin(c1*xlog) + sin(c2*xlog ))
}

t_asy <- function(x,beta){
  if(is.vector(x))
    x <- matrix(x,ncol=length(x))

  D <- ncol(x)

  i <- t(matrix(rep(1:D,nrow(x)),nrow =ncol(x)))
  beta_term <- beta * ((i-1)/(D-1))
  t1 <- x^(1+beta_term*(x^0.5))*as.integer(x>0)

  t2 <- x*as.integer(x<=0)

  t1[which(is.nan(t1))]  <- 0
  t <- t1+t2

  return(t)
}

cec_weight <- function(D){
  power_term <- 1*rnorm(D,0,1)
  return(10^power_term)
}

ill_lambda <- function(alpha,D){
  i <- 1:D
  power_term <- 0.5*(i-1)/(D-1)
  diagonal <- diag(D)%*%(alpha^power_term)
}

f_rastrigin <- function(x){
  if(is.vector(x))
    x <- matrix(x,ncol=length(x))
  if(nrow(x)>1){
    res <- rowSums(x*x - 10*cos(2*pi*x) + 10)
  }else{
    res <- sum(x*x - 10*cos(2*pi*x) + 10)
  }
  return(res)
}
