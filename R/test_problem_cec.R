#' CEC 2013 f1 function. Can accomodate different variable dimensionality m>0. The CEC2013 test problem uses m=1000.
#'
#' @title Shifted elliptic function.
#' @param x design variable. One row per point.
#' @param o the optimum
#' @examples
#' x <- runif(1000)
#' o <- double(1000)
#' f1cec(x,o)
#' @export
f1cec <- function(x,o=NULL){ # shifted elliptic
  if(is.vector(x))
    x <- matrix(x,ncol=length(x))

  if(is.null(o))
    o <- double(ncol(x))

  D <- ncol(x)
  z <- x-t(matrix( rep(o,nrow(x) ),nrow=D) )
  z <- t_osz(z)
  res <- f_elliptic(z)

  return(res)
}

#' CEC 2013 f2 function.
#'
#' @title Shifted Rastrigin's function with m>0 Dimensions variables. The CEC2013 test problem uses m=1000
#' @param x design variable
#' @param o the optimum
#' @examples
#' x <- runif(1000)
#' o <- double(1000)
#' f2cec(x,o)
#' @export
f2cec <- function(x,o=NULL){ # shifted rastrigin
  if(is.vector(x))
    x <- matrix(x,ncol=length(x))

  if(is.null(o))
    o <- double(ncol(x))

  z <- x-t(matrix( rep(o,nrow(x) ),nrow=ncol(x)) )
  z <- t_osz(z)
  z1 <- as.vector (ill_lambda(10,ncol(z)))

  z <- t(t_asy(z,0.2))*z1

  res <- f_rastrigin(t(z))

  return(res)
}

#' CEC 2010 f3 function.
#'
#' @title Shifted Ackley's function with m>0 variables. The CEC2013 test problem uses m=1000.
#' @param x design variable
#' @param o the optimum
#' @examples
#' x <- runif(1000)
#' o <- double(1000)
#' f3cec(x,o)
#' @export
f3cec <- function(x,o=NULL){ # shifted ackley
  if(is.vector(x))
    x <- matrix(x,ncol=length(x))

  if(is.null(o))
    o <- double(ncol(x))

  D <- ncol(x)

  z <- x-t(matrix( rep(o,nrow(x) ),nrow=D) )
  z <- t_osz(z)
  z1 <- as.vector (ill_lambda(10,ncol(z)))
  z <- t(t_asy(z,0.2))*z1
  z <- t(z)
  res <- f_ackley(z)

  return(res)
}

#' CEC 2013 f4 function. The permutation can be modified by setting the seed. The problem is fixed to have 1000D.
#'
#' @title Rotated elliptic function with 1000D variables.
#' @param x design variable
#' @param o the optimum
#' @param seed Seed for random permutation and rotation matrix
#' @examples
#' x <- runif(100)
#' o <- double(100)
#' f4cec(x,o)
#' @export
f4cec<- function(x,o=NULL,rotation_matrix,permutation,weight){ # rotated elliptic
  x <- matrix(x,ncol=1000)
  D <- ncol(x)
  if(is.null(o))
    o <- double(ncol(x))

  res <- NULL
  group_size <- c(50,25,25,100,50,25,25,700)
  group <- NULL
  groupEnd <- 0
  for(i in 1:8){
    group <- append(group,list(groupEnd+(1:group_size[i])))
    #  permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
    groupEnd <- groupEnd + group_size[i]
  }
  nPoint <- nrow(x)
  y <- x-t(matrix( rep(o,nrow(x) ),nrow=D) )
  y <- y[,permutation,drop=F]
  y <- t(rotation_matrix %*% t(y))

  y <- t_osz(y)

  group_sum <-numeric(nrow(x))
  for(i in 1:7){
    group_sum <- group_sum + weight[i]*f_elliptic(y[,group[[i]]])
  }
  group_sum <- group_sum + f_elliptic(y[,group[[8]]])
  res <-  (group_sum)

  return(res)
}



#' CEC 2013 f5 function.
#'
#' @title Rotated rastrigin function with 1000D variables.
#' @param x design variable
#' @param o the optimum
#' @param seed Seed for random permutation and rotation matrix
#' @examples
#' x <- runif(100)
#' o <- double(100)
#' f5cec(x,o)
#' @export
f5cec<- function(x,o=NULL, rotation_matrix,permutation,weight){ # rotated elliptic
  x <- matrix(x,ncol=1000)
  D <- ncol(x)
  if(is.null(o))
    o <- double(ncol(x))

  res <- NULL
  group_size <- c(25,10,10,55)
  group_size <- c(50,25,25,100,50,25,25,700)
  group <- NULL
  groupEnd <- 0
  for(i in 1:8){
    group <- append(group,list(groupEnd+(1:group_size[i])))
    #  permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
    groupEnd <- groupEnd + group_size[i]
  }

  y <- x-t(matrix( rep(o,nrow(x) ),nrow=D) )
  y <- y[,permutation,drop=F]
  y <- t(rotation_matrix %*% t(y))
  y <- t_osz(y)
  z1 <- as.vector(ill_lambda(10,ncol(y)))
  y <- t(t_asy(y,0.2))*z1
  y <- t(y)

  group_sum <-numeric(nrow(x))

  for(i in 1:7){
    group_sum <- group_sum + weight[i]*f_rastrigin(y[,group[[i]]])
  }
  group_sum <- group_sum + f_rastrigin(y[,group[[8]]])
  res <-  (group_sum)
  return(res)
}

#' CEC 2013 f6 function.
#'
#' @title Rotated Ackley function with 1000D variables.
#' @param x design variable
#' @param o the optimum
#' @param seed Seed for random permutation and rotation matrix
#' @examples
#' x <- runif(100)
#' o <- double(100)
#' f6cec(x,o)
#' @export
f6cec<- function(x,o=NULL, rotation_matrix,permutation,weight){ # rotated elliptic
  x <- matrix(x,ncol=1000)
  D <- ncol(x)
  if(is.null(o))
    o <- double(ncol(x))

  res <- NULL
  group_size <- c(50,25,25,100,50,25,25,700)
  group <- NULL
  groupEnd <- 0
  for(i in 1:8){
    group <- append(group,list(groupEnd+(1:group_size[i])))
    #  permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
    groupEnd <- groupEnd + group_size[i]
  }

  y <- x-t(matrix( rep(o,nrow(x) ),nrow=D) )
  y <- y[,permutation,drop=F]
  y <- t(rotation_matrix %*% t(y))
  y <- t_osz(y)
  z1 <- as.vector(ill_lambda(10,ncol(y)))
  y <- t(t_asy(y,0.2))*z1
  y <- t(y)

  group_sum <-numeric( nrow(x) )
  for(i in 1:7){
    group_sum <- group_sum + weight[i]*f_ackley(y[,group[[i]]])
  }
  group_sum <- group_sum + f_ackley(y[,group[[8]]])
  res <-  (group_sum)
  return(res)
}


#' CEC 2013 f7 function.
#'
#' @title Rotated Schwefel function with 1000D variables.
#' @param x design variable
#' @param o the optimum
#' @param seed Seed for random permutation and rotation matrix
#' @examples
#' x <- runif(100)
#' o <- double(100)
#' f7cec(x,o)
#' @export
f7cec<- function(x,o=NULL, rotation_matrix,permutation,weight){ # rotated Schwefel
  x <- matrix(x,ncol=1000)
  D <- ncol(x)
  if(is.null(o))
    o <- double(ncol(x))

  res <- NULL
  group_size <- c(50,25,25,100,50,25,25,700)
  group <- NULL
  groupEnd <- 0
  for(i in 1:8){
    group <- append(group,list(groupEnd+(1:group_size[i])))
    #  permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
    groupEnd <- groupEnd + group_size[i]
  }
  nPoint <- nrow(x)
  y <- x-t(matrix( rep(o,nrow(x) ),nrow=D) )
  y <- y[,permutation,drop=F]
  y <- t(rotation_matrix %*% t(y))
  y <- t_osz(y)
  y <- t_asy(y,0.2)
  #y <-t_asy(y,0.2)*ill_lambda(10,ncol(y))

  group_sum <-numeric(nrow(x))
  for(i in 1:7){
    group_sum <- group_sum + weight[i]*f_schewefel_1_2(y[,group[[i]],drop=F])
  }
  group_sum <- group_sum + f_schewefel_1_2(y[,group[[8]],drop=F])
  res <-  (group_sum)
  return(res)
}


#' CEC 2013 f8 function.
#'
#' @title Groups of nonseparable Rotated elliptic function with 1000D variables.
#' @param x design variable
#' @param o the optimum
#' @param seed Seed for random permutation and rotation matrix
#' @examples
#' x <- runif(100)
#' o <- double(100)
#' f8cec(x,o)
#' @export
f8cec<- function(x,o=NULL, rotation_matrix,permutation,weight){ # rotated elliptic
  x <- matrix(x,ncol=1000)
  D <- ncol(x)
  if(is.null(o))
    o <- double(ncol(x))

  res <- NULL
  group_size <- c(50,50,25,25,100,100,25,25,50,25,100,25,100,50,25,25,25,100,50,25)
  nGroup <- length(group_size)
  group <- NULL
  groupEnd <- 0

  for(i in 1:nGroup){
    group <- append(group,list(groupEnd+(1:group_size[i])))
    # permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
    groupEnd <- groupEnd + group_size[i]
  }
  nPoint <- nrow(x)
  y <- x-t(matrix( rep(o,nrow(x) ),nrow=D) )
  y <- y[,permutation,drop=F]
  y <- t(rotation_matrix %*% t(y))
  y <- t_osz(y)
  #y <- t(y)
  group_sum <-numeric(nrow(x))
  for(i in 1:nGroup){
    group_sum <- group_sum + weight[i]*f_elliptic(y[,group[[i]]])
  }

  res <-  (group_sum)
  return(res)
}

#' CEC 2013 f9 function.
#'
#' @title Rotated rastrigin function with 1000D variables. no separable variable
#' @param x design variable
#' @param o the optimum
#' @param seed Seed for random permutation and rotation matrix
#' @examples
#' x <- runif(100)
#' o <- double(100)
#' f9cec(x,o)
#' @export
f9cec<- function(x,o=NULL, rotation_matrix,permutation,weight){ # rotated rastrigin
  x <- matrix(x,ncol=1000)
  D <- ncol(x)
  if(is.null(o))
    o <- double(ncol(x))

  res <- NULL
  group_size <- c(50,50,25,25,100,100,25,25,50,25,100,25,100,50,25,25,25,100,50,25)
  nGroup <- length(group_size)
  group <- NULL
  groupEnd <- 0

  for(i in 1:nGroup){
    group <- append(group,list(groupEnd+(1:group_size[i])))
    # permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
    groupEnd <- groupEnd + group_size[i]
  }


  y <- x-t(matrix( rep(o,nrow(x) ),nrow=D) )
  y <- y[,permutation,drop=F]
  y <- t(rotation_matrix %*% t(y))
  y <- t_osz(y)
  z1 <- as.vector(ill_lambda(10,ncol(y)))
  y <- t(t_asy(y,0.2))*z1
  y <- t(y)

  group_sum <-numeric(nrow(x))
  for(i in 1:nGroup){
    group_sum <- group_sum + weight[i]*f_rastrigin(y[,group[[i]]])
  }
  res <-  (group_sum)
  return(res)
}


#' CEC 2013 f10 function.
#'
#' @title Rotated Ackley function with 1000D variables. no separable variable
#' @param x design variable
#' @param o the optimum
#' @param seed Seed for random permutation and rotation matrix
#' @examples
#' x <- runif(100)
#' o <- double(100)
#' f10cec(x,o)
#' @export
f10cec<- function(x,o=NULL, rotation_matrix,permutation,weight){ # rotated ackley
  x <- matrix(x,ncol=1000)
  D <- ncol(x)
  if(is.null(o))
    o <- double(ncol(x))
  res <- NULL
  group_size <- c(50,50,25,25,100,100,25,25,50,25,100,25,100,50,25,25,25,100,50,25)
  nGroup <- length(group_size)
  group <- NULL
  groupEnd <- 0

  for(i in 1:nGroup){
    group <- append(group,list(groupEnd+(1:group_size[i])))
    # permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
    groupEnd <- groupEnd + group_size[i]
  }


  y <- x-t(matrix( rep(o,nrow(x) ),nrow=D) )
  y <- y[,permutation,drop=F]
  y <- t(rotation_matrix %*% t(y))
  y <- t_osz(y)
  z1 <- as.vector(ill_lambda(10,ncol(y)))
  y <- t(t_asy(y,0.2))*z1
  y <- t(y)

  group_sum <-numeric(nrow(x))
  for(i in 1:nGroup){
    group_sum <- group_sum + weight[i]*f_ackley(y[,group[[i]]])
  }
  res <-  (group_sum)
  return(res)
}

#' CEC 2013 f11 function.
#'
#' @title Rotated Schwefel function with 1000D variables. no separable variable
#' @param x design variable
#' @param o the optimum
#' @param seed Seed for random permutation and rotation matrix
#' @examples
#' x <- runif(100)
#' o <- double(100)
#' f11cec(x,o)
#' @export
f11cec<- function(x,o=NULL, rotation_matrix,permutation,weight){ # rotated elliptic
  x <- matrix(x,ncol=1000)
  D <- ncol(x)
  if(is.null(o))
    o <- double(ncol(x))
  res <- NULL
  group_size <- c(50,50,25,25,100,100,25,25,50,25,100,25,100,50,25,25,25,100,50,25)
  nGroup <- length(group_size)
  group <- NULL
  groupEnd <- 0

  for(i in 1:nGroup){
    group <- append(group,list(groupEnd+(1:group_size[i])))
    # permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
    groupEnd <- groupEnd + group_size[i]
  }


  y <- x-t(matrix( rep(o,nrow(x) ),nrow=D) )
  y <- y[,permutation,drop=F]
  y <- t(rotation_matrix %*% t(y))
  y <- t_osz(y)
  y <- t_asy(y,0.2)
  #y <-t_asy(y,0.2)*ill_lambda(10,ncol(y))

  group_sum <-numeric(nrow(x))
  for(i in 1:nGroup){
    group_sum <- group_sum + weight[i]*f_schewefel_1_2(y[,group[[i]]])
  }
  res <-  (group_sum)
  return(res)
}

#' CEC 2013 f12 function.
#'
#' @title Shifted Rosenbrock function with 100D variables. no separable variable
#' @param x design variable
#' @param o the optimum
#' @param seed Seed for random permutation and rotation matrix
#' @examples
#' x <- runif(100)
#' o <- double(100)
#' f12cec(x,o)
#' @export
f12cec<- function(x,o=NULL){ # rotated Rosenbrock
  x <- matrix(x,ncol=1000)
  D <- ncol(x)
  if(is.null(o))
    o <- double(ncol(x))

  res <- NULL

  y <- x-t(matrix( rep(o,nrow(x) ),nrow=D) )
  #  y <- t(y)

  group_sum <- numeric(nrow(x))
  for(i in 1:(ncol(x)-1)){
    group_sum <- group_sum + 100*((y[,i]^2 - y[,i+1])^2) + (y[,i]-1)^2
  }
  res <-  (group_sum)

  return(res)
}

#' CEC 2013 f13 function.
#'
#' @title Schwefel function with Conforming Overlapping Subcomponents
#' @param x design variable
#' @param o the optimum
#' @param seed Seed for random permutation and rotation matrix
#' @examples
#' x <- runif(94)
#' o <- double(94)
#' f13cec(x,o)
#' @export
f13cec<- function(x,o=NULL, rotation_matrix,permutation,weight){ # rotated elliptic
  x <- matrix(x,ncol=905)
  D <- ncol(x)
  if(is.null(o))
    o <- double(ncol(x))
  res <- NULL
  group_size <- c(50,50,25,25,100,100,25,25,50,25,100,25,100,50,25,25,25,100,50,25)
  nGroup <- length(group_size)
  groupStart <- 1
  groupEnd <- 0
  m <- 5

  group <- NULL

  for(i in 1:nGroup){
    group <- append(group,list(groupEnd+(1:group_size[i])))
    #permutation <- append(permutation,sample.int(group_size[i])+groupStart-1)
    groupStart <- groupStart + group_size[i] - m
    groupEnd <- groupEnd + group_size[i]
  }

  y <- x-t(matrix( rep(o,nrow(x) ),nrow=D) )

  y <- y[,permutation,drop=F]

  y <- t(rotation_matrix %*% t(y))
  y <- t_osz(y)
  y <- t_asy(y,0.2)
  #y <-t_asy(y,0.2)*ill_lambda(10,ncol(y))
  group_sum <-numeric(nrow(x))
  for(i in 1:nGroup){
    group_sum <- group_sum + weight[i]*f_schewefel_1_2(y[,group[[i]]])
  }
  res <-  (group_sum)
  return(res)
}

#' CEC 2013 f14 function.
#'
#' @title Schwefel function with Conflicting Overlapping Subcomponents
#' @param x design variable
#' @param o the optimum
#' @param seed Seed for random permutation and rotation matrix
#' @examples
#' x <- runif(94)
#' o <- double(94)
#' f14cec(x,o)
#' @export
f14cec<- function(x,o=NULL, rotation_matrix,permutation,weight){
  x <- matrix(x,ncol=905)
  D <- ncol(x)
  if(is.null(o))
    o <- double(1000)
  res <- NULL
  group_size <- c(50,50,25,25,100,100,25,25,50,25,100,25,100,50,25,25,25,100,50,25)
  nGroup <- length(group_size)
  groupStart <- 1
  groupEnd <- 0
  m <- 5

  group <- NULL

  for(i in 1:nGroup){
    group <- append(group,list(groupEnd+(1:group_size[i])))
    #permutation <- append(permutation,sample.int(group_size[i])+groupStart-1)
    groupStart <- groupStart + group_size[i] - m
    groupEnd <- groupEnd + group_size[i]
  }

  y <- x[,permutation,drop=F]-t(matrix( rep(o,nrow(x) ),nrow=1000) )
  #y <- y[,permutation,drop=F]
  y <- t(rotation_matrix %*% t(y))
  y <- t_osz(y)
  y <- t_asy(y,0.2)

  #y <-t_asy(y,0.2)*ill_lambda(10,ncol(y))

  group_sum <-numeric(nrow(x))
  for(i in 1:nGroup){
    group_sum <- group_sum + weight[i]*f_schewefel_1_2(y[,group[[i]]])
  }
  res <-  (group_sum)
  return(res)
}

#' CEC 2013 f15 function.
#'
#' @title shifted Schwefel function
#' @param x design variable
#' @param o the optimum
#' @param seed Seed for random permutation and rotation matrix
#' @examples
#' x <- runif(100)
#' o <- double(100)
#' f15cec(x,o)
#' @export
f15cec<- function(x,o=NULL){ # rotated elliptic
  x <- matrix(x,ncol=1000)
  D <- ncol(x)
  if(is.null(o))
    o <- double(ncol(x))
  res <- NULL

  set.seed(seed)

  y <- x-t(matrix( rep(o,nrow(x) ),nrow=D) )
  y <- t_osz(y)
  y <- t_asy(y,0.2)

  #y <-t_asy(y,0.2)*ill_lambda(10,ncol(y))

  return <- f_schewefel_1_2(y)
}
