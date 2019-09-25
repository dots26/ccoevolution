#' CEC 2013 f4 function. The permutation can be modified by setting the seed. The problem is fixed to have 1000D.
#'
#' @title Rotated elliptic function with 1000D variables.
#' @param x design variable
#' @param o the optimum
#' @param seed Seed for random permutation and rotation matrix
#' @examples
#' x <- runif(100)
#' o <- double(100)
#' f4cec_s(x,o)
#' @export
f4cec_s <- function(x,o=NULL,rotation_matrix,permutation,weight){ # rotated elliptic
  x <- matrix(x,ncol=100)
  D <- ncol(x)
  if(is.null(o))
    o <- double(ncol(x))

  res <- NULL
  # group_size <- c(25,10,10,55)
  # rotation_matrix <- matrix(numeric(100*100),nrow=100)
  # set.seed(seed)
  # weight <- cec_weight(3)
  # weight <- append(weight,1)
  # group <- NULL
  # groupEnd <- 0
  # permutation <- NULL
  # for(i in 1:4){
  #   group <- append(group,list(groupEnd+(1:group_size[i])))
  #   permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
  #   groupEnd <- groupEnd + group_size[i]
  # }
  #
  # for(i in 1:3){
  #   rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
  # }
  # rotation_matrix[group[[4]],group[[4]]] <- diag(1,group_size[4])
  nPoint <- nrow(x)
  y <- x-t(matrix( rep(o,nrow(x) ),nrow=D) )
  y <- y[,permutation,drop=F]
  y <- t(rotation_matrix %*% t(y))

  y <- t_osz(y)

  group_sum <-numeric(nrow(x))
  for(i in 1:3){
    group_sum <- group_sum + weight[i]*f_elliptic(y[,group[[i]]])
  }
  group_sum <- group_sum + f_elliptic(y[,group[[4]]])
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
#' f5cec_s(x,o)
#' @export
f5cec_s <- function(x,o=NULL, rotation_matrix,permutation,weight){ # rotated elliptic
  x <- matrix(x,ncol=100)
  D <- ncol(x)
  if(is.null(o))
    o <- double(ncol(x))

  res <- NULL
  group_size <- c(25,10,10,55)
  # rotation_matrix <- matrix(numeric(100*100),nrow=100)
  # groupEnd <- 0
  #
  # set.seed(seed)
  # weight <- cec_weight(3)
  # group <- NULL
  # permutation <- NULL
  #
  # for(i in 1:4){
  #   group <- append(group,list(groupEnd+(1:group_size[i])))
  #   permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
  #   groupEnd <- groupEnd + group_size[i]
  # }
  # for(i in 1:3){
  #   rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
  # }
  # rotation_matrix[group[[4]],group[[4]]] <- diag(1,group_size[4])

  y <- x-t(matrix( rep(o,nrow(x) ),nrow=D) )
  y <- y[,permutation,drop=F]
  y <- t(rotation_matrix %*% t(y))
  y <- t_osz(y)
  z1 <- as.vector(ill_lambda(10,ncol(y)))
  y <- t(t_asy(y,0.2))*z1
  y <- t(y)

  group_sum <-numeric(nrow(x))

  for(i in 1:3){
    group_sum <- group_sum + weight[i]*f_rastrigin(y[,group[[i]]])
  }
  group_sum <- group_sum + f_rastrigin(y[,group[[4]]])
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
#' f6cec_s(x,o)
#' @export
f6cec_s <- function(x,o=NULL, rotation_matrix,permutation,weight){ # rotated elliptic
  x <- matrix(x,ncol=100)
  D <- ncol(x)
  if(is.null(o))
    o <- double(ncol(x))

  res <- NULL
  group_size <- c(25,10,10,55)
  # rotation_matrix <- matrix(numeric(100*100),nrow=100)
  # groupEnd <- 0
  #
  # set.seed(seed)
  # weight <- cec_weight(3)
  # #weight[1] <- 10000
  # group <- NULL
  # permutation <- NULL
  #
  # for(i in 1:4){
  #   group <- append(group,list(groupEnd+(1:group_size[i])))
  #   permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
  #   groupEnd <- groupEnd + group_size[i]
  # }
  # for(i in 1:3){
  #   rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
  # }
  # rotation_matrix[group[[4]],group[[4]]] <- diag(1,group_size[4])

  y <- x-t(matrix( rep(o,nrow(x) ),nrow=D) )
  y <- y[,permutation,drop=F]
  y <- t(rotation_matrix %*% t(y))
  y <- t_osz(y)
  z1 <- as.vector(ill_lambda(10,ncol(y)))
  y <- t(t_asy(y,0.2))*z1
  y <- t(y)

  group_sum <-numeric( nrow(x) )
  for(i in 1:3){
    group_sum <- group_sum + weight[i]*f_ackley(y[,group[[i]]])
  }
  group_sum <- group_sum + f_ackley(y[,group[[4]]])
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
f7cec_s <- function(x,o=NULL, rotation_matrix,permutation,weight){ # rotated Schwefel
  x <- matrix(x,ncol=100)
  D <- ncol(x)
  if(is.null(o))
    o <- double(ncol(x))

  res <- NULL
  group_size <- c(25,10,10,55)
  # rotation_matrix <- matrix(numeric(100*1000),nrow=100)
  # groupEnd <- 0
  #
  # set.seed(seed)
  # weight <- cec_weight(3)
  # group <- NULL
  # permutation <- NULL
  #
  # for(i in 1:4){
  #   group <- append(group,list(groupEnd+(1:group_size[i])))
  #   permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
  #   groupEnd <- groupEnd + group_size[i]
  # }
  # for(i in 1:3){
  #   rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
  # }
  # rotation_matrix[group[[4]],group[[4]]] <- diag(1,group_size[4])

  nPoint <- nrow(x)
  y <- x-t(matrix( rep(o,nrow(x) ),nrow=D) )
  y <- y[,permutation,drop=F]
  y <- t(rotation_matrix %*% t(y))
  y <- t_osz(y)
  y <- t_asy(y,0.2)
  #y <-t_asy(y,0.2)*ill_lambda(10,ncol(y))

  print(y)
  group_sum <-numeric(nrow(x))
  for(i in 1:3){
    group_sum <- group_sum + weight[i]*f_schewefel_1_2(y[,group[[i]],drop=F])
  }
  group_sum <- group_sum + f_schewefel_1_2(y[,group[[4]],drop=F])
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
#' f8cec_s(x,o)
#' @export
f8cec_s <- function(x,o=NULL, rotation_matrix,permutation,weight){ # rotated elliptic
  x <- matrix(x,ncol=100)
  D <- ncol(x)
  if(is.null(o))
    o <- double(ncol(x))

  res <- NULL
  group_size <- c(25,10,10,55)
  nGroup <- length(group_size)
  # rotation_matrix <- matrix(numeric(100*100),nrow=100)
  # set.seed(seed)
  # weight <- cec_weight(nGroup)
  # group <- NULL
  # groupEnd <- 0
  # permutation <- NULL
  #
  # for(i in 1:nGroup){
  #   group <- append(group,list(groupEnd+(1:group_size[i])))
  #   permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
  #   groupEnd <- groupEnd + group_size[i]
  # }
  # for(i in 1:nGroup){
  #   rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
  # }
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
#' f9cec_s(x,o)
#' @export
f9cec_s <- function(x,o=NULL, rotation_matrix,permutation,weight){ # rotated rastrigin
  x <- matrix(x,ncol=100)
  D <- ncol(x)
  if(is.null(o))
    o <- double(ncol(x))

  res <- NULL
  group_size <- c(25,10,10,55)
  nGroup <- length(group_size)
  # rotation_matrix <- matrix(numeric(100*100),nrow=100)
  # groupEnd <- 0
  #
  # set.seed(seed)
  # weight <- cec_weight(nGroup)
  # group <- NULL
  # permutation <- NULL
  #
  # for(i in 1:nGroup){
  #   group <- append(group,list(groupEnd+(1:group_size[i])))
  #   permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
  #   groupEnd <- groupEnd + group_size[i]
  # }
  # for(i in 1:nGroup){
  #   rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
  # }


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
#' x <- runif(1000)
#' o <- double(1000)
#' f10cec_s(x,o)
#' @export
f10cec_s <- function(x,o=NULL, rotation_matrix,permutation,weight){ # rotated ackley
  x <- matrix(x,ncol=100)
  D <- ncol(x)
  if(is.null(o))
    o <- double(ncol(x))
  res <- NULL
  group_size <- c(25,10,10,55)
  nGroup <- length(group_size)
  # rotation_matrix <- matrix(numeric(100*100),nrow=100)
  # groupEnd <- 0
  #
  # set.seed(seed)
  # weight <- cec_weight(nGroup)
  # group <- NULL
  # permutation <- NULL
  #
  # for(i in 1:nGroup){
  #   group <- append(group,list(groupEnd+(1:group_size[i])))
  #   permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
  #   groupEnd <- groupEnd + group_size[i]
  # }
  # for(i in 1:nGroup){
  #   rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
  # }


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
#' f11cec_s(x,o)
#' @export
f11cec_s <- function(x,o=NULL, rotation_matrix,permutation,weight){ # rotated elliptic
  x <- matrix(x,ncol=100)
  D <- ncol(x)
  if(is.null(o))
    o <- double(ncol(x))
  res <- NULL
  group_size <- c(25,10,10,55)
  nGroup <- length(group_size)
  # rotation_matrix <- matrix(numeric(100*100),nrow=100)
  # groupEnd <- 0
  #
  # set.seed(seed)
  # weight <- cec_weight(nGroup)
  # group <- NULL
  # permutation <- NULL
  #
  # for(i in 1:nGroup){
  #   group <- append(group,list(groupEnd+(1:group_size[i])))
  #   permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
  #   groupEnd <- groupEnd + group_size[i]
  # }
  # for(i in 1:nGroup){
  #   rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
  # }


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
#' @title Shifted Rosenbrock function with 1000D variables. no separable variable
#' @param x design variable
#' @param o the optimum
#' @param seed Seed for random permutation and rotation matrix
#' @examples
#' x <- runif(100)
#' o <- double(100)
#' f12cec_s(x,o)
#' @export
f12cec_s <- function(x,o=NULL){ # rotated Rosenbrock
  x <- matrix(x,ncol=100)
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
#' f13cec_s(x,o)
#' @export
f13cec_s <- function(x,o=NULL, rotation_matrix,permutation,weight){ # rotated elliptic
  x <- matrix(x,ncol=94)
  D <- ncol(x)
  if(is.null(o))
    o <- double(ncol(x))
  res <- NULL
  group_size <- c(25,10,10,55)
  nGroup <- length(group_size)

  # rotation_matrix <- matrix(numeric(100*100),nrow=100)
  #
  # groupStart <- 1
  # groupEnd <- 0
  # m <- 2
  #
  # set.seed(seed)
  # weight <- cec_weight(nGroup)
  # group <- NULL
  # permutation <- NULL
  #
  # for(i in 1:nGroup){
  #   group <- append(group,list(groupEnd+(1:group_size[i])))
  #   permutation <- append(permutation,sample.int(group_size[i])+groupStart-1)
  #   groupStart <- groupStart + group_size[i] - m
  #   groupEnd <- groupEnd + group_size[i]
  # }
  # for(i in 1:nGroup){
  #   rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
  # }

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
#' f14cec_s(x,o)
#' @export
f14cec_s <- function(x,o=NULL, rotation_matrix,permutation,weight){
  x <- matrix(x,ncol=94)
  D <- ncol(x)
  if(is.null(o))
    o <- double(ncol(x))
  res <- NULL
  group_size <- c(25,10,10,55)
  nGroup <- length(group_size)
  # rotation_matrix <- matrix(numeric(100*100),nrow=100)
  # groupStart <- 1
  # groupEnd <- 0
  # m <- 5
  #
  # set.seed(seed)
  # weight <- cec_weight(nGroup)
  # group <- NULL
  # permutation <- NULL
  #
  # for(i in 1:nGroup){
  #   group <- append(group,list(groupEnd+(1:group_size[i])))
  #   permutation <- append(permutation,sample.int(group_size[i])+groupStart-1)
  #   groupStart <- groupStart + group_size[i] - m
  #   groupEnd <- groupEnd + group_size[i]
  # }
  #
  # for(i in 1:nGroup){
  #   rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
  # }

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

#' CEC 2013 f15 function.
#'
#' @title shifted Schwefel function
#' @param x design variable
#' @param o the optimum
#' @param seed Seed for random permutation and rotation matrix
#' @examples
#' x <- runif(100)
#' o <- double(100)
#' f15cec_s(x,o)
#' @export
f15cec_s <- function(x,o=NULL){ # rotated elliptic
  x <- matrix(x,ncol=100)
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
