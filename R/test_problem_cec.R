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
f1cec <- function(x,o){ # shifted elliptic
  if(is.vector(x))
    x <- matrix(x,ncol=length(x))

  res <- NULL
  res <- foreach(rowIndex=1:nrow(x),.combine = 'c')%dopar%{
    z <- x[rowIndex,]-o
    z <- t_osz(z)
    f_elliptic(z)
  }
  return(res)
}

#' CEC 2013 f2 function.
#'
#' @title Shifted Rastrigin's function with m>0 Dimensions variables. The CEC2013 test problem uses m=1000
#' @param x design variable
#' @param o the optimum
#' @examples
#' x <- runif(100)
#' o <- double(100)
#' f2cec(x,o)
#' @export
f2cec <- function(x,o){ # shifted rastrigin
  if(is.vector(x))
    x <- matrix(x,ncol=length(x))

  res <- NULL
  res <- foreach(rowIndex = 1:nrow(x),.combine = 'c')%dopar%{
    z <- x[rowIndex,]-o
    z <- t_osz(z)
    z <- ill_lambda(10,length(z))*t_asy(z,0.2)
    cmaes::f_rastrigin(z)
  }
  return(res)
}

#' CEC 2010 f3 function.
#'
#' @title Shifted Ackley's function with m>0 variables. The CEC2013 test problem uses m=1000.
#' @param x design variable
#' @param o the optimum
#' @examples
#' x <- runif(100)
#' o <- double(100)
#' f3cec(x,o)
#' @export
f3cec <- function(x,o){ # shifted ackley
  if(is.vector(x))
    x <- matrix(x,ncol=length(x))
  res <- NULL
  res <- foreach(rowIndex = 1:nrow(x),.combine = 'c')%dopar%{
    z <- x[rowIndex,]-o
    z <- t_osz(z)
    z <- ill_lambda(10,length(z))*t_asy(z,0.2)
    f_ackley(z)
  }
  return(res)
}

#' CEC 2013 f4 function. The permutation can be modified by setting the seed. The problem is fixed to have 1000D.
#'
#' @title Rotated elliptic function with 1000D variables.
#' @param x design variable
#' @param o the optimum
#' @param seed Seed for random permutation and rotation matrix
#' @examples
#' x <- runif(1000)
#' o <- double(1000)
#' f4cec(x,o)
#' @export
f4cec <- function(x,o, seed=2000){ # rotated elliptic
  x <- matrix(x,ncol=1000)
  res <- NULL
  group_size <- c(50,25,25,100,50,25,25,700)
  rotation_matrix <- matrix(numeric(1000*1000),nrow=1000)
  set.seed(seed)
  weight <- cec_weight(7)
  group <- NULL
  groupEnd <- 0
  permutation <- NULL

  for(i in 1:8){
    group <- append(group,list(groupEnd+(1:group_size[i])))
    permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
    groupEnd <- groupEnd + group_size[i]
  }
  for(i in 1:7){
    rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
  }
  rotation_matrix[group[[8]],group[[8]]] <- diag(1,group_size[8])
  nPoint <- nrow(x)
  res <- foreach(rowIndex = 1:nrow(x),.combine = 'c')%dopar%{
    y <- x[rowIndex,]-o
    y <- y[permutation]
    y <- rotation_matrix %*% y
    y <- t_osz(y)

    group_sum <-0
    for(i in 1:7){
      group_sum <- group_sum + weight[i]*f_elliptic(y[group[[i]]])
    }
    group_sum <- group_sum + f_elliptic(y[group[[8]]])
    group_sum
  }
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
f5cec <- function(x,o, seed=1000){ # rotated elliptic
  x <- matrix(x,ncol=1000)
  res <- NULL
  group_size <- c(50,25,25,100,50,25,25,700)
  rotation_matrix <- matrix(numeric(1000*1000),nrow=1000)
  groupEnd <- 0

  set.seed(seed)
  weight <- cec_weight(7)
  group <- NULL
  permutation <- NULL

  for(i in 1:8){
    group <- append(group,list(groupEnd+(1:group_size[i])))
    permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
    groupEnd <- groupEnd + group_size[i]
  }
  for(i in 1:7){
    rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
  }
  rotation_matrix[group[[8]],group[[8]]] <- diag(1,group_size[8])


  res <- foreach(rowIndex = 1:nrow(x),.combine = 'c')%dopar%{
    y <- x[rowIndex,]-o
    y <- y[permutation]
    y <- rotation_matrix %*% y

    y <- t_osz(y)
    y <- t_asy(y,0.2)
    y <- ill_lambda(10,length(y)) %*% y
    group_sum <-0
    for(i in 1:7){
      group_sum <- group_sum + weight[i]*cmaes:::f_rastrigin(y[group[[i]]])
    }
    group_sum <- group_sum + cmaes:::f_rastrigin(y[group[[8]]])
    group_sum
  }
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
f6cec <- function(x,o, seed=1000){ # rotated elliptic
  x <- matrix(x,ncol=1000)
  res <- NULL
  group_size <- c(50,25,25,100,50,25,25,700)
  rotation_matrix <- matrix(numeric(1000*1000),nrow=1000)
  groupEnd <- 0

  set.seed(seed)
  weight <- cec_weight(7)
  #weight[1] <- 10000
  group <- NULL
  permutation <- NULL

  for(i in 1:8){
    group <- append(group,list(groupEnd+(1:group_size[i])))
    permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
    groupEnd <- groupEnd + group_size[i]
  }
  for(i in 1:7){
    rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
  }
  rotation_matrix[group[[8]],group[[8]]] <- diag(1,group_size[8])
  res <- foreach(rowIndex = 1:nrow(x),.combine = 'c')%dopar%{
    y <- x[rowIndex,]-o
    y <- y[permutation]
    y <- rotation_matrix %*% y
    y <- t_osz(y)
    y <- t_asy(y,0.2)

    y <- ill_lambda(10,length(y)) * y

    group_sum <-0

    for(i in 1:7){
      group_sum <- group_sum + weight[i]*f_ackley(y[group[[i]]])
    }
    group_sum <- group_sum + f_ackley(y[group[[8]]])
    group_sum
  }

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
f7cec <- function(x,o, seed=1000){ # rotated elliptic
  x <- matrix(x,ncol=1000)
  res <- NULL
  group_size <- c(50,25,25,100,50,25,25,700)
  rotation_matrix <- matrix(numeric(1000*1000),nrow=1000)
  groupEnd <- 0

  set.seed(seed)
  weight <- cec_weight(7)
  group <- NULL
  permutation <- NULL

  for(i in 1:8){
    group <- append(group,list(groupEnd+(1:group_size[i])))
    permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
    groupEnd <- groupEnd + group_size[i]
  }
  for(i in 1:7){
    rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
  }
  rotation_matrix[group[[8]],group[[8]]] <- diag(1,group_size[8])

  nPoint <- nrow(x)
  res <- foreach(rowIndex = 1:nrow(x),.combine = 'c')%dopar%{
    y <- x[rowIndex,]-o
    y <- y[permutation]
    y <- rotation_matrix %*% y

    y <- t_osz(y)
    y <- t_asy(y,0.2)

    group_sum <-0
    for(i in 1:7){
      group_sum <- group_sum + weight[i]*f_schewefel_1_2(y[group[[i]]])
    }
    group_sum <- group_sum + f_schewefel_1_2(y[group[[8]]])
    group_sum
  }
  return(res)
}
