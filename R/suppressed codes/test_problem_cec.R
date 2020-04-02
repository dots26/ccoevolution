#' #' CEC 2013 f1 function. Can accomodate different variable dimensionality m>0. The CEC2013 test problem uses m=1000.
#' #'
#' #' @title Shifted elliptic function.
#' #' @param x design variable. One row per point.
#' #' @param o the optimum
#' #' @examples
#' #' x <- runif(1000)
#' #' o <- double(1000)
#' #' f1cec(x,o)
#' #' @export
#' f1cec <- function(x,o=NULL){ # shifted elliptic
#'   if(is.vector(x))
#'     x <- matrix(x,ncol=length(x))
#'
#'   if(is.null(o))
#'     o <- double(ncol(x))
#'
#'   res <- NULL
#'   res <- foreach(rowIndex=1:nrow(x),.combine = 'c')%dopar%{
#'     z <- x[rowIndex,]-o
#'     z <- t_osz(z)
#'     f_elliptic(z)
#'   }
#'   return(res)
#' }
#'
#' #' CEC 2013 f2 function.
#' #'
#' #' @title Shifted Rastrigin's function with m>0 Dimensions variables. The CEC2013 test problem uses m=1000
#' #' @param x design variable
#' #' @param o the optimum
#' #' @examples
#' #' x <- runif(1000)
#' #' o <- double(1000)
#' #' f2cec(x,o)
#' #' @export
#' f2cec <- function(x,o=NULL){ # shifted rastrigin
#'   if(is.vector(x))
#'     x <- matrix(x,ncol=length(x))
#'
#'
#'   if(is.null(o))
#'     o <- double(ncol(x))
#'
#'   res <- NULL
#'   res <- foreach(rowIndex = 1:nrow(x),.combine = 'c')%do%{
#'     z <- x[rowIndex,]-o
#'     z <- t_osz(z)
#'     z <- as.vector(ill_lambda(10,length(z)))*t_asy(z,0.2)
#' if(rowIndex==1) print(z)
#'     cmaes::f_rastrigin(z)
#'   }
#'   return(res)
#' }
#'
#' #' CEC 2010 f3 function.
#' #'
#' #' @title Shifted Ackley's function with m>0 variables. The CEC2013 test problem uses m=1000.
#' #' @param x design variable
#' #' @param o the optimum
#' #' @examples
#' #' x <- runif(1000)
#' #' o <- double(1000)
#' #' f3cec(x,o)
#' #' @export
#' f3cec <- function(x,o=NULL){ # shifted ackley
#'   if(is.vector(x))
#'     x <- matrix(x,ncol=length(x))
#'
#'
#'   if(is.null(o))
#'     o <- double(ncol(x))
#'
#'   res <- NULL
#'   res <- foreach(rowIndex = 1:nrow(x),.combine = 'c')%dopar%{
#'     z <- x[rowIndex,]-o
#'     z <- t_osz(z)
#'     z <- as.vector(ill_lambda(10,length(z)))*t_asy(z,0.2)
#'     if(rowIndex==1)print(z)
#'     f_ackley(z)
#'   }
#'   return(res)
#' }
#'
#' #' CEC 2013 f4 function. The permutation can be modified by setting the seed. The problem is fixed to have 1000D.
#' #'
#' #' @title Rotated elliptic function with 1000D variables.
#' #' @param x design variable
#' #' @param o the optimum
#' #' @param seed Seed for random permutation and rotation matrix
#' #' @examples
#' #' x <- runif(1000)
#' #' o <- double(1000)
#' #' f4cec(x,o)
#' #' @export
#' f4cec <- function(x,o=NULL, seed=2000){ # rotated elliptic
#'   x <- matrix(x,ncol=1000)
#'
#'   if(is.null(o))
#'     o <- double(ncol(x))
#'
#'   res <- NULL
#'   group_size <- c(50,25,25,100,50,25,25,700)
#'   rotation_matrix <- matrix(numeric(1000*1000),nrow=1000)
#'   set.seed(seed)
#'   weight <- cec_weight(7)
#'   group <- NULL
#'   groupEnd <- 0
#'   permutation <- NULL
#'
#'   for(i in 1:8){
#'     group <- append(group,list(groupEnd+(1:group_size[i])))
#'     permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
#'     groupEnd <- groupEnd + group_size[i]
#'   }
#'   for(i in 1:7){
#'     rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
#'   }
#'   rotation_matrix[group[[8]],group[[8]]] <- diag(1,group_size[8])
#'   nPoint <- nrow(x)
#'   res <- foreach(rowIndex = 1:nrow(x),.combine = 'c')%do%{
#'     y <- x[rowIndex,]-o
#'     y <- y[permutation]
#'     y <- rotation_matrix %*% y
#'     y <- t_osz(y)
#'
#'     group_sum <-0
#'     for(i in 1:7){
#'       group_sum <- group_sum + weight[i]*f_elliptic(y[group[[i]]])
#'     }
#'     group_sum <- group_sum + f_elliptic(y[group[[8]]])
#'     group_sum
#'   }
#'   return(res)
#' }
#'
#'
#'
#' #' CEC 2013 f5 function.
#' #'
#' #' @title Rotated rastrigin function with 1000D variables.
#' #' @param x design variable
#' #' @param o the optimum
#' #' @param seed Seed for random permutation and rotation matrix
#' #' @examples
#' #' x <- runif(1000)
#' #' o <- double(1000)
#' #' f5cec(x,o)
#' #' @export
#' f5cec <- function(x,o=NULL, seed=1000){ # rotated elliptic
#'   x <- matrix(x,ncol=1000)
#'
#'   if(is.null(o))
#'     o <- double(ncol(x))
#'
#'   res <- NULL
#'   group_size <- c(50,25,25,100,50,25,25,700)
#'   rotation_matrix <- matrix(numeric(1000*1000),nrow=1000)
#'   groupEnd <- 0
#'
#'   set.seed(seed)
#'   weight <- cec_weight(7)
#'   group <- NULL
#'   permutation <- NULL
#'
#'   for(i in 1:8){
#'     group <- append(group,list(groupEnd+(1:group_size[i])))
#'     permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
#'     groupEnd <- groupEnd + group_size[i]
#'   }
#'   for(i in 1:7){
#'     rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
#'   }
#'   rotation_matrix[group[[8]],group[[8]]] <- diag(1,group_size[8])
#'
#'
#'   res <- foreach(rowIndex = 1:nrow(x),.combine = 'c')%dopar%{
#'     y <- x[rowIndex,]-o
#'     y <- y[permutation]
#'     y <- rotation_matrix %*% y
#'
#'     y <- t_osz(y)
#'     y <- t_asy(y,0.2)
#'     y <- as.vector(ill_lambda(10,length(z)))* y
#'     group_sum <-0
#'     for(i in 1:7){
#'       group_sum <- group_sum + weight[i]*cmaes::f_rastrigin(y[group[[i]]])
#'     }
#'     group_sum <- group_sum + cmaes::f_rastrigin(y[group[[8]]])
#'     group_sum
#'   }
#'   return(res)
#' }
#'
#' #' CEC 2013 f6 function.
#' #'
#' #' @title Rotated Ackley function with 1000D variables.
#' #' @param x design variable
#' #' @param o the optimum
#' #' @param seed Seed for random permutation and rotation matrix
#' #' @examples
#' #' x <- runif(1000)
#' #' o <- double(1000)
#' #' f6cec(x,o)
#' #' @export
#' f6cec <- function(x,o=NULL, seed=1000){ # rotated elliptic
#'   x <- matrix(x,ncol=1000)
#'
#'   if(is.null(o))
#'     o <- double(ncol(x))
#'
#'   res <- NULL
#'   group_size <- c(50,25,25,100,50,25,25,700)
#'   rotation_matrix <- matrix(numeric(1000*1000),nrow=1000)
#'   groupEnd <- 0
#'
#'   set.seed(seed)
#'   weight <- cec_weight(7)
#'   #weight[1] <- 10000
#'   group <- NULL
#'   permutation <- NULL
#'
#'   for(i in 1:8){
#'     group <- append(group,list(groupEnd+(1:group_size[i])))
#'     permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
#'     groupEnd <- groupEnd + group_size[i]
#'   }
#'   for(i in 1:7){
#'     rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
#'   }
#'   rotation_matrix[group[[8]],group[[8]]] <- diag(1,group_size[8])
#'   res <- foreach(rowIndex = 1:nrow(x),.combine = 'c')%dopar%{
#'     y <- x[rowIndex,]-o
#'     y <- y[permutation]
#'     y <- rotation_matrix %*% y
#'     y <- t_osz(y)
#'     y <- t_asy(y,0.2)
#'
#'     y <- as.vector(ill_lambda(10,length(z))) * y
#'
#'     group_sum <-0
#'
#'     for(i in 1:7){
#'       group_sum <- group_sum + weight[i]*f_ackley(y[group[[i]]])
#'     }
#'     group_sum <- group_sum + f_ackley(y[group[[8]]])
#'     group_sum
#'   }
#'
#'   return(res)
#' }
#'
#'
#' #' CEC 2013 f7 function.
#' #'
#' #' @title Rotated Schwefel function with 1000D variables.
#' #' @param x design variable
#' #' @param o the optimum
#' #' @param seed Seed for random permutation and rotation matrix
#' #' @examples
#' #' x <- runif(1000)
#' #' o <- double(1000)
#' #' f7cec(x,o)
#' #' @export
#' f7cec <- function(x,o=NULL, seed=1000){ # rotated elliptic
#'   x <- matrix(x,ncol=1000)
#'
#'   if(is.null(o))
#'     o <- double(ncol(x))
#'
#'   res <- NULL
#'   group_size <- c(50,25,25,100,50,25,25,700)
#'   rotation_matrix <- matrix(numeric(1000*1000),nrow=1000)
#'   groupEnd <- 0
#'
#'   set.seed(seed)
#'   weight <- cec_weight(7)
#'   group <- NULL
#'   permutation <- NULL
#'
#'   for(i in 1:8){
#'     group <- append(group,list(groupEnd+(1:group_size[i])))
#'     permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
#'     groupEnd <- groupEnd + group_size[i]
#'   }
#'   for(i in 1:7){
#'     rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
#'   }
#'   rotation_matrix[group[[8]],group[[8]]] <- diag(1,group_size[8])
#'
#'   nPoint <- nrow(x)
#'   res <- foreach(rowIndex = 1:nrow(x),.combine = 'c')%dopar%{
#'     y <- x[rowIndex,]-o
#'     y <- y[permutation]
#'     y <- rotation_matrix %*% y
#'
#'     y <- t_osz(y)
#'     y <- t_asy(y,0.2)
#'
#'     group_sum <-0
#'     for(i in 1:7){
#'       group_sum <- group_sum + weight[i]*f_schewefel_1_2(y[group[[i]]])
#'     }
#'     group_sum <- group_sum + f_schewefel_1_2(y[group[[8]]])
#'     group_sum
#'   }
#'   return(res)
#' }
#'
#'
#' #' CEC 2013 f8 function.
#' #'
#' #' @title Groups of nonseparable Rotated elliptic function with 1000D variables.
#' #' @param x design variable
#' #' @param o the optimum
#' #' @param seed Seed for random permutation and rotation matrix
#' #' @examples
#' #' x <- runif(1000)
#' #' o <- double(1000)
#' #' f8cec(x,o)
#' #' @export
#' f8cec <- function(x,o=NULL, seed=1000){ # rotated elliptic
#'   x <- matrix(x,ncol=1000)
#'
#'   if(is.null(o))
#'     o <- double(ncol(x))
#'
#'   res <- NULL
#'   group_size <- c(50,50,25,25,100,100,25,25,50,25,100,25,100,50,25,25,25,100,50,25)
#'   nGroup <- length(group_size)
#'   rotation_matrix <- matrix(numeric(1000*1000),nrow=1000)
#'   set.seed(seed)
#'   weight <- cec_weight(nGroup)
#'   group <- NULL
#'   groupEnd <- 0
#'   permutation <- NULL
#'
#'   for(i in 1:nGroup){
#'     group <- append(group,list(groupEnd+(1:group_size[i])))
#'     permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
#'     groupEnd <- groupEnd + group_size[i]
#'   }
#'   for(i in 1:nGroup){
#'     rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
#'   }
#'   nPoint <- nrow(x)
#'   res <- foreach(rowIndex = 1:nrow(x),.combine = 'c')%dopar%{
#'     y <- x[rowIndex,]-o
#'     y <- y[permutation]
#'     y <- rotation_matrix %*% y
#'     y <- t_osz(y)
#'
#'     group_sum <-0
#'     for(i in 1:nGroup){
#'       group_sum <- group_sum + weight[i]*f_elliptic(y[group[[i]]])
#'     }
#'     group_sum
#'   }
#'   return(res)
#' }
#'
#' #' CEC 2013 f9 function.
#' #'
#' #' @title Rotated rastrigin function with 1000D variables. no separable variable
#' #' @param x design variable
#' #' @param o the optimum
#' #' @param seed Seed for random permutation and rotation matrix
#' #' @examples
#' #' x <- runif(1000)
#' #' o <- double(1000)
#' #' f9cec(x,o)
#' #' @export
#' f9cec <- function(x,o=NULL, seed=1000){ # rotated elliptic
#'   x <- matrix(x,ncol=1000)
#'   if(is.null(o))
#'     o <- double(ncol(x))
#'
#'   res <- NULL
#'   group_size <- c(50,50,25,25,100,100,25,25,50,25,100,25,100,50,25,25,25,100,50,25)
#'   nGroup <- length(group_size)
#'   rotation_matrix <- matrix(numeric(1000*1000),nrow=1000)
#'   groupEnd <- 0
#'
#'   set.seed(seed)
#'   weight <- cec_weight(nGroup)
#'   group <- NULL
#'   permutation <- NULL
#'
#'   for(i in 1:nGroup){
#'     group <- append(group,list(groupEnd+(1:group_size[i])))
#'     permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
#'     groupEnd <- groupEnd + group_size[i]
#'   }
#'   for(i in 1:nGroup){
#'     rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
#'   }
#'
#'
#'   res <- foreach(rowIndex = 1:nrow(x),.combine = 'c')%dopar%{
#'     y <- x[rowIndex,]-o
#'     y <- y[permutation]
#'     y <- rotation_matrix %*% y
#'
#'     y <- t_osz(y)
#'     y <- t_asy(y,0.2)
#'     y <- as.vector(ill_lambda(10,length(z)))* y
#'     group_sum <-0
#'     for(i in 1:nGroup){
#'       group_sum <- group_sum + weight[i]*cmaes::f_rastrigin(y[group[[i]]])
#'     }
#'     group_sum
#'   }
#'   return(res)
#' }
#'
#'
#' #' CEC 2013 f10 function.
#' #'
#' #' @title Rotated Ackley function with 1000D variables. no separable variable
#' #' @param x design variable
#' #' @param o the optimum
#' #' @param seed Seed for random permutation and rotation matrix
#' #' @examples
#' #' x <- runif(1000)
#' #' o <- double(1000)
#' #' f10cec(x,o)
#' #' @export
#' f10cec <- function(x,o=NULL, seed=1000){ # rotated elliptic
#'   x <- matrix(x,ncol=1000)
#'   if(is.null(o))
#'     o <- double(ncol(x))
#'   res <- NULL
#'   group_size <- c(50,50,25,25,100,100,25,25,50,25,100,25,100,50,25,25,25,100,50,25)
#'   nGroup <- length(group_size)
#'   rotation_matrix <- matrix(numeric(1000*1000),nrow=1000)
#'   groupEnd <- 0
#'
#'   set.seed(seed)
#'   weight <- cec_weight(nGroup)
#'   group <- NULL
#'   permutation <- NULL
#'
#'   for(i in 1:nGroup){
#'     group <- append(group,list(groupEnd+(1:group_size[i])))
#'     permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
#'     groupEnd <- groupEnd + group_size[i]
#'   }
#'   for(i in 1:nGroup){
#'     rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
#'   }
#'
#'
#'   res <- foreach(rowIndex = 1:nrow(x),.combine = 'c')%dopar%{
#'     y <- x[rowIndex,]-o
#'     y <- y[permutation]
#'     y <- rotation_matrix %*% y
#'
#'     y <- t_osz(y)
#'     y <- t_asy(y,0.2)
#'     y <- as.vector(ill_lambda(10,length(z)))*y
#'     group_sum <-0
#'     for(i in 1:nGroup){
#'       group_sum <- group_sum + weight[i]*f_ackley(y[group[[i]]])
#'     }
#'     group_sum
#'   }
#'   return(res)
#' }
#'
#' #' CEC 2013 f11 function.
#' #'
#' #' @title Rotated Schwefel function with 1000D variables. no separable variable
#' #' @param x design variable
#' #' @param o the optimum
#' #' @param seed Seed for random permutation and rotation matrix
#' #' @examples
#' #' x <- runif(1000)
#' #' o <- double(1000)
#' #' f11cec(x,o)
#' #' @export
#' f11cec <- function(x,o=NULL, seed=1000){ # rotated elliptic
#'   x <- matrix(x,ncol=1000)
#'   if(is.null(o))
#'     o <- double(ncol(x))
#'   res <- NULL
#'   group_size <- c(50,50,25,25,100,100,25,25,50,25,100,25,100,50,25,25,25,100,50,25)
#'   nGroup <- length(group_size)
#'   rotation_matrix <- matrix(numeric(1000*1000),nrow=1000)
#'   groupEnd <- 0
#'
#'   set.seed(seed)
#'   weight <- cec_weight(nGroup)
#'   group <- NULL
#'   permutation <- NULL
#'
#'   for(i in 1:nGroup){
#'     group <- append(group,list(groupEnd+(1:group_size[i])))
#'     permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
#'     groupEnd <- groupEnd + group_size[i]
#'   }
#'   for(i in 1:nGroup){
#'     rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
#'   }
#'
#'
#'   res <- foreach(rowIndex = 1:nrow(x),.combine = 'c')%dopar%{
#'     y <- x[rowIndex,]-o
#'     y <- y[permutation]
#'     y <- rotation_matrix %*% y
#'
#'     y <- t_osz(y)
#'     y <- t_asy(y,0.2)
#'    # y <- as.vector(ill_lambda(10,length(z)))* y
#'     group_sum <-0
#'     for(i in 1:nGroup){
#'       group_sum <- group_sum + weight[i]*f_schewefel_1_2(y[group[[i]]])
#'     }
#'     group_sum
#'   }
#'   return(res)
#' }
#'
#' #' CEC 2013 f12 function.
#' #'
#' #' @title Shifted Rosenbrock function with 1000D variables. no separable variable
#' #' @param x design variable
#' #' @param o the optimum
#' #' @param seed Seed for random permutation and rotation matrix
#' #' @examples
#' #' x <- runif(1000)
#' #' o <- double(1000)
#' #' f12cec(x,o)
#' #' @export
#' f12cec <- function(x,o=NULL, seed=1000){ # rotated elliptic
#'   x <- matrix(x,ncol=1000)
#'   if(is.null(o))
#'     o <- double(ncol(x))
#'   res <- NULL
#'
#'   res <- foreach(rowIndex = 1:nrow(x),.combine = 'c')%dopar%{
#'     y <- x[rowIndex,]-o
#'
#'     group_sum <-0
#'     for(i in 1:ncol(x)){
#'       group_sum <- group_sum + 100*((y[i]^2 - y[i+1])^2) + (y[i]-1)^2
#'     }
#'     group_sum
#'   }
#'   return(res)
#' }
#'
#' #' CEC 2013 f13 function.
#' #'
#' #' @title Schwefel function with Conforming Overlapping Subcomponents
#' #' @param x design variable
#' #' @param o the optimum
#' #' @param seed Seed for random permutation and rotation matrix
#' #' @examples
#' #' x <- runif(905)
#' #' o <- double(905)
#' #' f13cec(x,o)
#' #' @export
#' f13cec <- function(x,o=NULL, seed=1000){ # rotated elliptic
#'   x <- matrix(x,ncol=905)
#'   if(is.null(o))
#'     o <- double(ncol(x))
#'   res <- NULL
#'   group_size <- c(50,50,25,25,100,100,25,25,50,25,100,25,100,50,25,25,25,100,50,25)
#'   nGroup <- length(group_size)
#'   rotation_matrix <- matrix(numeric(1000*1000),nrow=1000)
#'
#'   groupStart <- 1
#'   groupEnd <- 0
#'   m <- 5
#'
#'   set.seed(seed)
#'   weight <- cec_weight(nGroup)
#'   group <- NULL
#'   permutation <- NULL
#'
#'   for(i in 1:nGroup){
#'     group <- append(group,list(groupEnd+(1:group_size[i])))
#'     permutation <- append(permutation,sample.int(group_size[i])+groupStart-1)
#'     groupStart <- groupStart + group_size[i] - m
#'     groupEnd <- groupEnd + group_size[i]
#'   }
#'   for(i in 1:nGroup){
#'     rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
#'   }
#'
#'   res <- foreach(rowIndex = 1:nrow(x),.combine = 'c')%do%{
#'     y <- x[rowIndex,]-o
#'     y <- y[permutation]
#'     y <- rotation_matrix %*% y
#'     #print(y)
#'     y <- t_osz(y)
#'     y <- t_asy(y,0.2)
#'
#'     group_sum <-0
#'     for(i in 1:nGroup){
#'       group_sum <- group_sum + weight[i]*f_schewefel_1_2(y[group[[i]]])
#'     }
#'     group_sum
#'   }
#'   return(res)
#' }
#'
#' #' CEC 2013 f14 function.
#' #'
#' #' @title Schwefel function with Conflicting Overlapping Subcomponents
#' #' @param x design variable
#' #' @param o the optimum
#' #' @param seed Seed for random permutation and rotation matrix
#' #' @examples
#' #' x <- runif(905)
#' #' o <- double(905)
#' #' f14cec(x,o)
#' #' @export
#' f14cec <- function(x,o=NULL, seed=1000){
#'   x <- matrix(x,ncol=905)
#'   if(is.null(o))
#'     o <- double(ncol(x))
#'   res <- NULL
#'   group_size <- c(50,50,25,25,100,100,25,25,50,25,100,25,100,50,25,25,25,100,50,25)
#'   nGroup <- length(group_size)
#'   rotation_matrix <- matrix(numeric(1000*1000),nrow=1000)
#'   groupStart <- 1
#'   groupEnd <- 0
#'   m <- 5
#'
#'   set.seed(seed)
#'   weight <- cec_weight(nGroup)
#'   group <- NULL
#'   permutation <- NULL
#'
#'   for(i in 1:nGroup){
#'     group <- append(group,list(groupEnd+(1:group_size[i])))
#'     permutation <- append(permutation,sample.int(group_size[i])+groupStart-1)
#'     groupStart <- groupStart + group_size[i] - m
#'     groupEnd <- groupEnd + group_size[i]
#'   }
#'
#'   for(i in 1:nGroup){
#'     rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
#'   }
#'
#'   res <- foreach(rowIndex = 1:nrow(x),.combine = 'c')%dopar%{
#'     y <- x[rowIndex,permutation]-o[permutation]
#'     y <- rotation_matrix %*% y
#'     y <- t_osz(y)
#'     y <- t_asy(y,0.2)
#'
#'     group_sum <-0
#'     for(i in 1:nGroup){
#'       group_sum <- group_sum + weight[i]*f_schewefel_1_2(y[group[[i]]])
#'     }
#'     group_sum
#'   }
#'   return(res)
#' }
#'
#' #' CEC 2013 f15 function.
#' #'
#' #' @title shifted Schwefel function
#' #' @param x design variable
#' #' @param o the optimum
#' #' @param seed Seed for random permutation and rotation matrix
#' #' @examples
#' #' x <- runif(1000)
#' #' o <- double(1000)
#' #' f15cec(x,o)
#' #' @export
#' f15cec <- function(x,o=NULL, seed=1000){ # rotated elliptic
#'   x <- matrix(x,ncol=1000)
#'   if(is.null(o))
#'     o <- double(ncol(x))
#'   res <- NULL
#'
#'   set.seed(seed)
#'
#'   res <- foreach(rowIndex = 1:nrow(x),.combine = 'c')%dopar%{
#'     y <- x[rowIndex,]-o
#'     y <- t_osz(y)
#'     y <- t_asy(y,0.2)
#'
#'     f_schewefel_1_2(y)
#'   }
#'   return(res)
#' }
