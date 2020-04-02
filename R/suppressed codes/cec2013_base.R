# f_elliptic <- function(x){
#   D <- length(x)
#   i <- 1:D
#   conditioning <- 10^(6*(i-1)/(D-1))
#   return(sum(conditioning*x*x))
# }
#
# f_ackley <- function(z){
#   D <- length(z)
#   ackley <- -20 * exp(-0.2 * (1/D*sum(z*z))^0.5 ) - exp(1/D*sum(cos(2*pi*z)) ) + 20 + exp(1)
# }
#
# f_schewefel_1_2 <- function(x){
#   D <- length(x)
#
#   summand <- 0
#   for(i in 1:D){
#     summand <- summand + ( sum(x[1:i])*sum(x[1:i]))
#   }
#   return(summand)
# }
#
# t_osz <- function(x){
#   xlog <- log(abs(x))* as.integer(x!=0)
#
#   xlog[which(is.nan(xlog))] <- 0
#
#   c1 <- 10*as.integer(x>0) +  5.5*as.integer(x<=0)
#   c2 <- 7.9*as.integer(x>0) +  3.1*as.integer(x<=0)
#
#   t <- sign(x) * exp( xlog + 0.0049*sin(c1*xlog) + sin(c2*xlog ))
# }
#
# t_asy <- function(x,beta){
#   D <- length(x)
#
#   i <- 1:D
#   beta_term <- beta * ((i-1)/(D-1))
#   t <- NULL
#   for(iter in i){
#     if(x[iter]>0){
#       t <- append(t,x[iter]^(1+beta_term[iter]*(x[iter]^0.5)))
#     }else{
#       t <- append(t,x[iter])
#     }
#   }
#   return(t)
# }
#
# cec_weight <- function(D){
#   power_term <- 1*rnorm(D,0,1)
#   return(10^power_term)
# }
#
# ill_lambda <- function(alpha,D){
#   i <- 1:D
#   power_term <- 0.5*(i-1)/(D-1)
#   diagonal <- diag(D)%*%(alpha^power_term)
# }
