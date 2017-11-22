rm(list = ls())
setwd("~/GitHub/CoupledPF/rscript/")
# load package
library(CoupledPF)
library(doRNG)
library(ggthemes)
ncores <- 8
registerDoMC(cores = ncores)
# load custom theme for the plots
setmytheme()
# fix the random seed
set.seed(17)
dimension <- 5
ar1 <- get_ar(dimension)
alpha_star <- 0.95
A_star <- create_A(alpha_star,dimension)
datalength <- 10000
observations <- matrix(nrow = datalength, ncol = dimension)
x_t <- fast_rmvnorm(1,rep(0,dimension), diag(1,nrow = dimension,ncol = dimension))
for (time in 1:datalength) {
  x_t <- t(A_star %*% t(x_t)) + fast_rmvnorm(1,rep(0,dimension), diag(1,nrow = dimension,ncol = dimension))
  observations[time,] <- x_t + fast_rmvnorm(1,rep(0,dimension), diag(1,nrow = dimension,ncol = dimension))
}

save(observations,alpha_star,datalength,file='ar5data.Rdata')

