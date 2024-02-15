library(tidyverse)
library(ggplot2)
library(mvtnorm)


###################################
RD_ess <- function(muvec, covmat, 
                   rand_ratio=rand_ratio, correction=correction, grid_width=grid_width){
  p1 <- seq(grid_width,1-grid_width, grid_width)
  p0 <- seq(grid_width,1-grid_width, grid_width)
  rddt <- expand_grid(p1, p0) %>% 
    mutate(lp0=log(p0/(1-p0)), RD=p1-p0, 
           var_iu=p1*(1-p1)/rand_ratio[1] + p0*(1-p0)/rand_ratio[2])
  
  J <- with(rddt, 1/(p0*(1-p0)))
  tt <- with(rddt, cbind(lp0, RD))
  bvnden <- dmvnorm(tt, mean=muvec, sigma=covmat)
  
  rddt <- rddt %>% mutate(J=J, bvnden=bvnden, transden=J*bvnden)
  
  ess_iu <- (rddt$var_iu %*% rddt$transden) * grid_width^2 / s^2
  ess_pt <- ess_iu * sum(rand_ratio)
  
  # c(ess_iu, ess_pt)
  return(list(ess_iu, ess_pt, rddt))
}


theta0 <- 0.3
s <- 0.1
mu0 <- -1
m0 <- 2
corr <- -0.8
muvec <- c(mu0, theta0)
covmat <- matrix(c(m0^2, corr*m0*s, corr*m0*s, s^2), 2, 2)

iu_size <- c(2, 1)
iu_multiplier <- c(1)
grid_width <- 0.005

rand_ratio <- iu_size * iu_multiplier
prior_ESS <- RD_ess(muvec, covmat, rand_ratio=rand_ratio, correction=correction, grid_width=grid_width)
##########################
##########################
source("02-RD_Newton.R")
nsim <- 100
n1 <- 80; n0 <- 40
p1 <- 0.75; p0 <- 0.4

for(sim in 1:nsim){
  y1 <- rbinom(1, n1, p1); y0 <- rbinom(1, n0, p0)
  #n, y: control and treatment
  result <- rd_est(initial_guess = c(0,0), n = c(n0,n1), y = c(y0, y1))
  
  par_hat <- result$x
  par_covmat <- result$covariance_matrix
  
  Mmat <- solve(covmat) + solve(par_covmat)
  bvec <- solve(par_covmat) %*% par_hat + solve(covmat) %*% matrix(muvec, nrow = 2)
  
  post_meanvec <- solve(Mmat) %*% bvec
  post_covmat <- solve(Mmat)
  
  tt <- RD_ess(c(post_meanvec), post_covmat, rand_ratio=rand_ratio, correction=correction, grid_width=grid_width)
  if(sim==1) {out <- tt}
  else {out <- rbind(out, tt)}
}

prior_ESS
colMeans(out)
colMeans(out) - prior_ESS
