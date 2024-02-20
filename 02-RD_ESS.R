rm(list=ls())
library(tidyverse)
library(ggplot2)
library(mvtnorm)

source("00-Functions.R")

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

  ess_iu <- (rddt$var_iu %*% rddt$transden) * grid_width^2 / covmat[2,2]
  ess_pt <- ess_iu * sum(rand_ratio)
  
  return(c(ess_iu, ess_pt))
}


theta0 <- 0.4
s <- 0.1
mu0 <- -1
m0 <- 1
corr <- -0.8
prior_mean <- c(mu0, theta0)
prior_covmat <- matrix(c(m0^2, corr*m0*s, corr*m0*s, s^2), 2, 2)

iu_size <- c(2, 1)
# iu_multiplier <- c(1,2,5)
iu_multiplier <- c(1)

grid_width <- 0.005

plot_prior_bin(prior_mean, prior_covmat, emtype = "rd", IU = c(2,1), VarIU = F)



for(m in 1:length(iu_multiplier)){
  rand_ratio <- iu_size * iu_multiplier[m]
  rr <- RD_ess(prior_mean, prior_covmat, rand_ratio=rand_ratio, correction=correction, grid_width=grid_width)
  if(m==1){out_prior <- rr}
  else {out_prior <- rbind(out_prior, rr)}
}
out_prior

####################
####################
####################
set.seed(712)
source("00-Functions.R")

nsim <- 100
n1 <- 100; n0 <- 50
p1 <- 0.65; p0 <- 0.4

for(sim in 1:nsim){
  y1 <- rbinom(1, n1, p1); y0 <- rbinom(1, n0, p0)
  #n, y: control and treatment
  result <- reparfit_bin(n = c(n0,n1), y = c(y0, y1), emtype = "rd", initial_guess = c(0,0))
  
  post <- prior2post(prior_mean = prior_mean, prior_covmat = prior_covmat,
                     par_est = result$x, par_covmat = result$covariance_matrix)
  
  # par_hat <- result$x
  # par_covmat <- result$covariance_matrix
  # 
  # Mmat <- solve(covmat) + solve(par_covmat)
  # bvec <- solve(par_covmat) %*% par_hat + solve(covmat) %*% matrix(muvec, nrow = 2)
  # 
  # post_meanvec <- solve(Mmat) %*% bvec
  # post_covmat <- solve(Mmat)
  
  tt <- RD_ess(post$post_mean, post$post_covmat, rand_ratio=rand_ratio, correction=correction, grid_width=grid_width)
  if(sim==1) {out_post <- tt}
  else {out_post <- rbind(out_post, tt)}
}

out_prior
colMeans(out_post)
colMeans(out_post) - out_prior - n1 - n0

