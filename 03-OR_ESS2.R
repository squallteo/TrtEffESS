rm(list=ls())

library(tidyverse)
library(ggplot2)
library(mvtnorm)

source("00-Functions.R")

###################################
logOR_ess <- function(muvec, covmat, 
                      rand_ratio=rand_ratio, grid_width=grid_width){
  p1 <- seq(grid_width,1-grid_width, grid_width)
  p0 <- seq(grid_width,1-grid_width, grid_width)
  lordt <- expand_grid(p1, p0) %>% 
    mutate(lp0=log(p0/(1-p0)), lp1=log(p1/(1-p1)), logOR=lp1 - lp0,
           var_iu=1/(rand_ratio[1]*p1*(1-p1)) + 1/(rand_ratio[2]*p0*(1-p0))
    )
  
  J <- with(lordt, 1/(p0*(1-p0)*p1*(1-p1)))
  tt <- with(lordt, cbind(lp0, logOR))
  bvnden <- dmvnorm(tt, mean=muvec, sigma=covmat)
  
  lordt <- lordt %>% mutate(J=J, bvnden=bvnden, transden=J*bvnden)
  ess_iu <- (lordt$var_iu %*% lordt$transden) * grid_width^2 / covmat[2,2]
  ess_pt <- ess_iu * sum(rand_ratio)
  
  c(ess_iu, ess_pt)
}


theta0 <- 0.4
s <- 0.5
mu0 <- -1
m0 <- 1
corr <- -0.8
prior_mean <- c(mu0, theta0)
prior_covmat <- matrix(c(m0^2, corr*m0*s, corr*m0*s, s^2), 2, 2)

iu_size <- c(2, 1)
# iu_multiplier <- c(1, 5, 10, 20)
iu_multiplier <- c(1)

grid_width <- 0.005

plot_prior_bin(prior_mean, prior_covmat, emtype = "or", IU = c(2,1), VarIU = T)



for(m in 1:length(iu_multiplier)){
  rand_ratio <- iu_size * iu_multiplier[m]
  rr <- logOR_ess(muvec, covmat, rand_ratio=rand_ratio, grid_width=grid_width)
  if(m==1){out_prior <- rr}
  else {out_prior <- rbind(out_prior, rr)}
}
out_prior
# write.csv(out, "out.csv")

####################
####################
####################
source("00-Functions.R")


nsim <- 100
n1 <- 100; n0 <- 50
p1 <- 0.65; p0 <- 0.1

for(sim in 1:nsim){
  y1 <- rbinom(1, n1, p1); y0 <- rbinom(1, n0, p0)
  #n, y: control and treatment
  result <- reparfit_bin(n = c(n0,n1), y = c(y0, y1), emtype = "or", initial_guess = c(0,0))
  
  par_hat <- result$x
  par_covmat <- result$covariance_matrix
  
  Mmat <- solve(covmat) + solve(par_covmat)
  bvec <- solve(par_covmat) %*% par_hat + solve(covmat) %*% matrix(muvec, nrow = 2)
  
  post_meanvec <- solve(Mmat) %*% bvec
  post_covmat <- solve(Mmat)
  
  tt <- logOR_ess(c(post_meanvec), post_covmat, rand_ratio=rand_ratio, grid_width=grid_width)
  if(sim==1) {out_post <- tt}
  else {out_post <- rbind(out_post, tt)}
}

out_prior
colMeans(out_post)
colMeans(out_post) - out_prior - n1 - n0

