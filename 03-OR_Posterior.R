library(tidyverse)
library(ggplot2)
library(mvtnorm)


###################################
OR_ess <- function(muvec, covmat, 
                      rand_ratio=rand_ratio, correction=correction, grid_width=grid_width){
  p1 <- seq(grid_width,1-grid_width, grid_width)
  p0 <- seq(grid_width,1-grid_width, grid_width)
  lordt <- expand_grid(p1, p0) %>% 
    mutate(lp0=log(p0/(1-p0)), lp1=log(p1/(1-p1)), logOR=lp1 - lp0)
  
  J <- with(lordt, 1/(p0*(1-p0)*p1*(1-p1)))
  tt <- with(lordt, cbind(lp0, logOR))
  bvnden <- dmvnorm(tt, mean=muvec, sigma=covmat)
  
  lordt <- lordt %>% mutate(J=J, bvnden=bvnden, transden=J*bvnden)
  
  lor_if <- function(p1, p0, rand_ratio, correction){
    n11 <- rbinom(10000, rand_ratio[1], p1)
    n10 <- rand_ratio[1] - n11
    n01 <- rbinom(10000, rand_ratio[2], p0)
    n00 <- rand_ratio[2] - n01
    1/mean(1/(1/(n11+correction) + 1/(n01+correction) + 1/(n10+correction) + 1/(n00+correction)), na.rm = T)
  }
  if_iu <- lordt %>% select(p1, p0) %>% pmap_vec(lor_if, rand_ratio=rand_ratio, correction=correction)
  ess_iu <- (if_iu %*% lordt$transden) * grid_width^2 / s^2
  ess_pt <- ess_iu * sum(rand_ratio)
  
  c(ess_iu, ess_pt)
}


theta0 <- 0.4
s <- 0.5
mu0 <- -1
m0 <- 2
corr <- -0.8
muvec <- c(mu0, theta0)
covmat <- matrix(c(m0^2, corr*m0*s, corr*m0*s, s^2), 2, 2)

iu_size <- c(20, 10)
iu_multiplier <- c(50)

correction <- 0.5
grid_width <- 0.01

rand_ratio <- iu_size * iu_multiplier
prior_ESS <- OR_ess(muvec, covmat, rand_ratio=rand_ratio, correction=correction, grid_width=grid_width)

#####################################
#####################################
source("03-OR_Newton.R")
nsim <- 30
n1 <- 100; n0 <- 50
p1 <- 0.5; p0 <- 0.4

for(sim in 1:nsim){
  y1 <- rbinom(1, n1, p1); y0 <- rbinom(1, n0, p0)
  #n, y: control and treatment
  result <- or_est(initial_guess = c(0,0), n = c(n0,n1), y = c(y0, y1))
  
  par_hat <- result$x
  par_covmat <- result$covariance_matrix
  
  Mmat <- solve(covmat) + solve(par_covmat)
  bvec <- solve(par_covmat) %*% par_hat + solve(covmat) %*% matrix(muvec, nrow = 2)
  
  post_meanvec <- solve(Mmat) %*% bvec
  post_covmat <- solve(Mmat)
  
  tt <- OR_ess(c(post_meanvec), post_covmat, rand_ratio=rand_ratio, correction=correction, grid_width=grid_width)
  if(sim==1) {out <- tt}
  else {out <- rbind(out, tt)}
}

prior_ESS
colMeans(out)
colMeans(out) - prior_ESS
