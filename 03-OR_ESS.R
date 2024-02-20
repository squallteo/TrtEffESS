library(tidyverse)
library(ggplot2)
library(mvtnorm)


###################################
logOR_ess <- function(muvec, covmat, 
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
iu_multiplier <- c(40, 50, 60)

correction <- 0.5
grid_width <- 0.005



for(m in 1:length(iu_multiplier)){
  rand_ratio <- iu_size * iu_multiplier[m]
  rr <- logOR_ess(muvec, covmat, rand_ratio=rand_ratio, correction=correction, grid_width=grid_width)
  if(m==1){out <- rr}
  else {out <- rbind(out, rr)}
}
out
# write.csv(out, "out.csv")
