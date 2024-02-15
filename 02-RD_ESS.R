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
  
  c(ess_iu, ess_pt)
}


theta0 <- 0.2
s <- 0.1
mu0 <- -1
m0 <- 2
corr <- -0.8
muvec <- c(mu0, theta0)
covmat <- matrix(c(m0^2, corr*m0*s, corr*m0*s, s^2), 2, 2)

iu_size <- c(20, 10)
iu_multiplier <- c(1)

correction <- 0.5
grid_width <- 0.01



for(m in 1:length(iu_multiplier)){
  rand_ratio <- iu_size * iu_multiplier[m]
  rr <- RD_ess(muvec, covmat, rand_ratio=rand_ratio, correction=correction, grid_width=grid_width)
  if(m==1){out <- rr}
  else {out <- rbind(out, rr)}
}
out
# write.csv(out, "out.csv")
