library(tidyverse)
library(ggplot2)


###################################
logOR_ess <- function(theta0, s, 
                      rand_ratio=rand_ratio, correction=correction, grid_width=grid_width){
  p1 <- seq(grid_width,1-grid_width, grid_width)
  p0 <- seq(grid_width,1-grid_width, grid_width)
  lordt <- expand_grid(p1, p0) %>% mutate(logOR=log((p1*(1-p0))/(p0*(1-p1)))) %>% 
    mutate(lor_min=qnorm(0.0025,theta0,s), lor_max=qnorm(1-0.0025, theta0, s)) %>%
    mutate(surface=(logOR>=lor_min & logOR<=lor_max), density=dnorm(logOR, theta0, s))
  
  lordt_int <- lordt %>% filter(surface) 
  
  lor_if <- function(p1, p0, rand_ratio, correction){
    n11 <- rbinom(10000, rand_ratio[1], p1)
    n10 <- rand_ratio[1] - n11
    n01 <- rbinom(10000, rand_ratio[2], p0)
    n00 <- rand_ratio[2] - n01
    1/mean(1/(1/(n11+correction) + 1/(n01+correction) + 1/(n10+correction) + 1/(n00+correction)), na.rm = T)
  }
  if_iu <- lordt_int %>% select(p1, p0) %>% pmap_vec(lor_if, rand_ratio=rand_ratio, correction=correction)
  ess_iu <- (if_iu %*% lordt_int$density) * grid_width^2 / s^2
  ess_pt <- ess_iu * sum(rand_ratio)
  
  c(ess_iu, ess_pt)
}


theta0 <- 0
s <- 1
iu_size <- c(20, 10)
iu_multiplier <- c(1, 5, 10, 20, 40, 50, 60)

correction <- 0.2
grid_width <- 0.01

for(m in 1:length(iu_multiplier)){
  rand_ratio <- iu_size * iu_multiplier[m]
  rr <- logOR_ess(theta0, s, rand_ratio=rand_ratio, correction=correction, grid_width=grid_width)
  if(m==1){out <- rr}
  else {out <- rbind(out, rr)}
}
out
# write.csv(out, "out.csv")
