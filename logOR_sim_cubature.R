rm(list=ls())

set.seed(712)
library(ggplot2)
library(tidyverse)

# x <- 0.4
# y <- 0.2
# n0vec <- 1:150
n0 <- 100
correction <- 0.5
theta0 <- 0
s <- 0.1
rr <- 2




#x[1] is p1, x[2] is p0
# lOR_iF <- function(x, n0=n0, rr=rr, correction=correction){
lOR_iF <- function(x, n0, rr, correction){
  n11 <- rbinom(10000, n0*rr, x[1])
  n10 <- n0*rr - n11
  n01 <- rbinom(10000, n0, x[2])
  n00 <- n0 - n01
  
  mean(1/(1/(n11+correction) + 1/(n01+correction) + 1/(n10+correction) + 1/(n00+correction)), na.rm = T)
}


int_func <- function(x, n0=n0, rr=rr, correction=correction, theta0=theta0, s=s){
  dnorm(log((x[1]*(1-x[2])) / (x[2]*(1-x[1]))), mean=theta0, sd=s) * (1/lOR_iF(x, n0=n0, rr=rr, correction=correction))
}

library(cubature)

nudge <- 0.01
(
tt <- 
pcubature(int_func,
          lowerLimit = c(0,0)+nudge, upperLimit = c(1,1)-nudge,
          n0 = n0, rr = rr, correction = correction,
          theta0 = theta0, s = s
)
)

(s^-2 * tt$integral) * n0 * (rr+1)
