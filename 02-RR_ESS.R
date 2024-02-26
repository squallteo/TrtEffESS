rm(list=ls())

library(tidyverse)
library(ggplot2)
library(mvtnorm)

source("00-Functions.R")


theta0 <- 0.1
s <- 0.5
mu0 <- -0.25
m0 <- 1
corr <- -0.8
prior_mean <- c(mu0, theta0)
prior_covmat <- matrix(c(m0^2, corr*m0*s, corr*m0*s, s^2), 2, 2)

iu_size <- c(2, 1)
# iu_multiplier <- c(1, 5, 10, 20)
iu_multiplier <- c(1)

grid_width <- 0.001

plot_prior_bin(prior_mean, prior_covmat, emtype = "rr", IU = c(2,1), VarIU = T)

