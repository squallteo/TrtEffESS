rm(list=ls())
library(RBesT)

theta0 <- 1
s <- 0.5
normal_prior <- mixnorm(c(1, theta0, s), param="ms")
sigma1 <- 1; sigma0 <- 1
rand_ratio <- c(2,1)
var_iu <- sigma1^2/rand_ratio[1] + sigma0^2/rand_ratio[2]
sigma(normal_prior) <- sqrt(var_iu)


############
n1 <- 200; n0 <- 100
mean1 <- 6; mean0 <- 4.5
sample1 <- rnorm(n1, mean = mean1, sd = sigma1)
sample0 <- rnorm(n0, mean = mean0, sd = sigma0)

theta_hat <- mean(sample1) - mean(sample0)
sd_theta_hat <- sqrt(sigma1^2/n1 + sigma0^2/n0)

post_mean <- (theta0/s^2 + theta_hat/sd_theta_hat^2) / (1/s^2 + 1/sd_theta_hat^2)
post_sd <- sqrt((1/sd_theta_hat^2 + 1/s^2)^-1)


normal_posterior <- mixnorm(c(1, post_mean, post_sd), param="ms")
sigma(normal_posterior) <- sqrt(var_iu)


ess(normal_prior)
ess(normal_prior)*sum(rand_ratio)

ess(normal_posterior)
ess(normal_posterior)*sum(rand_ratio)
