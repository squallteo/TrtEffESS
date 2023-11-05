library(RBesT)

normal_prior <- mixnorm(c(1, 0, 0.5), param="ms")

sigma1 <- 1.5
sigma2 <- 1.5

# rand_ratio <- c(2,1)
rand_ratio <- c(4,2)
ref_scale <- sqrt(sigma1^2/rand_ratio[1] + sigma2^2/rand_ratio[2])

sigma(normal_prior) <- ref_scale

ess(normal_prior, method = "elir")
