library(RBesT)
normal_prior <- mixnorm(c(1, 0.03, 0.077), param="ms")
# normal_prior <- mixnorm(c(1, 0.03, 0.2), param="ms")

sigma1 <- 1.2
sigma0 <- 1.2

#2-to-1
# rand_ratio <- c(2,1)
# rand_ratio <- c(4,2)
# rand_ratio <- c(10,5)
#1-to-1
# rand_ratio <- c(1,1)
# rand_ratio <- c(2,2)
rand_ratio <- c(5,5)

(var_iu <- sigma1^2/rand_ratio[1] + sigma0^2/rand_ratio[2])
sigma(normal_prior) <- sqrt(var_iu)
ess(normal_prior, method = "elir")
ess(normal_prior, method = "elir") * sum(rand_ratio)


