library(RBesT)
# normal_prior <- mixnorm(c(1, 0, 0.5), param="ms")
normal_prior <- mixnorm(c(1, 1, 0.5), param="ms")
sigma1 <- 1
sigma0 <- 1
rand_ratio <- c(2,1)
# rand_ratio <- c(4,2)
# rand_ratio <- c(10,5)
(var_iu <- sigma1^2/rand_ratio[1] + sigma0^2/rand_ratio[2])
sigma(normal_prior) <- sqrt(var_iu)
ess(normal_prior, method = "elir")
ess(normal_prior, method = "elir") * sum(rand_ratio)

#mixture prior
w1 <- 0.5
mix_prior <- mixnorm(prior1 = c(w1, 1, 0.5), prior2 = c(1-w1, 0, 2), param="ms")
sigma(mix_prior) <- sqrt(var_iu)
ess(mix_prior, method = "elir")
ess(mix_prior, method = "elir") * sum(rand_ratio)
