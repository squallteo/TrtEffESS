library(RBesT)

#107/245
mu0 <- -0.254
m0 <- sqrt(0.0166)

tt <- rnorm(10000, mu0, m0)
tt1 <- inv_logit(tt)

fit <- automixfit(tt1, type = "beta", Nc=1:4)

ess(fit)
