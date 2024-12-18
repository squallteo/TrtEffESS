library(RBesT)

mu0 <- -0.254
m0 <- sqrt(0.0166)

tt <- rnorm(1000, mu0, m0)
tt1 <- inv_logit(tt)

fit <- automixfit(tt1, type = "beta", Nc=1:6)

ess(fit)
