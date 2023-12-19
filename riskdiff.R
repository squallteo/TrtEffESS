library(cubature)
rand_ratio <- c(2,1)
# rand_ratio <- c(4,2)
# rand_ratio <- c(10,5)
theta0 <- 0.3
s <- 0.1
#x[1] is p1, x[2] is p0
int_func <- function(x, rand_ratio=rand_ratio, theta0=theta0, s=s){
  (x[1]*(1-x[1])/rand_ratio[1] + x[2]*(1-x[2])/rand_ratio[2]) * 
    dnorm(x[1]-x[2], mean=theta0, sd=s)
}
Evar_iu <- pcubature(int_func, lowerLimit = c(0,0), upperLimit = c(1,1),
                     rand_ratio=rand_ratio, theta0=theta0, s=s)
Evar_iu$integral / s^2
(Evar_iu$integral / s^2) * sum(rand_ratio)