logit <- function(x){
  log(x/(1-x))
}
source("00-Functions.R")
nsim <- 100

n1 <- 100; n0 <- 50
p1 <- 0.7; p0 <- 0.4
trueval <- c(logit(p0), log(p1)-log(p0))
# trueval <- c(logit(p0), logit(p1)-logit(p0))
# trueval <- c(logit(p0), p1-p0)


for(sim in 1:nsim){
  y1 <- rbinom(1, n1, p1); y0 <- rbinom(1, n0, p0)
  p0hat <- y0/n0; p1hat <- y1/n1
  #n, y: control and treatment
  # fit <- reparfit_bin(n = c(n0,n1), y = c(y0, y1), emtype = "rr", 
                      # initial_guess = c(logit(p0hat), log(p1hat)-log(p0hat)))
  fit <- reparfit_bin(n = c(n0,n1), y = c(y0, y1), emtype = "rr", 
                      initial_guess = c(0,0))
  est_vec <- c(fit$x)
  se_vec <- sqrt(c(fit$covariance_matrix[1,1], fit$covariance_matrix[2,2]))
  lb_vec <- est_vec - qnorm(0.975) * se_vec
  ub_vec <- est_vec + qnorm(0.975) * se_vec
  
  tt <- as_tibble(t(rbind(est_vec, lb_vec, ub_vec)))
  tt1 <- tibble(par = c("par1", "par2"))
  ttt <- cbind(tt1, tt)
  if(sim==1){out=ttt}
  else{out=rbind(out, ttt)}
}

dt1 <- out %>% filter(par=="par1") %>% 
  mutate(truth=trueval[1], cover=(lb_vec<=truth & truth<=ub_vec)*1) 
mean(dt1$est_vec)
mean(dt1$cover)

dt2 <- out %>% filter(par=="par2") %>% 
  mutate(truth=trueval[2], cover=(lb_vec<=truth & truth<=ub_vec)*1) 
mean(dt2$est_vec)
mean(dt2$cover)

trueval
