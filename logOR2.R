set.seed(712)
library(RBesT)
library(ggplot2)
library(tidyverse)

ess_lOR2 <- function(prior, p1, p0, n0, rr, correction){
  IU <- c(rr*n0, n0)
  n11 <- rbinom(100000, n0*rr, p1)
  n10 <- n0*rr - n11
  n01 <- rbinom(100000, n0, p0)
  n00 <- n0 - n01

  iF <- mean(1/(1/(n11+correction) + 1/(n01+correction) + 1/(n10+correction) + 1/(n00+correction)))
  ref_scale <- sqrt(1/iF)
  sigma(prior) <- ref_scale
  tt <- ess(prior, method = "elir")
  
  c(correction, n0*(1+rr), tt, tt*sum(IU))
  
}

p1 <- 0.4
p0 <- 0.2
n0vec <- 1:150

normal_prior <- mixnorm(c(1, 0, 0.4), param="ms")

for(s in n0vec){
  tt1 <- ess_lOR2(prior = normal_prior, p1, p0, n0=n0vec[s], rr=2, correction = 0.5)
  tt2 <- ess_lOR2(prior = normal_prior, p1, p0, n0=n0vec[s], rr=2, correction = 0.75)
  tt3 <- ess_lOR2(prior = normal_prior, p1, p0, n0=n0vec[s], rr=2, correction = 1)
  ttt <- rbind(tt1, tt2, tt3)
  if(s==1) {out2=ttt}
  else {out2=rbind(out2,ttt)}
}
out2 <- data.frame(out2)
colnames(out2) <- c("Correction", "sizeIU", "nIU2", "nSubj2")
out2 <- out2 %>% mutate(Correction = as.factor(Correction))

# png("logOR2.png", width = 800, height = 400)

ggplot(out2, aes(x=sizeIU, y=nSubj2, group=Correction, color=Correction)) + geom_line(size = 1.2) + 
  xlab("Size of Information Unit") + ylab("ESS in Total Number of Subjects") + theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, vjust = -8, size = 15),
        legend.text = element_text(size = 15)
  )

# dev.off()
