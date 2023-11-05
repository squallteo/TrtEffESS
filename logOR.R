rm(list=ls())
library(RBesT)
library(ggplot2)
library(tidyverse)

my_exp <- function(n, p, correction = 0.5){
  tt1 <- dbinom(0:n, size = n, prob = p)
  tt2 <- 1/(0:n + correction)
  sum(tt1*tt2)
}

ess_lOR <- function(prior, p1, p0, n0, rr, correction){
  IU <- c(rr*n0, n0)
  
  e11 <- my_exp(n = IU[1], p = p1, correction = correction)
  e10 <- my_exp(n = IU[1], p = 1-p1, correction = correction)
  e01 <- my_exp(n = IU[2], p = p0, correction = correction)
  e00 <- my_exp(n = IU[2], p = 1-p0, correction = correction)
  
  ref_scale <- sqrt(e11 + e10 + e01 + e00)
  sigma(prior) <- ref_scale
  tt <- ess(prior, method = "elir")
  
  c(correction, n0*(1+rr), tt, tt*sum(IU))
  
}


normal_prior <- mixnorm(c(1, 0, 0.4), param="ms")

p1 <- 0.4
p0 <- 0.2
n0vec <- 1:150

for(s in n0vec){
  tt1 <- ess_lOR(prior = normal_prior, p1, p0, n0=n0vec[s], rr=2, correction = 0.5)
  tt2 <- ess_lOR(prior = normal_prior, p1, p0, n0=n0vec[s], rr=2, correction = 0.75)
  tt3 <- ess_lOR(prior = normal_prior, p1, p0, n0=n0vec[s], rr=2, correction = 1)
  ttt <- rbind(tt1, tt2, tt3)
  if(s==1) {out=ttt}
  else {out=rbind(out,ttt)}
}
out <- data.frame(out)
colnames(out) <- c("Correction", "sizeIU", "nIU", "nSubj")
out <- out %>% mutate(Correction = as.factor(Correction))


png("logOR.png", width = 800, height = 400)

ggplot(out, aes(x=sizeIU, y=nSubj, group=Correction, color=Correction)) + geom_line(size = 1.2) + 
  xlab("Size of Information Unit") + ylab("ESS in Total Number of Subjects") + theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, vjust = -8, size = 15),
        legend.text = element_text(size = 15)
  )

dev.off()
