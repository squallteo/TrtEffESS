library(tidyverse)
library(ggplot2)
library(ggpubr)

# theta0 <- 0
# s <- 1

logOR_prior <- function(theta0, s){
  p1 <- seq(0.01,0.99, 0.0025)
  p0 <- seq(0.01,0.99, 0.0025)
  lordt <- expand_grid(p1, p0) %>% mutate(logOR=log((p1*(1-p0))/(p0*(1-p1)))) %>% 
    mutate(lor_min=qnorm(0.0025,theta0,s), lor_max=qnorm(1-0.0025, theta0, s)) %>%
    mutate(surface=(logOR>=lor_min & logOR<=lor_max), density=dnorm(logOR, theta0, s))
  
  lordt_int <- lordt %>% filter(surface) 
  
  ggplot(lordt_int, aes(x=p0,y=p1,fill=density)) + geom_tile() +
    scale_fill_gradient(low="blue", high="orange") + 
    geom_abline(slope=1,intercept=0, color="black", linewidth=2, linetype="dashed") + 
    xlab("Control Rate p0") + ylab("Treatment Rate p1") + 
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)
    )
}


(plot1 <- logOR_prior(theta0=0.4, s=0.2))
(plot2 <- logOR_prior(theta0=0, s=1))

# png("logOR_prior.png", width = 800, height = 400)
ggarrange(plot1, plot2, labels=c("A","B"),ncol=2,nrow=1)
# dev.off()
