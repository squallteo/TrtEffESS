rm(list=ls())
library(tidyverse)
library(ggplot2)
library(mvtnorm)


theta0 <- 0.3
s <- 0.1
mu0 <- -1
m0 <- 1
corr <- -0.8
meanvec <- c(mu0, theta0)
covmat <- matrix(c(m0^2, corr*m0*s, corr*m0*s, s^2), 2, 2)
IU <- c(2, 1)
grid_width <- 0.001
############
p1 <- seq(grid_width,1-grid_width, grid_width)
p0 <- seq(grid_width,1-grid_width, grid_width)
plotdt <- expand_grid(p1, p0) %>% 
  mutate(lp0=log(p0/(1-p0)), RD=p1-p0, 
         var_iu=p1*(1-p1)/IU[1] + p0*(1-p0)/IU[2])
J <- with(plotdt, 1/(p0*(1-p0)))
tt <- with(plotdt, cbind(lp0, RD))
bvnden <- dmvnorm(tt, mean=meanvec, sigma=covmat)
plotdt <- plotdt %>% mutate(J=J, bvnden=bvnden, transden=J*bvnden)

plot1 <- 
  plotdt %>%
  ggplot(aes(x=p0,y=p1,z=transden)) + geom_raster(aes(fill=transden)) +
  geom_contour(color="white") +
  scale_fill_gradient(low="cornflowerblue", high="red") +
  geom_abline(slope=1,intercept=0, color="black", linewidth=2, linetype="dashed") +
  xlab("Control Rate p0") + ylab("Treatment Rate p1") +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,1,0.1)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,1,0.1)) +
  ggtitle("Induced Joint Density of (p0, p1)") + guides(fill=guide_legend(title="Density"))
theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)
  )

#############
theta0 <- 0.4
meanvec <- c(mu0, theta0)
covmat <- matrix(c(m0^2, corr*m0*s, corr*m0*s, s^2), 2, 2)

p1 <- seq(grid_width,1-grid_width, grid_width)
p0 <- seq(grid_width,1-grid_width, grid_width)
plotdt <- expand_grid(p1, p0) %>% 
  mutate(lp0=log(p0/(1-p0)), RD=p1-p0, 
         var_iu=p1*(1-p1)/IU[1] + p0*(1-p0)/IU[2])
J <- with(plotdt, 1/(p0*(1-p0)))
tt <- with(plotdt, cbind(lp0, RD))
bvnden <- dmvnorm(tt, mean=meanvec, sigma=covmat)
plotdt <- plotdt %>% mutate(J=J, bvnden=bvnden, transden=J*bvnden)

plot2 <- 
  plotdt %>%
  ggplot(aes(x=p0,y=p1,z=transden)) + geom_raster(aes(fill=transden)) +
  geom_contour(color="white") +
  scale_fill_gradient(low="cornflowerblue", high="red") +
  geom_abline(slope=1,intercept=0, color="black", linewidth=2, linetype="dashed") +
  xlab("Control Rate p0") + ylab("Treatment Rate p1") +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,1,0.1)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,1,0.1)) +
  ggtitle("Induced Joint Density of (p0, p1)") + guides(fill=guide_legend(title="Density"))
theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)
  )

rd <- plot1 + plot2
ggsave(filename = "RD_priors.eps",plot = rd, device = "eps",
       dpi = 72, width = 10, height = 5, units = "cm")
