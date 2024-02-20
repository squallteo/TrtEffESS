#visualize the bivariate normal prior
plot_prior_bin <- function(meanvec, covmat, emtype=c("rd","or"), IU, VarIU=FALSE, grid_width=0.005){
  p1 <- seq(grid_width,1-grid_width, grid_width)
  p0 <- seq(grid_width,1-grid_width, grid_width)
  
  if(emtype=="rd"){
    plotdt <- expand_grid(p1, p0) %>% 
      mutate(lp0=log(p0/(1-p0)), RD=p1-p0, 
             var_iu=p1*(1-p1)/IU[1] + p0*(1-p0)/IU[2])
    J <- with(plotdt, 1/(p0*(1-p0)))
    tt <- with(plotdt, cbind(lp0, RD))
    bvnden <- dmvnorm(tt, mean=meanvec, sigma=covmat)
    plotdt <- plotdt %>% mutate(J=J, bvnden=bvnden, transden=J*bvnden)
  }
  
  if(emtype=="or"){
    plotdt <- expand_grid(p1, p0) %>% 
      mutate( lp0=log(p0/(1-p0)), lp1=log(p1/(1-p1)), logOR=lp1 - lp0,
             var_iu=1/(IU[1]*p1*(1-p1)) + 1/(IU[2]*p0*(1-p0)) )
    J <- with(plotdt, 1/(p0*(1-p0)*p1*(1-p1)))
    tt <- with(plotdt, cbind(lp0, logOR))
    bvnden <- dmvnorm(tt, mean=meanvec, sigma=covmat)
    plotdt <- plotdt %>% mutate(J=J, bvnden=bvnden, transden=J*bvnden)
  }
  
  plot1 <- 
  plotdt %>% 
  ggplot(aes(x=p0,y=p1,fill=transden)) + geom_tile() +
    scale_fill_gradient(low="blue", high="orange") + 
    geom_abline(slope=1,intercept=0, color="black", linewidth=2, linetype="dashed") + 
    xlab("Control Rate p0") + ylab("Treatment Rate p1") + 
    ggtitle("Induced Joint Density of (p0, p1)") +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)
    )
  print(plot1)
  
  if(VarIU){
    plot2 <- 
    plotdt %>% 
      ggplot(aes(x=p0,y=p1,fill=var_iu)) + geom_tile() +
      scale_fill_gradient(low="blue", high="orange") + 
      geom_abline(slope=1,intercept=0, color="black", linewidth=2, linetype="dashed") + 
      xlab("Control Rate p0") + ylab("Treatment Rate p1") + 
      ggtitle("Variance of IU") +
      theme(axis.text = element_text(size = 15),
            axis.title = element_text(size = 12),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15)
      )
    print(plot2)
  }
  
}

#fit the reparameterized models
reparfit_bin <- function(n, y, emtype=c("rd","or"), initial_guess, tolerance = 1e-5, max_iter = 100){
  if(emtype=="rd"){
    grad_f <- function(par, n, y) {
      a <- exp(par[1])
      b <- 1 + a
      g1 <- y[1] - n[1]*a/b + y[2]*((a/b + par[2])^-1)*(a/b^2) - 
        (n[2]-y[2])*((1/b - par[2])^-1)*(a/b^2)
      g2 <- y[2]*((a/b + par[2])^-1) - (n[2]-y[2])*((1/b - par[2])^-1)
      c(g1, g2)
    }
    hessian_f <- function(par, n, y) {
      a <- exp(par[1])
      b <- 1 + a
      H11 <- -n[1]*(a/b^2) + 
        y[2]*(-((a/b + par[2])^-2)*(a/b^2)^2 + ((a/b + par[2])^-1)*(a*(1-a)/b^3)) -
        (n[2]-y[2])*(((1/b - par[2])^-2)*(a/b^2)^2 + ((1/b - par[2])^-1)*(a*(1-a)/b^3)) 
      H12 <- -y[2]*((a/b + par[2])^-2)*(a/b^2) - (n[2]-y[2])*((1/b - par[2])^-2)*(a/b^2)
      H22 <- -y[2]*((a/b + par[2])^-2) - (n[2]-y[2])*((1/b - par[2])^-2)
      matrix(c(H11, H12, H12, H22), nrow = 2)
    }
  }
  
  if(emtype=="or"){
    grad_f <- function(par, n, y) {
      g1 <- y[1] - n[1]*exp(par[1])/(1+exp(par[1])) + 
        y[2] - n[2]*exp(par[1]+par[2])/(1+exp(par[1]+par[2]))
      g2 <- y[2] - n[2]*exp(par[1]+par[2])/(1+exp(par[1]+par[2]))
      c(g1, g2)
    }
    hessian_f <- function(par, n, y) {
      H11 <- -n[1]*exp(par[1])/(1+exp(par[1]))^2 - n[2]*exp(par[1]+par[2])/(1+exp(par[1]+par[2]))^2
      H12 <- -n[2]*exp(par[1]+par[2])/(1+exp(par[1]+par[2]))^2
      H22 <- -n[2]*exp(par[1]+par[2])/(1+exp(par[1]+par[2]))^2
      matrix(c(H11, H12, H12, H22), nrow = 2)
    }
  }
  
  #begin Newton-Raphson
  x <- initial_guess
  iter <- 0
  
  while (iter < max_iter) {
    gradient <- grad_f(par=x, n=n, y=y)
    hessian <- hessian_f(par=x, n=n, y=y)
    
    # Check for singularity of the Hessian
    if (det(hessian) == 0) {
      stop("Hessian is singular. Unable to continue.")
    }
    
    # Update the parameter using Newton-Raphson update
    x <- x - solve(hessian) %*% gradient
    # Check for convergence
    if (sqrt(sum(gradient^2)) < tolerance) {
      break
    }
    iter <- iter + 1
  }
  
  # Calculate the covariance matrix using the inverse of the Hessian
  covariance_matrix <- solve(-hessian)
  
  return(list(x = x, covariance_matrix = covariance_matrix, iterations = iter))
}

#update the prior to get the posterior
prior2post <- function(prior_mean, prior_covmat, par_est, par_covmat){
  Mmat <- solve(prior_covmat) + solve(par_covmat)
  bvec <- solve(par_covmat) %*% par_est + solve(prior_covmat) %*% matrix(prior_mean, nrow = 2)
  
  post_mean <- solve(Mmat) %*% bvec
  post_covmat <- solve(Mmat)
  
  return(list(post_mean = post_mean, post_covmat = post_covmat))
}
