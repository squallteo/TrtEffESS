rd_est <- function(initial_guess, n, y, tolerance = 1e-5, max_iter = 100) {
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
  
  x <- initial_guess
  iter <- 0
  
  while (iter < max_iter) {
    # Calculate gradient and Hessian at the current point
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


#par: l0 and theta
#n, y: control and treatment
# n <- c(100, 200)
# y <- c(22, 70)
# initial_guess <- c(0,0)
# 
# result <- rd_est(initial_guess, n=n, y=y)
# result$covariance_matrix[1,2] / sqrt(result$covariance_matrix[1,1]*result$covariance_matrix[2,2])
# 
