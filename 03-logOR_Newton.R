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


#Newton Raphson
newton_raphson <- function(initial_guess, n, y, tolerance = 1e-5, max_iter = 500) {
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
n <- c(100, 200)
y <- c(20, 60)
p <- y/n
# initial_guess <- c(log(p[1]/(1-p[1])), log(p[2]/(1-p[2]))-log(p[1]/(1-p[1])))
initial_guess <- c(0,0)
result <- newton_raphson(initial_guess, n=n, y=y)
result$covariance_matrix[1,2] / sqrt(result$covariance_matrix[1,1]*result$covariance_matrix[2,2])

