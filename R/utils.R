chiprod <- function(dimension, degrees, iterations = 1) {
  A <- matrix(stats::rchisq(dimension * iterations, (degrees - dimension + 1):degrees),
    nrow = dimension, ncol = iterations
  )
  product <- apply(A, 2, prod)
  return(product)
}
