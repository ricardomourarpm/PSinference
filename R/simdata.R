#' Plug-in Sampling Single Synthetic Dataset Generation
#'
#' This function calculates
#'
#' @param part .
#' @param nsample .
#' @param pvariates .
#' @param iterations .
#'
#' @references
#'  ref
#' @importFrom stats rWishart
#' @examples
#'
#' @export





library(MASS)

simOrigData <- function(n_sample, n_var, mu=1:n_var, Sigma=diag(1,n_var,n_var) , seed = 1234){
  set.seed(seed)
  X <- mvrnorm(n_sample, mu, Sigma)
  return (X)
}

simSynthData <- function(X, n_imp = dim(X)[1]){
  mean_X <- colMeans(X)
  S_X <- t(X - mean_X) %*% (X - mean_X)
  V <- mvrnorm(n_imp, mean_X, var(X) * dim(X)[1]/ (dim(X)[1] - 1))
  return (V)
}

# Create original simulated dataset

df_o = simOrigData(100, 4, mu = c(0,0,0,0), Sigma = diag(1,4,4))

# Create singly imputed synthetic dataset

df_s = simSynthData(df_o)

#Estimators Original
mean_o <- colMeans(df_o)
S_o <- (t(df_o)- mean_o) %*% t(t(df_o)- mean_o) # careful about this computation
# mean_o is a column vector and if you are thinking as n X p matrices and row vectors
# you should be aware of this

#Estimators synthetic
mean_s <- colMeans(df_s)
S_s <- (t(df_s)- mean_s) %*% t(t(df_s)- mean_s) # careful about this computation

print(mean_o)
print(mean_s)
print(c(0,0,0,0))
print(S_o/(dim(df_o)[1]-1))
print(S_s/(dim(df_s)[1]-1))
print(diag(1,4,4))
