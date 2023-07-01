#' Plug-in Sampling Dataset Generation
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
#' Inddist(part = 2, nsample = 100, pvariates = 4, iterations = 2)
#'
#' @export













library(MASS)

n <- 1000
alpha <- 4
mu <- c(1, 2, 3, 4, 5)
Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
p <- length(mu)
set.seed(123)
X <- mvrnorm(n, mu, Sigma)
mean_X <- colMeans(X)
S_X <- t(X - mean_X) %*% (X - mean_X)
#Plug-in generation of dataset
V <- mvrnorm(n, mean_X, S_X / (n - 1))
#PPS generation of dataset
tilde_inverse_Sigma <- rwish(n + alpha - p - 2, solve(S_X))
tilde_Sigma <- solve(tilde_inverse_Sigma)
tilde_mu <- mvrnorm(1, mean_X, tilde_Sigma / n)
W <- mvrnorm(n, tilde_mu, tilde_Sigma)
#Estimators Plug-in
mean_V <- colMeans(V)
S_star <- t(V - mean_V) %*% (V - mean_V)
#Estimators Posterior Predictive Sampling
mean_W <- colMeans(W)
S_bullet <- t(W - mean_W) %*% (W - mean_W)
