#' Plug-in Sampling Single Synthetic Dataset Generation
#'
#' This function is used to generate a single synthetic version of the original data via Plug-in Sampling.
#'
#' Assume that \eqn{\mathbf{X}=\left(\mathbf{x}_1, \dots, \mathbf{x}_n\right)} is the original data, assumed to be normally distributed,
#' we compute \eqn{\bar{\mathbf{x}}} as the sample mean and \eqn{\hat{\boldsymbol{\Sigma}}=\mathbf{S}/(n-1)} as the sample covariance matrix,
#' where \eqn{\mathbf{S}} is the sample Wishart matrix.
#' We generate \eqn{\mathbf{V}=\left(\mathbf{v}_1, \dots, \mathbf{v}_n\right)}, by drawing
#'
#' \deqn{\mathbf{v}_i\stackrel{i.i.d.}{\sim}N_p(\bar{\mathbf{x}},\hat{\boldsymbol{\Sigma}}).}
#'
#' @param X matrix or dataframe
#' @param n_imp sample size
#'
#' @references
#'  Klein, M., Moura, R. and Sinha, B. (2021). Multivariate Normal Inference based on Singly Imputed Synthetic Data under Plug-in Sampling. Sankhya B 83, 273–287.
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats var
#' @examples
#' library(MASS)
#' n_sample = 100
#' mu=c(0,0,0,0)
#' Sigma=diag(1,4,4)
#' # Create original simulated dataset
#' df_o = mvrnorm(n_sample, mu, Sigma)
#' # Create singly imputed synthetic dataset
#' df_s = simSynthData(df_o)
#' #Estimators synthetic
#' mean_s <- colMeans(df_s)
#' S_s <- (t(df_s)- mean_s) %*% t(t(df_s)- mean_s)
#' # careful about this computation
#' # mean_o is a column vector and if you are thinking as n X p matrices and
#' # row vectors you should be aware of this.
#' print(mean_s)
#' print(S_s/(dim(df_s)[1]-1))
#' @export

simSynthData <- function(X, n_imp = dim(X)[1]){
  X <- as.matrix(X)
  mean_X <- colMeans(X)
  S_X <- t(X - mean_X) %*% (X - mean_X)
  V <- MASS::mvrnorm(n_imp, mean_X, var(X) * dim(X)[1]/ (dim(X)[1] - 1))
  return (V)
}
