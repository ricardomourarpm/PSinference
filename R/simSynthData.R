#' Plug-in Sampling Single Synthetic Dataset Generation
#'
#' This function calculates
#'
#' @param X
#' @param nsample .
#' @param pvariates .
#' @param iterations .
#'
#' @references
#'  ref
#' @importFrom stats rWishart
#' @examples
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
  mean_X <- colMeans(X)
  S_X <- t(X - mean_X) %*% (X - mean_X)
  V <- mvrnorm(n_imp, mean_X, var(X) * dim(X)[1]/ (dim(X)[1] - 1))
  return (V)
}
