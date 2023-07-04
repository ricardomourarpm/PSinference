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
#' # Create original simulated dataset
#' df_o = simOrigData(100, 4, mu = c(0,0,0,0), Sigma = diag(1,4,4))
#'
#' #Estimators Original
#' mean_o <- colMeans(df_o)
#' S_o <- (t(df_o)- mean_o) %*% t(t(df_o)- mean_o)
#' # careful about this computation
#' # mean_o is a column vector and if you are thinking as n X p matrices and
#' # row vectors you should be aware of this.
#' print(mean_o)
#' print(S_o/(dim(df_o)[1]-1))
#' @export
#'
#'

simOrigData <- function(n_sample, n_var, mu=1:n_var, Sigma=diag(1,n_var,n_var) , seed = 1234){
  set.seed(seed)
  X <- mvrnorm(n_sample, mu, Sigma)
  return (X)
}
