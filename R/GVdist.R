#' GVdist
#'
#' This function calculates the near-exact distribution of the pivotal random variable
#' \deqn{T_1^* = (n-1)\frac{|S^*|}{|\Sigma|}}
#' The generalized variance is defined as \eqn{\\Sigma\}
#'
#'
#'
#' @param nsample .
#' @param pvariates .
#' @param iterations .
#'
#' @references
#'  ref
#'
#' @examples
#' library(MASS)
library(matrixStats)
library(ggplot2)

# Original data created
# two different population variances (determinants)
Sigma1 <- matrix(c(1, 0.5, 0.5, 0.5,
                   0.5, 1, 0.5, 0.5,
                   0.5, 0.5, 1, 0.5,
                   0.5, 0.5, 0.5, 1), nrow = 4, ncol = 4, byrow = TRUE)

Sigma2 <- matrix(c(1, 0.5, 0, 0,
                   0.5, 2, 0, 0,
                   0, 0, 3, 0.2,
                   0, 0, 0.2, 4), nrow = 4, ncol = 4, byrow = TRUE)

df_o1 = simOrigData(100, 4, mu = c(1,0,5,0), Sigma = Sigma1)
df_o2 = simOrigData(100, 4, mu = c(1,0,5,0), Sigma = Sigma2)

# Sim data created

df_s1 = simSynthData(df_o1)
df_s2 = simSynthData(df_o2)


# Gather the 0.025 and 0.975 quantiles and construct confident interval for sigma^2
# Check that sigma^2 is inside in both cases
p1 = dim(df_s1)[2]
p2 = dim(df_s2)[2]

T1 <- GVdist(100, p1, 10000)
q975 <- quantile(T1, 0.975)
q025 <- quantile(T1, 0.025)

left <- (100-1)^p1*det(cov(df_s1)*99)/q975
right <- (100-1)^p1*det(cov(df_s1)*99)/q025

cat(left,right,'\n')
print(det(Sigma1))
# The synthetic value is inside the confidence interval of GV
#'
#' @export


GVdist <- function(nsample, pvariates, iterations) {
  q <- chiprod(pvariates, nsample - 1,iterations)
  p <- chiprod(pvariates, nsample - 1,iterations)
  return(p * q)
}

