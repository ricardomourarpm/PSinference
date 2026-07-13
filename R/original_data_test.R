#' @title Classical Generalized Variance Test (Original Data)
#'
#' @description
#' Tests \eqn{H_0 : |\Sigma| = |\Sigma_0|} using the original data matrix
#' \code{X}. The null value \code{Sigma0} is typically supplied by the
#' user.
#'
#' If \code{Sigma0 = cov(X)}, the null determinant is estimated from the
#' same data used to compute the test statistic. In that case, the
#' determinant ratio is one, the chi-square statistic is zero, and the
#' \eqn{p}-value is one by construction.
#'
#' Supplying a different \code{Sigma0} gives a likelihood-ratio
#' chi-square test with Bartlett correction.
#'
#' The test statistic is
#' \deqn{
#'   \chi^2 = -\Bigl[(n-1) - \tfrac{2p^2+p+2}{6p}\Bigr]
#'             \log\!\Bigl(\tfrac{|\hat\Sigma|}{|\Sigma_0|}\Bigr),
#' }
#' referred to a \eqn{\chi^2} distribution with
#' \eqn{\tfrac{1}{2}p(p+1) - 1} degrees of freedom.
#'
#' @param X Original data matrix with dimension \eqn{n \times p}.
#' @param Sigma0 A \eqn{p \times p} positive-definite null covariance
#'   matrix. The default is \code{cov(X)}.
#' @param alpha Significance level. The default is \code{0.05}.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{\code{statistic}}{Observed \eqn{\chi^2} statistic.}
#'   \item{\code{p.value}}{\eqn{p}-value.}
#'   \item{\code{df}}{Degrees of freedom, \eqn{p(p+1)/2 - 1}.}
#'   \item{\code{det.Sigma.hat}}{Value of \eqn{|\hat\Sigma|}.}
#'   \item{\code{decision}}{Character string: \code{"Reject H0"} or
#'     \code{"Fail to reject H0"}.}
#'   \item{\code{alpha}}{Significance level used.}
#'   \item{\code{n}, \code{p}}{Sample size and number of variables.}
#' }
#'
#' @seealso
#' \code{\link{gv_test}}
#'
#' @references
#' Anderson, T. W. (1984). An Introduction to Multivariate Statistical Analysis
#' (2nd ed.). John Wiley & Sons, New York.
#'
#' @export
#'
#' @examples
#' data(brittany_soil_ps)
#' X <- brittany_soil_ps
#'
#' ## Null = MLE  =>  p-value = 1 by construction
#' original_gv_test(X)
#'
#' ## Test against a specific null
#' Sigma0 <- diag(ncol(X))
#' original_gv_test(X, Sigma0 = Sigma0)
original_gv_test <- function(X, Sigma0 = NULL, alpha = 0.05) {
  X <- .validate_X(X)
  n <- nrow(X)
  p <- ncol(X)

  S_hat <- stats::cov(X) # p x p sample covariance matrix

  if (is.null(Sigma0)) Sigma0 <- S_hat

  Sigma0 <- as.matrix(Sigma0)
  if (!isTRUE(all.equal(dim(Sigma0), c(p, p)))) {
    stop(sprintf("'Sigma0' must be a %d x %d matrix.", p, p), call. = FALSE)
  }
  .check_pd(Sigma0, "'Sigma0'")

  det_S <- det(S_hat)
  det_S0 <- det(Sigma0)

  # Bartlett-corrected LR statistic
  # (equals 0 when Sigma0 = S_hat, i.e., the plug-in MLE)
  correction <- (n - 1) - (2 * p^2 + p + 2) / (6 * p)
  lambda_lr <- det_S / det_S0 # ratio of determinants
  chi2_stat <- -correction * log(lambda_lr)
  df <- p * (p + 1L) / 2L - 1L
  pval <- stats::pchisq(chi2_stat, df = df, lower.tail = FALSE)

  dec <- if (pval < alpha) "Reject H0" else "Fail to Reject H0"

  structure(
    list(
      statistic = chi2_stat,
      p.value = pval,
      df = df,
      det.Sigma.hat = det_S,
      decision = dec,
      alpha = alpha,
      n = n,
      p = p
    ),
    class = "original_test"
  )
}


#' @title Classical Sphericity Test (Original Data)
#'
#' @description
#' Tests \eqn{H_0 : \Sigma = \sigma^2 I_p} using the Bartlett-Box
#' chi-square approximation applied to the original data matrix \code{X}.
#'
#' The Mauchly statistic is
#' \deqn{
#'   W = \frac{|\hat\Sigma|}{(\mathrm{tr}(\hat\Sigma)/p)^p},
#' }
#' and the test statistic is
#' \deqn{
#'   \chi^2 = -\Bigl[(n-1) - \tfrac{2p^2+p+2}{6p}\Bigr]\log W,
#' }
#' referred to a \eqn{\chi^2} distribution with
#' \eqn{\tfrac{1}{2}p(p+1)-1} degrees of freedom.
#'
#' @param X Original data matrix with dimension \eqn{n \times p}.
#' @param alpha Significance level. The default is \code{0.05}.
#'
#' @return
#' A list with components \code{statistic} (\eqn{\chi^2}),
#' \code{p.value}, \code{df}, \code{W} (Mauchly statistic),
#' \code{sigma2.hat} (plug-in estimate
#' \eqn{\hat\sigma^2 = \mathrm{tr}(\hat\Sigma)/p}), \code{decision},
#' \code{alpha}, \code{n}, and \code{p}.
#'
#' @seealso
#' \code{\link{sphericity_test}}
#'
#' @references
#' Anderson, T. W. (1984). An Introduction to Multivariate Statistical Analysis
#' (2nd ed.). John Wiley & Sons, New York.
#'
#' Bartlett, M. S. (1954). A note on the multiplying factors for various
#' chi-square approximations. \emph{Journal of the Royal Statistical
#' Society: Series B}, \strong{16}, 296--298.
#'
#'
#'
#' @export
#'
#' @examples
#' data(brittany_soil_ps)
#' original_sphericity_test(brittany_soil_ps)
original_sphericity_test <- function(X, alpha = 0.05) {
  X <- .validate_X(X)
  n <- nrow(X)
  p <- ncol(X)

  S_hat <- stats::cov(X)

  sigma2_hat <- sum(diag(S_hat)) / p
  W <- det(S_hat) / sigma2_hat^p
  correction <- (n - 1) - (2 * p^2 + p + 2) / (6 * p)
  chi2_stat <- -correction * log(W)
  df <- p * (p + 1L) / 2L - 1L
  pval <- stats::pchisq(chi2_stat, df = df, lower.tail = FALSE)

  dec <- if (pval < alpha) "Reject H0" else "Fail to Reject H0"

  structure(
    list(
      statistic  = chi2_stat,
      p.value    = pval,
      df         = df,
      W          = W,
      sigma2.hat = sigma2_hat,
      decision   = dec,
      alpha      = alpha,
      n          = n,
      p          = p
    ),
    class = "original_test"
  )
}

#' @title Classical Independence Test (Original Data)
#'
#' @description
#' Tests \eqn{H_0 : \Sigma_{12} = \mathbf{0}}, corresponding to block
#' independence, using Bartlett's factored-likelihood chi-square
#' approximation applied to the original data matrix \code{X}.
#'
#' The Wilks statistic is
#' \deqn{
#'   \Lambda =
#'   \frac{|\hat\Sigma|}
#'        {|\hat\Sigma_{11}| |\hat\Sigma_{22}|},
#' }
#' and the test statistic is
#' \deqn{
#'   \chi^2 = -\Bigl[(n-1) - \tfrac{p+3}{2}\Bigr]\log\Lambda,
#' }
#' referred to a \eqn{\chi^2} distribution with \eqn{p_1 p_2} degrees
#' of freedom.
#'
#' @param X Original data matrix with dimension \eqn{n \times p}.
#' @param part Integer scalar. The first \code{part} columns form Block 1.
#'   Ignored when \code{group_a} and \code{group_b} are supplied.
#' @param group_a Integer indices or column names identifying Block 1.
#' @param group_b Integer indices or column names identifying Block 2.
#' @param alpha Significance level. The default is \code{0.05}.
#'
#' @return
#' A list with components \code{statistic} (\eqn{\chi^2}),
#' \code{p.value}, \code{df}, \code{Lambda} (Wilks statistic),
#' \code{decision}, \code{alpha}, \code{n}, \code{p}, \code{p1},
#' \code{p2}, \code{lbl1}, and \code{lbl2}.
#'
#' @seealso
#' \code{\link{independence_test}}
#'
#' @references
#' Anderson, T. W. (1984). An Introduction to Multivariate Statistical Analysis
#' (2nd ed.). John Wiley & Sons, New York.
#'
#' @export
#'
#' @examples
#' data(brittany_soil_ps)
#' original_independence_test(brittany_soil_ps,
#'   group_a = c("pH_water", "pH_KCl"),
#'   group_b = c(
#'     "log_CEC_Metson", "log_Organic_C",
#'     "log_Total_N", "log_P_Olsen"
#'   )
#' )
original_independence_test <- function(X,
                                       part = NULL,
                                       group_a = NULL,
                                       group_b = NULL,
                                       alpha = 0.05) {
  X <- .validate_X(X)
  n <- nrow(X)
  p <- ncol(X)

  res <- .resolve_blocks(
    X, p, part, group_a, group_b,
    "original_independence_test"
  )
  X <- res$V # possibly reordered
  p1 <- res$p1
  p2 <- p - p1
  lbl1 <- res$lbl1
  lbl2 <- res$lbl2

  S_hat <- stats::cov(X)
  idx1 <- seq_len(p1)
  idx2 <- seq_len(p2) + p1
  S_11 <- S_hat[idx1, idx1, drop = FALSE]
  S_22 <- S_hat[idx2, idx2, drop = FALSE]

  Lambda <- det(S_hat) / (det(S_11) * det(S_22))
  correction <- (n - 1) - (p + 3) / 2
  chi2_stat <- -correction * log(Lambda)
  df <- p1 * p2
  pval <- stats::pchisq(chi2_stat, df = df, lower.tail = FALSE)

  dec <- if (pval < alpha) "Reject H0" else "Fail to Reject H0"

  structure(
    list(
      statistic = chi2_stat,
      p.value   = pval,
      df        = df,
      Lambda    = Lambda,
      decision  = dec,
      alpha     = alpha,
      n         = n,
      p         = p,
      p1        = p1,
      p2        = p2,
      lbl1      = lbl1,
      lbl2      = lbl2
    ),
    class = "original_test"
  )
}

#' @title Classical Regression Test (Original Data)
#'
#' @description
#' Tests \eqn{H_0 : \Delta = \Delta_0} for the population regression
#' matrix \eqn{\Delta = \Sigma_{12}\Sigma_{22}^{-1}}, using Rao's
#' \eqn{F}-approximation to Wilks' \eqn{\Lambda} statistic applied to
#' the original data matrix \code{X}.
#'
#' The Wilks statistic is
#' \deqn{
#'   \Lambda =
#'   \frac{|\hat\Sigma_{11.2(\Delta_0)}|}
#'        {|\hat\Sigma_{11}|},
#' }
#' where
#' \deqn{
#'   \hat\Sigma_{11.2(\Delta_0)}
#'   =
#'   \hat\Sigma_{11}
#'   -
#'   (\hat\Delta - \Delta_0)
#'   \hat\Sigma_{22}
#'   (\hat\Delta - \Delta_0)^\top
#' }
#' is the residual Schur complement under \eqn{H_0}. This equals the
#' ordinary Schur complement \eqn{\hat\Sigma_{11.2}} when
#' \eqn{\Delta_0 = \hat\Delta}, giving \eqn{\Lambda = 1} and
#' \eqn{F = 0} by construction. The default \code{Delta0 = NULL} tests
#' \eqn{H_0 : \Delta = 0}, corresponding to zero regression.
#'
#'
#' @param X Original data matrix with dimension \eqn{n \times p}.
#' @param part Integer scalar giving the size of the response block
#'   (Block 1). The first \code{part} columns form the response block.
#'   Ignored when \code{response} and \code{predictors} are supplied.
#'   Must satisfy \eqn{p_1 \leq p_2}.
#' @param Delta0 A \eqn{p_1 \times p_2} null regression matrix. The
#'   default \code{NULL} sets \eqn{\Delta_0 = 0}.
#' @param response Integer or character vector identifying the response
#'   block.
#' @param predictors Integer or character vector identifying the predictor
#'   block.
#' @param alpha Significance level. The default is \code{0.05}.
#'
#' @return
#' A list with components \code{statistic} (\eqn{F}), \code{p.value},
#' \code{df1}, \code{df2}, \code{Lambda}, \code{Delta.hat}
#' (\eqn{\hat\Delta = \hat\Sigma_{12}\hat\Sigma_{22}^{-1}}),
#' \code{decision}, \code{alpha}, \code{n}, \code{p}, \code{p1},
#' \code{p2}, \code{lbl1}, and \code{lbl2}.
#'
#' @seealso
#' \code{\link{regression_test}}
#'
#' @references
#' Anderson, T. W. (1984). An Introduction to Multivariate Statistical Analysis
#' (2nd ed.). John Wiley & Sons, New York.
#'
#' Rao, C. R. (1951). An asymptotic expansion of the distribution of
#' Wilks' criterion. \emph{Bulletin of the International Statistical
#' Institute}, \strong{33}, 177--180.
#'
#' @export
#'
#' @examples
#' data(brittany_soil_ps)
#' X <- brittany_soil_ps
#'
#' ## Test H0: Delta = 0 (default)
#' original_regression_test(X,
#'   response = c("pH_water", "pH_KCl"),
#'   predictors = c(
#'     "log_CEC_Metson", "log_Organic_C",
#'     "log_Total_N", "log_P_Olsen"
#'   )
#' )
#'
#' ## Test H0: Delta = Delta_hat (MLE) => F = 0, p = 1 by construction
#' blk <- partition(cov(X), part1 = c("pH_water", "pH_KCl"))
#' Delta0 <- blk$B %*% solve(blk$D)
#' original_regression_test(X,
#'   response = c("pH_water", "pH_KCl"),
#'   predictors = c(
#'     "log_CEC_Metson", "log_Organic_C",
#'     "log_Total_N", "log_P_Olsen"
#'   ),
#'   Delta0 = Delta0
#' )
original_regression_test <- function(X,
                                     part = NULL,
                                     Delta0 = NULL,
                                     response = NULL,
                                     predictors = NULL,
                                     alpha = 0.05) {
  X <- .validate_X(X)
  n <- nrow(X)
  p <- ncol(X)
  nu <- n - 1L # Wishart degrees of freedom

  res <- .resolve_blocks(
    X, p, part, response, predictors,
    "original_regression_test"
  )
  X <- res$V
  p1 <- res$p1
  p2 <- p - p1
  lbl1 <- res$lbl1
  lbl2 <- res$lbl2

  if (p1 > p2) {
    stop(
      sprintf(
        "Regression test requires p1 <= p2. Got p1 = %d, p2 = %d.", p1, p2
      ),
      call. = FALSE
    )
  }

  if (is.null(Delta0)) {
    Delta0 <- matrix(0.0, nrow = p1, ncol = p2)
  } else {
    Delta0 <- as.matrix(Delta0)
    if (!isTRUE(all.equal(dim(Delta0), c(p1, p2)))) {
      stop(sprintf("'Delta0' must be a %d x %d matrix.", p1, p2),
        call. = FALSE
      )
    }
  }

  idx1 <- seq_len(p1)
  idx2 <- seq_len(p2) + p1
  S_hat <- stats::cov(X)
  S_11 <- S_hat[idx1, idx1, drop = FALSE]
  S_12 <- S_hat[idx1, idx2, drop = FALSE]
  S_22 <- S_hat[idx2, idx2, drop = FALSE]
  S22_inv <- solve(S_22)
  Delta_hat <- S_12 %*% S22_inv

  # Residual Schur complement under H0: Delta = Delta0
  diff <- Delta_hat - Delta0
  S_11_2_H0 <- S_11 - diff %*% S_22 %*% t(diff)

  # Wilks Lambda
  Lambda <- det(S_11_2_H0) / det(S_11)

  # Rao F-approximation
  fp <- .wilks_f_params(p1, p2, nu)
  s <- fp$s
  df1 <- fp$df1
  df2 <- fp$df2

  Lam_s <- Lambda^(1.0 / s)
  F_stat <- ((1 - Lam_s) / Lam_s) * (df2 / df1)
  pval <- stats::pf(F_stat, df1 = df1, df2 = df2, lower.tail = FALSE)

  dec <- if (pval < alpha) "Reject H0" else "Fail to Reject H0"

  structure(
    list(
      statistic = F_stat,
      p.value   = pval,
      df1       = df1,
      df2       = df2,
      Lambda    = Lambda,
      Delta.hat = Delta_hat,
      decision  = dec,
      alpha     = alpha,
      n         = n,
      p         = p,
      p1        = p1,
      p2        = p2,
      lbl1      = lbl1,
      lbl2      = lbl2
    ),
    class = "original_test"
  )
}

#' @title Print an \code{original_test} Object
#'
#' @description
#' Prints a concise summary of a classical test result.
#'
#' @param x An object of class \code{original_test}.
#' @param ... Further arguments, currently ignored.
#'
#' @return
#' Invisibly returns \code{x}.
#'
#' @exportS3Method print original_test
print.original_test <- function(x, ...) {
  bar <- strrep("-", 52)

  # Detect test type from available fields
  type <- if (!is.null(x$W) && !is.null(x$sigma2.hat)) {
    "Sphericity (Bartlett-Box)"
  } else if (!is.null(x$Lambda) && !is.null(x$Delta.hat)) {
    "Regression (Wilks Lambda F)"
  } else if (!is.null(x$Lambda)) {
    "Independence (Bartlett)"
  } else {
    "Generalised Variance (Bartlett-Box)"
  }

  cat("\n")
  cat("Classical LR Test on Original Data:", type, "\n")
  cat(bar, "\n")
  cat(sprintf("  n = %d   p = %d\n", x$n, x$p))
  cat(bar, "\n")

  if (!is.null(x$df)) {
    cat(sprintf("  Test statistic (chi2) : %.4g\n", x$statistic))
    cat(sprintf("  df                    : %d\n", x$df))
  } else {
    cat(sprintf("  Test statistic (F)    : %.4g\n", x$statistic))
    cat(sprintf("  df1 / df2             : %d / %.1f\n", x$df1, x$df2))
  }

  pval_str <- if (!is.na(x$p.value) && x$p.value < 1e-4) {
    sprintf("< %.0e", signif(x$p.value, 1))
  } else {
    sprintf("%.4f", x$p.value)
  }
  cat(sprintf("  p-value               : %s\n", pval_str))
  cat(sprintf("  alpha                 : %.2f\n", x$alpha))
  cat(sprintf("  Decision              : %s\n", x$decision))

  if (!is.null(x$sigma2.hat)) {
    cat(sprintf("  sigma2_hat            : %.4f\n", x$sigma2.hat))
  }
  if (!is.null(x$det.Sigma.hat)) {
    cat(sprintf("  |Sigma_hat|           : %.4e\n", x$det.Sigma.hat))
  }
  if (!is.null(x$Delta.hat)) {
    cat("  Delta_hat (plug-in slope):\n")
    print(round(x$Delta.hat, 4))
  }
  if (!is.null(x$lbl1) && !is.null(x$lbl2)) {
    cat(sprintf("  Blocks: {%s}  vs  {%s}\n", x$lbl1, x$lbl2))
  }

  cat(bar, "\n\n")
  invisible(x)
}
