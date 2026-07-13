#' @title Generalized Variance Test and Confidence Interval
#'
#' @description
#' Tests \eqn{H_0 : |\Sigma| = |\Sigma_0|} and computes a
#' \eqn{(1-\alpha)}-level confidence interval for the generalized
#' variance \eqn{|\Sigma|}, based on \eqn{M} released plug-in sampling
#' synthetic data sets stacked into \code{V}. Setting \code{M = 1}
#' recovers the single-release procedure of Klein et al. (2021).
#'
#' @param V Stacked synthetic data set, given as an \eqn{Mn \times p}
#'   numeric matrix, as returned by \code{\link{simSynthData}}.
#' @param M Positive integer giving the number of synthetic releases.
#'   The default is \code{1L}.
#' @param Sigma A \eqn{p \times p} positive-definite matrix specifying
#'   the null covariance matrix \eqn{\Sigma_0}. Typically, this is
#'   \code{cov(X)}, where \code{X} is the original data matrix.
#' @param alpha Significance level. The default is \code{0.05}.
#' @param iterations Monte Carlo sample size used to approximate the
#'   null distribution. The default is \code{10000L}.
#' @param null_dist Optional numeric vector containing a precomputed null
#'   distribution. If supplied, \code{iterations} is ignored.
#'
#' @return
#' An object of class \code{\link{ps_test}} with component
#' \code{conf.int} giving the exact \eqn{(1-\alpha)} confidence interval for
#' \eqn{|\Sigma|}. The usual S3 methods, including \code{print},
#' \code{summary}, and \code{plot}, are available.
#'
#' @seealso
#' \code{\link{GVdist}},
#' \code{\link{ps_test}},
#' \code{\link{simSynthData}}
#'
#' @references
#' Klein, M., Moura, R., and Sinha, B. (2021). Multivariate normal
#' inference based on singly imputed synthetic data under plug-in
#' sampling. \emph{Sankhya B}, \strong{83}, 273--287.
#' \doi{10.1007/s13571-019-00215-9}
#'
#' @export
#'
#' @examples
#' data(attitude)
#'
#' set.seed(1)
#' V1 <- simSynthData(attitude)
#'
#' \donttest{
#' res <- gv_test(V1,
#'   M = 1, Sigma = cov(attitude),
#'   iterations = 1000L
#' )
#' print(res)
#' plot(res)
#'
#' set.seed(1)
#' V5 <- simSynthData(attitude, M = 5)
#' res5 <- gv_test(V5,
#'   M = 5, Sigma = cov(attitude),
#'   iterations = 1000L
#' )
#'
#' print(res5)
#' plot(res5)
#' }
gv_test <- function(V, M = 1L, Sigma,
                    alpha = 0.05,
                    iterations = 10000L,
                    null_dist = NULL) {
  V <- .validate_X(V)
  M <- .validate_M(M)
  alpha <- .validate_alpha(alpha)
  iterations <- .validate_iterations(iterations)

  dims <- .resolve_ps_dimensions(V, M)
  N <- dims$N
  n <- dims$n
  p <- dims$p

  if (missing(Sigma)) {
    stop(
      paste0(
        "'Sigma' must be supplied: the null covariance matrix Sigma_0 ",
        "for H0: |Sigma| = |Sigma_0|. Typically use cov(X)."
      ),
      call. = FALSE
    )
  }

  Sigma <- as.matrix(Sigma)

  if (!isTRUE(all.equal(dim(Sigma), c(p, p)))) {
    stop(
      sprintf("'Sigma' must be a %d x %d matrix.", p, p),
      call. = FALSE
    )
  }

  .check_pd(Sigma, "'Sigma'")

  S_star <- .compute_S_star(V)

  null_d <- if (!is.null(null_dist)) {
    .validate_null_dist(null_dist)
  } else {
    GVdist(
      nsample = n,
      pvariates = p,
      M = M,
      iterations = iterations
    )
  }

  q_lo <- as.numeric(stats::quantile(null_d, probs = alpha / 2))
  q_hi <- as.numeric(stats::quantile(null_d, probs = 1 - alpha / 2))

  scale <- (n - 1L)^p * det(S_star)
  T_obs <- scale / det(Sigma)

  pval <- min(
    2 * min(
      mean(null_d <= T_obs),
      mean(null_d >= T_obs)
    ),
    1
  )

  dec <- if (T_obs < q_lo || T_obs > q_hi) {
    "Reject H0"
  } else {
    "Fail to reject H0"
  }

  ci <- c(
    lower = scale / q_hi,
    upper = scale / q_lo
  )

  new_ps_test(
    statistic = T_obs,
    p.value   = pval,
    alpha     = alpha,
    decision  = dec,
    null.dist = null_d,
    test      = "gv",
    n         = n,
    M         = M,
    p         = p,
    conf.int  = ci
  )
}

#' @title Generalized Variance Confidence Interval
#'
#' @description
#' Backward-compatible alias for \code{\link{gv_test}}.
#'
#' @inheritParams gv_test
#'
#' @export
gv_ci <- function(V, M = 1L, Sigma,
                  alpha = 0.05,
                  iterations = 10000L,
                  null_dist = NULL) {
  if (missing(Sigma)) {
    return(
      gv_test(
        V = V,
        M = M,
        alpha = alpha,
        iterations = iterations,
        null_dist = null_dist
      )
    )
  }

  gv_test(
    V = V,
    M = M,
    Sigma = Sigma,
    alpha = alpha,
    iterations = iterations,
    null_dist = null_dist
  )
}

#' @title Sphericity Test
#'
#' @description
#' Tests \eqn{H_0 : \Sigma = \sigma^2 I_p}, that is, all variables uncorrelated
#' with equal variance. The test is based on \eqn{M} released plug-in sampling
#' synthetic data sets stacked into \code{V}. The test is left-tailed. Setting
#' \code{M = 1} recovers the single-release procedure of Klein et al. (2021).
#'
#' @param V Stacked synthetic data set, given as an \eqn{Mn \times p}
#'   numeric matrix.
#' @param M Positive integer giving the number of synthetic releases.
#'   The default is \code{1L}.
#' @param alpha Significance level. The default is \code{0.05}.
#' @param iterations Monte Carlo sample size used to approximate the
#'   null distribution. The default is \code{10000L}.
#' @param null_dist Optional numeric vector containing a precomputed null
#'   distribution. If supplied, \code{iterations} is ignored.
#'
#' @return
#' An object of class \code{\link{ps_test}}. Component
#' \code{sigma2.hat} gives the plug-in estimator
#' \eqn{\hat\sigma^2 = \mathrm{tr}(S^\star)/(p(N-1))} under \eqn{H_0}.
#'
#' @seealso
#' \code{\link{Sphdist}},
#' \code{\link{ps_test}}
#'
#' @references
#' Klein, M., Moura, R., and Sinha, B. (2021). Multivariate normal
#' inference based on singly imputed synthetic data under plug-in
#' sampling. \emph{Sankhya B}, \strong{83}, 273--287.
#' \doi{10.1007/s13571-019-00215-9}
#'
#' @export
#'
#' @examples
#' data(attitude)
#'
#' set.seed(1)
#' V5 <- simSynthData(attitude, M = 5)
#'
#' \donttest{
#' res <- sphericity_test(V5, M = 5, iterations = 1000L)
#' print(res)
#' plot(res)
#' }
sphericity_test <- function(V, M = 1L,
                            alpha = 0.05,
                            iterations = 10000L,
                            null_dist = NULL) {
  V <- .validate_X(V)
  M <- .validate_M(M)
  alpha <- .validate_alpha(alpha)
  iterations <- .validate_iterations(iterations)

  dims <- .resolve_ps_dimensions(V, M)
  N <- dims$N
  n <- dims$n
  p <- dims$p

  S_star <- .compute_S_star(V)

  T_obs <- det(S_star)^(1.0 / p) / (sum(diag(S_star)) / p)

  null_d <- if (!is.null(null_dist)) {
    .validate_null_dist(null_dist)
  } else {
    Sphdist(
      nsample = n,
      pvariates = p,
      M = M,
      iterations = iterations
    )
  }

  crit <- as.numeric(stats::quantile(null_d, probs = alpha))
  pval <- mean(null_d <= T_obs)

  dec <- if (T_obs < crit) {
    "Reject H0"
  } else {
    "Fail to reject H0"
  }

  sigma2hat <- sum(diag(S_star)) / (p * (N - 1L))

  new_ps_test(
    statistic  = T_obs,
    p.value    = pval,
    alpha      = alpha,
    decision   = dec,
    null.dist  = null_d,
    test       = "sphericity",
    n          = n,
    M          = M,
    p          = p,
    sigma2.hat = sigma2hat
  )
}

#' @title Independence Test
#'
#' @description
#' Tests \eqn{H_0 : \Sigma_{12} = \mathbf{0}}, that is, independence
#' between two subsets of variables, based on \eqn{M} released plug-in
#' sampling synthetic data sets stacked into \code{V}. Setting
#' \code{M = 1} recovers the single-release procedure of Klein et al.
#' (2021).
#'
#' The two variable blocks can be specified in exactly one of two ways:
#' \describe{
#'   \item{\code{part}}{
#'     An integer scalar. The first \code{part} columns form Block 1, and
#'     the remaining columns form Block 2. This is the original
#'     backward-compatible interface.
#'   }
#'   \item{\code{group_a} and \code{group_b}}{
#'     Integer indices or column names identifying the two blocks. Together,
#'     they must cover all columns of \code{V} exactly once. If names are
#'     used, \code{V} must have column names.
#'   }
#' }
#'
#' @param V Stacked synthetic data set, given as an \eqn{Mn \times p}
#'   numeric matrix.
#' @param M Positive integer giving the number of synthetic releases.
#'   The default is \code{1L}.
#' @param part Integer scalar giving the size of Block 1. The first
#'   \code{part} columns form Block 1, and the remaining columns form
#'   Block 2. Ignored when \code{group_a} and \code{group_b} are supplied.
#' @param group_a Integer indices or column names identifying Block 1.
#' @param group_b Integer indices or column names identifying Block 2.
#'   Together with \code{group_a}, these must cover all columns.
#' @param alpha Significance level. The default is \code{0.05}.
#' @param iterations Monte Carlo sample size used to approximate the
#'   null distribution. The default is \code{10000L}.
#' @param null_dist Optional numeric vector containing a precomputed null
#'   distribution. If supplied, \code{iterations} is ignored.
#'
#'
#' @return
#' An object of class \code{\link{ps_test}}. The null hypothesis string
#' in \code{$null.value} names the two blocks explicitly.
#'
#' @seealso
#' \code{\link{Inddist}},
#' \code{\link{ps_test}}
#'
#' @references
#' Klein, M., Moura, R., and Sinha, B. (2021). Multivariate normal
#' inference based on singly imputed synthetic data under plug-in
#' sampling. \emph{Sankhya B}, \strong{83}, 273--287.
#' \doi{10.1007/s13571-019-00215-9}
#'
#' @export
#'
#' @examples
#' data(attitude)
#'
#' set.seed(1)
#' V5 <- simSynthData(attitude, M = 5)
#'
#' \donttest{
#' # Integer interface
#' independence_test(V5, M = 5, part = 2L, iterations = 1000L)
#'
#' # Named interface
#' independence_test(
#'   V5,
#'   M = 5,
#'   group_a = c("rating", "complaints", "advance"),
#'   group_b = c("privileges", "learning", "raises", "critical"),
#'   iterations = 1000L
#' )
#' }
independence_test <- function(V, M = 1L,
                              part = NULL,
                              group_a = NULL,
                              group_b = NULL,
                              alpha = 0.05,
                              iterations = 10000L,
                              null_dist = NULL) {
  V <- .validate_X(V)
  M <- .validate_M(M)
  alpha <- .validate_alpha(alpha)
  iterations <- .validate_iterations(iterations)

  dims <- .resolve_ps_dimensions(V, M)
  n <- dims$n
  p <- dims$p

  res <- .resolve_blocks(
    V = V,
    p = p,
    part = part,
    group_a = group_a,
    group_b = group_b,
    fun_name = "independence_test"
  )

  V <- res$V
  p1 <- res$p1
  p2 <- p - p1
  lbl1 <- res$lbl1
  lbl2 <- res$lbl2

  S_star <- .compute_S_star(V)

  idx1 <- seq_len(p1)
  idx2 <- seq_len(p2) + p1

  S_11 <- S_star[idx1, idx1, drop = FALSE]
  S_22 <- S_star[idx2, idx2, drop = FALSE]

  T_obs <- det(S_star) / (det(S_11) * det(S_22))

  null_d <- if (!is.null(null_dist)) {
    .validate_null_dist(null_dist)
  } else {
    Inddist(
      part = p1,
      nsample = n,
      pvariates = p,
      M = M,
      iterations = iterations
    )
  }

  crit <- as.numeric(stats::quantile(null_d, probs = alpha))
  pval <- mean(null_d <= T_obs)

  dec <- if (T_obs < crit) {
    "Reject H0"
  } else {
    "Fail to reject H0"
  }

  new_ps_test(
    statistic = T_obs,
    p.value   = pval,
    alpha     = alpha,
    decision  = dec,
    null.dist = null_d,
    test      = "independence",
    n         = n,
    M         = M,
    p         = p,
    lbl1      = lbl1,
    lbl2      = lbl2
  )
}


#' @title Regression Test
#'
#' @description
#' Tests \eqn{H_0 : \Delta = \Delta_0} for the population regression
#' matrix \eqn{\Delta = \Sigma_{12}\Sigma_{22}^{-1}}, based on \eqn{M}
#' released plug-in sampling synthetic data sets stacked into \code{V}.
#' The test requires \eqn{p_1 \leq p_2}. Setting \code{M = 1} recovers
#' the single-release procedure of Klein et al. (2021).
#'
#' The two variable blocks can be specified in exactly one of two ways:
#' \describe{
#'   \item{\code{part}}{
#'     An integer scalar. The first \code{part} columns form the response
#'     block, and the remaining columns form the predictor block. This is
#'     the original backward-compatible interface.
#'   }
#'   \item{\code{response} and \code{predictors}}{
#'     Integer indices or column names identifying the response and
#'     predictor blocks. Together, they must cover all columns of
#'     \code{V} exactly once.
#'   }
#' }
#'
#'
#' @param V Stacked synthetic data set, given as an \eqn{Mn \times p}
#'   numeric matrix.
#' @param M Positive integer giving the number of synthetic releases.
#'   The default is \code{1L}.
#' @param part Integer scalar giving the size of the response block.
#'   The first \code{part} columns form the response block, and the
#'   remaining columns form the predictor block. Must satisfy
#'   \eqn{p_1 \leq p_2}. Ignored when \code{response} and
#'   \code{predictors} are supplied.
#' @param Delta0 A \eqn{p_1 \times p_2} matrix giving the null value
#'   \eqn{\Delta_0}. The default is the zero matrix, corresponding to a
#'   test of zero regression.
#' @param response Integer or character vector identifying the response
#'   block.
#' @param predictors Integer or character vector identifying the predictor
#'   block. Together with \code{response}, these must cover all columns.
#' @param alpha Significance level. The default is \code{0.05}.
#' @param iterations Monte Carlo sample size used to approximate the
#'   null distribution. The default is \code{10000L}.
#' @param null_dist Optional numeric vector containing a precomputed null
#'   distribution. If supplied, \code{iterations} is ignored.
#'
#' @return
#' An object of class \code{\link{ps_test}}. Component \code{Delta.hat}
#' gives the plug-in slope estimator
#' \eqn{\hat\Delta = S_{12}^\star (S_{22}^\star)^{-1}}.
#' The null hypothesis string in \code{$null.value} names both blocks.
#'
#' @seealso
#' \code{\link{canodist}},
#' \code{\link{ps_test}}
#'
#' @references
#' Klein, M., Moura, R., and Sinha, B. (2021). Multivariate normal
#' inference based on singly imputed synthetic data under plug-in
#' sampling. \emph{Sankhya B}, \strong{83}, 273--287.
#' \doi{10.1007/s13571-019-00215-9}
#'
#'
#' @export
#'
#' @examples
#' data(attitude)
#'
#' set.seed(1)
#' V5 <- simSynthData(attitude, M = 5)
#'
#' \donttest{
#' # Integer interface: zero regression
#' regression_test(V5, M = 5, part = 2L, iterations = 1000L)
#'
#' # Named interface with Delta0 estimated from the original data
#' S0 <- cov(attitude)
#' response <- c("rating", "complaints", "advance")
#' predictors <- c("privileges", "learning", "raises", "critical")
#' b <- partition(S0,
#'   part1 = response,
#'   part2 = predictors
#' )
#' Delta0 <- b$B %*% solve(b$D)
#'
#' regression_test(
#'   V5,
#'   M = 5,
#'   response = response,
#'   predictors = predictors,
#'   Delta0 = Delta0,
#'   iterations = 1000L
#' )
#' }
regression_test <- function(V, M = 1L,
                            part = NULL,
                            Delta0 = NULL,
                            response = NULL,
                            predictors = NULL,
                            alpha = 0.05,
                            iterations = 10000L,
                            null_dist = NULL) {
  V <- .validate_X(V)
  M <- .validate_M(M)
  alpha <- .validate_alpha(alpha)
  iterations <- .validate_iterations(iterations)

  dims <- .resolve_ps_dimensions(V, M)
  n <- dims$n
  p <- dims$p

  res <- .resolve_blocks(
    V = V,
    p = p,
    part = part,
    group_a = response,
    group_b = predictors,
    fun_name = "regression_test"
  )

  V <- res$V
  p1 <- res$p1
  p2 <- p - p1
  lbl1 <- res$lbl1
  lbl2 <- res$lbl2

  if (p1 > p2) {
    stop(
      sprintf(
        "Regression test requires p1 <= p2. Got p1 = %d and p2 = %d.",
        p1, p2
      ),
      call. = FALSE
    )
  }

  if (is.null(Delta0)) {
    Delta0 <- matrix(0.0, nrow = p1, ncol = p2)
  } else {
    Delta0 <- as.matrix(Delta0)

    if (!isTRUE(all.equal(dim(Delta0), c(p1, p2)))) {
      stop(
        sprintf("'Delta0' must be a %d x %d matrix.", p1, p2),
        call. = FALSE
      )
    }

    if (any(!is.finite(Delta0))) {
      stop(
        "'Delta0' contains non-finite values: NA, NaN, Inf, or -Inf.",
        call. = FALSE
      )
    }
  }

  idx1 <- seq_len(p1)
  idx2 <- seq_len(p2) + p1

  S_star <- .compute_S_star(V)

  S_11 <- S_star[idx1, idx1, drop = FALSE]
  S_12 <- S_star[idx1, idx2, drop = FALSE]
  S_22 <- S_star[idx2, idx2, drop = FALSE]

  .check_pd(S_22, "S_22")

  S22_inv <- solve(S_22)

  Delta_hat <- S_12 %*% S22_inv
  diff <- Delta_hat - Delta0

  num <- diff %*% S_22 %*% t(diff)
  den <- S_11 - S_12 %*% S22_inv %*% t(S_12)

  T_obs <- det(num) / det(den)

  null_d <- if (!is.null(null_dist)) {
    .validate_null_dist(null_dist)
  } else {
    canodist(
      part = p1,
      nsample = n,
      pvariates = p,
      M = M,
      iterations = iterations
    )
  }

  crit <- as.numeric(stats::quantile(null_d, probs = 1 - alpha))
  pval <- mean(null_d >= T_obs)

  dec <- if (T_obs > crit) {
    "Reject H0"
  } else {
    "Fail to reject H0"
  }

  new_ps_test(
    statistic = T_obs,
    p.value   = pval,
    alpha     = alpha,
    decision  = dec,
    null.dist = null_d,
    test      = "regression",
    n         = n,
    M         = M,
    p         = p,
    Delta.hat = Delta_hat,
    lbl1      = lbl1,
    lbl2      = lbl2
  )
}
