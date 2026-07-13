#' @title Simulate the Generalized Variance Null Distribution
#'
#' @description
#' Simulates the null distribution of the generalized variance pivotal
#' statistic \eqn{T_1^\star} under plug-in sampling.
#'
#' Under the multiple-release stacking result,
#' \deqn{
#'   \mathbf{(n-1)S}^\star_{\mathrm{M}} \mid \mathbf{S}
#'   \sim
#'   \mathcal{W}_p\left(
#'     Mn - 1,\;
#'     \mathbf{S}
#'   \right),
#' }
#' the Bartlett decomposition gives
#' \deqn{
#'   T_1^\star
#'   =
#'   (n-1)^p \frac{|\mathbf{S}^\star|}{|\Sigma|}
#'   \;\overset{d}{=}\;
#'   \left(\prod_{j=1}^p A_j\right)
#'   \left(\prod_{j=1}^p B_j\right),
#' }
#' where \eqn{A_j \sim \chi^2_{Mn-j}} comes from the synthetic Wishart
#' distribution with degrees of freedom \eqn{Mn - 1}, and
#' \eqn{B_j \sim \chi^2_{n-j}} comes from the original-sample Wishart
#' distribution with degrees of freedom \eqn{n - 1}. All \eqn{2p}
#' variables are mutually independent.
#'
#' For \eqn{M = 1}, \eqn{A_j} and \eqn{B_j} have the same
#' \eqn{\chi^2_{n-j}} distribution, recovering the single-release result
#' of Klein et al. (2021).
#'
#' @param nsample Original sample size \eqn{n}.
#' @param pvariates Number of variables \eqn{p}.
#' @param iterations Number of Monte Carlo draws. The default is
#'   \code{10000L}.
#' @param M Number of synthetic releases. The default is \code{1L}.
#'   The effective sample size is \eqn{N = Mn}.
#'
#' @return
#' A numeric vector of length \code{iterations} containing draws from the
#' null distribution of \eqn{T_1^\star}.
#'
#' @seealso
#' \code{\link{gv_test}}
#'
#' @references
#' Klein, M., Moura, R., and Sinha, B. (2021). Multivariate normal
#' inference based on singly imputed synthetic data under plug-in
#' sampling. \emph{Sankhya B}, \strong{83}, 273--287.
#' \doi{10.1007/s13571-019-00215-9}
#'
#' @export
#'
#'
#' @examples
#' set.seed(1)
#'
#' \donttest{
#' # Single release
#' nd1 <- GVdist(nsample = 50, pvariates = 4, M = 1, iterations = 1000L)
#' stats::quantile(nd1, probs = c(0.025, 0.975))
#'
#' # Five releases
#' nd5 <- GVdist(nsample = 50, pvariates = 4, M = 5, iterations = 1000L)
#' stats::quantile(nd5, probs = c(0.025, 0.975))
#' }
GVdist <- function(nsample, pvariates, iterations = 10000L, M = 1L) {
  args <- .validate_distribution_args(
    nsample = nsample,
    pvariates = pvariates,
    iterations = iterations,
    M = M
  )

  n <- args$n
  p <- args$p
  it <- args$iterations
  N <- args$N

  df_A <- (N - 1L) - seq_len(p) + 1L
  df_B <- (n - 1L) - seq_len(p) + 1L

  replicate(it, {
    A <- prod(stats::rchisq(p, df = df_A))
    B <- prod(stats::rchisq(p, df = df_B))
    A * B
  })
}


#' @title Simulate the Sphericity Null Distribution
#'
#' @description
#' Simulates the null distribution of the sphericity pivotal statistic
#' \eqn{T_2^\star} under plug-in sampling for \eqn{M \geq 1} releases.
#'
#' Under the stacking result, the compound Wishart representation gives
#' \deqn{
#'   T_2^\star
#'   \overset{d}{=}
#'   \frac{|W_1 W_2|^{1/p}}{\mathrm{tr}(W_1 W_2)/p},
#' }
#' where
#' \eqn{W_1 \sim \mathcal{W}_p((n - 1)^{-1} I_p,\, n - 1)} and
#' \eqn{W_2 \sim \mathcal{W}_p(I_p,\, Mn - 1)} independently.
#'
#' For \eqn{M = 1}, both Wishart distributions have degrees of freedom
#' \eqn{n - 1}, recovering the single-release result of Klein et al.
#' (2021).
#'
#' @param nsample Original sample size \eqn{n}.
#' @param pvariates Number of variables \eqn{p}.
#' @param iterations Number of Monte Carlo draws. The default is
#'   \code{10000L}.
#' @param M Number of synthetic releases. The default is \code{1L}.
#'
#' @return
#' A numeric vector of length \code{iterations}.
#'
#' @seealso
#' \code{\link{sphericity_test}}
#'
#' @export
#'
#'
#' @examples
#' set.seed(1)
#'
#' \donttest{
#' nd1 <- Sphdist(nsample = 50, pvariates = 4, M = 1, iterations = 1000L)
#' stats::quantile(nd1, probs = 0.05)
#'
#' nd5 <- Sphdist(nsample = 50, pvariates = 4, M = 5, iterations = 1000L)
#' stats::quantile(nd5, probs = 0.05)
#' }
#'
Sphdist <- function(nsample, pvariates, iterations = 10000L, M = 1L) {
  args <- .validate_distribution_args(
    nsample = nsample,
    pvariates = pvariates,
    iterations = iterations,
    M = M
  )

  n <- args$n
  p <- args$p
  it <- args$iterations
  N <- args$N

  Ip <- diag(p)

  replicate(it, {
    W1 <- stats::rWishart(
      n = 1L,
      df = n - 1L,
      Sigma = Ip / (n - 1L)
    )[, , 1L]

    W2 <- stats::rWishart(
      n = 1L,
      df = N - 1L,
      Sigma = Ip
    )[, , 1L]

    WW <- W1 %*% W2

    det(WW)^(1.0 / p) / (sum(diag(WW)) / p)
  })
}

#' @title Simulate the Independence Null Distribution
#'
#' @description
#' Simulates the null distribution of the independence pivotal statistic
#' \eqn{T_3^\star} under plug-in sampling for \eqn{M \geq 1} releases.
#'
#' Under the stacking result, the compound Wishart representation gives
#' \deqn{
#'   T_3^\star
#'   \overset{d}{=}
#'   \frac{|\Omega_2|}
#'        {|\Omega_{2,11}| |\Omega_{2,22}|},
#' }
#' where
#' \eqn{\Omega_1 \sim \mathcal{W}_p(n - 1, I_p)} and
#' \eqn{(n - 1)\Omega_2 \mid \Omega_1
#' \sim
#' \mathcal{W}_p(Mn - 1,\,\Omega_1)}.
#'
#'
#' @param part Size of the first variable block, \eqn{p_1}.
#' @param nsample Original sample size \eqn{n}.
#' @param pvariates Total number of variables, \eqn{p}.
#' @param iterations Number of Monte Carlo draws. The default is
#'   \code{10000L}.
#' @param M Number of synthetic releases. The default is \code{1L}.
#'
#' @return
#' A numeric vector of length \code{iterations}.
#'
#' @seealso
#' \code{\link{independence_test}}
#'
#' @export
#'
#' @examples
#' set.seed(1)
#'
#' \donttest{
#' nd1 <- Inddist(
#'   part = 2, nsample = 50, pvariates = 4,
#'   M = 1, iterations = 1000L
#' )
#' stats::quantile(nd1, probs = 0.05)
#'
#' nd5 <- Inddist(
#'   part = 2, nsample = 50, pvariates = 4,
#'   M = 5, iterations = 1000L
#' )
#' stats::quantile(nd5, probs = 0.05)
#' }
Inddist <- function(part, nsample, pvariates,
                    iterations = 10000L,
                    M = 1L) {
  args <- .validate_distribution_args(
    nsample = nsample,
    pvariates = pvariates,
    iterations = iterations,
    M = M
  )

  n <- args$n
  p <- args$p
  it <- args$iterations
  N <- args$N

  block <- .validate_distribution_part(part, p)
  p1 <- block$p1
  p2 <- block$p2

  idx1 <- seq_len(p1)
  idx2 <- seq_len(p2) + p1

  replicate(it, {
    Om1 <- stats::rWishart(
      n = 1L,
      df = n - 1L,
      Sigma = diag(p)
    )[, , 1L]

    Om2 <- stats::rWishart(
      n = 1L,
      df = N - 1L,
      Sigma = Om1 / (n - 1L)
    )[, , 1L]

    det(Om2) /
      (
        det(Om2[idx1, idx1, drop = FALSE]) *
          det(Om2[idx2, idx2, drop = FALSE])
      )
  })
}



#' @title Simulate the Canonical Regression Null Distribution
#'
#' @description
#' Simulates the null distribution of the regression pivotal statistic
#' \eqn{T_4^\star} under plug-in sampling for \eqn{M \geq 1} releases.
#'
#' The simulation uses the same compound Wishart structure as
#' \code{\link{Inddist}}:
#' \deqn{
#'   \Omega_1 \sim \mathcal{W}_p(n - 1, I_p),
#'   \qquad
#'   (n - 1)\Omega_2 \mid \Omega_1
#'   \sim
#'   \mathcal{W}_p(Mn - 1,\,\Omega_1).
#' }
#'
#' @param part Size of the first variable block, \eqn{p_1}. Must satisfy
#'   \eqn{p_1 \leq p_2}.
#' @param nsample Original sample size \eqn{n}.
#' @param pvariates Total number of variables, \eqn{p}.
#' @param iterations Number of Monte Carlo draws. The default is
#'   \code{10000L}.
#' @param M Number of synthetic releases. The default is \code{1L}.
#'
#' @return
#' A numeric vector of length \code{iterations}.
#'
#' @seealso
#' \code{\link{regression_test}}
#'
#' @export
#'
#' @examples
#' set.seed(1)
#'
#' \donttest{
#' nd1 <- canodist(
#'   part = 2, nsample = 50, pvariates = 4,
#'   M = 1, iterations = 1000L
#' )
#' stats::quantile(nd1, probs = 0.95)
#'
#' nd5 <- canodist(
#'   part = 2, nsample = 50, pvariates = 4,
#'   M = 5, iterations = 1000L
#' )
#' stats::quantile(nd5, probs = 0.95)
#' }

#' canodist(part = 2, nsample = 50, pvariates = 4, M = 5) |> quantile(0.95)
canodist <- function(part, nsample, pvariates,
                     iterations = 10000L,
                     M = 1L) {
  args <- .validate_distribution_args(
    nsample = nsample,
    pvariates = pvariates,
    iterations = iterations,
    M = M
  )

  n <- args$n
  p <- args$p
  it <- args$iterations
  N <- args$N

  block <- .validate_distribution_part(
    part = part,
    p = p,
    require_p1_le_p2 = TRUE
  )

  p1 <- block$p1
  p2 <- block$p2

  idx1 <- seq_len(p1)
  idx2 <- seq_len(p2) + p1

  replicate(it, {
    Om1 <- stats::rWishart(
      n = 1L,
      df = n - 1L,
      Sigma = diag(p)
    )[, , 1L]

    Om2 <- stats::rWishart(
      n = 1L,
      df = N - 1L,
      Sigma = Om1 / (n - 1L)
    )[, , 1L]

    Om_11 <- Om2[idx1, idx1, drop = FALSE]
    Om_12 <- Om2[idx1, idx2, drop = FALSE]
    Om_21 <- Om2[idx2, idx1, drop = FALSE]
    Om_22 <- Om2[idx2, idx2, drop = FALSE]

    Om22i <- solve(Om_22)

    num <- Om_12 %*% Om22i %*% Om_21
    den <- Om_11 - num

    det(num) / det(den)
  })
}
