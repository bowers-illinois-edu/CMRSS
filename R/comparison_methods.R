#' Comparison Methods for Benchmarking
#'
#' This file contains comparison methods from the literature for benchmarking
#' against the CMRSS combined methods. These methods are wrappers around the
#' RIQITE package implementation.
#'
#' @name comparison_methods
NULL


#' Rank score function for comparison methods
#'
#' Computes rank scores for Wilcoxon or Stephenson statistics.
#' This version normalizes scores by the maximum value, matching
#' the RIQITE package convention.
#'
#' @param n Number of units.
#' @param method.list A list specifying the scoring method:
#'   \itemize{
#'     \item \code{name}: "Wilcoxon" or "Stephenson"
#'     \item \code{s}: (for Stephenson) the parameter s
#'   }
#'
#' @return A numeric vector of length n containing the normalized rank scores.
#'
#' @keywords internal
rank_score_riqite <- function(n, method.list = list(name = "Wilcoxon")) {

  if (method.list$name == "DIM") {
    stop("Can't calculate rank scores for DIM")
  }

  if (method.list$name == "Wilcoxon") {
    score <- c(1:n)
    score <- score / max(score)
    return(score)
  }

  if (method.list$name == "Stephenson") {
    score <- choose(c(1:n) - 1, method.list$s - 1)
    score <- score / max(score)
    return(score)
  }
}


#' Confidence intervals for effect quantiles (experimental)
#'
#' Computes (1-alpha) simultaneous confidence/prediction intervals (lower bounds)
#' for effect quantiles among treated, control, or all units using the RIQITE
#' package.
#'
#' @param Z An n-dimensional binary treatment assignment vector (1 = treated, 0 = control).
#' @param Y An n-dimensional observed outcome vector.
#' @param treat.method.list Method specification for treated units.
#' @param control.method.list Method specification for control units.
#' @param score Optional pre-computed score vector.
#' @param stat.null Optional pre-computed null distribution.
#' @param nperm Number of permutations for null distribution.
#' @param Z.perm Optional permutation matrix.
#' @param alpha Significance level.
#' @param set Set of quantiles: "treat", "control", or "all".
#' @param alpha.ratio.treat For set="all", proportion of alpha allocated to treated.
#' @param tol Tolerance for root-finding.
#'
#' @return A numeric vector of lower confidence limits.
#'
#' @keywords internal
ci_lower_quantile_exp <- function(Z, Y,
                                  treat.method.list = list(name = "Stephenson", s = 6),
                                  control.method.list = list(name = "Stephenson", s = 6),
                                  score = NULL, stat.null = NULL, nperm = 10^4, Z.perm = NULL,
                                  alpha = 0.05, set = "treat", alpha.ratio.treat = 0.5, tol = 10^(-3)) {
  n <- length(Z)

  if (set == "treat") {
    ci.treat <- RIQITE::ci_quantile(Z = Z, Y = Y, k.vec = (n - sum(Z) + 1):n,
                                    alternative = "greater",
                                    method.list = treat.method.list,
                                    score = score, stat.null = stat.null,
                                    nperm = nperm, Z.perm = Z.perm,
                                    alpha = alpha, tol = tol, switch = FALSE)
    ci.treat <- as.numeric(ci.treat$lower)
    return(ci.treat)
  }

  if (set == "control") {
    if (is.null(stat.null)) {
      stat.null.control <- NULL
    } else {
      if (is.null(score)) {
        score <- rank_score_riqite(n, method.list = control.method.list)
      }
      stat.null.control <- sum(score) - stat.null
    }

    if (is.null(Z.perm)) {
      Z.perm.control <- NULL
    } else {
      Z.perm.control <- 1 - Z.perm
    }
    ci.control <- RIQITE::ci_quantile(Z = 1 - Z, Y = -Y, k.vec = (n - sum(1 - Z) + 1):n,
                                      alternative = "greater",
                                      method.list = control.method.list,
                                      score = score, stat.null = stat.null.control,
                                      nperm = nperm, Z.perm = Z.perm.control,
                                      alpha = alpha, tol = tol, switch = FALSE)
    ci.control <- as.numeric(ci.control$lower)
    return(ci.control)
  }

  if (set == "all") {
    ci.treat <- RIQITE::ci_quantile(Z = Z, Y = Y, k.vec = (n - sum(Z) + 1):n,
                                    alternative = "greater",
                                    method.list = treat.method.list,
                                    score = score, stat.null = stat.null,
                                    nperm = nperm, Z.perm = Z.perm,
                                    alpha = alpha * alpha.ratio.treat,
                                    tol = tol, switch = FALSE)
    ci.treat <- as.numeric(ci.treat$lower)

    if (is.null(stat.null)) {
      stat.null.control <- NULL
    } else {
      if (is.null(score)) {
        score <- rank_score_riqite(n, method.list = control.method.list)
      }
      stat.null.control <- sum(score) - stat.null
    }

    if (is.null(Z.perm)) {
      Z.perm.control <- NULL
    } else {
      Z.perm.control <- 1 - Z.perm
    }
    ci.control <- RIQITE::ci_quantile(Z = 1 - Z, Y = -Y, k.vec = (n - sum(1 - Z) + 1):n,
                                      alternative = "greater",
                                      method.list = control.method.list,
                                      score = score, stat.null = stat.null.control,
                                      nperm = nperm, Z.perm = Z.perm.control,
                                      alpha = alpha * (1 - alpha.ratio.treat),
                                      tol = tol, switch = FALSE)
    ci.control <- as.numeric(ci.control$lower)

    ci.all <- sort(c(ci.treat, ci.control))
    return(ci.all)
  }
}


#' Multivariate hypergeometric correction
#'
#' Computes the correction term based on the multivariate hypergeometric
#' distribution for generalizing inference to larger populations.
#'
#' @param N Population size.
#' @param n Sample size.
#' @param K.vec Vector of quantile ranks in the population.
#' @param kk.vec Vector of corresponding sample quantile ranks.
#' @param ndraw Number of draws for Monte Carlo approximation.
#' @param hg.draw Optional pre-computed hypergeometric draws.
#'
#' @return Probability value for the correction.
#'
#' @keywords internal
correct_hypergeom <- function(N = 100, n = 50, K.vec = c(60, 80),
                              kk.vec = c(30, 40), ndraw = 10^4, hg.draw = NULL) {
  if (length(K.vec) > 1) {
    if (is.null(hg.draw)) {
      hg.draw <- extraDistr::rmvhyper(ndraw, diff(c(0, K.vec, N)), n)
      hg.draw <- hg.draw[, -1, drop = FALSE]
    }
    H.sum <- apply(hg.draw[, ncol(hg.draw):1, drop = FALSE], 1, cumsum)
    H.sum <- matrix(H.sum, nrow = ncol(hg.draw))
    H.sum <- H.sum[nrow(H.sum):1, , drop = FALSE]
    prob <- mean(apply(H.sum - (n - kk.vec) > 0, 2, max))
    return(prob)
  }

  if (length(K.vec) == 1) {
    prob <- 1 - stats::phyper(n - kk.vec, N - K.vec, K.vec, n)
    return(prob)
  }
}


#' Threshold for hypergeometric correction
#'
#' Calculates the correction term Delta using the proposed choice of kk.vec
#' for generalizing inference to larger populations.
#'
#' @inheritParams correct_hypergeom
#' @param alpha Significance level.
#' @param tol Tolerance for bisection search.
#'
#' @return A list with:
#'   \item{kappa}{The correction factor}
#'   \item{prob}{The resulting probability}
#'   \item{kk.vec}{The adjusted sample quantile ranks}
#'
#' @keywords internal
threshold_correct_hypergeom <- function(N = 100, n = 50, K.vec = c(60, 80),
                                        alpha = 0.05, ndraw = 10^4,
                                        hg.draw = NULL, tol = 10^(-2)) {

  if (length(K.vec) > 1) {
    if (is.null(hg.draw)) {
      hg.draw <- extraDistr::rmvhyper(ndraw, diff(c(0, K.vec, N)), n)
      hg.draw <- hg.draw[, -1, drop = FALSE]
    }

    kappa <- NULL

    kappa.lower <- 1 / length(K.vec)
    kappa.upper <- 1
    prob.lower <- correct_hypergeom(N = N, n = n, K.vec = K.vec,
                                    kk.vec = n - stats::qhyper(1 - kappa.lower * alpha, N - K.vec, K.vec, n),
                                    hg.draw = hg.draw)
    prob.upper <- correct_hypergeom(N = N, n = n, K.vec = K.vec,
                                    kk.vec = n - stats::qhyper(1 - kappa.upper * alpha, N - K.vec, K.vec, n),
                                    hg.draw = hg.draw)

    if (prob.upper <= alpha) {
      kappa <- kappa.upper
    }

    if (prob.lower > alpha) {
      kappa <- kappa.lower
    }

    while (is.null(kappa)) {
      kappa.mid <- (kappa.upper + kappa.lower) / 2
      prob.mid <- correct_hypergeom(N = N, n = n, K.vec = K.vec,
                                    kk.vec = n - stats::qhyper(1 - kappa.mid * alpha, N - K.vec, K.vec, n),
                                    hg.draw = hg.draw)

      if (prob.mid <= alpha) {
        kappa.lower <- kappa.mid
      } else {
        kappa.upper <- kappa.mid
      }

      if (kappa.upper - kappa.lower <= tol) {
        kappa <- kappa.lower
      }
    }

    kk.vec <- n - stats::qhyper(1 - kappa * alpha, N - K.vec, K.vec, n)
    prob <- correct_hypergeom(N = N, n = n, K.vec = K.vec, kk.vec = kk.vec, hg.draw = hg.draw)
    return(list(kappa = kappa, prob = prob, kk.vec = kk.vec))
  }

  if (length(K.vec) == 1) {
    kappa <- 1
    kk.vec <- n - stats::qhyper(1 - alpha, N - K.vec, K.vec, n)
    prob <- correct_hypergeom(N = N, n = n, K.vec = K.vec, kk.vec = kk.vec)
    return(list(kappa = kappa, prob = prob, kk.vec = kk.vec))
  }
}


#' Generalized confidence intervals with hypergeometric correction
#'
#' Internal function for computing confidence intervals with correction
#' for generalizing to a larger population.
#'
#' @inheritParams ci_lower_quantile_exp
#' @param N Population size.
#' @param K.vec Vector of quantile ranks in the population.
#' @param gamma Proportion of alpha for hypergeometric correction.
#' @param ndraw Number of Monte Carlo draws.
#'
#' @return A data frame with columns k, lower, and upper.
#'
#' @keywords internal
ci_lower_quantile_gen <- function(Z, Y, N = 2 * length(Z),
                                  K.vec = ceiling(c(0.6, 0.7, 0.8, 0.9) * N),
                                  alpha = 0.05, gamma = 0.5, ndraw = 10^4,
                                  treat.method.list = list(name = "Stephenson", s = 6),
                                  control.method.list = list(name = "Stephenson", s = 6),
                                  score = NULL, stat.null = NULL, nperm = 10^4, Z.perm = NULL,
                                  set = "all", alpha.ratio.treat = 0.5, tol = 10^(-3)) {

  if (set == "all") {
    n.star <- length(Z)
  }

  if (set == "control") {
    n.star <- sum(1 - Z)
  }

  if (set == "treat") {
    n.star <- sum(Z)
  }

  if (N > n.star) {
    step1 <- threshold_correct_hypergeom(N = N, n = n.star, K.vec = K.vec,
                                         alpha = gamma * alpha, ndraw = ndraw)
  }

  if (N == n.star) {
    step1 <- list(prob = 0, kk.vec = K.vec)
  }

  step2 <- ci_lower_quantile_exp(Z, Y, treat.method.list = treat.method.list,
                                 control.method.list = control.method.list,
                                 score = score, stat.null = stat.null,
                                 nperm = nperm, Z.perm = Z.perm,
                                 alpha = alpha - step1$prob, set = set,
                                 alpha.ratio.treat = alpha.ratio.treat, tol = tol)
  step2 <- c(-Inf, step2)  # add lower confidence limit "-Inf" for case kk.vec = 0

  ci <- step2[step1$kk.vec + 1]  # index +1 for case kk.vec = 0

  K.vec <- sort(K.vec)
  conf.int <- data.frame(k = K.vec, lower = ci, upper = Inf)
  return(conf.int)
}


#' Generalized confidence intervals for effect quantiles
#'
#' Computes confidence intervals for effect quantiles that can be generalized
#' to a larger population using the hypergeometric correction approach.
#'
#' @inheritParams ci_lower_quantile_gen
#' @param k_vec Vector of quantile ranks of interest.
#' @param simul Logical; if TRUE, compute simultaneous intervals.
#'
#' @return A data frame with columns k, lower, and upper.
#'
#' @keywords internal
ci_lower_quantile_generalize <- function(Z, Y, N, k_vec = ceiling(c(0.6, 0.7, 0.8, 0.9) * N),
                                         alpha = 0.05, gamma = 0.5, ndraw = 10^4,
                                         treat.method.list = list(name = "Stephenson", s = 6),
                                         control.method.list = list(name = "Stephenson", s = 6),
                                         score = NULL, stat.null = NULL, nperm = 10^4, Z.perm = NULL,
                                         set = "all", alpha.ratio.treat = 0.5, tol = 10^(-3), simul = TRUE) {

  k_vec <- sort(k_vec)
  if (simul) {
    res <- ci_lower_quantile_gen(Z = Z, Y = Y, N = N, K.vec = k_vec,
                                 alpha = alpha, gamma = gamma, ndraw = ndraw,
                                 treat.method.list = treat.method.list,
                                 control.method.list = control.method.list,
                                 score = score, stat.null = stat.null,
                                 nperm = nperm, Z.perm = Z.perm,
                                 set = set, alpha.ratio.treat = alpha.ratio.treat, tol = tol)
  } else {
    res <- NULL
    for (k in k_vec) {
      res <- rbind(res,
                   ci_lower_quantile_gen(Z = Z, Y = Y, N = N, K.vec = c(k),
                                         alpha = alpha, gamma = gamma, ndraw = ndraw,
                                         treat.method.list = treat.method.list,
                                         control.method.list = control.method.list,
                                         score = score, stat.null = stat.null,
                                         nperm = nperm, Z.perm = Z.perm,
                                         set = set, alpha.ratio.treat = alpha.ratio.treat, tol = tol))
    }
  }
  return(res)
}


#' Caughey et al. (2023) Method
#'
#' Wrapper for the method from Caughey et al. (2023) using the RIQITE package.
#' This provides confidence intervals for effect quantiles using a single
#' rank statistic.
#'
#' @param Z An n-dimensional binary treatment assignment vector (1 = treated, 0 = control).
#' @param Y An n-dimensional observed outcome vector.
#' @param k_vec Vector of quantile ranks to compute intervals for.
#' @param method.list A list specifying the rank statistic method:
#'   \itemize{
#'     \item \code{name}: "Wilcoxon" or "Stephenson"
#'     \item \code{s}: (for Stephenson) the parameter s
#'   }
#' @param nperm Number of permutations for the null distribution.
#' @param alpha Significance level.
#'
#' @return A data frame with columns k, lower, and upper.
#'
#' @details
#' This method is based on Caughey, Dafoe, Li, and Miratrix (2023) and uses
#' the RIQITE package implementation. It provides simultaneous confidence
#' intervals for specified quantiles using a single rank statistic.
#'
#' @references
#' Caughey, D., Dafoe, A., Li, X., and Miratrix, L. (2023). Randomization
#' Inference for Treatment Effect Quantiles. \emph{Journal of the American
#' Statistical Association}.
#'
#' @examples
#' \dontrun{
#' data(electric_teachers)
#' Z <- electric_teachers$TxAny
#' Y <- electric_teachers$gain
#'
#' # Using Stephenson statistic with s = 6
#' result <- method_caughey(Z, Y,
#'                          k_vec = floor(c(0.7, 0.8, 0.9) * length(Z)),
#'                          method.list = list(name = "Stephenson", s = 6),
#'                          nperm = 10000,
#'                          alpha = 0.05)
#' print(result)
#' }
#'
#' @seealso \code{\link{method_chen_li}} for the Chen and Li combination method,
#'   \code{\link{com_conf_quant_larger_cre}} for the CMRSS combined method
#' @export
method_caughey <- function(Z, Y, k_vec, method.list = list(name = "Wilcoxon"),
                           nperm = 10^4, alpha = 0.05) {
  return(RIQITE::ci_quantile(Z = Z, Y = Y, k.vec = k_vec,
                             alternative = "greater",
                             method.list = method.list,
                             nperm = nperm, alpha = alpha))
}


#' Chen and Li Combination Method
#'
#' Method combining treated and control inference using Bonferroni correction,
#' based on Chen and Li (2024).
#'
#' @param Z An n-dimensional binary treatment assignment vector (1 = treated, 0 = control).
#' @param Y An n-dimensional observed outcome vector.
#' @param N Population size for generalization. Set to length(Z) for sample inference.
#' @param k_vec Vector of quantile ranks to compute intervals for.
#' @param control.method.list Method specification for control units.
#' @param treat.method.list Method specification for treated units.
#' @param simul Logical; if TRUE (default), compute simultaneous intervals.
#' @param stat.null Optional pre-computed null distribution.
#' @param score Optional pre-computed score vector.
#' @param Z.perm Optional permutation matrix.
#' @param alpha Significance level.
#' @param gamma Proportion of alpha for hypergeometric correction.
#' @param alpha.ratio.treat Proportion of alpha allocated to treated units.
#' @param ndraw Number of Monte Carlo draws for hypergeometric correction.
#' @param nperm Number of permutations for null distribution.
#' @param tol Tolerance for root-finding.
#'
#' @return A data frame with columns k, lower, and upper.
#'
#' @details
#' This method combines inference from both treated and control perspectives
#' using a Bonferroni-style correction. It provides confidence intervals for
#' all effect quantiles by separately computing prediction intervals for
#' treated and control units, then combining them.
#'
#' When N > n (the sample size), the method incorporates a hypergeometric
#' correction to generalize inference to the larger population.
#'
#' @references
#' Chen, H. and Li, X. (2024). Randomization Inference for Quantile Treatment
#' Effects with Varying Numbers of Treated and Control Units.
#'
#' @examples
#' \dontrun{
#' data(electric_teachers)
#' Z <- electric_teachers$TxAny
#' Y <- electric_teachers$gain
#' n <- length(Z)
#'
#' # Confidence intervals for all effect quantiles
#' result <- method_chen_li(Z, Y,
#'                          N = n,  # sample inference
#'                          k_vec = 1:n,
#'                          treat.method.list = list(name = "Stephenson", s = 6),
#'                          control.method.list = list(name = "Stephenson", s = 6),
#'                          nperm = 10000,
#'                          alpha = 0.10)
#' print(head(result))
#' }
#'
#' @seealso \code{\link{method_caughey}} for single-direction inference,
#'   \code{\link{com_conf_quant_larger_cre}} for the CMRSS combined method
#' @export
method_chen_li <- function(Z, Y, N, k_vec,
                           control.method.list = list(name = "Stephenson", s = 6),
                           treat.method.list = list(name = "Stephenson", s = 6),
                           simul = TRUE,
                           stat.null = NULL, score = NULL, Z.perm = NULL,
                           alpha = 0.05, gamma = 0.5, alpha.ratio.treat = 0.5,
                           ndraw = 10^5, nperm = 10^4, tol = 10^(-3)) {

  return(ci_lower_quantile_generalize(Z = Z, Y = Y, N = N, k_vec = k_vec, simul = simul,
                                      alpha = alpha, gamma = gamma, ndraw = ndraw,
                                      control.method.list = control.method.list,
                                      treat.method.list = treat.method.list,
                                      score = score, stat.null = stat.null,
                                      nperm = nperm, Z.perm = Z.perm,
                                      set = "all", alpha.ratio.treat = alpha.ratio.treat,
                                      tol = tol))
}


#' Berger and Boos (1994) Method
#'
#' Method using the approach from Berger and Boos (1994) for confidence
#' intervals that combine treated and control inference.
#'
#' @inheritParams method_chen_li
#'
#' @return A data frame with columns k, lower, and upper, containing the
#'   element-wise maximum of treated and control confidence limits.
#'
#' @details
#' This method separately computes confidence intervals from the treated and
#' control perspectives, then takes the element-wise maximum. This provides
#' a more conservative but potentially more robust interval.
#'
#' @references
#' Berger, R. L. and Boos, D. D. (1994). P Values Maximized Over a Confidence
#' Set for the Nuisance Parameter. \emph{Journal of the American Statistical
#' Association}, 89(427), 1012-1016.
#'
#' @examples
#' \dontrun{
#' data(electric_teachers)
#' Z <- electric_teachers$TxAny
#' Y <- electric_teachers$gain
#' n <- length(Z)
#'
#' # Using Berger-Boos approach
#' result <- method_berger_boos(Z, Y,
#'                              N = n,
#'                              k_vec = floor(c(0.6, 0.7, 0.8, 0.9) * n),
#'                              treat.method.list = list(name = "Stephenson", s = 6),
#'                              control.method.list = list(name = "Stephenson", s = 6),
#'                              nperm = 10000,
#'                              alpha = 0.05)
#' print(result)
#' }
#'
#' @seealso \code{\link{method_chen_li}} for the Chen and Li combination method
#' @export
method_berger_boos <- function(Z, Y, N, k_vec,
                               treat.method.list = list(name = "Stephenson", s = 6),
                               control.method.list = list(name = "Stephenson", s = 6),
                               simul = TRUE,
                               stat.null = NULL, score = NULL, Z.perm = NULL,
                               alpha = 0.05, gamma = 0.5, alpha.ratio.treat = 0.5,
                               ndraw = 10^5, nperm = 10^4, tol = 10^(-3)) {

  ci.gen1 <- ci_lower_quantile_generalize(Z, Y, N, k_vec,
                                          treat.method.list = treat.method.list,
                                          stat.null = stat.null, set = "treat", simul = simul,
                                          alpha = alpha / 2, gamma = gamma,
                                          alpha.ratio.treat = alpha.ratio.treat,
                                          score = score, Z.perm = Z.perm,
                                          ndraw = ndraw, nperm = nperm, tol = tol)

  ci.gen2 <- ci_lower_quantile_generalize(Z, Y, N, k_vec,
                                          control.method.list = control.method.list,
                                          stat.null = stat.null, set = "control", simul = simul,
                                          alpha = alpha / 2, gamma = gamma,
                                          alpha.ratio.treat = alpha.ratio.treat,
                                          score = score, Z.perm = Z.perm,
                                          ndraw = ndraw, nperm = nperm, tol = tol)

  return(pmax(ci.gen1, ci.gen2))
}
