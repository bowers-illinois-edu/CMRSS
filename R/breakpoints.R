# Helpers that avoid recomputing what `min_stat` already determines.
#
# `min_stat` (R/CMRSS_CRE.R:253-281) is called in tight loops: the matrix
# builders in R/CMRSS_SRE.R need its value at every exempt count from zero
# up to the block's treated count, for each block and each rank score
# function, at every threshold a search visits.  Profiling in
# memos/memo_cmrss_reuse_2026-09-07.qmd put `min_stat` at 61 percent of a
# confidence-interval run and `rank`, called from inside it, at 48.
#
# Most of that work is repeated.  The functions here compute the same
# quantities from fewer rankings.

#' Thresholds at which the combined statistic can change
#'
#' Returns the sorted, unique set of within-block treated-minus-control
#' outcome differences.  These are the only values of the threshold \code{c}
#' at which \code{\link{min_stat}}, and therefore the combined statistic
#' built from it, can change value.
#'
#' @param Z An \eqn{n} dimensional treatment assignment vector.
#' @param Y An \eqn{n} dimensional observed outcome vector.
#' @param block An \eqn{n} dimensional vector giving each unit's block.
#' @param block.sum Optional precomputed output of \code{summary_block}.
#'
#' @return A sorted numeric vector of candidate thresholds, with at most
#'   \eqn{\sum_b m_b (n_b - m_b)} elements.
#'
#' @keywords internal
min_stat_breakpoints <- function(Z, Y, block, block.sum = NULL) {
  if (is.null(block.sum)) {
    block.sum <- summary_block(Z, block)
  }
  block <- block.sum$block

  # Differences only matter within a block, because ranks are computed
  # within a block.  A treated unit in one block never crosses a control
  # unit in another.
  diffs <- lapply(block.sum$units.block, function(units) {
    z <- Z[units]
    if (all(z == 1) || all(z == 0)) {
      return(numeric(0))
    }
    as.vector(outer(Y[units][z == 1], Y[units][z == 0], "-"))
  })

  sort(unique(unlist(diffs, use.names = FALSE)))
}


#' Smallest index at which a monotone predicate first holds
#'
#' Binary search over \code{1:n} for the smallest index where
#' \code{is_true} holds, assuming it is false up to some point and true
#' thereafter.  Used so that locating a confidence limit costs
#' \eqn{O(\log n)} evaluations of the test statistic instead of the linear
#' walk the previous search performed.
#'
#' @param is_true A function of an integer index returning a single logical,
#'   false then true as the index increases.
#' @param n Length of the index range.
#'
#' @return The smallest index where \code{is_true} holds, or
#'   \code{NA_integer_} if it never does.
#'
#' @keywords internal
binary_search_first <- function(is_true, n) {
  if (n == 0L) {
    return(NA_integer_)
  }
  if (!is_true(n)) {
    return(NA_integer_)
  }

  lo <- 1L
  hi <- n
  # Invariant: is_true(hi) holds, and is_true(lo - 1) does not (vacuous at
  # lo == 1).  Each pass halves hi - lo, so this ends in ceiling(log2(n))
  # evaluations.
  while (lo < hi) {
    mid <- lo + (hi - lo) %/% 2L
    if (is_true(mid)) {
      hi <- mid
    } else {
      lo <- mid + 1L
    }
  }

  lo
}


#' Cache of coefficient matrices keyed on the threshold
#'
#' The coefficient matrices that \code{comb_matrix_block} and
#' \code{max_comb_matrix_block_stratum} build depend on the threshold
#' \code{c} but not on the quantile index \code{k}, while the optimization
#' that consumes them depends on \code{k} through the coverage bound
#' \code{p}.  The confidence-interval loop runs over every \code{k} and
#' revisits the same thresholds, so without a cache the same matrices are
#' rebuilt many times.  Profiling in
#' \code{memos/memo_cmrss_reuse_2026-09-07.qmd} recorded 520 rebuilds in a
#' single interval.
#'
#' Keys are formatted with 17 significant digits, which round-trips a double
#' exactly, so two thresholds share a cache entry only when they are the
#' same number.
#'
#' @param builder A function of one numeric argument returning the
#'   coefficient object for that threshold.
#'
#' @return A function of one numeric argument with the same values as
#'   \code{builder}, backed by a cache private to this call.
#'
#' @keywords internal
memoize_on_threshold <- function(builder) {
  cache <- new.env(parent = emptyenv())

  function(c) {
    key <- sprintf("%.17g", c)
    hit <- cache[[key]]
    if (!is.null(hit)) {
      return(hit)
    }
    val <- builder(c)
    assign(key, val, envir = cache)
    val
  }
}


#' One probe point strictly inside each constancy interval
#'
#' The breakpoints cut the real line into intervals on which the combined
#' statistic is constant.  This returns one point from each, chosen so it
#' does not depend on whatever bracket a particular search happens to use.
#' That matters because \code{com_block_conf_quant_larger_trt} runs a
#' separate search for every quantile index with a different upper bracket:
#' bracket-dependent probes would make each search ask for thresholds no
#' other search asks for, and the coefficient cache would rarely hit.
#'
#' @param breakpoints Sorted candidate thresholds from
#'   \code{\link{min_stat_breakpoints}}.
#'
#' @return A numeric vector the same length as \code{breakpoints}, whose
#'   \eqn{i}th entry lies strictly between the \eqn{i}th breakpoint and the
#'   next one, and above the last breakpoint for the final entry.
#'
#' @keywords internal
interval_probes <- function(breakpoints) {
  M <- length(breakpoints)
  if (M == 0L) {
    return(numeric(0))
  }
  if (M == 1L) {
    return(breakpoints + 1)
  }

  # Above the last breakpoint the statistic never changes again, so any
  # point there serves; one unit above keeps the choice canonical.
  c((breakpoints[-M] + breakpoints[-1L]) / 2, breakpoints[M] + 1)
}


#' Confidence limit as the exact infimum of the non-rejected thresholds
#'
#' Given a bracket the caller has already checked --- \code{f} positive at
#' \code{c.min} and not positive at \code{c.max} --- returns the infimum of
#' the set of thresholds at which the test does not reject.  Because
#' \code{f} is constant strictly between consecutive breakpoints, that
#' infimum is a breakpoint, and the interval it starts is found by a binary
#' search rather than by stepping.
#'
#' @section Why the probe points sit between breakpoints:
#' \code{min_stat} ranks with \code{ties.method = "first"}
#' (\code{R/CMRSS_CRE.R:278}).  Exactly at a breakpoint \eqn{Y_i - c} equals
#' some \eqn{Y_j}, and which of the two the ranking puts first depends on
#' their order in the input vector, not on the data.  So the value of
#' \code{f} at a breakpoint agrees sometimes with the interval below it and
#' sometimes with the interval above: on the 12-unit block used in
#' \code{tests/testthat/test-min_stat-breakpoints.R} the split is 18 and 18.
#' Evaluating there would therefore make the answer depend on unit ordering.
#' Probing the open intervals instead avoids the ties entirely, and leaves
#' this function unaffected by any later change to the tie rule.
#'
#' Reporting the infimum is the conservative choice.  Where the non-rejected
#' set is open at its left end, no smallest non-rejected threshold exists;
#' returning the infimum makes the reported interval contain every
#' non-rejected threshold, which is what an inverted test has to do.
#'
#' @param f Non-increasing function of the threshold, positive where the
#'   test rejects.
#' @param breakpoints Sorted candidate thresholds from
#'   \code{\link{min_stat_breakpoints}}.
#' @param probes Output of \code{\link{interval_probes}} for the same
#'   breakpoints.
#' @param c.min,c.max The bracket the caller has already evaluated.
#'
#' @return The confidence limit, exactly.
#'
#' @keywords internal
search_limit <- function(f, breakpoints, probes, c.min, c.max) {
  # Intervals whose left endpoint lies inside the bracket.  A probe may sit
  # above c.max when the interval straddles it, which is harmless: the
  # statistic is constant across the whole interval, so its value there is
  # the value below c.max too.
  idx <- which(breakpoints > c.min & breakpoints < c.max)

  if (length(idx) == 0L) {
    # Nothing inside the bracket changes the statistic, so the non-rejected
    # set begins at the top of the bracket.
    return(c.max)
  }

  j <- binary_search_first(function(j) f(probes[idx[j]]) <= 0, length(idx))

  if (is.na(j)) {
    # The test rejects on every interval below c.max, so the non-rejected
    # set begins at c.max itself.
    return(c.max)
  }

  breakpoints[idx[j]]
}


#' Every exempt count's statistic from a single ranking
#'
#' Returns \code{min_stat} for a block at exempt counts
#' \eqn{0, 1, \ldots, m_b}, the whole sequence \code{comb_matrix_block}
#' needs, computed from one call to \code{rank} instead of \eqn{m_b + 1} of
#' them.
#'
#' @details
#' Write \eqn{r_0} for the ranks of the vector at exempt count zero, where
#' treated units sit at \eqn{Y_i - c} and controls at \eqn{Y_j}.  The null
#' exempts treated units in order of decreasing outcome, and every treated
#' unit shifts by the same \eqn{c}, so the treated units' order under
#' \eqn{r_0} is their outcome order and the \eqn{ii} exempted units are the
#' \eqn{ii} treated units with the largest \eqn{r_0}.
#'
#' Two things follow.  No remaining unit has an exempted unit below it, so
#' each remaining rank is its \eqn{r_0} shifted up by \eqn{ii}.  And the
#' exempted units occupy ranks \eqn{1, \ldots, ii}, all of them treated, so
#' they contribute \code{sum(score[1:ii])} however they are ordered among
#' themselves.  Hence
#'
#'   \deqn{\mathrm{min\_stat}(ii) = \sum_{l \le ii} s_l +
#'         \sum_{j > ii} s_{ii + r_0(u_j)},}
#'
#' where \eqn{s} is the score vector and \eqn{u_1, \ldots, u_{m_b}} are the
#' treated units ordered by \eqn{r_0}.  The identity is checked against
#' \code{min_stat} on random designs, with ties and at exact breakpoints, in
#' \code{tests/testthat/test-min_stat-path.R}.
#'
#' @param Zb Treatment assignment within one block.
#' @param Yb Observed outcomes within the same block.
#' @param c The threshold.
#' @param score Rank score vector of length \code{length(Yb)}.
#'
#' @return A numeric vector of length \code{sum(Zb) + 1} giving the
#'   statistic at exempt counts \code{0:sum(Zb)}.
#'
#' @keywords internal
min_stat_path <- function(Zb, Yb, c, score) {
  # The one ranking everything else is read off. ties.method matches
  # min_stat at R/CMRSS_CRE.R:278 so tie-breaking carries over unchanged.
  r0 <- rank(Yb - Zb * c, ties.method = "first")

  # Treated ranks in increasing order: exempting the ii largest outcomes
  # removes the last ii of these.
  treated_ranks <- sort(r0[Zb == 1])
  mb <- length(treated_ranks)

  # Scores of the exempt block, whose size grows by one at each step.
  exempt_total <- c(0, cumsum(score[seq_len(mb)]))

  vapply(0:mb, function(ii) {
    kept <- treated_ranks[seq_len(mb - ii)]
    exempt_total[ii + 1L] + sum(score[ii + kept])
  }, numeric(1))
}
