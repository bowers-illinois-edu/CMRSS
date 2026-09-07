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
