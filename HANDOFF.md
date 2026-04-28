# HANDOFF: CMRSS paper-vs-code audit

Originally written 2026-04-26 after a paper-vs-code audit. Updated
2026-04-28 to reflect a repo move, David Kim's reply on item 1A, and
a follow-on finding that complicates David's reply. Read top to bottom.

## Current state (2026-04-28)

Work is **paused indefinitely**, awaiting David Kim's reply on a
follow-up question. Jake's framing in the Slack draft is "tackle this
after the dissertation is submitted," so the pause may run weeks or
months. A fresh Claude picking this up should not assume the pause is
short.

The two commits from this session are pushed to
`bowers-illinois-edu/CMRSS` (a fork of `davidk91919/CMRSS`):

```
53001b0 Validate k in pval_comb_block; record David's reply on item 1A
a3d09fa Add renv-based dependency management for R 4.6.0
```

The Slack follow-up message to David is drafted at
`/tmp/slack_to_david_followup.md` (not in the repo). Jake has edited
it; the version on disk is the one to send. As of this write Jake had
not yet sent it. When the next session starts, check whether the file
still exists, whether David has replied (Slack is the channel, not
email), and update the file or this handoff accordingly.

## Repo layout (post 2026-04-28 move)

The package's canonical home is now `davidk91919/CMRSS`. Prior to
2026-04-28 it lived at `bowers-illinois-edu/CMRSS`; that org repo was
mirror-pushed to David's account, renamed to
`bowers-illinois-edu/CMRSS_archive`, and archived (read-only). An
`archive` tag and `archive` branch on the archive repo mark the
pre-move HEAD `f2703fb` for posterity.

`bowers-illinois-edu/CMRSS` was then re-created as a fresh GitHub fork
of `davidk91919/CMRSS`. The local clone at
`/Users/jwbowers/repos/CMRSS_jake` has:

- `origin` -> `git@github.com:bowers-illinois-edu/CMRSS.git` (the org fork)
- `upstream` -> `git@github.com:davidk91919/CMRSS.git` (David's canonical)
- `main` tracks `origin/main`

Standard fork workflow: branch off `main`, push to `origin`, open PRs
against `upstream/main`. Pull David's work with
`git fetch upstream && git merge upstream/main`.

## The open question for David (item 1A, second pass)

David replied earlier on 2026-04-28 that line `R/CMRSS_SRE.R:1034`
(`p <- m - k`, introduced in his commit `762c4d08` on 2025-12-18) is
correct as code. He said the bug is in the documentation -- the
docstring claiming "k between 1 and n" -- and that he would fix the
docs upstream. Item 1A in `PLAN.md` was marked resolved at that point.

Subsequent investigation in the same session turned up a stronger
signal that complicates David's reply:

- `tests/testthat/test-pval-cre.R:43-56` cross-validates
  `CMRSS::pval_comb_block` against `RIQITE::pval_quantile` at the
  same `(k, c)` with `k = floor(0.8 * n)`. RIQITE uses the all-units
  convention. Pre-`762c4d08` this test passed. Post-, the call has
  `k = 40 > sum(Z) = 25`, the LP is infeasible, and the function
  silently returns `p.value = 0` -- the cross-validation no longer
  cross-validates.
- The same all-units `k = floor(0.8 * N)` (or `floor(0.9 * n)`)
  pattern appears in **ten** existing test sites:
  - `tests/testthat/test-CMRSS_SRE.R:171, 216, 344, 386, 442`
  - `tests/testthat/test-pval_scre.R:42, 90`
  - `tests/testthat/test-pval-cre.R:49`
  - `tests/testthat/test-solvers.R:95, 233`

That is not a documentation problem. It is the function's behavior
having shifted from `H_{k,c}` (all units) to `H_{k,c}^treat`
(treated only). The Slack follow-up to David asks him to choose
between:

(a) **Intentional shift to treated-only.** The function now tests
    `H_{k,c}^treat`. Then we need to update the docstring, examples,
    the ten test sites, and audit anything downstream -- in
    particular `com_block_conf_quant_larger`'s `set = "all"` path and
    the paper's experiment scripts in
    `/Users/jwbowers/repos/combined_stephenson_tests/code/` -- for
    assumptions about all-units `k`. (Per the original 2026-04-26
    audit, the paper's experiment scripts call only the wrapper
    `com_block_conf_quant_larger` and `com_conf_quant_larger_cre`,
    not `pval_comb_block` directly, so paper-side audit may be
    quick.)

(b) **`762c4d08` is a regression.** Line 1034 should revert to
    `p <- N - k`. The test suite goes back to green and there is no
    downstream audit.

**Both options are consistent with David's "the docs are wrong" reply
in some reading**, which is why the follow-up is needed before more
code moves.

## What was done in this session that survives

All committed and pushed to `bowers-illinois-edu/CMRSS`:

1. **renv setup** (commit `a3d09fa`). R 4.6.0 project library,
   `renv.lock`, `.Rprofile`, `renv/activate.R`, `renv/settings.json`,
   `renv/.gitignore`. `.Rbuildignore` updated to exclude renv from
   package build. ~305 packages installed; `gurobi` excluded
   (commercial, not on CRAN). `renv::status()` reports clean.

2. **Input validation in `pval_comb_block`** (commit `53001b0`).
   Added near `R/CMRSS_SRE.R:1036` (just below the `p <- m - k`
   line):
   ```r
   if (!is.numeric(k) || length(k) != 1L || k < 1 || k > m) {
     stop(sprintf(
       "k = %s is outside the feasible range 1..sum(Z) = %d for the current LP coverage constraint p = m - k at R/CMRSS_SRE.R:1034.",
       format(k), m
     ))
   }
   ```
   The error message is deliberately neutral about the hypothesis --
   it states the constraint without claiming `H_{k,c}` or
   `H_{k,c}^treat`. Correct under either resolution of the David
   question. Bumps `DESCRIPTION` to 0.2.6 and adds a `NEWS.md` entry.

3. **Test rewrite at
   `tests/testthat/test-pval-comb-block-p-convention.R`** (commit
   `53001b0`). Replaces the original three skipped all-units tests
   with three executable treated-only invariants:
   (i) no silent `p.value = 0 + test.stat = Inf` on out-of-range k
   (passes once the input validation guard is in place);
   (ii) valid k in `1..m_total` yields finite stat and p in `[0, 1]`;
   (iii) boundary `k = m_total` (so `p = 0`) is feasible.
   All three pass under the current code with the input-validation
   guard. They are neutral between (a) and (b).

4. **HANDOFF.md and PLAN.md** updated to record David's first reply
   and (in this rewrite) the follow-on finding.

Test edits I had **deliberately reverted** before committing:
`tests/testthat/test-CMRSS_SRE.R:344-345, 386-388, 442-444` (changing
`floor(0.8 * N)` to `floor(0.8 * m_total)`). These were assuming
interpretation (a). Reverting makes the failure pattern uniform
across all ten sites, which is easier for David to reproduce.

## How David should reproduce

The Slack message tells him:

```bash
git clone git@github.com:bowers-illinois-edu/CMRSS.git CMRSS_jake_audit
cd CMRSS_jake_audit
Rscript -e "renv::restore()"
Rscript -e "devtools::load_all('.', quiet = TRUE); testthat::test_dir('tests/testthat')"
```

He will see ten errors, all of the same shape:

```
Error in pval_comb_block(...):
k = 40 is outside the feasible range 1..sum(Z) = 25 for the current LP
coverage constraint p = m - k at R/CMRSS_SRE.R:1034.
```

The new test file `test-pval-comb-block-p-convention.R` passes 9/9 --
neutral signal, just confirms the guard works.

## Restarting after David replies

Likely scenario (a) -- intentional treated-only shift. Steps:

1. Update `pval_comb_block` docstring and example to use treated-only
   `k` (k in `1..sum(Z)`, example `k = floor(0.9 * sum(Z))` or
   similar). Match wording to David's preferred framing (`H_{k,c}^treat`).
2. Update the ten test sites listed above. Most just need
   `floor(0.8 * N)` -> `floor(0.8 * sum(Z))` or
   `floor(0.9 * n)` -> `floor(0.9 * sum(Z))`.
3. **Audit `com_block_conf_quant_larger`'s `set = "all"` path.**
   The wrapper at `R/CMRSS_SRE.R:~1422` swaps `Z <- 1 - Z` for the
   "control" and "all" branches, then calls into the same machinery.
   After the swap, `m` and `n - m` switch roles. Whether the wrapper
   correctly translates the user's intended (all-units) `k` into a
   treated-only `k` for the post-swap call needs to be checked
   carefully -- this is the most likely place for a subtle bug to
   hide.
4. **Decide what to do about `com_block_conf_quant_larger_trt`** at
   `R/CMRSS_SRE.R:1171, 1234`. Both lines still have `p <- n - k`
   (all-units). David did not touch them in `762c4d08`. If
   `pval_comb_block` is treated-only and the inversion sibling is
   all-units, they are testing different hypotheses despite the
   `_trt` suffix. Likely needs reconciliation.
5. **Audit paper experiment scripts in
   `/Users/jwbowers/repos/combined_stephenson_tests/code/`** for any
   call to `com_block_conf_quant_larger` that might depend on the
   pre-2025-12-18 behavior. The paper's reference code at
   `code/codes_20251026.R:1426` still has `p = N - k`. Re-running
   would probably be needed if the wrapper's behavior changed.
6. Update `tests/testthat/test-pval-comb-block-p-convention.R`'s
   header to remove "the inversion sibling uses `p <- n - k`" line
   if step 4 changes it.
7. `devtools::document()`, `devtools::check()`, and run the full
   test suite; expect 0 failures.

If scenario (b) -- regression -- steps:

1. Revert `R/CMRSS_SRE.R:1034` to `p <- N - k`. Keep the input
   validation guard but update its error message and feasibility
   range to `k <= n` (matching `n - k >= 0`).
2. The ten test sites and the new
   `test-pval-comb-block-p-convention.R` all pass without further
   edit. Update `test-pval-comb-block-p-convention.R` header to
   match.
3. `devtools::document()`, `devtools::check()`.

If David replies with something Jake hasn't anticipated: stop and
check with Jake, do not improvise.

## Other items still in PLAN.md (not touched this session)

Priority 1:
- **1B**. Verify column range `0:nb[i]` in `comb_matrix_block_stratum`
  (`R/CMRSS_SRE.R:704`) against paper's `eq:comb_per_stratum`. Possibly
  an over-enumeration. Independent of 1A; can move on this if 1A is
  blocked.

Priority 2:
- **2A**. Tie handling. `ties.method = "first"` in `rank()` calls at
  `R/CMRSS_CRE.R:13, 278` and `R/CMRSS_SRE.R:231, 799`. Paper may
  want `"random"` for finite-sample validity with discrete outcomes
  (`electric_teachers` has many ties). Fully independent of 1A.
- **2B**. Default `comb.method = 1` should flip to `2` per paper's
  Theorem 5. Depends on 1A being settled (the `p`/`k` semantics need
  to be fixed first).
- **2C**. `method_berger_boos` (`R/comparison_methods.R:539`) is
  mislabeled; it does a Bonferroni `pmax`, not Berger-Boos. Rename.

Priority 3:
- **3A**, **3B**, **3C**, **3D**, **3E**, **3F**. Cosmetic / docstring
  fixes. See `PLAN.md` for detail.

## Important context to preserve

- **Jake's role**: paper's first author and an applied statistician at
  UIUC. Treat the paper at
  `/Users/jwbowers/repos/combined_stephenson_tests/main.tex` as the
  authoritative description of intended behavior. Hypotheses
  `H_{k,c}` (eq:H_kc, ~line 511) and `H_{k,c}^treat` (Section 1.2)
  are both legitimate; which one this function tests is the open
  question.
- **David Kim's role**: co-author, current `cre` of the package, and
  active contributor. Do not assume his commits are wrong without
  evidence. Ask before reverting his work.
- **Paper's experiment scripts** in
  `/Users/jwbowers/repos/combined_stephenson_tests/code/` (excluding
  the embedded `codes_20251026.R` and `old_codes/`) produced the
  published numbers. They call `com_block_conf_quant_larger` (the
  wrapper) and `com_conf_quant_larger_cre`, not `pval_comb_block`
  directly. Per the original audit, paper-side numbers are likely
  unaffected by `762c4d08`, but if scenario (a) wins, the wrapper's
  `set = "all"` translation needs verification before claiming that.
- **Solvers**: HiGHS (open-source) is the default and what
  `R CMD check` will use unless Gurobi is installed locally. Tests
  use `skip_if_not(solver_available(...))`. Gurobi is in `Suggests`
  but is not in the renv lockfile.
- **renv**: project library is at
  `~/Library/Caches/org.R-project.R/R/renv/library/CMRSS_jake-2747cf55/...`
  (renv default cache layout). To restore on a fresh machine:
  `Rscript -e "renv::restore()"`. R 4.6.0 is the locked version.
- **ASCII only**: no unicode em dashes, en dashes, arrows, fancy
  quotes, ellipses, decorative bullets. Use `---`, `--`, `->`,
  straight quotes, `...`. Jake's global rule.
- **Workflow rules** at
  `/Users/jwbowers/repos/ai_workflow/CLAUDE_CODING.md`: tests-first;
  pause for review at three checkpoints (after tests, after
  implementation, at any unresolved design decision); WHY comments
  only; bump `DESCRIPTION` patch when an exported symbol or
  user-visible default changes. Project memory at
  `/Users/jwbowers/.claude/projects/-Users-jwbowers-repos-CMRSS-jake/memory/feedback_coding_workflow.md`
  points to the same file.
- **Bowers-illinois-edu org**: Jake's affiliation. The org repo is
  the working fork; David's repo is canonical.

## Files outside the repo to remember

- `/tmp/slack_to_david_followup.md` -- the active Slack draft to
  send. Edited by Jake; do not regenerate without checking with him.
- `/tmp/message_to_david.md`, `/tmp/slack_to_david.md` -- earlier
  drafts from 2026-04-26 (pre-David-reply). Probably stale; check
  before using.
- `/Users/jwbowers/repos/combined_stephenson_tests/` -- the paper
  repo. Main entry: `main.tex`. Reference R code:
  `code/codes_20251026.R`. Experiment scripts: rest of `code/`.

## Suggested next session opening

1. Read `HANDOFF.md` (this file) top to bottom.
2. Check whether `/tmp/slack_to_david_followup.md` still exists. If
   yes, ask Jake whether it has been sent and whether David has
   replied. If no, ask whether the conversation moved elsewhere.
3. If David has replied: follow the "Restarting after David replies"
   section above. Do not start coding before reading his reply.
4. If David has not replied: ask Jake whether to move on to an
   independent item (1B or 2A) or hold. Do not start work on 1A
   variations or anything that would conflict with either resolution.
5. Re-read the workflow memory at
   `/Users/jwbowers/.claude/projects/-Users-jwbowers-repos-CMRSS-jake/memory/feedback_coding_workflow.md`
   before any implementation.
