# HANDOFF: CMRSS paper-vs-code audit (2026-04-26)

A fresh Claude instance picking this up should be able to continue from
this file alone. Read it top to bottom.

## Context and goal

Jake Bowers (paper author, applied statistician at UIUC) recently submitted
the paper that develops the methods implemented in this R package. The paper
lives at `/Users/jwbowers/repos/combined_stephenson_tests/`. Jake asked
Claude to audit whether the package faithfully implements the paper's
methods and to flag any divergence.

Two parallel agent passes did the audit (paper-first and code-first). Their
findings agreed on "the core math is implemented correctly" and surfaced a
small number of items needing attention. Those items are tracked in
`PLAN.md` (in this repo).

The user's coding-workflow rules (`~/repos/ai_workflow/CLAUDE_CODING.md`)
require: tests before implementation; pause for review at three
checkpoints (after tests, after implementation, at any unresolved design
decision); WHY comments only; no unicode; bump DESCRIPTION patch version
when an exported symbol or user-visible default changes. There is also
project-level memory at
`/Users/jwbowers/.claude/projects/-Users-jwbowers-repos-CMRSS-jake/memory/`
that points to those rules.

## Files changed and why

- `PLAN.md` (created): full prioritized punch list from the audit.
  Investigation, tests-first, implementation, and Jake review checkpoints
  per item. Keep in repo or move to a tracker --- Jake's choice.
- `tests/testthat/test-pval-comb-block-p-convention.R` (created): three
  tests covering item 1A (the `p` convention bug). All `skip()`'d pending
  Jake's conversation with David Kim. Header explains the full situation.
  When David replies, drop `skip()` and the tests will pin the agreed
  behavior.
- `/tmp/message_to_david.md` (created, not in repo): draft email-style
  message Jake will send to David. Keep this one out of the repo.
- `/tmp/slack_to_david.md` (created, not in repo): Slack-friendly variant
  of the same message --- shorter, conversational, code-fenced. Jake asked
  for both. Pick whichever channel fits.

No source code in `R/` has been edited. No DESCRIPTION bump yet.

## The blocker (item 1A): `p` convention in `pval_comb_block`

The single open question right now. Until Jake hears from David, do not
edit `R/CMRSS_SRE.R:1034`.

**What we know:**

- `R/CMRSS_SRE.R:1033-1034`:
  ```
  # p <- N - k
  p <- m - k
  ```
- The change was made by David Kim in commit `762c4d08` on 2025-12-18.
  Commit message: "Update CMRSS_SRE.R / updated pval function in SRE"
  (uninformative). Most of the diff is whitespace cleanup; the substantive
  change is the one line above.
- Sibling function `com_block_conf_quant_larger_trt` at lines 1171 and 1234
  still uses `p <- n - k`. David did not touch it.
- Paper's reference code at
  `/Users/jwbowers/repos/combined_stephenson_tests/code/codes_20251026.R:1426`
  still has `p = N - k`. The paper's experiment scripts call only
  `com_block_conf_quant_larger` (the wrapper) and `com_conf_quant_larger_cre`,
  not `pval_comb_block` directly. So the **paper's published results are
  very likely unaffected** by David's change.
- Reproducer (Jake ran this with the existing SRE test inputs):
  `pval_comb_block(Z, Y, k=19, c=0, ...)` with N=24, m=12 yields `p = -7`,
  HiGHS reports infeasible, test statistic returns `Inf`, and the function
  silently returns `p.value = 0` (false rejection).
- Paper's hypotheses (main.tex Section 1.2): `H_{k,c}` over all n units uses
  `p = n - k` (eq:H_kc, line 511); `H_{k,c}^treat` over treated only uses
  `p = m - k` with k in 1..m. Both formulas can be "correct" depending on
  which hypothesis the function tests.
- Evidence the function is supposed to test the all-units hypothesis:
  docstring says "k between 1 and n", docstring example uses
  `k = floor(0.9 * n)`, the SRE test uses `k = floor(0.8 * N)`, the
  inversion sibling uses `p = n - k`. All four point to all-units.

**What's blocked on:** Jake will Slack/email David asking what bug he was
fixing in `762c4d08`. Two scenarios are open:

1. Regression: revert line 1034 to `p <- N - k`. Most likely.
2. Deliberate redirect to treated-only: keep `m - k`, but update the
   docstring, example, test, and wrapper accordingly, and audit any paper
   table that might depend on the original meaning.

When David replies, drop the `skip()` calls in
`tests/testthat/test-pval-comb-block-p-convention.R`, run the tests
(they should fail under current code), apply the fix, and bump the
DESCRIPTION patch version.

## Done

- Two-agent audit completed. Findings reported to Jake (in conversation,
  not a file).
- `PLAN.md` written: items 1A, 1B, 2A, 2B, 2C, 3A--3F, with priorities,
  independence map, suggested worktree split, and per-item investigation
  / tests-first / implementation / checkpoint structure.
- Item 1A reproduced and the bug confirmed: `pval_comb_block` returns
  `p.value = 0` silently for any `k > sum(Z)` due to infeasible LP.
- Git archaeology: the change is David's, not a refactor artifact.
- Paper-side check: paper's experiment scripts don't call
  `pval_comb_block` directly, so paper results likely unaffected.
- Test file written for item 1A, all skipped, header explains situation.
- Draft message to David written in `/tmp/message_to_david.md`.
- Project memory updated to point at `~/repos/ai_workflow/CLAUDE_CODING.md`.

## Remaining (in priority order; full detail in PLAN.md)

Priority 1 (possible bugs):

- **1A**: Resolve the `p` convention. Blocked on David. Tests already
  written and skipped. Once unblocked: drop skips, apply fix, bump
  DESCRIPTION patch.
- **1B**: Verify column range `0:nb[i]` in `comb_matrix_block_stratum`
  (`R/CMRSS_SRE.R:704`) vs paper's eq:comb_per_stratum which appears to
  index `0:mb[i]`. Brute-force test on a tiny SRE example. Possibly an
  over-enumeration; could change numerical results.

Priority 2 (paper-code consistency):

- **2A**: Tie handling. Code uses `ties.method = "first"` in `rank()`
  calls (`R/CMRSS_CRE.R:13, 278`; `R/CMRSS_SRE.R:231, 799`). Paper's
  finite-sample validity may require `"random"` for discrete outcomes
  (`electric_teachers` has many ties). Decide between (i) change default
  to `"random"`, (ii) document the choice, (iii) add a `ties.method`
  argument.
- **2B**: Default `comb.method = 1` in `pval_comb_block`
  (`R/CMRSS_SRE.R:1023`) and `com_block_conf_quant_larger`
  (`R/CMRSS_SRE.R:1422`) contradicts the paper's recommendation
  (`comb.method = 2`, Theorem 5). Flip default. Bump version. NEWS entry.
  This task depends on 1A being resolved first (the `p`/`k` semantics
  must be settled).
- **2C**: `method_berger_boos` (`R/comparison_methods.R:539`) is mislabeled.
  It does a Bonferroni `pmax` of two one-sided CIs at `alpha/2`, not the
  Berger-Boos nuisance-parameter procedure. Rename (suggest
  `method_bonferroni_two_sided`); add deprecated wrapper.

Priority 3 (cosmetic but worth doing):

- **3A**: Comment the `j = m_s - k_s` index relabeling in
  `solve_optimization` (`R/solvers.R`).
- **3B**: Comment that the LP coverage constraint is `<= p` not `= p` as
  in the paper, and explain monotonicity argument.
- **3C**: Docstring fix for `set = "all"` `alpha` semantics (joint level
  via `alpha/2` Bonferroni split, not paper's `1 - 2*alpha` framing).
- **3D**: Harmonize defaults across CRE and SRE (`null.max`: 10^4 vs 10^5;
  `set`: "treat" vs "all").
- **3E**: Docstring note on `Polynomial(std = FALSE)` vs paper
  eq:polynomial.
- **3F**: Add missing test invariants identified by the audit if not
  already covered by 1A/1B/2B/3C tests.

## Important context to preserve

- Jake is the paper's first author. Treat the paper at
  `/Users/jwbowers/repos/combined_stephenson_tests/main.tex` as the
  authoritative description of intended behavior.
- David Kim is a co-author and active package contributor. Do not assume
  David's commits are wrong without evidence. Ask before reverting.
- The paper's experiment scripts in
  `combined_stephenson_tests/code/` (excluding the embedded
  `codes_20251026.R` and `old_codes/`) are what produced the published
  numbers. Run paths through these, not through `pval_comb_block` direct.
- The package supports both HiGHS (open-source, default) and Gurobi
  (commercial). HiGHS is what `R CMD check` will use unless Gurobi is
  installed locally. Tests use `skip_if_not(solver_available(...))`.
- Plain ASCII only in any file. No unicode em dashes, arrows, fancy
  quotes, ellipses, etc. (Jake's global rule.)
- Memory: feedback memory at
  `/Users/jwbowers/.claude/projects/-Users-jwbowers-repos-CMRSS-jake/memory/feedback_coding_workflow.md`
  encodes the workflow rules from `~/repos/ai_workflow/CLAUDE_CODING.md`.
  Re-read both before doing implementation work.

## Suggested next session opening

If David has replied:
- Read his reply, decide which scenario applies (regression vs deliberate).
- Drop the `skip()` calls in
  `tests/testthat/test-pval-comb-block-p-convention.R`.
- Run the tests; they should fail.
- Apply the fix at `R/CMRSS_SRE.R:1034` (and any docstring/example/wrapper
  updates if scenario 2).
- Bump DESCRIPTION patch (0.2.5 -> 0.2.6).
- Run `devtools::document()` and `devtools::check()`.
- Pause for Jake to review (CHECKPOINT 1A.i).

If David has not replied:
- Move to item 1B (independent of 1A, same area of code). Per `PLAN.md`,
  start with the brute-force enumeration helper for `comb_matrix_block_stratum`.
- Or move to item 2A (tie handling, fully independent).
