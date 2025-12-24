# Repository Guidelines

## Project Structure & Module Organization

- `R/`: Package source code (exported and internal functions).
- `man/`: Generated `.Rd` documentation (from roxygen2 comments in
  `R/`).
- `tests/testthat/`: Unit tests (testthat).
- `tests/testthat.R`: testthat entrypoint used by `R CMD check`.
- `data/`: Package data assets (if present/added later).
- `_pkgdown.yml`: pkgdown site configuration.

## Build, Test, and Development Commands

Use the Makefile targets (preferred) or run the underlying R commands
directly:

- `make dependencies`: Install package dependencies via devtools.
- `make test`: Run testthat tests.
- `make check`: Run package checks (similar to `R CMD check`).
- `make document`: Regenerate docs (`NAMESPACE` + `man/`) from roxygen
  comments.
- `make build`: Build the package tarball.
- `make interactive`: Start an interactive R session with `load.R`
  profile tweaks.

## Coding Style & Naming Conventions

- Language: R. Prefer 2-space indentation, `<-` assignment, and explicit
  argument names in public APIs.
- Naming: `snake_case` for functions/variables (e.g.,
  [`summary_block()`](https://bowers-illinois-edu.github.io/CMRSS/reference/summary_block.md),
  [`assign_CRE()`](https://bowers-illinois-edu.github.io/CMRSS/reference/assign_CRE.md)).
- Documentation: Use roxygen2 blocks (`#'`) above functions; do not
  hand-edit `NAMESPACE`.
- Linting: `lintr` is configured via `.lintr` (run
  `lintr::lint_package()` when making non-trivial changes).

## Testing Guidelines

- Framework: testthat v3 (`tests/testthat/`).
- Conventions: keep tests in `tests/testthat/test-*.R` and cover bug
  fixes with a regression test.
- Run locally: `make test` (fast) and `make check` before opening a PR.

## Commit & Pull Request Guidelines

- Commits: history favors short, imperative messages (e.g., “Add lintr”,
  “Update description”, “Fix pkgdown workflow”).
- PRs: include a clear summary, link any relevant issue, and note
  user-facing changes (API/docs). Ensure `make check` passes and update
  docs/tests when behavior changes.

## Dependencies & Solver Notes

- Optional solvers: `highs` and `gurobi` are suggested; code should
  degrade gracefully when they’re not installed.
- Avoid adding hard dependencies unless required for core functionality.
