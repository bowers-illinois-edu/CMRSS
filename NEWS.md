# CMRSS 0.2.5

## Performance

- Avoid large `n x nperm` permutation matrices by default:
  - CRE null generation (`null_dist()`, `null_dist_multiple()`) now streams permutations in chunks when `Z.perm` is not supplied (new `chunk_size` argument).
  - SRE null generation (`com_null_dist_block()`, `com_null_dist_block_stratum()`) now generates permutations in chunks when `Z.perm` is not supplied (new `chunk_size` argument).
- Speed up stratified runs with many blocks:
  - `summary_block()` now builds `units.block` via `split()` (reduces `O(n·B)` scanning to `O(n)`).
  - Block-randomized chunk generation uses fast paths for common small-block cases (e.g., `mb == 1`, `mb == nb - 1`) and an internal cache for small `combn(nb, mb)` patterns.
- Speed up `comb.method = 2` null distribution:
  - `com_null_dist_block_stratum()` now computes per-block statistics with matrix operations instead of per-permutation nested loops.
- Reduce allocation/copying overhead in solver setup:
  - Solver constraint triplets are now assembled via list accumulation rather than repeated `c()` growth.
- Minor speed improvement in generalized CI code:
  - `ci_lower_quantile_generalize()` no longer grows results with repeated `rbind()` inside a loop.

## Bug Fixes

- `rank_score()` consistently returns a numeric vector (avoids 1-column matrix output when `scale = TRUE`).

## Testing

- Added unit tests for internal block-permutation chunk generation and for equivalence of the vectorized `com_null_dist_block_stratum()` implementation to a naive reference computation.
