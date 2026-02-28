# RepsampleR

`RepsampleR` is both:
- a port of Stata's `repsample` (Kontopantelis, 2013), and
- an expansion of it with additional performance and workflow features in R.

## Ported Functionality

- Population sampling mode (`cont`, `bincat`)
- Theoretical sampling mode (`mean`, `sd`, `perc`)
- `retain` (batched/continued sampling)
- `force` overwrite behavior for existing `repsample` column
- `randomperc` initial random fraction
- `srule` early stop by combined p-value
- `rrule` early stop by random candidate subsampling
- `wght` variable weights (must sum to 100; cont first, then bincat)
- `exact` behavior matching Stata scope
- Stata-like returned scalar bundle (`$r`)

## RepsampleR Extensions

- `n_cores` intra-iteration parallel candidate scoring (Unix-like systems)
- optional CUDA GPU fast path (`backend = "cuda"`) for supported theoretical continuous targets
- outer seed-search wrappers: `repsample_search()` and `repsample_search_auto()`
- simplified front-end: `repsample_easy()`
- one-command high-level API: `repsample_fit()`
- optional method dispatch: `method = "auto" | "greedy" | "importance"`
- nearest-neighbor theoretical matching path with multivariate continuous support:
  `method = "nearest"` with `nearest_replace = FALSE/TRUE`
- nearest-neighbor metric controls: `nearest_distance = "euclidean" | "mahalanobis" | "weighted"`
- nearest-neighbor backend controls: `nearest_backend = "auto" | "exact" | "hnsw"` (HNSW requires optional `RcppHNSW`)
- nearest-neighbor assignment controls: `nearest_match = "greedy" | "optimal"` (optimal requires optional `clue`)
- nearest calipers and optional strict caliper enforcement
- tolerance-based early stop (`mean_tol`, `sd_tol`, `ks_tol`, `perc_tol`)
- adaptive method search (`method = "adaptive"`) with pilot-then-exploit budget allocation
- staged search checkpoint/resume (`checkpoint_file`, `resume`) in `repsample_search_auto()`
- diagnostics helpers: `repsample_quality()`, `summary()`, and `plot()`

## Install from GitHub (public)

```r
install.packages("remotes")
remotes::install_github("ljlasker/RepsampleR")
```

## Usage

```r
library(RepsampleR)

out <- repsample(
  data = dat,
  size = 40,
  cont = c("age", "bmi"),
  bincat = c("sex", "site"),
  seednum = 10,
  randomperc = 10,
  srule = 0.8,
  rrule = 30,
  wght = c(35, 25, 20, 20),
  n_cores = max(1L, parallel::detectCores(logical = FALSE) - 1L),
  backend = "cpu"
)

# Selected sample indicator (0/1)
head(out$data$repsample)

# Stata-like r() scalars
out$r
```

## Simple Defaults + Distribution Targets

```r
library(RepsampleR)

# Single run with preset tuning and target distribution form
out <- repsample_easy(
  data = baseline,
  size = 1000,
  cont = "x",
  mean = 1,
  sd = 1,
  dist = "normal",      # normal | lognormal | poisson
  quality = "balanced", # fast | balanced | strict | manual
  method = "auto",      # auto | greedy | importance | nearest
  nearest_replace = FALSE,
  mode = "single",
  seed = 7
)
```

## One-Command API (`repsample_fit`)

```r
library(RepsampleR)

fit <- repsample_fit(
  data = baseline,
  size = 1000,
  cont = c("x1", "x2"),
  mean = c(1, -0.5),
  sd = c(1, 0.8),
  dist = c("normal", "normal"),
  method = "adaptive",   # auto | adaptive | greedy | importance | nearest
  mode = "search",
  control = list(
    n_seeds = 24,
    n_outer_workers = 4,
    nearest_distance = "mahalanobis",
    nearest_backend = "auto",
    nearest_match = "greedy",
    mean_tol = 0.05,
    sd_tol = 0.05,
    ks_tol = 0.03
  )
)

summary(fit$best)
plot(fit$best, var = "x1")
```

## Outer Search Preset

```r
library(RepsampleR)

search <- repsample_search(
  data = baseline,
  size = nrow(baseline) %/% 10,
  cont = "x",
  mean = 1,
  sd = 1,
  method = "auto",      # auto | greedy | importance | nearest
  nearest_replace = FALSE,
  n_seeds = 24,
  seed_start = 1001,
  n_outer_workers = max(1L, parallel::detectCores(logical = FALSE) - 1L),
  backend = "cpu",
  randomperc = 10,
  rrule = 800,
  srule = 0.85
)

best_fit <- search$best
best_fit$meta
```

## Multi-Stage Search with Resume

```r
library(RepsampleR)

search <- repsample_search_auto(
  data = baseline,
  size = nrow(baseline) %/% 10,
  cont = "x",
  mean = 1,
  sd = 1,
  n_stages = 3,
  n_seeds = c(24, 24, 24),
  checkpoint_file = "repsample_stage_checkpoint.rds",
  resume = TRUE
)
```

## Nearest Advanced Controls

```r
library(RepsampleR)

fit_nn <- repsample_fit(
  data = baseline,
  size = 1000,
  cont = c("x1", "x2"),
  mean = c(1, -0.5),
  sd = c(1, 0.8),
  method = "nearest",
  mode = "search",
  control = list(
    n_seeds = 12,
    nearest_replace = FALSE,
    nearest_distance = "weighted",
    nearest_feature_weights = c(2, 0.5),
    nearest_backend = "hnsw",       # auto fallback to exact if unavailable
    nearest_match = "optimal",      # auto fallback to greedy if `clue` unavailable
    nearest_caliper = c(3, 3),
    nearest_caliper_strict = FALSE
  )
)
```

## One-Command Multi-Stage Search

```r
library(RepsampleR)

search <- repsample_search_auto(
  data = baseline,
  size = nrow(baseline) %/% 10,
  cont = "x",
  mean = 1,
  sd = 1,
  n_stages = 3,
  n_seeds = c(24, 24, 24),
  top_k = 4,
  refine_radius = 128,
  stop_loss = 0.02,
  n_outer_workers = max(1L, parallel::detectCores(logical = FALSE) - 1L),
  outer_parallel = "auto",
  backend = "cpu",
  randomperc = 10,
  rrule = 800,
  srule = 0.85
)

best_fit <- search$best
search$stage_summary
```

For CUDA, pass `backend = "cuda"`. The outer search will use one worker by
default to avoid oversubscribing a single GPU.

For theoretical continuous targets (`cont` only, no `bincat`, `exact = FALSE`),
`method = "auto"` can dispatch to an importance-sampling fast path.

When `method = "nearest"` and `nearest_replace = TRUE`, duplicate draws are
tracked in `out$data$repsample_n` (counts per row), and draw order is preserved
in `out$selected_rows`.

## Speed notes

- Default mode uses C++ kernels for the KS-statistic bottleneck.
- Set `n_cores > 1` on macOS/Linux to parallelize candidate scoring within each greedy iteration.
- In RStudio on macOS, `RepsampleR` uses a PSOCK backend instead of forking,
  so `n_cores > 1` can still be used from the IDE.
- PSOCK parallelism has overhead; it is only engaged for sufficiently large
  candidate pools where it is likely to help.
- `srule` and `rrule` can substantially reduce runtime on large pools.
- `exact = TRUE` is intentionally slower (as in Stata).
- CUDA mode currently targets theoretical sampling with continuous variables
  and `exact = FALSE`; other models fall back to CPU scoring.
- By default, `backend = "cuda"` now falls back to CPU automatically when
  CUDA dependencies are unavailable (`cuda_fallback = TRUE`).
- CUDA mode requires Python `cupy` and R package `reticulate` on the target machine.
- Optional HNSW nearest backend requires R package `RcppHNSW`.
- Optional global optimal nearest assignment requires R package `clue`.
- On Windows with CuPy 12, install `nvidia-cuda-nvrtc-cu12`,
  `nvidia-cuda-runtime-cu12`, and `nvidia-cublas-cu12` in the same Python
  environment. Python 3.8+ may also require registering CUDA DLL paths via
  `os.add_dll_directory()` when launching from R.

## Benchmark Harness

Run local benchmark comparisons with:

```bash
Rscript scripts/benchmark_repsample_methods.R --n=10000 --size=1000 --n_seeds=8 --outfile=benchmark.csv
```

## CI and Commit Checks

This repository now includes GitHub Actions checks:
- `R-CMD-check`: cross-platform package checks (Ubuntu, macOS, Windows, plus R-devel on Ubuntu).
- `lint`: static linting with `lintr`.
- `coverage`: test coverage report (Cobertura + Codecov upload).
- `benchmark-smoke`: small benchmark run on pull requests.

# Citations

Kontopantelis, E. (2013). A Greedy Algorithm for Representative Sampling: repsample in Stata. *Journal of Statistical Software*, Code Snippets, 55(1), 1–19. https://doi.org/10.18637/jss.v055.c01
