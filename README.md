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
- optional method dispatch: `method = "auto" | "greedy" | "importance"`
- nearest-neighbor theoretical matching path with multivariate continuous support:
  `method = "nearest"` with `nearest_replace = FALSE/TRUE`

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
- CUDA mode requires Python `cupy` and R package `reticulate` on the target machine.
- On Windows with CuPy 12, install `nvidia-cuda-nvrtc-cu12`,
  `nvidia-cuda-runtime-cu12`, and `nvidia-cublas-cu12` in the same Python
  environment. Python 3.8+ may also require registering CUDA DLL paths via
  `os.add_dll_directory()` when launching from R.

# Citations

Kontopantelis, E. (2013). A Greedy Algorithm for Representative Sampling: repsample in Stata. *Journal of Statistical Software*, Code Snippets, 55(1), 1–19. https://doi.org/10.18637/jss.v055.c01
