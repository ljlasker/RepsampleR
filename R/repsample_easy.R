#' Simplified Front-End for Common RepsampleR Runs
#'
#' A convenience wrapper for quickly targeting mean/SD (and optional
#' distribution form) with simple quality presets, without manually setting
#' `randomperc`, `rrule`, and `srule`.
#'
#' @param data Data frame.
#' @param size Desired sample size.
#' @param cont Continuous variable names.
#' @param mean Target means for `cont`.
#' @param sd Target standard deviations for `cont`.
#' @param dist Target distribution form(s): `"normal"` (default),
#'   `"lognormal"`, `"poisson"`.
#' @param quality Tuning preset passed to `repsample()`: `"fast"`,
#'   `"balanced"` (default), `"strict"`, or `"manual"`.
#' @param mode `"single"` runs one seed; `"search"` runs multi-seed outer
#'   search with `repsample_search()`.
#' @param seed Seed for `"single"` mode.
#' @param n_seeds Number of seeds in `"search"` mode.
#' @param seed_start Starting seed in `"search"` mode.
#' @param n_outer_workers Outer workers in `"search"` mode.
#' @param n_cores Inner per-iteration CPU cores (passed to `repsample()`).
#' @param backend Compute backend (`"cpu"` or `"cuda"`).
#' @param method Sampling method:
#'   `"auto"` (default), `"greedy"`, `"importance"`, or `"nearest"`.
#' @param nearest_replace Logical flag used when `method = "nearest"`:
#'   `FALSE` (default) matches without replacement; `TRUE` allows replacement
#'   (duplicate draws are recorded in `data$repsample_n`).
#' @param nearest_distance Distance metric for `method = "nearest"`:
#'   `"euclidean"` (default), `"mahalanobis"`, or `"weighted"`.
#' @param nearest_backend Nearest-neighbor search backend:
#'   `"auto"` (default), `"exact"`, or `"hnsw"` (requires optional
#'   package `RcppHNSW`).
#' @param nearest_match Nearest-neighbor assignment mode:
#'   `"greedy"` (default) or `"optimal"` (global assignment; requires
#'   optional package `clue`).
#' @param nearest_feature_weights Optional positive weights for
#'   `nearest_distance = "weighted"` (length 1 or length `cont`).
#' @param nearest_caliper Optional positive caliper(s) on standardized
#'   per-variable absolute distances for `method = "nearest"`.
#' @param nearest_caliper_strict If `TRUE`, error when no candidate satisfies
#'   the caliper for a draw; otherwise fallback to nearest candidate.
#' @param nearest_optimal_max_size Maximum allowed `size` for
#'   `nearest_match = "optimal"` before automatically falling back to greedy.
#' @param nearest_hnsw_k Candidate count queried per draw when using HNSW.
#' @param nearest_hnsw_M HNSW graph degree parameter.
#' @param nearest_hnsw_ef_construction HNSW build-time search width.
#' @param nearest_hnsw_ef_search HNSW query-time search width.
#' @param mean_tol Optional tolerance on absolute mean error for early stopping
#'   in `mode = "search"`.
#' @param sd_tol Optional tolerance on absolute SD error for early stopping in
#'   `mode = "search"`.
#' @param ks_tol Optional tolerance on KS statistic for early stopping in
#'   `mode = "search"`.
#' @param perc_tol Optional tolerance on absolute percentage-point error for
#'   binary targets in `mode = "search"`.
#' @param cuda_fallback If `TRUE` (default), automatically fall back to CPU
#'   when CUDA is unavailable or a CUDA fast-path cannot be used.
#' @param ... Additional arguments forwarded to `repsample()` (for `"single"`)
#'   or to `repsample_search()` (for `"search"`), such as `rrule`, `srule`,
#'   `randomperc`, `objective`, etc.
#'
#' @return A `repsample_result` (`mode = "single"`) or
#'   `repsample_search_result` (`mode = "search"`).
#'
#' @export
repsample_easy <- function(data,
                           size,
                           cont,
                           mean,
                           sd,
                           dist = "normal",
                           quality = c("balanced", "fast", "strict", "manual"),
                           mode = c("single", "search"),
                           seed = 7,
                           n_seeds = 16,
                           seed_start = 1,
                           n_outer_workers = 1,
                           n_cores = 1,
                           backend = c("cpu", "cuda"),
                           method = c("auto", "greedy", "importance", "nearest"),
                           nearest_replace = FALSE,
                           nearest_distance = c("euclidean", "mahalanobis", "weighted"),
                           nearest_backend = c("auto", "exact", "hnsw"),
                           nearest_match = c("greedy", "optimal"),
                           nearest_feature_weights = NULL,
                           nearest_caliper = NULL,
                           nearest_caliper_strict = FALSE,
                           nearest_optimal_max_size = 2000L,
                           nearest_hnsw_k = 64L,
                           nearest_hnsw_M = 16L,
                           nearest_hnsw_ef_construction = 200L,
                           nearest_hnsw_ef_search = 64L,
                           mean_tol = NULL,
                           sd_tol = NULL,
                           ks_tol = NULL,
                           perc_tol = NULL,
                           cuda_fallback = TRUE,
                           ...) {
  quality <- match.arg(quality)
  mode <- match.arg(mode)
  backend <- match.arg(backend)
  method <- match.arg(method)
  nearest_distance <- match.arg(nearest_distance)
  nearest_backend <- match.arg(nearest_backend)
  nearest_match <- match.arg(nearest_match)
  check_scalar_flag(nearest_replace, "nearest_replace")
  check_scalar_flag(nearest_caliper_strict, "nearest_caliper_strict")
  check_scalar_flag(cuda_fallback, "cuda_fallback")

  extra_args <- list(...)
  exact_arg <- if ("exact" %in% names(extra_args)) isTRUE(extra_args$exact) else FALSE
  bincat_arg <- if ("bincat" %in% names(extra_args)) extra_args$bincat else NULL
  unsupported_is_args <- c("retain", "subset", "in_rows", "wght", "perc")
  has_unsupported_is <- any(unsupported_is_args %in% names(extra_args))

  can_use_is <- .use_importance_sampling(
    cont = cont,
    bincat = bincat_arg,
    mean = mean,
    sd = sd,
    dist = dist,
    exact = exact_arg
  ) && !has_unsupported_is
  can_use_nn <- .use_nearest_neighbor_sampling(
    cont = cont,
    bincat = bincat_arg,
    mean = mean,
    sd = sd,
    dist = dist,
    exact = exact_arg
  ) && !has_unsupported_is

  if (method == "importance" && !can_use_is) {
    stop(
      "`method = \"importance\"` cannot be used with the current arguments. ",
      "Use `method = \"greedy\"` or remove unsupported options.",
      call. = FALSE
    )
  }
  if (method == "nearest" && !can_use_nn) {
    stop(
      "`method = \"nearest\"` cannot be used with the current arguments. ",
      "Use `method = \"greedy\"` or remove unsupported options.",
      call. = FALSE
    )
  }

  base_args <- c(
    list(
      data = data,
      size = size,
      cont = cont,
      mean = mean,
      sd = sd,
      dist = dist,
      quality = quality,
      n_cores = n_cores,
      backend = backend,
      cuda_fallback = cuda_fallback
    ),
    extra_args
  )

  if (mode == "single") {
    if (method == "nearest" && can_use_nn) {
      if (backend != "cpu") {
        warning(
          "`backend` is ignored when nearest-neighbor sampling is selected.",
          call. = FALSE
        )
      }
      return(.nearest_sample_one(
        data = data,
        size = size,
        cont = cont,
        mean = mean,
        sd = sd,
        dist = dist,
        seed = seed,
        replace = nearest_replace,
        distance = nearest_distance,
        backend = nearest_backend,
        match = nearest_match,
        feature_weights = nearest_feature_weights,
        caliper = nearest_caliper,
        caliper_strict = nearest_caliper_strict,
        optimal_max_size = nearest_optimal_max_size,
        hnsw_k = nearest_hnsw_k,
        hnsw_M = nearest_hnsw_M,
        hnsw_ef_construction = nearest_hnsw_ef_construction,
        hnsw_ef_search = nearest_hnsw_ef_search
      ))
    }

    if (method != "greedy" && method != "nearest" && can_use_is) {
      if (backend != "cpu") {
        warning(
          "`backend` is ignored when importance sampling is selected.",
          call. = FALSE
        )
      }

      out_is <- .importance_sample_one(
        data = data,
        size = size,
        cont = cont,
        mean = mean,
        sd = sd,
        dist = dist,
        seed = seed,
        fallback_to_greedy = (method == "auto")
      )
      if (!is.null(out_is)) {
        return(out_is)
      }
    }

    return(do.call(repsample, c(base_args, list(seednum = seed))))
  }

  search_args <- c(
    base_args,
    list(
      n_seeds = n_seeds,
      seed_start = seed_start,
      n_outer_workers = n_outer_workers,
      method = method,
      nearest_replace = nearest_replace,
      nearest_distance = nearest_distance,
      nearest_backend = nearest_backend,
      nearest_match = nearest_match,
      nearest_feature_weights = nearest_feature_weights,
      nearest_caliper = nearest_caliper,
      nearest_caliper_strict = nearest_caliper_strict,
      nearest_optimal_max_size = nearest_optimal_max_size,
      nearest_hnsw_k = nearest_hnsw_k,
      nearest_hnsw_M = nearest_hnsw_M,
      nearest_hnsw_ef_construction = nearest_hnsw_ef_construction,
      nearest_hnsw_ef_search = nearest_hnsw_ef_search,
      mean_tol = mean_tol,
      sd_tol = sd_tol,
      ks_tol = ks_tol,
      perc_tol = perc_tol
    )
  )
  do.call(repsample_search, search_args)
}
