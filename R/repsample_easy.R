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
                           ...) {
  quality <- match.arg(quality)
  mode <- match.arg(mode)
  backend <- match.arg(backend)
  method <- match.arg(method)

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
      backend = backend
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
        seed = seed
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
      method = method
    )
  )
  do.call(repsample_search, search_args)
}
