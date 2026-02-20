default_repsample_loss <- function(out, cont, bincat, target_mean, target_sd, target_perc, theoretical) {
  chosen <- out$data[out$data$repsample == 1L, , drop = FALSE]
  if (nrow(chosen) == 0L) {
    return(Inf)
  }

  if (!theoretical) {
    if ("chi2" %in% names(out$r)) {
      return(as.numeric(out$r[["chi2"]]))
    }
    if ("p" %in% names(out$r)) {
      return(-as.numeric(out$r[["p"]]))
    }
    return(Inf)
  }

  loss <- 0
  if (length(cont) > 0L) {
    for (i in seq_along(cont)) {
      x <- chosen[[cont[i]]]
      dm <- base::mean(x) - target_mean[i]
      dsd <- stats::sd(x) - target_sd[i]
      loss <- loss + dm * dm + dsd * dsd
    }
  }

  if (length(bincat) > 0L) {
    for (i in seq_along(bincat)) {
      p_obs <- base::mean(chosen[[bincat[i]]] == 1)
      p_tar <- target_perc[i] / 100
      dp <- p_obs - p_tar
      loss <- loss + dp * dp
    }
  }

  loss
}

check_positive_int_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      x < 1 || abs(x - round(x)) > 0) {
    stop(sprintf("`%s` must be a positive integer.", name), call. = FALSE)
  }
  as.integer(x)
}

sanitize_seed_vector <- function(seeds, name = "`seeds`") {
  if (!is.numeric(seeds) || anyNA(seeds) || length(seeds) < 1L) {
    stop(sprintf("%s must be a non-empty numeric vector.", name), call. = FALSE)
  }
  if (any(abs(seeds - round(seeds)) > 0)) {
    stop(sprintf("%s must contain integer values.", name), call. = FALSE)
  }
  max_seed <- .Machine$integer.max
  if (any(seeds < 0 | seeds > max_seed)) {
    stop(sprintf("%s must be in [0, %d].", name, max_seed), call. = FALSE)
  }
  as.integer(unique(seeds))
}

normalize_stage_seed_counts <- function(n_stages, n_seeds) {
  if (!is.numeric(n_seeds) || anyNA(n_seeds) || length(n_seeds) < 1L) {
    stop("`n_seeds` must be a positive integer or a non-empty integer vector.", call. = FALSE)
  }
  if (any(n_seeds < 1 | abs(n_seeds - round(n_seeds)) > 0)) {
    stop("`n_seeds` values must be positive integers.", call. = FALSE)
  }
  n_seeds <- as.integer(n_seeds)

  if (length(n_seeds) == 1L) {
    return(rep.int(n_seeds, n_stages))
  }
  if (length(n_seeds) == n_stages) {
    return(n_seeds)
  }
  stop("When `n_seeds` has length > 1, it must have exactly one value per stage.", call. = FALSE)
}

normalize_refine_radius <- function(n_stages, stage_n_seeds, refine_radius, radius_shrink) {
  if (n_stages <= 1L) {
    return(integer(0))
  }

  if (is.null(refine_radius)) {
    refine_radius <- max(8L, as.integer(ceiling(stage_n_seeds[1L] / 2)))
  }

  if (!is.numeric(radius_shrink) || length(radius_shrink) != 1L || is.na(radius_shrink) ||
      radius_shrink <= 0 || radius_shrink > 1) {
    stop("`radius_shrink` must be in (0, 1].", call. = FALSE)
  }

  if (!is.numeric(refine_radius) || anyNA(refine_radius) || length(refine_radius) < 1L) {
    stop("`refine_radius` must be NULL, a positive integer, or an integer vector.", call. = FALSE)
  }
  if (any(refine_radius < 1 | abs(refine_radius - round(refine_radius)) > 0)) {
    stop("`refine_radius` values must be positive integers.", call. = FALSE)
  }
  refine_radius <- as.integer(refine_radius)

  if (length(refine_radius) == 1L) {
    out <- as.integer(round(refine_radius[1L] * (radius_shrink ^ (seq_len(n_stages - 1L) - 1L))))
    out[out < 1L] <- 1L
    return(out)
  }
  if (length(refine_radius) == (n_stages - 1L)) {
    return(refine_radius)
  }
  stop(
    "When `refine_radius` has length > 1, it must have one value for each refinement stage.",
    call. = FALSE
  )
}

fill_seed_budget <- function(center_seed, seed_pool, target_n, seen_seeds, allow_repeats) {
  out <- as.integer(unique(seed_pool))
  max_seed <- .Machine$integer.max
  offset <- 0L

  while (length(out) < target_n) {
    candidates <- if (offset == 0L) {
      center_seed
    } else {
      c(center_seed - offset, center_seed + offset)
    }

    for (cand in candidates) {
      if (cand < 0L || cand > max_seed) {
        next
      }
      if (!allow_repeats && (cand %in% seen_seeds)) {
        next
      }
      if (!(cand %in% out)) {
        out <- c(out, as.integer(cand))
      }
      if (length(out) >= target_n) {
        break
      }
    }

    offset <- offset + 1L
    if (offset > max_seed) {
      break
    }
  }

  if (length(out) < target_n) {
    stop("Unable to generate enough valid seeds for a refinement stage.", call. = FALSE)
  }

  out
}

build_refine_seed_set <- function(top_seeds,
                                  center_seed,
                                  radius,
                                  target_n,
                                  seen_seeds,
                                  allow_repeats) {
  out <- integer(0)
  max_seed <- .Machine$integer.max

  for (offset in seq.int(0L, radius)) {
    for (seed in top_seeds) {
      candidates <- if (offset == 0L) seed else c(seed - offset, seed + offset)
      for (cand in candidates) {
        if (cand < 0L || cand > max_seed) {
          next
        }
        if (!allow_repeats && (cand %in% seen_seeds)) {
          next
        }
        if (!(cand %in% out)) {
          out <- c(out, as.integer(cand))
        }
      }
    }
    if (length(out) >= target_n) {
      break
    }
  }

  if (length(out) > target_n) {
    ord <- order(abs(out - center_seed), out)
    out <- out[ord][seq_len(target_n)]
  } else if (length(out) < target_n) {
    out <- fill_seed_budget(center_seed, out, target_n, seen_seeds, allow_repeats)
  }

  as.integer(out)
}

#' Outer Seed Search Preset for `repsample()`
#'
#' Runs multiple `repsample()` calls over different seeds and returns the
#' best fit under a scalar loss function. This provides a concise preset for
#' outer parallelization.
#'
#' @param data Data frame.
#' @param size Desired sample size.
#' @param cont Continuous variable names.
#' @param bincat Binary/categorical variable names.
#' @param mean Means for continuous variables (theoretical sampling only).
#' @param sd Standard deviations for continuous variables (theoretical sampling only).
#' @param perc Percentages for binary variables (theoretical sampling only).
#' @param seeds Optional explicit seed vector.
#' @param n_seeds Number of seeds to search when `seeds` is not provided.
#' @param seed_start Starting seed when generating seeds.
#' @param n_outer_workers Number of workers for outer seed search.
#' @param outer_parallel Outer parallel backend: `"auto"`, `"serial"`,
#'   `"multicore"`, or `"psock"`.
#' @param method Sampling method:
#'   `"auto"` (default), `"greedy"`, `"importance"`, or `"nearest"`.
#' @param nearest_replace Logical flag used when `method = "nearest"`:
#'   `FALSE` (default) matches without replacement; `TRUE` allows replacement
#'   (duplicate draws are recorded in `best$data$repsample_n`).
#' @param objective Optional custom objective function that takes a
#'   `repsample_result` and returns a scalar loss (smaller is better).
#' @param keep_all If `TRUE`, return all per-seed fit objects.
#' @param ... Additional arguments passed to `repsample()`.
#'
#' @return Object of class `repsample_search_result` with:
#' - `best`: best `repsample_result`
#' - `best_seed`: seed with smallest loss
#' - `best_loss`: smallest loss
#' - `summary`: data frame of all searched seeds and losses
#' - `all`: optional full per-seed results when `keep_all = TRUE`
#' - `meta`: search metadata
#'
#' @export
repsample_search <- function(data,
                             size,
                             cont = NULL,
                             bincat = NULL,
                             mean = NULL,
                             sd = NULL,
                             perc = NULL,
                             seeds = NULL,
                             n_seeds = 16,
                             seed_start = 1,
                             n_outer_workers = 1,
                             outer_parallel = c("auto", "serial", "multicore", "psock"),
                             method = c("auto", "greedy", "importance", "nearest"),
                             nearest_replace = FALSE,
                             objective = NULL,
                             keep_all = FALSE,
                             ...) {
  outer_parallel <- match.arg(outer_parallel)
  method <- match.arg(method)
  check_scalar_flag(nearest_replace, "nearest_replace")

  n_seeds <- check_positive_int_scalar(n_seeds, "n_seeds")

  if (!is.numeric(seed_start) || length(seed_start) != 1L || is.na(seed_start)) {
    stop("`seed_start` must be a single numeric value.", call. = FALSE)
  }

  if (is.null(seeds)) {
    seeds <- sanitize_seed_vector(
      as.integer(seed_start) + seq_len(n_seeds) - 1L,
      name = "Generated `seeds`"
    )
  } else {
    seeds <- sanitize_seed_vector(seeds, name = "`seeds`")
  }

  n_outer_workers <- check_positive_int_scalar(n_outer_workers, "n_outer_workers")

  if (!is.null(objective) && !is.function(objective)) {
    stop("`objective` must be NULL or a function.", call. = FALSE)
  }
  if (!is.logical(keep_all) || length(keep_all) != 1L || is.na(keep_all)) {
    stop("`keep_all` must be TRUE or FALSE.", call. = FALSE)
  }

  extra_args <- list(...)
  if ("seednum" %in% names(extra_args)) {
    warning("Ignoring `seednum` in `...`; seeds are controlled by `seeds` / `n_seeds`.", call. = FALSE)
    extra_args$seednum <- NULL
  }

  exact_arg <- if ("exact" %in% names(extra_args)) isTRUE(extra_args$exact) else FALSE
  dist_arg <- if ("dist" %in% names(extra_args)) extra_args$dist else NULL
  unsupported_is_args <- c("retain", "subset", "in_rows", "wght", "perc")
  has_unsupported_is <- any(unsupported_is_args %in% names(extra_args))

  can_use_is <- .use_importance_sampling(
    cont = cont,
    bincat = bincat,
    mean = mean,
    sd = sd,
    dist = dist_arg,
    exact = exact_arg
  ) && !has_unsupported_is
  can_use_nn <- .use_nearest_neighbor_sampling(
    cont = cont,
    bincat = bincat,
    mean = mean,
    sd = sd,
    dist = dist_arg,
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

  if ((method %in% c("importance", "nearest")) &&
      !is.null(extra_args$backend) && extra_args$backend != "cpu") {
    label <- if (method == "nearest") "nearest-neighbor sampling" else "importance sampling"
    warning(
      sprintf("`backend` is ignored when %s is selected.", label),
      call. = FALSE
    )
  }

  backend_choice <- extra_args$backend
  if (is.null(backend_choice)) {
    backend_choice <- "cpu"
  }
  backend_choice <- match.arg(backend_choice, c("cpu", "cuda"))

  mode <- switch(
    outer_parallel,
    serial = "serial",
    multicore = "multicore",
    psock = "psock",
    auto = {
      if (n_outer_workers <= 1L || length(seeds) <= 1L) {
        "serial"
      } else if (.Platform$OS.type == "windows") {
        "psock"
      } else if (nzchar(Sys.getenv("RSTUDIO")) && identical(Sys.info()[["sysname"]], "Darwin")) {
        "psock"
      } else {
        "multicore"
      }
    }
  )

  if (mode != "serial" && n_outer_workers <= 1L) {
    mode <- "serial"
  }
  if (mode == "multicore" && .Platform$OS.type == "windows") {
    warning("Multicore outer mode is not available on Windows; using PSOCK.", call. = FALSE)
    mode <- "psock"
  }

  if (method == "nearest" && can_use_nn) {
    out_nn <- .nearest_sample_search(
      data = data,
      size = size,
      cont = cont,
      mean = mean,
      sd = sd,
      dist = dist_arg,
      seeds = seeds,
      n_outer_workers = n_outer_workers,
      outer_mode = mode,
      keep_all = keep_all,
      objective = objective,
      replace = nearest_replace
    )
    if (!is.null(out_nn)) {
      return(out_nn)
    }
  }

  if (method != "greedy" && method != "nearest" && can_use_is) {
    out_is <- .importance_sample_search(
      data = data,
      size = size,
      cont = cont,
      mean = mean,
      sd = sd,
      dist = dist_arg,
      seeds = seeds,
      n_outer_workers = n_outer_workers,
      outer_mode = mode,
      keep_all = keep_all,
      objective = objective,
      fallback_to_greedy = (method == "auto")
    )
    if (!is.null(out_is)) {
      return(out_is)
    }
  }

  if (backend_choice == "cuda" && n_outer_workers > 1L) {
    warning(
      "Using `backend = \"cuda\"` with multiple outer workers can oversubscribe one GPU; forcing one outer worker.",
      call. = FALSE
    )
    n_outer_workers <- 1L
    if (mode != "serial") {
      mode <- "serial"
    }
  }

  inner_cores <- extra_args$n_cores
  if (!is.null(inner_cores) &&
      backend_choice == "cpu" &&
      mode != "serial" &&
      n_outer_workers > 1L) {
    inner_cores_num <- suppressWarnings(as.numeric(inner_cores[1L]))
    if (is.finite(inner_cores_num) &&
        inner_cores_num >= 1 &&
        abs(inner_cores_num - round(inner_cores_num)) < 1e-9 &&
        as.integer(round(inner_cores_num)) > 1L) {
      warning(
        paste0(
          "Nested parallelism detected (`n_outer_workers > 1` and `n_cores > 1`). ",
          "Forcing `n_cores = 1` inside each outer worker to avoid CPU oversubscription. ",
          "Use `n_outer_workers = 1` if you want intra-iteration parallelism instead."
        ),
        call. = FALSE
      )
      extra_args$n_cores <- 1L
    }
  }

  theoretical <- (length(normalize_numlist(mean)) > 0L ||
                  length(normalize_numlist(sd)) > 0L ||
                  length(normalize_numlist(perc)) > 0L)
  cont_norm <- normalize_varlist(cont)
  bincat_norm <- normalize_varlist(bincat)
  mean_norm <- normalize_numlist(mean)
  sd_norm <- normalize_numlist(sd)
  perc_norm <- normalize_numlist(perc)

  eval_seed <- function(seed) {
    args <- c(
      list(
        data = data,
        size = size,
        cont = cont,
        bincat = bincat,
        mean = mean,
        sd = sd,
        perc = perc,
        seednum = seed
      ),
      extra_args
    )
    out <- do.call(repsample, args)
    loss <- if (is.null(objective)) {
      default_repsample_loss(out, cont_norm, bincat_norm, mean_norm, sd_norm, perc_norm, theoretical)
    } else {
      as.numeric(objective(out))
    }
    list(seed = seed, loss = loss, fit = out)
  }

  results <- NULL
  if (mode == "serial") {
    results <- lapply(seeds, eval_seed)
  } else if (mode == "multicore") {
    worker_n <- min(n_outer_workers, length(seeds))
    results <- parallel::mclapply(
      seeds,
      eval_seed,
      mc.cores = worker_n,
      mc.preschedule = TRUE,
      mc.set.seed = FALSE,
      mc.cleanup = FALSE,
      mc.allow.recursive = FALSE
    )
  } else {
    worker_n <- min(n_outer_workers, length(seeds))
    cl <- parallel::makeCluster(worker_n)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    libs <- .libPaths()
    parallel::clusterCall(cl, function(x) .libPaths(x), libs)
    parallel::clusterEvalQ(cl, library(RepsampleR, quietly = TRUE, warn.conflicts = FALSE))
    parallel::clusterExport(
      cl,
      varlist = c(
        "data", "size", "cont", "bincat", "mean", "sd", "perc", "extra_args",
        "objective", "theoretical", "cont_norm", "bincat_norm",
        "mean_norm", "sd_norm", "perc_norm", "eval_seed"
      ),
      envir = environment()
    )
    results <- parallel::parLapply(cl, seeds, eval_seed)
  }

  losses <- vapply(results, function(x) as.numeric(x$loss), numeric(1))
  best_idx <- which.min(losses)
  best <- results[[best_idx]]

  summary_df <- data.frame(
    seed = vapply(results, function(x) as.integer(x$seed), integer(1)),
    loss = losses,
    stringsAsFactors = FALSE
  )
  summary_df <- summary_df[order(summary_df$loss), , drop = FALSE]
  rownames(summary_df) <- NULL

  structure(
    list(
      best = best$fit,
      best_seed = as.integer(best$seed),
      best_loss = as.numeric(best$loss),
      summary = summary_df,
      all = if (keep_all) results else NULL,
      meta = list(
        n_seeds = length(seeds),
        n_outer_workers = n_outer_workers,
        outer_mode = mode,
        backend = backend_choice,
        theoretical = theoretical
      )
    ),
    class = "repsample_search_result"
  )
}

#' One-Command Multi-Stage Search Wrapper for `repsample()`
#'
#' Runs a coarse-to-fine outer seed search in one call by chaining
#' `repsample_search()` stages. Stage 1 performs a broad scan; each later
#' stage refines around the best seeds from the previous stage.
#'
#' @param data Data frame.
#' @param size Desired sample size.
#' @param cont Continuous variable names.
#' @param bincat Binary/categorical variable names.
#' @param mean Means for continuous variables (theoretical sampling only).
#' @param sd Standard deviations for continuous variables (theoretical sampling only).
#' @param perc Percentages for binary variables (theoretical sampling only).
#' @param n_stages Number of search stages.
#' @param n_seeds Seed budget per stage. Supply one integer to reuse across
#'   all stages, or one value per stage.
#' @param seed_start Starting seed for stage 1 when `seeds_stage1` is not provided.
#' @param seeds_stage1 Optional explicit seed vector for stage 1.
#' @param top_k Number of best seeds from each stage used to build the next
#'   stage's refinement neighborhood.
#' @param refine_radius Refinement radius around each retained seed for
#'   later stages. Provide one value to decay by `radius_shrink`, or one
#'   value per refinement stage.
#' @param radius_shrink Multiplicative decay for scalar `refine_radius`.
#' @param allow_repeats If `TRUE`, allow the same seed to be re-evaluated in
#'   multiple stages.
#' @param n_outer_workers Number of workers for each stage's outer search.
#' @param outer_parallel Outer parallel backend for each stage:
#'   `"auto"`, `"serial"`, `"multicore"`, or `"psock"`.
#' @param method Sampling method:
#'   `"auto"` (default), `"greedy"`, `"importance"`, or `"nearest"`.
#' @param objective Optional custom objective function passed to
#'   `repsample_search()`.
#' @param stop_loss Optional scalar threshold for early exit. If the best
#'   stage loss is less than or equal to `stop_loss`, the remaining stages
#'   are skipped.
#' @param keep_all If `TRUE`, retain all per-seed fit objects across all stages.
#' @param ... Additional arguments passed to `repsample()`.
#'
#' @return Object of class `repsample_search_result`. It follows the same
#' structure as `repsample_search()`, with `summary` including a `stage`
#' column and `meta$n_stages` describing the multi-stage run.
#'
#' @export
repsample_search_auto <- function(data,
                                  size,
                                  cont = NULL,
                                  bincat = NULL,
                                  mean = NULL,
                                  sd = NULL,
                                  perc = NULL,
                                  n_stages = 2,
                                  n_seeds = 24,
                                  seed_start = 1,
                                  seeds_stage1 = NULL,
                                  top_k = 4,
                                  refine_radius = NULL,
                                  radius_shrink = 0.5,
                                  allow_repeats = FALSE,
                                  n_outer_workers = 1,
                                  outer_parallel = c("auto", "serial", "multicore", "psock"),
                                  method = c("auto", "greedy", "importance", "nearest"),
                                  objective = NULL,
                                  stop_loss = NULL,
                                  keep_all = FALSE,
                                  ...) {
  outer_parallel <- match.arg(outer_parallel)
  method <- match.arg(method)
  n_stages <- check_positive_int_scalar(n_stages, "n_stages")
  stage_n_seeds <- normalize_stage_seed_counts(n_stages, n_seeds)
  n_outer_workers <- check_positive_int_scalar(n_outer_workers, "n_outer_workers")
  top_k <- check_positive_int_scalar(top_k, "top_k")

  if (!is.numeric(seed_start) || length(seed_start) != 1L || is.na(seed_start) ||
      abs(seed_start - round(seed_start)) > 0 || seed_start < 0 || seed_start > .Machine$integer.max) {
    stop(sprintf("`seed_start` must be an integer in [0, %d].", .Machine$integer.max), call. = FALSE)
  }
  seed_start <- as.integer(seed_start)

  if (!is.logical(allow_repeats) || length(allow_repeats) != 1L || is.na(allow_repeats)) {
    stop("`allow_repeats` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(objective) && !is.function(objective)) {
    stop("`objective` must be NULL or a function.", call. = FALSE)
  }
  if (!is.null(stop_loss)) {
    if (!is.numeric(stop_loss) || length(stop_loss) != 1L || is.na(stop_loss) || !is.finite(stop_loss)) {
      stop("`stop_loss` must be NULL or a single finite numeric value.", call. = FALSE)
    }
    stop_loss <- as.numeric(stop_loss)
  }
  if (!is.logical(keep_all) || length(keep_all) != 1L || is.na(keep_all)) {
    stop("`keep_all` must be TRUE or FALSE.", call. = FALSE)
  }

  stage_radius <- normalize_refine_radius(n_stages, stage_n_seeds, refine_radius, radius_shrink)
  stage_results <- vector("list", n_stages)
  seen_seeds <- integer(0)
  all_results <- if (keep_all) list() else NULL

  current_seeds <- if (is.null(seeds_stage1)) {
    sanitize_seed_vector(seed_start + seq_len(stage_n_seeds[1L]) - 1L, name = "Generated stage-1 `seeds`")
  } else {
    sanitize_seed_vector(seeds_stage1, name = "`seeds_stage1`")
  }
  if (!allow_repeats) {
    current_seeds <- current_seeds[!(current_seeds %in% seen_seeds)]
  }
  if (length(current_seeds) == 0L) {
    stop("Stage 1 has no seeds to evaluate.", call. = FALSE)
  }
  if (length(current_seeds) > stage_n_seeds[1L]) {
    current_seeds <- current_seeds[seq_len(stage_n_seeds[1L])]
  } else if (length(current_seeds) < stage_n_seeds[1L]) {
    current_seeds <- fill_seed_budget(
      center_seed = current_seeds[1L],
      seed_pool = current_seeds,
      target_n = stage_n_seeds[1L],
      seen_seeds = seen_seeds,
      allow_repeats = allow_repeats
    )
  }

  extra_args <- list(...)
  if ("seednum" %in% names(extra_args)) {
    warning(
      "Ignoring `seednum` in `...`; seeds are controlled by the multi-stage search.",
      call. = FALSE
    )
    extra_args$seednum <- NULL
  }

  backend_choice <- extra_args$backend
  if (is.null(backend_choice)) {
    backend_choice <- "cpu"
  }
  backend_choice <- match.arg(backend_choice, c("cpu", "cuda"))

  completed_stages <- 0L
  stopped_early <- FALSE
  stopped_stage <- NA_integer_

  for (stage_idx in seq_len(n_stages)) {
    call_args <- c(
      list(
        data = data,
        size = size,
        cont = cont,
        bincat = bincat,
        mean = mean,
        sd = sd,
        perc = perc,
        seeds = current_seeds,
        n_outer_workers = n_outer_workers,
        outer_parallel = outer_parallel,
        method = method,
        objective = objective,
        keep_all = keep_all
      ),
      extra_args
    )
    stage_out <- do.call(repsample_search, call_args)
    stage_results[[stage_idx]] <- stage_out
    completed_stages <- stage_idx

    if (keep_all && !is.null(stage_out$all)) {
      stage_all <- stage_out$all
      for (j in seq_along(stage_all)) {
        stage_all[[j]]$stage <- as.integer(stage_idx)
      }
      all_results <- c(all_results, stage_all)
    }

    seen_seeds <- unique(c(seen_seeds, as.integer(stage_out$summary$seed)))

    if (!is.null(stop_loss) && as.numeric(stage_out$best_loss) <= stop_loss) {
      stopped_early <- TRUE
      stopped_stage <- as.integer(stage_idx)
      break
    }

    if (stage_idx >= n_stages) {
      next
    }

    prev_summary <- stage_out$summary
    keep_top <- min(top_k, nrow(prev_summary))
    top_seeds <- as.integer(prev_summary$seed[seq_len(keep_top)])
    center_seed <- as.integer(prev_summary$seed[1L])

    current_seeds <- build_refine_seed_set(
      top_seeds = top_seeds,
      center_seed = center_seed,
      radius = stage_radius[stage_idx],
      target_n = stage_n_seeds[stage_idx + 1L],
      seen_seeds = seen_seeds,
      allow_repeats = allow_repeats
    )
  }

  if (completed_stages < 1L) {
    stop("No search stages were completed.", call. = FALSE)
  }
  stage_results <- stage_results[seq_len(completed_stages)]

  stage_best_loss <- vapply(stage_results, function(x) as.numeric(x$best_loss), numeric(1))
  stage_best_seed <- vapply(stage_results, function(x) as.integer(x$best_seed), integer(1))
  stage_best_idx <- which.min(stage_best_loss)

  summary_parts <- lapply(seq_along(stage_results), function(i) {
    x <- stage_results[[i]]$summary
    data.frame(
      stage = as.integer(i),
      seed = as.integer(x$seed),
      loss = as.numeric(x$loss),
      stringsAsFactors = FALSE
    )
  })
  combined_summary <- do.call(rbind, summary_parts)
  combined_summary <- combined_summary[order(combined_summary$loss), , drop = FALSE]
  rownames(combined_summary) <- NULL

  stage_summary <- data.frame(
    stage = seq_len(completed_stages),
    best_seed = stage_best_seed,
    best_loss = stage_best_loss,
    n_seeds = vapply(stage_results, function(x) nrow(x$summary), integer(1)),
    outer_mode = vapply(stage_results, function(x) as.character(x$meta$outer_mode), character(1)),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      best = stage_results[[stage_best_idx]]$best,
      best_seed = stage_best_seed[stage_best_idx],
      best_loss = stage_best_loss[stage_best_idx],
      summary = combined_summary,
      all = if (keep_all) all_results else NULL,
      stages = stage_results,
      stage_summary = stage_summary,
      meta = list(
        n_seeds = nrow(combined_summary),
        n_outer_workers = n_outer_workers,
        outer_mode = "multistage",
        backend = backend_choice,
        theoretical = as.logical(stage_results[[1L]]$meta$theoretical),
        n_stages = completed_stages,
        n_stages_requested = n_stages,
        top_k = top_k,
        stage_n_seeds = stage_n_seeds[seq_len(completed_stages)],
        stage_refine_radius = if (completed_stages > 1L) stage_radius[seq_len(completed_stages - 1L)] else integer(0),
        allow_repeats = allow_repeats,
        stop_loss = stop_loss,
        stopped_early = stopped_early,
        stopped_stage = stopped_stage
      )
    ),
    class = "repsample_search_result"
  )
}

#' Print a `repsample_search_result`
#'
#' @param x A `repsample_search_result` object.
#' @param ... Unused.
#'
#' @export
print.repsample_search_result <- function(x, ...) {
  cat("repsample_search_result\n")
  cat(sprintf("  Seeds searched: %d\n", x$meta$n_seeds))
  if (!is.null(x$meta$n_stages)) {
    cat(sprintf("  Stages: %d\n", x$meta$n_stages))
  }
  if (isTRUE(x$meta$stopped_early)) {
    cat(sprintf("  Early stop: yes (stage %d)\n", x$meta$stopped_stage))
  }
  cat(sprintf("  Outer mode: %s (%d workers)\n", x$meta$outer_mode, x$meta$n_outer_workers))
  cat(sprintf("  Backend: %s\n", x$meta$backend))
  cat(sprintf("  Best seed: %d\n", x$best_seed))
  cat(sprintf("  Best loss: %.6g\n", x$best_loss))
  invisible(x)
}
