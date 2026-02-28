.repsample_control_defaults <- function() {
  list(
    seed = 7L,
    n_seeds = 16L,
    seed_start = 1L,
    n_outer_workers = 1L,
    n_cores = 1L,
    backend = "cpu",
    cuda_fallback = TRUE,
    nearest_replace = FALSE,
    nearest_distance = "euclidean",
    nearest_backend = "auto",
    nearest_match = "greedy",
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
    objective = NULL,
    ks_weight = 1,
    n_stages = 2L,
    top_k = 4L,
    refine_radius = NULL,
    radius_shrink = 0.5,
    allow_repeats = FALSE,
    stop_loss = NULL,
    checkpoint_file = NULL,
    resume = FALSE,
    adaptive_methods = c("importance", "nearest", "greedy"),
    adaptive_pilot_n = 2L
  )
}

.merge_repsample_control <- function(control) {
  defaults <- .repsample_control_defaults()
  if (is.null(control)) {
    return(defaults)
  }
  if (!is.list(control)) {
    stop("`control` must be NULL or a list.", call. = FALSE)
  }
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown) > 0L) {
    stop(
      sprintf("Unknown `control` fields: %s", paste(unknown, collapse = ", ")),
      call. = FALSE
    )
  }
  defaults[names(control)] <- control
  defaults
}

.adaptive_seed_search <- function(data,
                                  size,
                                  cont,
                                  bincat,
                                  mean,
                                  sd,
                                  perc,
                                  dist,
                                  quality,
                                  ctrl,
                                  ...) {
  methods <- as.character(ctrl$adaptive_methods)
  methods <- unique(methods[methods %in% c("greedy", "importance", "nearest", "auto")])
  if (length(methods) < 1L) {
    methods <- c("importance", "nearest", "greedy")
  }

  pilot_n <- check_positive_int_scalar(ctrl$adaptive_pilot_n, "control$adaptive_pilot_n")
  total_budget <- check_positive_int_scalar(ctrl$n_seeds, "control$n_seeds")
  objective <- ctrl$objective
  if (is.null(objective)) {
    objective <- function(out) {
      repsample_quality(
        out = out,
        cont = cont,
        bincat = bincat,
        mean = mean,
        sd = sd,
        perc = perc,
        ks_weight = ctrl$ks_weight
      )$loss
    }
  }

  pilot_out <- list()
  pilot_rows <- list()
  seed_cursor <- as.integer(ctrl$seed_start)

  for (m in methods) {
    out <- tryCatch(
      repsample_search(
        data = data,
        size = size,
        cont = cont,
        bincat = bincat,
        mean = mean,
        sd = sd,
        perc = perc,
        n_seeds = pilot_n,
        seed_start = seed_cursor,
        n_outer_workers = ctrl$n_outer_workers,
        outer_parallel = "auto",
        method = m,
        nearest_replace = ctrl$nearest_replace,
        nearest_distance = ctrl$nearest_distance,
        nearest_backend = ctrl$nearest_backend,
        nearest_match = ctrl$nearest_match,
        nearest_feature_weights = ctrl$nearest_feature_weights,
        nearest_caliper = ctrl$nearest_caliper,
        nearest_caliper_strict = ctrl$nearest_caliper_strict,
        nearest_optimal_max_size = ctrl$nearest_optimal_max_size,
        nearest_hnsw_k = ctrl$nearest_hnsw_k,
        nearest_hnsw_M = ctrl$nearest_hnsw_M,
        nearest_hnsw_ef_construction = ctrl$nearest_hnsw_ef_construction,
        nearest_hnsw_ef_search = ctrl$nearest_hnsw_ef_search,
        mean_tol = ctrl$mean_tol,
        sd_tol = ctrl$sd_tol,
        ks_tol = ctrl$ks_tol,
        perc_tol = ctrl$perc_tol,
        objective = objective,
        keep_all = TRUE,
        dist = dist,
        quality = quality,
        n_cores = ctrl$n_cores,
        backend = ctrl$backend,
        cuda_fallback = ctrl$cuda_fallback,
        ...
      ),
      error = function(e) e
    )

    if (inherits(out, "error")) {
      seed_cursor <- seed_cursor + pilot_n
      next
    }

    pilot_out[[m]] <- out
    rows <- out$summary
    rows$method <- m
    rows$phase <- "pilot"
    pilot_rows[[length(pilot_rows) + 1L]] <- rows
    seed_cursor <- seed_cursor + pilot_n
  }

  if (length(pilot_out) < 1L) {
    stop("Adaptive search could not initialize any pilot method.", call. = FALSE)
  }

  pilot_best_losses <- vapply(pilot_out, function(x) as.numeric(x$best_loss), numeric(1))
  best_method <- names(which.min(pilot_best_losses))[1L]
  pilot_best <- pilot_out[[best_method]]

  spent <- pilot_n * length(pilot_out)
  remain <- max(0L, as.integer(total_budget) - spent)

  final_out <- NULL
  if (remain > 0L) {
    final_out <- repsample_search(
      data = data,
      size = size,
      cont = cont,
      bincat = bincat,
      mean = mean,
      sd = sd,
      perc = perc,
      n_seeds = remain,
      seed_start = seed_cursor,
      n_outer_workers = ctrl$n_outer_workers,
      outer_parallel = "auto",
      method = best_method,
      nearest_replace = ctrl$nearest_replace,
      nearest_distance = ctrl$nearest_distance,
      nearest_backend = ctrl$nearest_backend,
      nearest_match = ctrl$nearest_match,
      nearest_feature_weights = ctrl$nearest_feature_weights,
      nearest_caliper = ctrl$nearest_caliper,
      nearest_caliper_strict = ctrl$nearest_caliper_strict,
      nearest_optimal_max_size = ctrl$nearest_optimal_max_size,
      nearest_hnsw_k = ctrl$nearest_hnsw_k,
      nearest_hnsw_M = ctrl$nearest_hnsw_M,
      nearest_hnsw_ef_construction = ctrl$nearest_hnsw_ef_construction,
      nearest_hnsw_ef_search = ctrl$nearest_hnsw_ef_search,
      mean_tol = ctrl$mean_tol,
      sd_tol = ctrl$sd_tol,
      ks_tol = ctrl$ks_tol,
      perc_tol = ctrl$perc_tol,
      objective = objective,
      keep_all = TRUE,
      dist = dist,
      quality = quality,
      n_cores = ctrl$n_cores,
      backend = ctrl$backend,
      cuda_fallback = ctrl$cuda_fallback,
      ...
    )
  }

  all_entries <- list()
  summary_rows <- pilot_rows
  for (m in names(pilot_out)) {
    all_entries <- c(all_entries, pilot_out[[m]]$all)
  }
  if (!is.null(final_out)) {
    rows <- final_out$summary
    rows$method <- best_method
    rows$phase <- "exploit"
    summary_rows[[length(summary_rows) + 1L]] <- rows
    all_entries <- c(all_entries, final_out$all)
  }

  if (length(all_entries) < 1L) {
    stop("Adaptive search returned no evaluated seeds.", call. = FALSE)
  }

  all_losses <- vapply(all_entries, function(x) as.numeric(x$loss), numeric(1))
  all_seeds <- vapply(all_entries, function(x) as.integer(x$seed), integer(1))
  best_idx <- which.min(all_losses)
  summary_df <- do.call(rbind, summary_rows)
  summary_df <- summary_df[order(summary_df$loss), , drop = FALSE]
  rownames(summary_df) <- NULL

  pilot_table <- data.frame(
    method = names(pilot_out),
    best_loss = as.numeric(pilot_best_losses[names(pilot_out)]),
    stringsAsFactors = FALSE
  )
  pilot_table <- pilot_table[order(pilot_table$best_loss), , drop = FALSE]
  rownames(pilot_table) <- NULL

  structure(
    list(
      best = all_entries[[best_idx]]$fit,
      best_seed = all_seeds[best_idx],
      best_loss = all_losses[best_idx],
      summary = summary_df,
      all = all_entries,
      meta = list(
        n_seeds = nrow(summary_df),
        n_outer_workers = ctrl$n_outer_workers,
        outer_mode = "adaptive",
        backend = ctrl$backend,
        theoretical = TRUE,
        method = "adaptive",
        selected_method = best_method,
        pilot = pilot_table,
        pilot_n = pilot_n,
        total_budget = total_budget,
        exploit_budget = remain
      )
    ),
    class = "repsample_search_result"
  )
}

#' One-Command RepsampleR API with Simplified Defaults
#'
#' `repsample_fit()` is a high-level entry point that wraps
#' `repsample_easy()`, `repsample_search()`, and `repsample_search_auto()`
#' behind a single interface and a `control` list.
#'
#' @param data Data frame.
#' @param size Desired sample size.
#' @param cont Continuous variable names.
#' @param bincat Binary/categorical variable names.
#' @param mean Target means for `cont`.
#' @param sd Target SDs for `cont`.
#' @param perc Target percentages for `bincat`.
#' @param dist Target distribution form(s): `"normal"`, `"lognormal"`, `"poisson"`.
#' @param quality Tuning preset for greedy mode: `"fast"`, `"balanced"`,
#'   `"strict"`, or `"manual"`.
#' @param method Sampling method: `"auto"`, `"adaptive"`, `"greedy"`,
#'   `"importance"`, or `"nearest"`.
#' @param mode Execution mode: `"auto"`, `"single"`, `"search"`, or `"search_auto"`.
#' @param control Named list of advanced controls (seeds, parallelism,
#'   nearest options, tolerances, checkpointing, etc.).
#' @param ... Additional arguments forwarded to lower-level engines.
#'
#' @return A `repsample_result` (`single`) or `repsample_search_result`
#'   (`search`, `search_auto`, `adaptive`).
#'
#' @export
repsample_fit <- function(data,
                          size,
                          cont = NULL,
                          bincat = NULL,
                          mean = NULL,
                          sd = NULL,
                          perc = NULL,
                          dist = "normal",
                          quality = c("balanced", "fast", "strict", "manual"),
                          method = c("auto", "adaptive", "greedy", "importance", "nearest"),
                          mode = c("auto", "single", "search", "search_auto"),
                          control = NULL,
                          ...) {
  quality <- match.arg(quality)
  method <- match.arg(method)
  mode <- match.arg(mode)
  ctrl <- .merge_repsample_control(control)

  run_mode <- mode
  if (mode == "auto") {
    if (as.integer(ctrl$n_stages) > 1L && as.integer(ctrl$n_seeds) > 1L) {
      run_mode <- "search_auto"
    } else if (as.integer(ctrl$n_seeds) > 1L) {
      run_mode <- "search"
    } else {
      run_mode <- "single"
    }
  }

  if (method == "adaptive" && run_mode == "single") {
    method <- "auto"
  }

  if (method == "adaptive") {
    return(.adaptive_seed_search(
      data = data,
      size = size,
      cont = cont,
      bincat = bincat,
      mean = mean,
      sd = sd,
      perc = perc,
      dist = dist,
      quality = quality,
      ctrl = ctrl,
      ...
    ))
  }

  if (run_mode == "single") {
    return(repsample_easy(
      data = data,
      size = size,
      cont = cont,
      mean = mean,
      sd = sd,
      dist = dist,
      quality = quality,
      mode = "single",
      seed = ctrl$seed,
      n_cores = ctrl$n_cores,
      backend = ctrl$backend,
      method = method,
      nearest_replace = ctrl$nearest_replace,
      nearest_distance = ctrl$nearest_distance,
      nearest_backend = ctrl$nearest_backend,
      nearest_match = ctrl$nearest_match,
      nearest_feature_weights = ctrl$nearest_feature_weights,
      nearest_caliper = ctrl$nearest_caliper,
      nearest_caliper_strict = ctrl$nearest_caliper_strict,
      nearest_optimal_max_size = ctrl$nearest_optimal_max_size,
      nearest_hnsw_k = ctrl$nearest_hnsw_k,
      nearest_hnsw_M = ctrl$nearest_hnsw_M,
      nearest_hnsw_ef_construction = ctrl$nearest_hnsw_ef_construction,
      nearest_hnsw_ef_search = ctrl$nearest_hnsw_ef_search,
      cuda_fallback = ctrl$cuda_fallback,
      bincat = bincat,
      perc = perc,
      ...
    ))
  }

  if (run_mode == "search") {
    return(repsample_search(
      data = data,
      size = size,
      cont = cont,
      bincat = bincat,
      mean = mean,
      sd = sd,
      perc = perc,
      n_seeds = ctrl$n_seeds,
      seed_start = ctrl$seed_start,
      n_outer_workers = ctrl$n_outer_workers,
      outer_parallel = "auto",
      method = method,
      nearest_replace = ctrl$nearest_replace,
      nearest_distance = ctrl$nearest_distance,
      nearest_backend = ctrl$nearest_backend,
      nearest_match = ctrl$nearest_match,
      nearest_feature_weights = ctrl$nearest_feature_weights,
      nearest_caliper = ctrl$nearest_caliper,
      nearest_caliper_strict = ctrl$nearest_caliper_strict,
      nearest_optimal_max_size = ctrl$nearest_optimal_max_size,
      nearest_hnsw_k = ctrl$nearest_hnsw_k,
      nearest_hnsw_M = ctrl$nearest_hnsw_M,
      nearest_hnsw_ef_construction = ctrl$nearest_hnsw_ef_construction,
      nearest_hnsw_ef_search = ctrl$nearest_hnsw_ef_search,
      mean_tol = ctrl$mean_tol,
      sd_tol = ctrl$sd_tol,
      ks_tol = ctrl$ks_tol,
      perc_tol = ctrl$perc_tol,
      objective = ctrl$objective,
      dist = dist,
      quality = quality,
      n_cores = ctrl$n_cores,
      backend = ctrl$backend,
      cuda_fallback = ctrl$cuda_fallback,
      ...
    ))
  }

  repsample_search_auto(
    data = data,
    size = size,
    cont = cont,
    bincat = bincat,
    mean = mean,
    sd = sd,
    perc = perc,
    n_stages = ctrl$n_stages,
    n_seeds = ctrl$n_seeds,
    seed_start = ctrl$seed_start,
    top_k = ctrl$top_k,
    refine_radius = ctrl$refine_radius,
    radius_shrink = ctrl$radius_shrink,
    allow_repeats = ctrl$allow_repeats,
    n_outer_workers = ctrl$n_outer_workers,
    outer_parallel = "auto",
    method = method,
    nearest_replace = ctrl$nearest_replace,
    nearest_distance = ctrl$nearest_distance,
    nearest_backend = ctrl$nearest_backend,
    nearest_match = ctrl$nearest_match,
    nearest_feature_weights = ctrl$nearest_feature_weights,
    nearest_caliper = ctrl$nearest_caliper,
    nearest_caliper_strict = ctrl$nearest_caliper_strict,
    nearest_optimal_max_size = ctrl$nearest_optimal_max_size,
    nearest_hnsw_k = ctrl$nearest_hnsw_k,
    nearest_hnsw_M = ctrl$nearest_hnsw_M,
    nearest_hnsw_ef_construction = ctrl$nearest_hnsw_ef_construction,
    nearest_hnsw_ef_search = ctrl$nearest_hnsw_ef_search,
    objective = ctrl$objective,
    stop_loss = ctrl$stop_loss,
    mean_tol = ctrl$mean_tol,
    sd_tol = ctrl$sd_tol,
    ks_tol = ctrl$ks_tol,
    perc_tol = ctrl$perc_tol,
    checkpoint_file = ctrl$checkpoint_file,
    resume = ctrl$resume,
    dist = dist,
    quality = quality,
    n_cores = ctrl$n_cores,
    backend = ctrl$backend,
    cuda_fallback = ctrl$cuda_fallback,
    ...
  )
}
