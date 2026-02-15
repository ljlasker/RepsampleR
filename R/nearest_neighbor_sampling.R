.use_nearest_neighbor_sampling <- function(cont, bincat, mean, sd, dist, exact) {
  cont_norm <- normalize_varlist(cont)
  bincat_norm <- normalize_varlist(bincat)
  mean_norm <- normalize_numlist(mean)
  sd_norm <- normalize_numlist(sd)

  has_cont <- length(cont_norm) >= 1L
  no_bincat <- length(bincat_norm) == 0L
  has_target <- length(mean_norm) > 0L && length(sd_norm) > 0L
  not_exact <- !isTRUE(exact)

  if (!(has_cont && no_bincat && has_target && not_exact)) {
    return(FALSE)
  }
  if (length(mean_norm) != length(cont_norm) || length(sd_norm) != length(cont_norm)) {
    return(FALSE)
  }

  dist_norm <- tryCatch(
    {
      if (is.null(dist)) {
        rep.int("normal", length(cont_norm))
      } else {
        normalize_cont_dist(dist, cont_norm)
      }
    },
    error = function(e) character(0)
  )

  length(dist_norm) == length(cont_norm) &&
    all(dist_norm %in% c("normal", "lognormal", "poisson"))
}

.nearest_target_draw <- function(n, mean, sd, dist) {
  if (dist == "normal") {
    return(stats::rnorm(n, mean = mean, sd = sd))
  }

  if (dist == "lognormal") {
    sigma2 <- log1p((sd * sd) / (mean * mean))
    sdlog <- sqrt(sigma2)
    meanlog <- log(mean) - 0.5 * sigma2
    return(stats::rlnorm(n, meanlog = meanlog, sdlog = sdlog))
  }

  if (dist == "poisson") {
    return(as.numeric(stats::rpois(n, lambda = mean)))
  }

  stop(sprintf("Unsupported distribution `%s`.", dist), call. = FALSE)
}

.nearest_prepare <- function(data, size, cont, mean, sd, dist) {
  cont_norm <- normalize_varlist(cont)
  mean_norm <- normalize_numlist(mean)
  sd_norm <- normalize_numlist(sd)
  dist_norm <- if (is.null(dist)) rep.int("normal", length(cont_norm)) else normalize_cont_dist(dist, cont_norm)

  if (length(cont_norm) < 1L) {
    stop("Nearest-neighbor sampling requires at least one continuous variable.", call. = FALSE)
  }
  if (length(mean_norm) != length(cont_norm) || length(sd_norm) != length(cont_norm)) {
    stop("`mean`/`sd` lengths must match `cont`.", call. = FALSE)
  }

  keep <- stats::complete.cases(data[, cont_norm, drop = FALSE])
  pool_idx <- which(keep)
  if (length(pool_idx) < size) {
    stop("Not enough complete cases for the requested sample size.", call. = FALSE)
  }

  pool <- data[pool_idx, , drop = FALSE]
  x_mat <- as.matrix(pool[, cont_norm, drop = FALSE])
  storage.mode(x_mat) <- "double"

  scale_vec <- pmax(as.numeric(sd_norm), 1e-12)
  x_scaled <- sweep(x_mat, 2L, scale_vec, "/", check.margin = FALSE)

  list(
    pool = pool,
    pool_idx = pool_idx,
    x_mat = x_mat,
    x_scaled = x_scaled,
    cont = cont_norm,
    mean = mean_norm,
    sd = sd_norm,
    dist = dist_norm,
    scale = scale_vec
  )
}

.nearest_match_indices <- function(x_scaled, target_scaled) {
  n_pool <- nrow(x_scaled)
  n_take <- nrow(target_scaled)
  available <- seq_len(n_pool)
  selected <- integer(n_take)

  for (i in seq_len(n_take)) {
    x_avail <- x_scaled[available, , drop = FALSE]
    t_row <- target_scaled[i, ]
    d <- rowSums((x_avail - matrix(t_row, nrow = nrow(x_avail), ncol = ncol(x_avail), byrow = TRUE)) ^ 2)
    k <- which.min(d)
    selected[i] <- available[k]
    available <- available[-k]
  }

  selected
}

.nearest_eval_seed <- function(data, size, seed, prep, objective = NULL) {
  p <- length(prep$cont)

  target_mat <- matrix(0, nrow = size, ncol = p)
  for (j in seq_len(p)) {
    target_mat[, j] <- .nearest_target_draw(
      n = size,
      mean = prep$mean[j],
      sd = prep$sd[j],
      dist = prep$dist[j]
    )
  }
  target_scaled <- sweep(target_mat, 2L, prep$scale, "/", check.margin = FALSE)

  pick_local <- .nearest_match_indices(prep$x_scaled, target_scaled)
  selected_rows <- prep$pool_idx[pick_local]
  chosen <- data[selected_rows, , drop = FALSE]

  ks <- .importance_ks_score(chosen, prep$cont, prep$mean, prep$sd, prep$dist)

  out <- data
  out$repsample <- 0L
  out$repsample[selected_rows] <- 1L

  fit <- structure(
    list(
      data = out,
      r = numeric(0),
      selected_rows = which(out$repsample == 1L),
      meta = list(
        size = size,
        mode = "theoretical",
        seednum = as.integer(seed),
        randomperc = NA_real_,
        randsel = NA_integer_,
        srule = NULL,
        rrule = NULL,
        n_cores = 1L,
        parallel_enabled = FALSE,
        parallel_mode = "none",
        psock_stateful_fastpath = FALSE,
        backend = "nearest_neighbor",
        cuda_fastpath = FALSE,
        quality = "nearest",
        cont_dist = if (length(prep$cont) > 0L) stats::setNames(as.character(prep$dist), prep$cont) else character(0),
        exact = FALSE,
        weighted = FALSE,
        retain = NULL,
        method = "nearest_neighbor",
        ks = as.numeric(ks)
      )
    ),
    class = "repsample_result"
  )

  loss <- if (is.null(objective)) as.numeric(ks) else as.numeric(objective(fit))
  list(seed = as.integer(seed), loss = as.numeric(loss), fit = fit)
}

.nearest_sample_one <- function(data,
                                size,
                                cont,
                                mean,
                                sd,
                                dist = NULL,
                                seed = 7L) {
  prep <- .nearest_prepare(
    data = data,
    size = size,
    cont = cont,
    mean = mean,
    sd = sd,
    dist = dist
  )

  .importance_with_seed(seed, .nearest_eval_seed(data, size, seed, prep, objective = NULL)$fit)
}

.nearest_sample_search <- function(data,
                                   size,
                                   cont,
                                   mean,
                                   sd,
                                   dist = NULL,
                                   seeds,
                                   n_outer_workers = 1L,
                                   outer_mode = c("auto", "serial", "multicore", "psock"),
                                   keep_all = FALSE,
                                   objective = NULL) {
  outer_mode <- match.arg(outer_mode)
  if (length(seeds) < 1L) {
    stop("`seeds` must be non-empty.", call. = FALSE)
  }

  prep <- .nearest_prepare(
    data = data,
    size = size,
    cont = cont,
    mean = mean,
    sd = sd,
    dist = dist
  )

  eval_seed <- function(seed) {
    .importance_with_seed(seed, .nearest_eval_seed(data, size, seed, prep, objective = objective))
  }

  results <- NULL
  if (outer_mode == "serial" || length(seeds) == 1L || n_outer_workers <= 1L) {
    results <- lapply(seeds, eval_seed)
  } else if (outer_mode == "multicore") {
    worker_n <- min(as.integer(n_outer_workers), length(seeds))
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
    worker_n <- min(as.integer(n_outer_workers), length(seeds))
    cl <- parallel::makeCluster(worker_n)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    libs <- .libPaths()
    parallel::clusterCall(cl, function(x) .libPaths(x), libs)
    parallel::clusterExport(cl, varlist = c("eval_seed"), envir = environment())
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
      all = if (isTRUE(keep_all)) results else NULL,
      meta = list(
        n_seeds = length(seeds),
        n_outer_workers = as.integer(n_outer_workers),
        outer_mode = outer_mode,
        backend = "nearest_neighbor",
        theoretical = TRUE,
        method = "nearest_neighbor"
      )
    ),
    class = "repsample_search_result"
  )
}
