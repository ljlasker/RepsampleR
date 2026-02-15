.use_importance_sampling <- function(cont, bincat, mean, sd, dist, exact) {
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

.importance_target_density <- function(x, mean, sd, dist) {
  if (dist == "normal") {
    return(stats::dnorm(x, mean = mean, sd = sd))
  }

  if (dist == "lognormal") {
    sigma2 <- log1p((sd * sd) / (mean * mean))
    sdlog <- sqrt(sigma2)
    meanlog <- log(mean) - 0.5 * sigma2
    return(stats::dlnorm(x, meanlog = meanlog, sdlog = sdlog))
  }

  if (dist == "poisson") {
    return(stats::dpois(floor(x), lambda = mean))
  }

  stop(sprintf("Unsupported distribution `%s`.", dist), call. = FALSE)
}

.importance_target_cdf <- function(mean, sd, dist) {
  if (dist == "normal") {
    return(function(q) stats::pnorm(q, mean = mean, sd = sd))
  }

  if (dist == "lognormal") {
    sigma2 <- log1p((sd * sd) / (mean * mean))
    sdlog <- sqrt(sigma2)
    meanlog <- log(mean) - 0.5 * sigma2
    return(function(q) stats::plnorm(q, meanlog = meanlog, sdlog = sdlog))
  }

  if (dist == "poisson") {
    return(function(q) stats::ppois(floor(q), lambda = mean))
  }

  stop(sprintf("Unsupported distribution `%s`.", dist), call. = FALSE)
}

.importance_baseline_density <- function(x) {
  out <- rep.int(.Machine$double.eps, length(x))
  finite <- is.finite(x)
  if (sum(finite) < 2L) {
    return(out)
  }

  x_fin <- as.numeric(x[finite])
  if (length(unique(x_fin)) < 2L) {
    out[finite] <- 1
    return(out)
  }

  dens <- tryCatch(
    stats::density(
      x_fin,
      n = 2048L,
      from = min(x_fin),
      to = max(x_fin)
    ),
    error = function(e) NULL
  )

  if (is.null(dens)) {
    return(out)
  }

  interp <- stats::approx(dens$x, dens$y, xout = x_fin, rule = 2)$y
  interp[!is.finite(interp) | interp <= 0] <- .Machine$double.eps
  out[finite] <- interp
  out
}

.importance_weights_multi <- function(data, cont, means, sds, dist) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }

  cont <- normalize_varlist(cont)
  means <- normalize_numlist(means)
  sds <- normalize_numlist(sds)
  dist_vec <- if (is.null(dist)) rep.int("normal", length(cont)) else normalize_cont_dist(dist, cont)

  if (length(cont) < 1L) {
    stop("`cont` must include at least one variable.", call. = FALSE)
  }
  if (length(means) != length(cont) || length(sds) != length(cont)) {
    stop("`means` and `sds` must match `cont` length.", call. = FALSE)
  }
  if (length(dist_vec) != length(cont)) {
    stop("`dist` must be length 1 or one per `cont` variable.", call. = FALSE)
  }

  w <- rep.int(1, nrow(data))

  for (j in seq_along(cont)) {
    varname <- cont[j]
    xj <- as.numeric(data[[varname]])
    dtarget <- .importance_target_density(xj, means[j], sds[j], dist_vec[j])
    dbase <- .importance_baseline_density(xj)
    ratio <- dtarget / pmax(dbase, .Machine$double.eps)
    ratio[!is.finite(ratio) | ratio < 0] <- 0
    w <- w * ratio
  }

  w[!is.finite(w) | w < 0] <- 0
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    return(rep.int(0, nrow(data)))
  }
  w / sw
}

.ess <- function(w) {
  w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    return(0)
  }
  wn <- w / sw
  1 / sum(wn ^ 2)
}

.check_overlap <- function(w, size, fallback = TRUE) {
  ess <- .ess(w)
  n_pos <- sum(is.finite(w) & (w > 0))

  if (n_pos < size || ess < size) {
    if (fallback) {
      warning(
        sprintf(
          "ESS (%.0f) < requested size (%d). Poor target/baseline overlap. Falling back to greedy algorithm.",
          ess, size
        ),
        call. = FALSE
      )
      return(FALSE)
    }

    stop(
      sprintf(
        "ESS (%.0f) < requested size (%d). Baseline has poor overlap with target.",
        ess, size
      ),
      call. = FALSE
    )
  }

  if (ess < size * 1.5) {
    warning(
      sprintf(
        "ESS (%.0f) is close to requested size (%d). Results may be noisy.",
        ess, size
      ),
      call. = FALSE
    )
  }

  TRUE
}

.importance_prepare <- function(data, size, cont, mean, sd, dist, fallback_to_greedy) {
  cont_norm <- normalize_varlist(cont)
  mean_norm <- normalize_numlist(mean)
  sd_norm <- normalize_numlist(sd)
  dist_norm <- if (is.null(dist)) rep.int("normal", length(cont_norm)) else normalize_cont_dist(dist, cont_norm)

  if (length(cont_norm) < 1L) {
    stop("Importance sampling requires at least one continuous variable.", call. = FALSE)
  }

  keep <- stats::complete.cases(data[, cont_norm, drop = FALSE])
  pool_idx <- which(keep)
  if (length(pool_idx) < size) {
    stop("Not enough complete cases for the requested sample size.", call. = FALSE)
  }

  pool <- data[pool_idx, , drop = FALSE]
  w <- .importance_weights_multi(pool, cont_norm, mean_norm, sd_norm, dist_norm)
  ok <- .check_overlap(w, size = size, fallback = fallback_to_greedy)
  if (!ok) {
    return(NULL)
  }

  list(
    pool = pool,
    pool_idx = pool_idx,
    cont = cont_norm,
    mean = mean_norm,
    sd = sd_norm,
    dist = dist_norm,
    weights = w,
    ess = .ess(w)
  )
}

.importance_with_seed <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(as.integer(seed))
  force(expr)
}

.importance_ks_score <- function(chosen_df, cont, mean, sd, dist) {
  ks_vals <- numeric(length(cont))

  for (i in seq_along(cont)) {
    cdf_fun <- .importance_target_cdf(mean[i], sd[i], dist[i])
    stat <- tryCatch(
      as.numeric(suppressWarnings(stats::ks.test(chosen_df[[cont[i]]], cdf_fun)$statistic)),
      error = function(e) Inf
    )
    ks_vals[i] <- stat
  }

  if (any(!is.finite(ks_vals))) {
    return(Inf)
  }
  mean(ks_vals)
}

.importance_build_fit <- function(data, selected_rows, size, seed, prep, ks, include_weights = FALSE) {
  out <- data
  out$repsample <- 0L
  out$repsample[selected_rows] <- 1L

  meta <- list(
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
    backend = "importance_sampling",
    cuda_fastpath = FALSE,
    quality = "importance",
    cont_dist = if (length(prep$cont) > 0L) stats::setNames(as.character(prep$dist), prep$cont) else character(0),
    exact = FALSE,
    weighted = FALSE,
    retain = NULL,
    method = "importance_sampling",
    ks = as.numeric(ks),
    ess = as.numeric(prep$ess),
    weights = if (isTRUE(include_weights)) prep$weights else NULL
  )

  structure(
    list(
      data = out,
      r = numeric(0),
      selected_rows = which(out$repsample == 1L),
      meta = meta
    ),
    class = "repsample_result"
  )
}

.importance_eval_seed <- function(data, size, seed, prep, include_weights = FALSE, objective = NULL) {
  pick_local <- .importance_with_seed(seed, sample.int(nrow(prep$pool), size = size, replace = FALSE, prob = prep$weights))
  selected_rows <- prep$pool_idx[pick_local]
  chosen <- data[selected_rows, , drop = FALSE]
  ks <- .importance_ks_score(chosen, prep$cont, prep$mean, prep$sd, prep$dist)
  fit <- .importance_build_fit(data, selected_rows, size, seed, prep, ks, include_weights = include_weights)
  loss <- if (is.null(objective)) ks else as.numeric(objective(fit))

  list(
    seed = as.integer(seed),
    loss = as.numeric(loss),
    fit = fit
  )
}

.importance_sample_one <- function(data,
                                   size,
                                   cont,
                                   mean,
                                   sd,
                                   dist = NULL,
                                   seed = 7L,
                                   fallback_to_greedy = TRUE,
                                   include_weights = FALSE) {
  prep <- .importance_prepare(
    data = data,
    size = size,
    cont = cont,
    mean = mean,
    sd = sd,
    dist = dist,
    fallback_to_greedy = fallback_to_greedy
  )
  if (is.null(prep)) {
    return(NULL)
  }

  .importance_eval_seed(
    data = data,
    size = size,
    seed = seed,
    prep = prep,
    include_weights = include_weights,
    objective = NULL
  )$fit
}

.importance_sample_search <- function(data,
                                      size,
                                      cont,
                                      mean,
                                      sd,
                                      dist = NULL,
                                      seeds,
                                      n_outer_workers = 1L,
                                      outer_mode = c("auto", "serial", "multicore", "psock"),
                                      keep_all = FALSE,
                                      objective = NULL,
                                      fallback_to_greedy = TRUE,
                                      include_weights = FALSE) {
  outer_mode <- match.arg(outer_mode)
  if (length(seeds) < 1L) {
    stop("`seeds` must be non-empty.", call. = FALSE)
  }

  prep <- .importance_prepare(
    data = data,
    size = size,
    cont = cont,
    mean = mean,
    sd = sd,
    dist = dist,
    fallback_to_greedy = fallback_to_greedy
  )
  if (is.null(prep)) {
    return(NULL)
  }

  eval_seed <- function(seed) {
    .importance_eval_seed(
      data = data,
      size = size,
      seed = seed,
      prep = prep,
      include_weights = include_weights,
      objective = objective
    )
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
        backend = "importance_sampling",
        theoretical = TRUE,
        method = "importance_sampling"
      )
    ),
    class = "repsample_search_result"
  )
}
