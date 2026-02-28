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

.normalize_nearest_feature_weights <- function(feature_weights, p) {
  if (is.null(feature_weights)) {
    return(rep.int(1, p))
  }
  if (!is.numeric(feature_weights) || anyNA(feature_weights)) {
    stop("`nearest_feature_weights` must be NULL or a numeric vector.", call. = FALSE)
  }
  if (length(feature_weights) == 1L) {
    feature_weights <- rep.int(as.numeric(feature_weights), p)
  } else if (length(feature_weights) != p) {
    stop(
      sprintf("`nearest_feature_weights` must have length 1 or %d.", p),
      call. = FALSE
    )
  }
  if (any(feature_weights <= 0)) {
    stop("`nearest_feature_weights` values must be positive.", call. = FALSE)
  }
  as.numeric(feature_weights)
}

.normalize_nearest_caliper <- function(caliper, p) {
  if (is.null(caliper)) {
    return(NULL)
  }
  if (!is.numeric(caliper) || anyNA(caliper)) {
    stop("`nearest_caliper` must be NULL or a numeric scalar/vector.", call. = FALSE)
  }
  if (length(caliper) == 1L) {
    caliper <- rep.int(as.numeric(caliper), p)
  } else if (length(caliper) != p) {
    stop(
      sprintf("`nearest_caliper` must have length 1 or %d.", p),
      call. = FALSE
    )
  }
  if (any(caliper <= 0)) {
    stop("`nearest_caliper` values must be positive.", call. = FALSE)
  }
  as.numeric(caliper)
}

.stable_inverse_cov <- function(x_scaled) {
  p <- ncol(x_scaled)
  if (p <= 1L) {
    return(matrix(1, nrow = p, ncol = p))
  }

  cov_mat <- stats::cov(x_scaled)
  if (!is.matrix(cov_mat) || any(!is.finite(cov_mat))) {
    cov_mat <- diag(1, p)
  }

  eig <- eigen(cov_mat, symmetric = TRUE)
  vals <- pmax(as.numeric(eig$values), 1e-8)
  vecs <- eig$vectors
  vecs %*% (diag(1 / vals, nrow = p)) %*% t(vecs)
}

.nearest_distance_transform <- function(x_mat, distance, feature_weights, inv_cov) {
  if (distance == "weighted") {
    return(x_mat %*% diag(sqrt(feature_weights), nrow = length(feature_weights)))
  }
  if (distance == "mahalanobis") {
    if (is.null(inv_cov)) {
      return(x_mat)
    }
    return(x_mat %*% t(chol(inv_cov)))
  }
  x_mat
}

.nearest_hnsw_available <- function() {
  requireNamespace("RcppHNSW", quietly = TRUE)
}

.nearest_hnsw_build <- function(x, M = 16L, ef_construction = 200L) {
  ns <- asNamespace("RcppHNSW")

  if (exists("hnsw_build", envir = ns, inherits = FALSE)) {
    f <- get("hnsw_build", envir = ns, inherits = FALSE)
    calls <- list(
      list(X = x, distance = "l2", nlinks = as.integer(M), efConstruction = as.integer(ef_construction)),
      list(X = x, distance = "l2"),
      list(x = x, distance = "l2"),
      list(x = x)
    )
    for (args in calls) {
      out <- try(do.call(f, args), silent = TRUE)
      if (!inherits(out, "try-error")) {
        return(out)
      }
    }
  }

  if (exists("hnsw_build_index", envir = ns, inherits = FALSE)) {
    f <- get("hnsw_build_index", envir = ns, inherits = FALSE)
    calls <- list(
      list(X = x, distance = "l2", M = as.integer(M), ef = as.integer(ef_construction)),
      list(X = x, distance = "l2"),
      list(x = x, distance = "l2"),
      list(x = x)
    )
    for (args in calls) {
      out <- try(do.call(f, args), silent = TRUE)
      if (!inherits(out, "try-error")) {
        return(out)
      }
    }
  }

  stop("Unable to build HNSW index with installed `RcppHNSW` API.", call. = FALSE)
}

.nearest_hnsw_query <- function(index, query, k = 32L, ef_search = 64L) {
  ns <- asNamespace("RcppHNSW")
  out <- NULL

  if (exists("hnsw_search", envir = ns, inherits = FALSE)) {
    f <- get("hnsw_search", envir = ns, inherits = FALSE)
    calls <- list(
      list(index = index, query = query, k = as.integer(k), ef = as.integer(ef_search)),
      list(index = index, query = query, k = as.integer(k)),
      list(model = index, query = query, k = as.integer(k), ef = as.integer(ef_search)),
      list(model = index, query = query, k = as.integer(k))
    )
    for (args in calls) {
      out_try <- try(do.call(f, args), silent = TRUE)
      if (!inherits(out_try, "try-error")) {
        out <- out_try
        break
      }
    }
  }

  if (is.null(out) && exists("hnsw_knn", envir = ns, inherits = FALSE)) {
    f <- get("hnsw_knn", envir = ns, inherits = FALSE)
    calls <- list(
      list(index = index, query = query, k = as.integer(k), ef = as.integer(ef_search)),
      list(index = index, query = query, k = as.integer(k)),
      list(model = index, query = query, k = as.integer(k))
    )
    for (args in calls) {
      out_try <- try(do.call(f, args), silent = TRUE)
      if (!inherits(out_try, "try-error")) {
        out <- out_try
        break
      }
    }
  }

  if (is.null(out)) {
    stop("Unable to query HNSW index with installed `RcppHNSW` API.", call. = FALSE)
  }

  if (is.list(out)) {
    idx <- NULL
    for (nm in c("idx", "index", "indices", "id")) {
      if (!is.null(out[[nm]])) {
        idx <- out[[nm]]
        break
      }
    }
    if (is.null(idx) && length(out) >= 1L) {
      idx <- out[[1L]]
    }
    idx <- as.matrix(idx)
  } else {
    idx <- as.matrix(out)
  }

  storage.mode(idx) <- "integer"
  idx
}

.nearest_apply_caliper <- function(x_mat, t_row, caliper, caliper_strict) {
  if (is.null(caliper)) {
    return(list(mask = rep(TRUE, nrow(x_mat)), violated = FALSE))
  }
  delta <- abs(sweep(x_mat, 2L, t_row, "-", check.margin = FALSE))
  mask <- rowSums(delta <= matrix(caliper, nrow = nrow(delta), ncol = ncol(delta), byrow = TRUE)) == ncol(delta)
  if (!any(mask) && isTRUE(caliper_strict)) {
    stop("No candidate satisfied `nearest_caliper`.", call. = FALSE)
  }
  list(mask = mask, violated = !any(mask))
}

.nearest_optimal_assign <- function(cost_mat) {
  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("`nearest_match = \"optimal\"` requires the optional `clue` package.", call. = FALSE)
  }
  assign <- clue::solve_LSAP(cost_mat)
  as.integer(assign)
}

.resolve_nearest_backend <- function(backend,
                                     n_pool,
                                     n_take,
                                     hnsw_available) {
  if (backend == "exact") {
    return("exact")
  }
  if (backend == "hnsw") {
    if (!hnsw_available) {
      warning("HNSW backend requested but `RcppHNSW` is not installed; using exact backend.", call. = FALSE)
      return("exact")
    }
    return("hnsw")
  }

  # auto
  if (hnsw_available && (as.numeric(n_pool) * as.numeric(n_take) >= 1e6)) {
    return("hnsw")
  }
  "exact"
}

.nearest_prepare <- function(data,
                             size,
                             cont,
                             mean,
                             sd,
                             dist,
                             replace = FALSE,
                             distance = c("euclidean", "mahalanobis", "weighted"),
                             backend = c("auto", "exact", "hnsw"),
                             match = c("greedy", "optimal"),
                             feature_weights = NULL,
                             caliper = NULL,
                             caliper_strict = FALSE,
                             optimal_max_size = 2000L,
                             hnsw_k = 64L,
                             hnsw_M = 16L,
                             hnsw_ef_construction = 200L,
                             hnsw_ef_search = 64L) {
  cont_norm <- normalize_varlist(cont)
  mean_norm <- normalize_numlist(mean)
  sd_norm <- normalize_numlist(sd)
  dist_norm <- if (is.null(dist)) rep.int("normal", length(cont_norm)) else normalize_cont_dist(dist, cont_norm)
  distance <- match.arg(distance)
  backend <- match.arg(backend)
  match <- match.arg(match)
  check_scalar_flag(replace, "nearest_replace")
  check_scalar_flag(caliper_strict, "nearest_caliper_strict")
  optimal_max_size <- check_positive_int_scalar(optimal_max_size, "nearest_optimal_max_size")
  hnsw_k <- check_positive_int_scalar(hnsw_k, "nearest_hnsw_k")
  hnsw_M <- check_positive_int_scalar(hnsw_M, "nearest_hnsw_M")
  hnsw_ef_construction <- check_positive_int_scalar(hnsw_ef_construction, "nearest_hnsw_ef_construction")
  hnsw_ef_search <- check_positive_int_scalar(hnsw_ef_search, "nearest_hnsw_ef_search")

  if (length(cont_norm) < 1L) {
    stop("Nearest-neighbor sampling requires at least one continuous variable.", call. = FALSE)
  }
  if (length(mean_norm) != length(cont_norm) || length(sd_norm) != length(cont_norm)) {
    stop("`mean`/`sd` lengths must match `cont`.", call. = FALSE)
  }

  keep <- stats::complete.cases(data[, cont_norm, drop = FALSE])
  pool_idx <- which(keep)
  if (!isTRUE(replace) && length(pool_idx) < size) {
    stop("Not enough complete cases for the requested sample size.", call. = FALSE)
  }

  pool <- data[pool_idx, , drop = FALSE]
  x_mat <- as.matrix(pool[, cont_norm, drop = FALSE])
  storage.mode(x_mat) <- "double"

  scale_vec <- pmax(as.numeric(sd_norm), 1e-12)
  x_scaled <- sweep(x_mat, 2L, scale_vec, "/", check.margin = FALSE)
  p <- ncol(x_scaled)
  fw <- .normalize_nearest_feature_weights(feature_weights, p)
  caliper_norm <- .normalize_nearest_caliper(caliper, p)
  inv_cov <- if (distance == "mahalanobis") .stable_inverse_cov(x_scaled) else NULL
  hnsw_available <- .nearest_hnsw_available()
  backend_resolved <- .resolve_nearest_backend(
    backend = backend,
    n_pool = nrow(x_scaled),
    n_take = size,
    hnsw_available = hnsw_available
  )
  if (match == "optimal" && isTRUE(replace)) {
    warning("`nearest_match = \"optimal\"` is only relevant without replacement; using greedy matching.", call. = FALSE)
    match <- "greedy"
  }
  if (match == "optimal" && size > optimal_max_size) {
    warning(
      sprintf(
        "Requested optimal matching size (%d) exceeds `nearest_optimal_max_size` (%d); using greedy matching.",
        size,
        optimal_max_size
      ),
      call. = FALSE
    )
    match <- "greedy"
  }

  list(
    pool = pool,
    pool_idx = pool_idx,
    x_mat = x_mat,
    x_scaled = x_scaled,
    cont = cont_norm,
    mean = mean_norm,
    sd = sd_norm,
    dist = dist_norm,
    scale = scale_vec,
    replace = isTRUE(replace),
    distance = distance,
    nearest_backend = backend_resolved,
    nearest_match = match,
    feature_weights = fw,
    caliper = caliper_norm,
    caliper_strict = isTRUE(caliper_strict),
    inv_cov = inv_cov,
    nearest_backend_requested = backend,
    nearest_optimal_max_size = as.integer(optimal_max_size),
    nearest_hnsw_k = as.integer(hnsw_k),
    nearest_hnsw_M = as.integer(hnsw_M),
    nearest_hnsw_ef_construction = as.integer(hnsw_ef_construction),
    nearest_hnsw_ef_search = as.integer(hnsw_ef_search)
  )
}

.nearest_distance <- function(x_mat, t_row, distance, feature_weights, inv_cov) {
  delta <- sweep(x_mat, 2L, t_row, "-", check.margin = FALSE)
  if (distance == "weighted") {
    return(rowSums((delta * delta) * matrix(feature_weights, nrow = nrow(delta), ncol = ncol(delta), byrow = TRUE)))
  }
  if (distance == "mahalanobis") {
    proj <- delta %*% inv_cov
    return(rowSums(proj * delta))
  }
  rowSums(delta * delta)
}

.nearest_match_indices <- function(x_scaled,
                                   target_scaled,
                                   replace = FALSE,
                                   distance = c("euclidean", "mahalanobis", "weighted"),
                                   backend = c("auto", "exact", "hnsw"),
                                   match = c("greedy", "optimal"),
                                   feature_weights = NULL,
                                   inv_cov = NULL,
                                   caliper = NULL,
                                   caliper_strict = FALSE,
                                   hnsw_k = 64L,
                                   hnsw_M = 16L,
                                   hnsw_ef_construction = 200L,
                                   hnsw_ef_search = 64L) {
  distance <- match.arg(distance)
  backend <- match.arg(backend)
  match <- match.arg(match)
  check_scalar_flag(replace, "nearest_replace")
  check_scalar_flag(caliper_strict, "nearest_caliper_strict")
  hnsw_k <- check_positive_int_scalar(hnsw_k, "nearest_hnsw_k")
  hnsw_M <- check_positive_int_scalar(hnsw_M, "nearest_hnsw_M")
  hnsw_ef_construction <- check_positive_int_scalar(hnsw_ef_construction, "nearest_hnsw_ef_construction")
  hnsw_ef_search <- check_positive_int_scalar(hnsw_ef_search, "nearest_hnsw_ef_search")
  n_pool <- nrow(x_scaled)
  n_take <- nrow(target_scaled)
  if (!isTRUE(replace) && n_take > n_pool) {
    stop("`size` cannot exceed available rows when `nearest_replace = FALSE`.", call. = FALSE)
  }

  caliper_violations <- 0L
  backend_fallbacks <- 0L
  used_backend <- backend

  if (isTRUE(replace)) {
    selected <- integer(n_take)

    hnsw_index <- NULL
    x_backend <- NULL
    if (backend == "hnsw") {
      x_backend <- .nearest_distance_transform(x_scaled, distance, feature_weights, inv_cov)
      hnsw_index_try <- try(
        .nearest_hnsw_build(
          x_backend,
          M = hnsw_M,
          ef_construction = hnsw_ef_construction
        ),
        silent = TRUE
      )
      if (!inherits(hnsw_index_try, "try-error")) {
        hnsw_index <- hnsw_index_try
      } else {
        used_backend <- "exact"
        backend_fallbacks <- backend_fallbacks + 1L
      }
    }

    for (i in seq_len(n_take)) {
      t_row <- target_scaled[i, ]

      if (!is.null(hnsw_index)) {
        q_backend <- matrix(.nearest_distance_transform(matrix(t_row, nrow = 1L), distance, feature_weights, inv_cov), nrow = 1L)
        k_eff <- min(n_pool, max(1L, as.integer(hnsw_k)))
        idx_try <- try(
          .nearest_hnsw_query(hnsw_index, q_backend, k = k_eff, ef_search = hnsw_ef_search),
          silent = TRUE
        )
        if (!inherits(idx_try, "try-error")) {
          cand <- unique(as.integer(idx_try[1L, ]))
          cand <- cand[cand >= 1L & cand <= n_pool]
          if (length(cand) > 0L) {
            x_cand <- x_scaled[cand, , drop = FALSE]
            d_cand <- .nearest_distance(x_cand, t_row, distance, feature_weights, inv_cov)
            cal_sub <- .nearest_apply_caliper(x_cand, t_row, caliper, caliper_strict = FALSE)
            if (isTRUE(cal_sub$violated)) {
              caliper_violations <- caliper_violations + 1L
            } else {
              d_cand[!cal_sub$mask] <- Inf
            }
            if (any(is.finite(d_cand))) {
              selected[i] <- cand[which.min(d_cand)]
              next
            }
          }
        }
        backend_fallbacks <- backend_fallbacks + 1L
      }

      d <- .nearest_distance(x_scaled, t_row, distance, feature_weights, inv_cov)
      cal <- .nearest_apply_caliper(x_scaled, t_row, caliper, caliper_strict = caliper_strict)
      if (isTRUE(cal$violated)) {
        caliper_violations <- caliper_violations + 1L
      } else {
        d[!cal$mask] <- Inf
      }
      selected[i] <- which.min(d)
    }
    return(list(
      selected = selected,
      caliper_violations = caliper_violations,
      backend_fallbacks = backend_fallbacks,
      backend_used = used_backend,
      match_used = "greedy"
    ))
  }

  if (match == "optimal") {
    if (!requireNamespace("clue", quietly = TRUE)) {
      warning("`nearest_match = \"optimal\"` requested but `clue` is not installed; using greedy matching.", call. = FALSE)
      match <- "greedy"
    }
  }

  if (match == "optimal") {
    cost <- matrix(0, nrow = n_take, ncol = n_pool)
    for (i in seq_len(n_take)) {
      t_row <- target_scaled[i, ]
      d <- .nearest_distance(x_scaled, t_row, distance, feature_weights, inv_cov)
      cal <- .nearest_apply_caliper(x_scaled, t_row, caliper, caliper_strict = caliper_strict)
      if (isTRUE(cal$violated)) {
        caliper_violations <- caliper_violations + 1L
      } else {
        d[!cal$mask] <- Inf
      }
      cost[i, ] <- d
    }

    finite_vals <- cost[is.finite(cost)]
    if (length(finite_vals) < 1L) {
      stop("No finite distances available for optimal matching.", call. = FALSE)
    }
    penalty <- max(finite_vals) * 1e6 + 1
    cost[!is.finite(cost)] <- penalty

    selected <- .nearest_optimal_assign(cost)
    return(list(
      selected = selected,
      caliper_violations = caliper_violations,
      backend_fallbacks = backend_fallbacks,
      backend_used = "exact",
      match_used = "optimal"
    ))
  }

  hnsw_index <- NULL
  x_backend <- NULL
  if (backend == "hnsw") {
    x_backend <- .nearest_distance_transform(x_scaled, distance, feature_weights, inv_cov)
    hnsw_index_try <- try(
      .nearest_hnsw_build(
        x_backend,
        M = hnsw_M,
        ef_construction = hnsw_ef_construction
      ),
      silent = TRUE
    )
    if (!inherits(hnsw_index_try, "try-error")) {
      hnsw_index <- hnsw_index_try
    } else {
      used_backend <- "exact"
      backend_fallbacks <- backend_fallbacks + 1L
    }
  }

  available <- seq_len(n_pool)
  selected <- integer(n_take)

  for (i in seq_len(n_take)) {
    t_row <- target_scaled[i, ]

    if (!is.null(hnsw_index)) {
      q_backend <- matrix(.nearest_distance_transform(matrix(t_row, nrow = 1L), distance, feature_weights, inv_cov), nrow = 1L)
      k_eff <- min(n_pool, max(1L, as.integer(hnsw_k)))
      idx_try <- try(
        .nearest_hnsw_query(hnsw_index, q_backend, k = k_eff, ef_search = hnsw_ef_search),
        silent = TRUE
      )
      if (!inherits(idx_try, "try-error")) {
        cand <- unique(as.integer(idx_try[1L, ]))
        cand <- cand[cand >= 1L & cand <= n_pool]
        cand <- cand[cand %in% available]

        if (length(cand) > 0L) {
          x_cand <- x_scaled[cand, , drop = FALSE]
          d_cand <- .nearest_distance(x_cand, t_row, distance, feature_weights, inv_cov)
          cal_sub <- .nearest_apply_caliper(x_cand, t_row, caliper, caliper_strict = FALSE)
          if (isTRUE(cal_sub$violated)) {
            caliper_violations <- caliper_violations + 1L
          } else {
            d_cand[!cal_sub$mask] <- Inf
          }
          if (any(is.finite(d_cand))) {
            pick <- cand[which.min(d_cand)]
            selected[i] <- pick
            available <- available[available != pick]
            next
          }
        }
      }
      backend_fallbacks <- backend_fallbacks + 1L
    }

    x_avail <- x_scaled[available, , drop = FALSE]
    d <- .nearest_distance(x_avail, t_row, distance, feature_weights, inv_cov)
    cal <- .nearest_apply_caliper(x_avail, t_row, caliper, caliper_strict = caliper_strict)
    if (isTRUE(cal$violated)) {
      caliper_violations <- caliper_violations + 1L
    } else {
      d[!cal$mask] <- Inf
    }
    k <- which.min(d)
    selected[i] <- available[k]
    available <- available[-k]
  }

  list(
    selected = selected,
    caliper_violations = caliper_violations,
    backend_fallbacks = backend_fallbacks,
    backend_used = used_backend,
    match_used = "greedy"
  )
}

.nearest_eval_seed <- function(data, size, seed, prep, objective = NULL, replace = FALSE) {
  check_scalar_flag(replace, "nearest_replace")
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

  match_out <- .nearest_match_indices(
    x_scaled = prep$x_scaled,
    target_scaled = target_scaled,
    replace = replace,
    distance = prep$distance,
    backend = prep$nearest_backend,
    match = prep$nearest_match,
    feature_weights = prep$feature_weights,
    inv_cov = prep$inv_cov,
    caliper = prep$caliper,
    caliper_strict = prep$caliper_strict,
    hnsw_k = prep$nearest_hnsw_k,
    hnsw_M = prep$nearest_hnsw_M,
    hnsw_ef_construction = prep$nearest_hnsw_ef_construction,
    hnsw_ef_search = prep$nearest_hnsw_ef_search
  )
  pick_local <- match_out$selected
  selected_rows <- prep$pool_idx[pick_local]
  chosen <- data[selected_rows, , drop = FALSE]

  ks <- .importance_ks_score(chosen, prep$cont, prep$mean, prep$sd, prep$dist)

  out <- data
  draw_counts <- tabulate(selected_rows, nbins = nrow(data))
  out$repsample <- as.integer(draw_counts > 0L)
  if (isTRUE(replace)) {
    out$repsample_n <- as.integer(draw_counts)
  }

  fit <- structure(
    list(
      data = out,
      r = numeric(0),
      selected_rows = as.integer(selected_rows),
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
        cont_vars = prep$cont,
        bincat_vars = character(0),
        target_mean = as.numeric(prep$mean),
        target_sd = as.numeric(prep$sd),
        target_perc = numeric(0),
        cont_dist = if (length(prep$cont) > 0L) stats::setNames(as.character(prep$dist), prep$cont) else character(0),
        exact = FALSE,
        weighted = FALSE,
        retain = NULL,
        method = "nearest_neighbor",
        ks = as.numeric(ks),
        nearest_replace = isTRUE(replace),
        selected_unique = as.integer(sum(draw_counts > 0L)),
        nearest_distance = as.character(prep$distance),
        nearest_backend = as.character(match_out$backend_used),
        nearest_backend_requested = as.character(prep$nearest_backend_requested),
        nearest_match = as.character(match_out$match_used),
        nearest_feature_weights = as.numeric(prep$feature_weights),
        nearest_caliper = prep$caliper,
        nearest_caliper_strict = isTRUE(prep$caliper_strict),
        nearest_caliper_violations = as.integer(match_out$caliper_violations),
        nearest_backend_fallbacks = as.integer(match_out$backend_fallbacks),
        nearest_hnsw_k = as.integer(prep$nearest_hnsw_k),
        nearest_hnsw_ef_search = as.integer(prep$nearest_hnsw_ef_search)
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
                                seed = 7L,
                                replace = FALSE,
                                distance = c("euclidean", "mahalanobis", "weighted"),
                                backend = c("auto", "exact", "hnsw"),
                                match = c("greedy", "optimal"),
                                feature_weights = NULL,
                                caliper = NULL,
                                caliper_strict = FALSE,
                                optimal_max_size = 2000L,
                                hnsw_k = 64L,
                                hnsw_M = 16L,
                                hnsw_ef_construction = 200L,
                                hnsw_ef_search = 64L) {
  prep <- .nearest_prepare(
    data = data,
    size = size,
    cont = cont,
    mean = mean,
    sd = sd,
    dist = dist,
    replace = replace,
    distance = distance,
    backend = backend,
    match = match,
    feature_weights = feature_weights,
    caliper = caliper,
    caliper_strict = caliper_strict,
    optimal_max_size = optimal_max_size,
    hnsw_k = hnsw_k,
    hnsw_M = hnsw_M,
    hnsw_ef_construction = hnsw_ef_construction,
    hnsw_ef_search = hnsw_ef_search
  )

  .importance_with_seed(
    seed,
    .nearest_eval_seed(data, size, seed, prep, objective = NULL, replace = replace)$fit
  )
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
                                   objective = NULL,
                                   replace = FALSE,
                                   distance = c("euclidean", "mahalanobis", "weighted"),
                                   backend = c("auto", "exact", "hnsw"),
                                   match = c("greedy", "optimal"),
                                   feature_weights = NULL,
                                   caliper = NULL,
                                   caliper_strict = FALSE,
                                   optimal_max_size = 2000L,
                                   hnsw_k = 64L,
                                   hnsw_M = 16L,
                                   hnsw_ef_construction = 200L,
                                   hnsw_ef_search = 64L,
                                   mean_tol = NULL,
                                   sd_tol = NULL,
                                   ks_tol = NULL,
                                   perc_tol = NULL) {
  outer_mode <- match.arg(outer_mode)
  distance <- match.arg(distance)
  backend <- match.arg(backend)
  match <- match.arg(match)
  check_scalar_flag(replace, "nearest_replace")
  if (length(seeds) < 1L) {
    stop("`seeds` must be non-empty.", call. = FALSE)
  }

  prep <- .nearest_prepare(
    data = data,
    size = size,
    cont = cont,
    mean = mean,
    sd = sd,
    dist = dist,
    replace = replace,
    distance = distance,
    backend = backend,
    match = match,
    feature_weights = feature_weights,
    caliper = caliper,
    caliper_strict = caliper_strict,
    optimal_max_size = optimal_max_size,
    hnsw_k = hnsw_k,
    hnsw_M = hnsw_M,
    hnsw_ef_construction = hnsw_ef_construction,
    hnsw_ef_search = hnsw_ef_search
  )

  eval_seed <- function(seed) {
    res <- .importance_with_seed(
      seed,
      .nearest_eval_seed(
        data,
        size,
        seed,
        prep,
        objective = objective,
        replace = replace
      )
    )
    tol <- .fit_tolerance_status(
      out = res$fit,
      cont = prep$cont,
      mean = prep$mean,
      sd = prep$sd,
      mean_tol = mean_tol,
      sd_tol = sd_tol,
      ks_tol = ks_tol,
      perc_tol = perc_tol
    )
    res$meets_tolerance <- isTRUE(tol$enabled) && isTRUE(tol$ok)
    res$tolerance <- tol
    res
  }

  has_tol <- any(!vapply(list(mean_tol, sd_tol, ks_tol, perc_tol), is.null, logical(1)))
  results <- NULL
  if (outer_mode == "serial" || length(seeds) == 1L || n_outer_workers <= 1L) {
    if (!has_tol) {
      results <- lapply(seeds, eval_seed)
    } else {
      results <- list()
      for (seed in seeds) {
        item <- eval_seed(seed)
        results[[length(results) + 1L]] <- item
        if (isTRUE(item$meets_tolerance)) {
          break
        }
      }
    }
  } else if (outer_mode == "multicore") {
    worker_n <- min(as.integer(n_outer_workers), length(seeds))
    if (!has_tol) {
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
      results <- list()
      for (start_i in seq.int(1L, length(seeds), by = worker_n)) {
        idx <- start_i:min(length(seeds), start_i + worker_n - 1L)
        batch <- parallel::mclapply(
          seeds[idx],
          eval_seed,
          mc.cores = min(worker_n, length(idx)),
          mc.preschedule = TRUE,
          mc.set.seed = FALSE,
          mc.cleanup = FALSE,
          mc.allow.recursive = FALSE
        )
        results <- c(results, batch)
        if (any(vapply(batch, function(x) isTRUE(x$meets_tolerance), logical(1)))) {
          break
        }
      }
    }
  } else {
    worker_n <- min(as.integer(n_outer_workers), length(seeds))
    cl <- parallel::makeCluster(worker_n)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    libs <- .libPaths()
    parallel::clusterCall(cl, function(x) .libPaths(x), libs)
    parallel::clusterExport(cl, varlist = c("eval_seed"), envir = environment())
    if (!has_tol) {
      results <- parallel::parLapply(cl, seeds, eval_seed)
    } else {
      results <- list()
      for (start_i in seq.int(1L, length(seeds), by = worker_n)) {
        idx <- start_i:min(length(seeds), start_i + worker_n - 1L)
        batch <- parallel::parLapply(cl, seeds[idx], eval_seed)
        results <- c(results, batch)
        if (any(vapply(batch, function(x) isTRUE(x$meets_tolerance), logical(1)))) {
          break
        }
      }
    }
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
        n_seeds = nrow(summary_df),
        n_outer_workers = as.integer(n_outer_workers),
        outer_mode = outer_mode,
        backend = "nearest_neighbor",
        theoretical = TRUE,
        method = "nearest_neighbor",
        nearest_replace = isTRUE(replace),
        nearest_distance = as.character(distance),
        nearest_backend = if (!is.null(best$fit$meta$nearest_backend)) as.character(best$fit$meta$nearest_backend) else as.character(prep$nearest_backend),
        nearest_backend_requested = as.character(prep$nearest_backend_requested),
        nearest_match = if (!is.null(best$fit$meta$nearest_match)) as.character(best$fit$meta$nearest_match) else as.character(prep$nearest_match),
        nearest_feature_weights = as.numeric(prep$feature_weights),
        nearest_caliper = prep$caliper,
        nearest_caliper_strict = isTRUE(prep$caliper_strict),
        nearest_optimal_max_size = as.integer(prep$nearest_optimal_max_size),
        nearest_hnsw_k = as.integer(prep$nearest_hnsw_k),
        nearest_hnsw_M = as.integer(prep$nearest_hnsw_M),
        nearest_hnsw_ef_construction = as.integer(prep$nearest_hnsw_ef_construction),
        nearest_hnsw_ef_search = as.integer(prep$nearest_hnsw_ef_search),
        stopped_early_by_tolerance = isTRUE(has_tol) &&
          any(vapply(results, function(x) isTRUE(x$meets_tolerance), logical(1)))
      )
    ),
    class = "repsample_search_result"
  )
}
