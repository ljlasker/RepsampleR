.normalize_tol <- function(x, n, name) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.numeric(x) || anyNA(x) || length(x) < 1L) {
    stop(sprintf("`%s` must be NULL or a numeric scalar/vector.", name), call. = FALSE)
  }
  if (any(x < 0)) {
    stop(sprintf("`%s` values must be non-negative.", name), call. = FALSE)
  }
  if (length(x) == 1L) {
    return(rep.int(as.numeric(x), n))
  }
  if (length(x) == n) {
    return(as.numeric(x))
  }
  stop(sprintf("`%s` must have length 1 or length %d.", name, n), call. = FALSE)
}

.repsample_selected_data <- function(out, include_duplicates = TRUE) {
  if (!is.list(out) || is.null(out$data) || !is.data.frame(out$data)) {
    stop("`out` must be a valid `repsample_result` object.", call. = FALSE)
  }

  d <- out$data
  has_rep_n <- ("repsample_n" %in% names(d))
  with_replacement <- !is.null(out$meta$nearest_replace) && isTRUE(out$meta$nearest_replace)

  if (isTRUE(include_duplicates) && has_rep_n && with_replacement) {
    counts <- as.integer(d$repsample_n)
    counts[!is.finite(counts) | counts < 0L] <- 0L
    idx <- rep.int(seq_len(nrow(d)), counts)
    if (length(idx) > 0L) {
      return(d[idx, , drop = FALSE])
    }
  }

  if (isTRUE(include_duplicates) &&
      with_replacement &&
      !is.null(out$selected_rows) &&
      length(out$selected_rows) > 0L) {
    idx <- as.integer(out$selected_rows)
    idx <- idx[idx >= 1L & idx <= nrow(d)]
    if (length(idx) > 0L) {
      return(d[idx, , drop = FALSE])
    }
  }

  if (!("repsample" %in% names(d))) {
    return(d[0, , drop = FALSE])
  }
  d[d$repsample == 1L, , drop = FALSE]
}

.extract_ks_metric <- function(out, cont_norm) {
  if (!is.null(out$meta$ks) && is.finite(as.numeric(out$meta$ks))) {
    return(as.numeric(out$meta$ks))
  }
  if (length(cont_norm) == 0L || length(out$r) == 0L) {
    return(NA_real_)
  }
  d_names <- paste0(cont_norm, "_D")
  d_vals <- suppressWarnings(as.numeric(out$r[d_names]))
  d_vals <- d_vals[is.finite(d_vals)]
  if (length(d_vals) == 0L) {
    return(NA_real_)
  }
  mean(d_vals)
}

.repsample_metrics <- function(out,
                               cont = NULL,
                               bincat = NULL,
                               mean = NULL,
                               sd = NULL,
                               perc = NULL) {
  cont_norm <- normalize_varlist(cont)
  bincat_norm <- normalize_varlist(bincat)
  mean_norm <- normalize_numlist(mean)
  sd_norm <- normalize_numlist(sd)
  perc_norm <- normalize_numlist(perc)

  chosen <- .repsample_selected_data(out, include_duplicates = TRUE)
  if (nrow(chosen) == 0L) {
    return(list(
      n_selected = 0L,
      cont = data.frame(
        variable = cont_norm,
        observed_mean = rep(NA_real_, length(cont_norm)),
        target_mean = if (length(mean_norm) == length(cont_norm)) mean_norm else rep(NA_real_, length(cont_norm)),
        mean_abs_err = rep(NA_real_, length(cont_norm)),
        observed_sd = rep(NA_real_, length(cont_norm)),
        target_sd = if (length(sd_norm) == length(cont_norm)) sd_norm else rep(NA_real_, length(cont_norm)),
        sd_abs_err = rep(NA_real_, length(cont_norm)),
        stringsAsFactors = FALSE
      ),
      bincat = data.frame(
        variable = bincat_norm,
        observed_perc = rep(NA_real_, length(bincat_norm)),
        target_perc = if (length(perc_norm) == length(bincat_norm)) perc_norm else rep(NA_real_, length(bincat_norm)),
        perc_abs_err = rep(NA_real_, length(bincat_norm)),
        stringsAsFactors = FALSE
      ),
      ks = NA_real_
    ))
  }

  cont_df <- data.frame(
    variable = cont_norm,
    observed_mean = rep(NA_real_, length(cont_norm)),
    target_mean = rep(NA_real_, length(cont_norm)),
    mean_abs_err = rep(NA_real_, length(cont_norm)),
    observed_sd = rep(NA_real_, length(cont_norm)),
    target_sd = rep(NA_real_, length(cont_norm)),
    sd_abs_err = rep(NA_real_, length(cont_norm)),
    stringsAsFactors = FALSE
  )
  if (length(cont_norm) > 0L) {
    for (i in seq_along(cont_norm)) {
      x <- as.numeric(chosen[[cont_norm[i]]])
      cont_df$observed_mean[i] <- mean(x)
      cont_df$observed_sd[i] <- stats::sd(x)
      if (length(mean_norm) == length(cont_norm)) {
        cont_df$target_mean[i] <- mean_norm[i]
        cont_df$mean_abs_err[i] <- abs(cont_df$observed_mean[i] - mean_norm[i])
      }
      if (length(sd_norm) == length(cont_norm)) {
        cont_df$target_sd[i] <- sd_norm[i]
        cont_df$sd_abs_err[i] <- abs(cont_df$observed_sd[i] - sd_norm[i])
      }
    }
  }

  bin_df <- data.frame(
    variable = bincat_norm,
    observed_perc = rep(NA_real_, length(bincat_norm)),
    target_perc = rep(NA_real_, length(bincat_norm)),
    perc_abs_err = rep(NA_real_, length(bincat_norm)),
    stringsAsFactors = FALSE
  )
  if (length(bincat_norm) > 0L) {
    for (i in seq_along(bincat_norm)) {
      vals <- chosen[[bincat_norm[i]]]
      obs <- 100 * mean(vals == 1, na.rm = TRUE)
      bin_df$observed_perc[i] <- obs
      if (length(perc_norm) == length(bincat_norm)) {
        bin_df$target_perc[i] <- perc_norm[i]
        bin_df$perc_abs_err[i] <- abs(obs - perc_norm[i])
      }
    }
  }

  list(
    n_selected = nrow(chosen),
    cont = cont_df,
    bincat = bin_df,
    ks = .extract_ks_metric(out, cont_norm)
  )
}

#' Compute Match Quality Metrics for a RepsampleR Fit
#'
#' Returns a compact set of fit-quality diagnostics and a scalar composite loss
#' that is comparable across methods (`greedy`, `importance`, `nearest`).
#'
#' @param out A `repsample_result` object.
#' @param cont Continuous variable names used for matching.
#' @param bincat Binary variable names used for matching.
#' @param mean Target means for `cont`.
#' @param sd Target standard deviations for `cont`.
#' @param perc Target percentages for `bincat`.
#' @param ks_weight Weight on the KS component in the composite loss.
#'
#' @return A named list with components `loss`, `ks`, `cont`, and `bincat`.
#'
#' @export
repsample_quality <- function(out,
                              cont = NULL,
                              bincat = NULL,
                              mean = NULL,
                              sd = NULL,
                              perc = NULL,
                              ks_weight = 1) {
  if (!is.numeric(ks_weight) || length(ks_weight) != 1L || is.na(ks_weight) || ks_weight < 0) {
    stop("`ks_weight` must be a single non-negative numeric value.", call. = FALSE)
  }

  m <- .repsample_metrics(
    out = out,
    cont = cont,
    bincat = bincat,
    mean = mean,
    sd = sd,
    perc = perc
  )

  cont_loss <- 0
  if (nrow(m$cont) > 0L) {
    mean_err <- m$cont$mean_abs_err
    sd_err <- m$cont$sd_abs_err
    mean_err[!is.finite(mean_err)] <- 0
    sd_err[!is.finite(sd_err)] <- 0
    cont_loss <- sum(mean_err ^ 2 + sd_err ^ 2)
  }

  bin_loss <- 0
  if (nrow(m$bincat) > 0L) {
    p_err <- m$bincat$perc_abs_err
    p_err[!is.finite(p_err)] <- 0
    bin_loss <- sum((p_err / 100) ^ 2)
  }

  ks_term <- if (is.finite(m$ks)) (as.numeric(ks_weight) * (m$ks ^ 2)) else 0

  list(
    loss = as.numeric(cont_loss + bin_loss + ks_term),
    ks = as.numeric(m$ks),
    cont = m$cont,
    bincat = m$bincat,
    n_selected = m$n_selected
  )
}

.fit_tolerance_status <- function(out,
                                  cont = NULL,
                                  bincat = NULL,
                                  mean = NULL,
                                  sd = NULL,
                                  perc = NULL,
                                  mean_tol = NULL,
                                  sd_tol = NULL,
                                  ks_tol = NULL,
                                  perc_tol = NULL) {
  cont_norm <- normalize_varlist(cont)
  bincat_norm <- normalize_varlist(bincat)

  mean_tol_norm <- .normalize_tol(mean_tol, length(cont_norm), "mean_tol")
  sd_tol_norm <- .normalize_tol(sd_tol, length(cont_norm), "sd_tol")
  perc_tol_norm <- .normalize_tol(perc_tol, length(bincat_norm), "perc_tol")

  ks_tol_norm <- NULL
  if (!is.null(ks_tol)) {
    if (!is.numeric(ks_tol) || length(ks_tol) != 1L || is.na(ks_tol) || ks_tol < 0) {
      stop("`ks_tol` must be NULL or a single non-negative numeric value.", call. = FALSE)
    }
    ks_tol_norm <- as.numeric(ks_tol)
  }

  has_tol <- !is.null(mean_tol_norm) || !is.null(sd_tol_norm) || !is.null(ks_tol_norm) || !is.null(perc_tol_norm)
  metrics <- .repsample_metrics(out, cont = cont, bincat = bincat, mean = mean, sd = sd, perc = perc)
  if (!has_tol) {
    return(list(enabled = FALSE, ok = FALSE, metrics = metrics))
  }

  ok <- TRUE

  if (!is.null(mean_tol_norm) && length(cont_norm) > 0L) {
    errs <- metrics$cont$mean_abs_err
    errs[!is.finite(errs)] <- Inf
    ok <- ok && all(errs <= mean_tol_norm)
  }

  if (!is.null(sd_tol_norm) && length(cont_norm) > 0L) {
    errs <- metrics$cont$sd_abs_err
    errs[!is.finite(errs)] <- Inf
    ok <- ok && all(errs <= sd_tol_norm)
  }

  if (!is.null(perc_tol_norm) && length(bincat_norm) > 0L) {
    errs <- metrics$bincat$perc_abs_err
    errs[!is.finite(errs)] <- Inf
    ok <- ok && all(errs <= perc_tol_norm)
  }

  if (!is.null(ks_tol_norm)) {
    ks <- metrics$ks
    if (!is.finite(ks)) {
      ok <- FALSE
    } else {
      ok <- ok && (ks <= ks_tol_norm)
    }
  }

  list(
    enabled = TRUE,
    ok = isTRUE(ok),
    metrics = metrics,
    mean_tol = mean_tol_norm,
    sd_tol = sd_tol_norm,
    ks_tol = ks_tol_norm,
    perc_tol = perc_tol_norm
  )
}

#' Summarize a RepsampleR Fit
#'
#' @param object A `repsample_result` object.
#' @param ... Unused.
#'
#' @return A list with selected-size and per-variable matching diagnostics.
#'
#' @export
summary.repsample_result <- function(object, ...) {
  cont <- if (!is.null(object$meta$cont_vars)) object$meta$cont_vars else NULL
  bincat <- if (!is.null(object$meta$bincat_vars)) object$meta$bincat_vars else NULL
  mean <- if (!is.null(object$meta$target_mean)) object$meta$target_mean else NULL
  sd <- if (!is.null(object$meta$target_sd)) object$meta$target_sd else NULL
  perc <- if (!is.null(object$meta$target_perc)) object$meta$target_perc else NULL

  m <- .repsample_metrics(
    out = object,
    cont = cont,
    bincat = bincat,
    mean = mean,
    sd = sd,
    perc = perc
  )

  structure(
    list(
      n_selected = m$n_selected,
      cont = m$cont,
      bincat = m$bincat,
      ks = m$ks
    ),
    class = "summary.repsample_result"
  )
}

#' @export
print.summary.repsample_result <- function(x, ...) {
  cat("summary.repsample_result\n")
  cat(sprintf("  Selected rows: %d\n", as.integer(x$n_selected)))
  if (nrow(x$cont) > 0L) {
    cat("  Continuous diagnostics:\n")
    print(x$cont, row.names = FALSE)
  }
  if (nrow(x$bincat) > 0L) {
    cat("  Binary diagnostics:\n")
    print(x$bincat, row.names = FALSE)
  }
  if (is.finite(x$ks)) {
    cat(sprintf("  Mean KS: %.6f\n", x$ks))
  }
  invisible(x)
}

#' Summarize a RepsampleR Search Result
#'
#' @param object A `repsample_search_result` object.
#' @param ... Unused.
#'
#' @return A list with top-level search metadata and best-fit summary.
#'
#' @export
summary.repsample_search_result <- function(object, ...) {
  best_summary <- summary(object$best)
  structure(
    list(
      best_seed = object$best_seed,
      best_loss = object$best_loss,
      seeds_searched = object$meta$n_seeds,
      outer_mode = object$meta$outer_mode,
      best = best_summary
    ),
    class = "summary.repsample_search_result"
  )
}

#' @export
print.summary.repsample_search_result <- function(x, ...) {
  cat("summary.repsample_search_result\n")
  cat(sprintf("  Best seed: %d\n", as.integer(x$best_seed)))
  cat(sprintf("  Best loss: %.6g\n", as.numeric(x$best_loss)))
  cat(sprintf("  Seeds searched: %d\n", as.integer(x$seeds_searched)))
  cat(sprintf("  Outer mode: %s\n", as.character(x$outer_mode)))
  print(x$best)
  invisible(x)
}

#' Plot a RepsampleR Fit
#'
#' Produces a simple density comparison for one continuous variable.
#'
#' @param x A `repsample_result`.
#' @param var Optional continuous variable to plot. Defaults to the first
#'   variable in `x$meta$cont_vars`.
#' @param ... Graphical parameters passed to `plot()`.
#'
#' @export
plot.repsample_result <- function(x, var = NULL, ...) {
  cont <- if (!is.null(x$meta$cont_vars)) x$meta$cont_vars else character(0)
  if (is.null(var)) {
    if (length(cont) < 1L) {
      stop("No continuous variable available to plot.", call. = FALSE)
    }
    var <- cont[1L]
  }
  if (!is.character(var) || length(var) != 1L || !(var %in% names(x$data))) {
    stop("`var` must name a column in `x$data`.", call. = FALSE)
  }

  selected <- .repsample_selected_data(x, include_duplicates = TRUE)
  if (nrow(selected) < 2L) {
    stop("Not enough selected rows to plot a density.", call. = FALSE)
  }

  x_all <- as.numeric(x$data[[var]])
  x_sel <- as.numeric(selected[[var]])

  d_all <- stats::density(x_all)
  d_sel <- stats::density(x_sel)

  plot(d_all,
       col = "grey40",
       lwd = 2,
       main = sprintf("RepsampleR Density: %s", var),
       xlab = var,
       ylab = "Density",
       ...)
  lines(d_sel, col = "#1f77b4", lwd = 2)

  if (!is.null(x$meta$target_mean) &&
      !is.null(x$meta$target_sd) &&
      length(cont) > 0L &&
      var %in% cont) {
    idx <- match(var, cont)
    mu <- as.numeric(x$meta$target_mean[idx])
    sigma <- as.numeric(x$meta$target_sd[idx])
    if (is.finite(mu) && is.finite(sigma) && sigma > 0) {
      curve(stats::dnorm(x, mean = mu, sd = sigma),
            add = TRUE,
            col = "#d62728",
            lwd = 2,
            lty = 2)
      legend("topright",
             legend = c("baseline", "selected", "target"),
             col = c("grey40", "#1f77b4", "#d62728"),
             lwd = c(2, 2, 2),
             lty = c(1, 1, 2),
             bty = "n")
      return(invisible(NULL))
    }
  }

  legend("topright",
         legend = c("baseline", "selected"),
         col = c("grey40", "#1f77b4"),
         lwd = c(2, 2),
         bty = "n")
  invisible(NULL)
}
