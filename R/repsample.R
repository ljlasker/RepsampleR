stata_round <- function(x) {
  ifelse(x >= 0, floor(x + 0.5), ceiling(x - 0.5))
}

normalize_varlist <- function(x) {
  if (is.null(x)) {
    return(character(0))
  }

  if (is.character(x) && length(x) == 1L) {
    x <- trimws(x)
    if (!nzchar(x)) {
      return(character(0))
    }
    x <- strsplit(x, "\\s+")[[1]]
  }

  out <- as.character(x)
  out <- out[nzchar(out)]
  out
}

normalize_numlist <- function(x) {
  if (is.null(x)) {
    return(numeric(0))
  }

  if (is.character(x) && length(x) == 1L) {
    x <- trimws(x)
    if (!nzchar(x)) {
      return(numeric(0))
    }
    x <- strsplit(x, "\\s+")[[1]]
  }

  out <- as.numeric(x)
  if (anyNA(out)) {
    stop("Numeric lists contain non-numeric values.", call. = FALSE)
  }

  out
}

normalize_cont_dist <- function(dist, cont) {
  ccnt <- length(cont)
  if (ccnt == 0L) {
    return(character(0))
  }

  if (is.null(dist)) {
    return(rep.int("normal", ccnt))
  }

  normalize_one <- function(x) {
    x <- tolower(trimws(as.character(x)))
    if (!nzchar(x)) {
      stop("`dist` values must be non-empty strings.", call. = FALSE)
    }
    aliases <- c(
      normal = "normal",
      gaussian = "normal",
      gauss = "normal",
      norm = "normal",
      lognormal = "lognormal",
      "log-normal" = "lognormal",
      lnorm = "lognormal",
      poisson = "poisson",
      pois = "poisson"
    )
    if (!x %in% names(aliases)) {
      stop(
        sprintf(
          "Unsupported distribution `%s`. Supported distributions: normal, lognormal, poisson.",
          x
        ),
        call. = FALSE
      )
    }
    aliases[[x]]
  }

  out <- NULL
  if (is.character(dist) || is.factor(dist)) {
    d <- as.character(dist)
    if (length(d) == 1L) {
      out <- rep.int(normalize_one(d), ccnt)
    } else if (length(d) == ccnt && is.null(names(d))) {
      out <- vapply(d, normalize_one, character(1))
    } else if (!is.null(names(d))) {
      out <- rep.int("normal", ccnt)
      names(out) <- cont
      unknown <- setdiff(names(d), cont)
      if (length(unknown) > 0L) {
        stop(
          sprintf("`dist` has unknown named variables: %s", paste(unknown, collapse = ", ")),
          call. = FALSE
        )
      }
      for (nm in names(d)) {
        out[nm] <- normalize_one(d[[nm]])
      }
      out <- unname(out)
    } else {
      stop("`dist` must be length 1, length(cont), or a named vector.", call. = FALSE)
    }
  } else if (is.list(dist)) {
    if (length(dist) == 1L) {
      out <- rep.int(normalize_one(dist[[1L]]), ccnt)
    } else if (length(dist) == ccnt && is.null(names(dist))) {
      out <- vapply(dist, normalize_one, character(1))
    } else if (!is.null(names(dist))) {
      out <- rep.int("normal", ccnt)
      names(out) <- cont
      unknown <- setdiff(names(dist), cont)
      if (length(unknown) > 0L) {
        stop(
          sprintf("`dist` has unknown named variables: %s", paste(unknown, collapse = ", ")),
          call. = FALSE
        )
      }
      for (nm in names(dist)) {
        out[nm] <- normalize_one(dist[[nm]])
      }
      out <- unname(out)
    } else {
      stop("`dist` list must be length 1, length(cont), or a named list.", call. = FALSE)
    }
  } else {
    stop("`dist` must be NULL, character/factor, or list.", call. = FALSE)
  }

  unname(out)
}

quality_defaults <- function(quality, size) {
  quality <- match.arg(quality, c("manual", "fast", "balanced", "strict"))

  if (quality == "manual") {
    return(list(randomperc = NULL, rrule = NULL, srule = NULL))
  }

  if (quality == "fast") {
    return(list(
      randomperc = 70,
      rrule = as.integer(max(50L, min(300L, round(size * 0.15)))),
      srule = 0.70
    ))
  }

  if (quality == "balanced") {
    return(list(
      randomperc = 35,
      rrule = as.integer(max(150L, min(900L, round(size * 0.45)))),
      srule = 0.82
    ))
  }

  list(
    randomperc = 5,
    rrule = as.integer(max(500L, min(4000L, round(size * 1.25)))),
    srule = NULL
  )
}

ks_p_asymptotic <- function(d, n_eff) {
  if (!is.finite(d) || d <= 0) {
    return(1)
  }
  if (d >= 1) {
    return(0)
  }

  lambda <- sqrt(n_eff) * d
  j <- seq_len(100L)
  p <- 2 * sum(((-1) ^ (j - 1L)) * exp(-2 * (j ^ 2) * (lambda ^ 2)))
  p <- max(min(p, 1), 0)
  p
}

chisq_from_counts <- function(sample_counts, pop_counts) {
  if (length(sample_counts) <= 1L) {
    return(list(chi2 = 0, p = 1))
  }

  n1 <- sum(sample_counts)
  n0 <- sum(pop_counts)
  grand <- n1 + n0

  row_totals <- sample_counts + pop_counts
  exp1 <- row_totals * (n1 / grand)
  exp0 <- row_totals * (n0 / grand)

  term1 <- ifelse(exp1 > 0, (sample_counts - exp1) ^ 2 / exp1, 0)
  term0 <- ifelse(exp0 > 0, (pop_counts - exp0) ^ 2 / exp0, 0)

  chi2 <- sum(term1 + term0)
  df <- length(sample_counts) - 1L
  p <- stats::pchisq(chi2, df = df, lower.tail = FALSE)

  list(chi2 = chi2, p = p)
}

check_scalar_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
}

cuda_prepare_windows_dll_search <- function(py_run_string) {
  if (.Platform$OS.type != "windows") {
    return(invisible(NULL))
  }

  py_run_string(
    paste(
      "import os",
      "import sys",
      "",
      "if os.name == 'nt' and hasattr(os, 'add_dll_directory'):",
      "    _added = set()",
      "",
      "    def _add_dir(path):",
      "        if not path:",
      "            return",
      "        if os.path.isdir(path) and path not in _added:",
      "            try:",
      "                os.add_dll_directory(path)",
      "                _added.add(path)",
      "            except (FileNotFoundError, OSError):",
      "                pass",
      "",
      "    for _env in ('CUDA_PATH', 'CUDA_HOME', 'CONDA_PREFIX'):",
      "        _root = os.environ.get(_env, '')",
      "        if _root:",
      "            _add_dir(os.path.join(_root, 'bin'))",
      "",
      "    for _base in list(sys.path):",
      "        if not _base:",
      "            continue",
      "        _nvidia = os.path.join(_base, 'nvidia')",
      "        _add_dir(os.path.join(_nvidia, 'cuda_runtime', 'bin'))",
      "        _add_dir(os.path.join(_nvidia, 'cuda_nvrtc', 'bin'))",
      "        _add_dir(os.path.join(_nvidia, 'cublas', 'bin'))",
      sep = "\n"
    ),
    local = FALSE
  )

  invisible(NULL)
}

cuda_backend_available <- function() {
  reticulate_path <- suppressWarnings(
    tryCatch(base::find.package("reticulate", quiet = TRUE), error = function(e) character(0))
  )
  has_reticulate <- is.character(reticulate_path) &&
    length(reticulate_path) >= 1L &&
    any(nzchar(reticulate_path))
  if (!isTRUE(has_reticulate)) {
    return(FALSE)
  }
  rns <- asNamespace("reticulate")
  py_module_available <- get("py_module_available", envir = rns)
  py_import <- get("import", envir = rns)
  py_run_string <- get("py_run_string", envir = rns)

  try(cuda_prepare_windows_dll_search(py_run_string), silent = TRUE)

  if (!py_module_available("cupy")) {
    return(FALSE)
  }

  ok <- FALSE
  try({
    cp <- py_import("cupy", delay_load = FALSE)
    n_dev <- as.integer(cp$cuda$runtime$getDeviceCount())
    ok <- is.finite(n_dev) && n_dev > 0L
  }, silent = TRUE)

  ok
}

cuda_ks1_batch <- local({
  py_fun <- NULL

  function(selected_sorted, candidates, mean, sd) {
    reticulate_path <- suppressWarnings(
      tryCatch(base::find.package("reticulate", quiet = TRUE), error = function(e) character(0))
    )
    has_reticulate <- is.character(reticulate_path) &&
      length(reticulate_path) >= 1L &&
      any(nzchar(reticulate_path))
    if (!isTRUE(has_reticulate)) {
      stop("CUDA backend requires the `reticulate` package.", call. = FALSE)
    }
    rns <- asNamespace("reticulate")
    py_module_available <- get("py_module_available", envir = rns)
    py_run_string <- get("py_run_string", envir = rns)
    py_to_r <- get("py_to_r", envir = rns)

    try(cuda_prepare_windows_dll_search(py_run_string), silent = TRUE)

    if (is.null(py_fun)) {
      if (!py_module_available("cupy")) {
        stop("CUDA backend requires Python package `cupy`.", call. = FALSE)
      }

      tryCatch(
        py_run_string(
          paste(
            "import cupy as cp",
            "",
            "try:",
            "    from cupyx.scipy import special as _cpx_special",
            "    _gpu_erf = _cpx_special.erf",
            "except Exception:",
            "    if hasattr(cp, 'special') and hasattr(cp.special, 'erf'):",
            "        _gpu_erf = cp.special.erf",
            "    elif hasattr(cp, 'erf'):",
            "        _gpu_erf = cp.erf",
            "    else:",
            "        raise RuntimeError('No erf implementation found in CuPy installation.')",
            "",
            "def _cummax(x):",
            "    if int(x.size) == 0:",
            "        return x",
            "    if hasattr(cp, 'cummax'):",
            "        return cp.cummax(x)",
            "    out = cp.empty_like(x)",
            "    out[0] = x[0]",
            "    for idx in range(1, int(x.size)):",
            "        out[idx] = cp.maximum(out[idx - 1], x[idx])",
            "    return out",
            "",
            "def repsampler_ks1_batch(selected_sorted, candidates, mean, sd):",
            "    s = cp.asarray(selected_sorted, dtype=cp.float64)",
            "    c = cp.asarray(candidates, dtype=cp.float64)",
            "    n = int(s.size)",
            "    n1 = float(n + 1)",
            "    sqrt2 = cp.sqrt(2.0)",
            "",
            "    if n > 0:",
            "        F_s = 0.5 * (1.0 + _gpu_erf((s - mean) / (sd * sqrt2)))",
            "        i = cp.arange(1, n + 1, dtype=cp.float64)",
            "        A = i / n1 - F_s",
            "        B = F_s - (i - 1.0) / n1",
            "        step = 1.0 / n1",
            "        A_shift = A + step",
            "        B_shift = B - step",
            "",
            "        prefixA = cp.concatenate((cp.array([-cp.inf], dtype=cp.float64), _cummax(A)))",
            "        prefixB = cp.concatenate((cp.array([-cp.inf], dtype=cp.float64), _cummax(B)))",
            "        suffixAshift = cp.concatenate((_cummax(A_shift[::-1])[::-1], cp.array([-cp.inf], dtype=cp.float64)))",
            "        suffixBshift = cp.concatenate((_cummax(B_shift[::-1])[::-1], cp.array([-cp.inf], dtype=cp.float64)))",
            "    else:",
            "        prefixA = cp.array([-cp.inf], dtype=cp.float64)",
            "        prefixB = cp.array([-cp.inf], dtype=cp.float64)",
            "        suffixAshift = cp.array([-cp.inf], dtype=cp.float64)",
            "        suffixBshift = cp.array([-cp.inf], dtype=cp.float64)",
            "",
            "    k = cp.searchsorted(s, c, side='right').astype(cp.int64)",
            "    F_c = 0.5 * (1.0 + _gpu_erf((c - mean) / (sd * sqrt2)))",
            "    dplus_c = (k + 1.0) / n1 - F_c",
            "    dminus_c = F_c - k / n1",
            "",
            "    dplus = cp.maximum(cp.maximum(prefixA[k], suffixAshift[k]), dplus_c)",
            "    dminus = cp.maximum(cp.maximum(prefixB[k], suffixBshift[k]), dminus_c)",
            "    D = cp.maximum(dplus, dminus)",
            "",
            "    lam2 = (cp.sqrt(n1) * D) ** 2",
            "    p = cp.zeros_like(D)",
            "    for j in range(1, 101):",
            "        coeff = 2.0 if (j % 2 == 1) else -2.0",
            "        p = p + coeff * cp.exp(-2.0 * (j ** 2) * lam2)",
            "    p = cp.clip(p, 0.0, 1.0)",
            "",
            "    return cp.asnumpy(D), cp.asnumpy(p)",
            sep = "\n"
          ),
          local = FALSE
        ),
        error = function(e) {
          msg <- conditionMessage(e)
          if (.Platform$OS.type == "windows") {
            stop(
              paste0(
                "Failed to initialize CUDA backend. On Windows, ensure CUDA DLL paths are registered ",
                "(Python 3.8+ may require `os.add_dll_directory`) and install Python packages ",
                "`nvidia-cuda-nvrtc-cu12`, `nvidia-cuda-runtime-cu12`, and `nvidia-cublas-cu12` ",
                "alongside `cupy-cuda12x`.\nOriginal error: ", msg
              ),
              call. = FALSE
            )
          }
          stop(paste0("Failed to initialize CUDA backend: ", msg), call. = FALSE)
        }
      )

      py_fun <<- get("py", envir = rns)$repsampler_ks1_batch
    }

    res <- py_fun(selected_sorted, candidates, mean, sd)
    list(
      D = as.numeric(py_to_r(res[[1]])),
      p = as.numeric(py_to_r(res[[2]]))
    )
  }
})

#' Greedy Representative Sampling
#'
#' `repsample()` performs greedy representative sampling with support for
#' population and theoretical targets, retain/force behavior, exact/asymptotic
#' testing, variable weights, and early stopping rules.
#'
#' @param data Data frame.
#' @param size Desired sample size.
#' @param cont Continuous variable names.
#' @param bincat Binary/categorical variable names.
#' @param mean Means for continuous variables (theoretical sampling only).
#' @param sd Standard deviations for continuous variables (theoretical sampling only).
#' @param dist Target distribution form(s) for `cont` variables in theoretical
#'   mode. Supported: `"normal"` (default), `"lognormal"`, `"poisson"`.
#'   You can pass one value, one per `cont`, or a named vector/list.
#' @param perc Percentages for binary variables (theoretical sampling only), in `(0, 100)`.
#' @param seednum Random seed. Default is `7`.
#' @param randomperc Percentage of random picks at start. Default is `10`.
#' @param srule Early stopping threshold on Fisher combined p-value in `[0.5, 1)`.
#' @param rrule Early stopping candidate subsample size per iteration (minimum `10`).
#' @param quality Optional tuning preset for faster setup:
#'   `"manual"` (default), `"fast"`, `"balanced"`, `"strict"`.
#'   Presets only fill in `randomperc` / `rrule` / `srule` when those are not
#'   explicitly provided.
#' @param wght Variable weights that must sum to `100`; continuous first, then binary/categorical.
#' @param retain Optional existing binary 0/1 column to carry forward selected cases.
#' @param exact Use exact tests where Stata does (`TRUE/FALSE`).
#' @param force Force overwrite of existing `repsample` column.
#' @param subset Optional logical expression evaluated in `data` (Stata `if` analog).
#' @param in_rows Optional integer row indices to keep (Stata `in` analog).
#' @param n_cores Number of CPU cores for candidate scoring in each greedy
#'   iteration. `1` (default) is serial. Values greater than `1` use
#'   fork-based parallelism in terminal R on Unix-like systems and PSOCK
#'   clusters in environments where forking is unstable (for example,
#'   RStudio on macOS).
#' @param backend Compute backend. `"cpu"` uses the standard C++/CPU path.
#'   `"cuda"` enables an additional CUDA fast-path for theoretical sampling
#'   with continuous variables and `exact = FALSE` (uses Python `cupy`
#'   through `reticulate`).
#' @param cuda_fallback If `TRUE` (default), automatically falls back to CPU
#'   when CUDA dependencies are unavailable or a CUDA fast-path cannot be used.
#'
#' @return An object of class `repsample_result` with:
#' - `data`: input data plus `repsample` (`0/1`)
#' - `r`: named numeric vector equivalent to Stata `r()` scalars (empty if fully random)
#' - `selected_rows`: row indices selected in the returned data
#' - `meta`: run metadata
#'
#' @export
repsample <- function(data,
                      size,
                      cont = NULL,
                      bincat = NULL,
                      mean = NULL,
                      sd = NULL,
                      dist = NULL,
                      perc = NULL,
                      seednum = 7,
                      randomperc = 10,
                      srule = NULL,
                      rrule = NULL,
                      quality = c("manual", "fast", "balanced", "strict"),
                      wght = NULL,
                      retain = NULL,
                      exact = FALSE,
                      force = FALSE,
                      subset = NULL,
                      in_rows = NULL,
                      n_cores = 1,
                      backend = c("cpu", "cuda"),
                      cuda_fallback = TRUE) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }

  if (length(size) != 1L || is.na(size) || size < 1 || abs(size - round(size)) > 0) {
    stop("A single positive integer is required as the desired sample size.", call. = FALSE)
  }
  size <- as.integer(size)
  quality <- match.arg(quality)

  randomperc_missing <- missing(randomperc)
  srule_missing <- missing(srule)
  rrule_missing <- missing(rrule)

  if (quality != "manual") {
    q_defaults <- quality_defaults(quality, size)
    if (randomperc_missing) {
      randomperc <- q_defaults$randomperc
    }
    if (rrule_missing) {
      rrule <- q_defaults$rrule
    }
    if (srule_missing) {
      srule <- q_defaults$srule
    }
  }

  check_scalar_flag(exact, "exact")
  check_scalar_flag(force, "force")

  if (!is.numeric(seednum) || length(seednum) != 1L || is.na(seednum)) {
    stop("`seednum` must be a single numeric value.", call. = FALSE)
  }

  if (!is.numeric(randomperc) || length(randomperc) != 1L || is.na(randomperc) ||
      randomperc < 0 || randomperc > 100) {
    stop("Random percentage parameter must be in the [0, 100] range.", call. = FALSE)
  }

  backend <- match.arg(backend)
  cuda_requested <- identical(backend, "cuda")
  requested_backend <- backend
  check_scalar_flag(cuda_fallback, "cuda_fallback")

  if (!is.numeric(n_cores) || length(n_cores) != 1L || is.na(n_cores) ||
      n_cores < 1 || abs(n_cores - round(n_cores)) > 0) {
    stop("`n_cores` must be a positive integer.", call. = FALSE)
  }
  n_cores <- as.integer(n_cores)

  parallel_mode <- "none"
  in_rstudio <- nzchar(Sys.getenv("RSTUDIO"))
  if (n_cores > 1L) {
    if (.Platform$OS.type == "windows") {
      parallel_mode <- "psock"
    } else if (in_rstudio && identical(Sys.info()[["sysname"]], "Darwin")) {
      parallel_mode <- "psock"
      warning(
        paste0(
          "Detected RStudio on macOS. Using PSOCK parallel backend for n_cores > 1."
        ),
        call. = FALSE
      )
    } else {
      parallel_mode <- "fork"
    }
  }
  parallel_enabled <- (parallel_mode != "none")

  if (in_rstudio && identical(Sys.info()[["sysname"]], "Darwin") && parallel_mode == "fork") {
    warning(
      paste0(
        "Detected RStudio on macOS. Fork backend is unstable there; switching to PSOCK."
      ),
      call. = FALSE
    )
    parallel_mode <- "psock"
  }

  psock_parallel_threshold <- max(5000L, n_cores * 1024L)

  if (cuda_requested && !cuda_backend_available()) {
    if (!isTRUE(cuda_fallback)) {
      stop(
        paste0(
          "CUDA backend requires `reticulate`, Python `cupy`, and a visible CUDA GPU. ",
          "On Windows, also ensure CUDA DLL search paths are registered and install ",
          "`nvidia-cuda-nvrtc-cu12`, `nvidia-cuda-runtime-cu12`, and `nvidia-cublas-cu12`. ",
          "Install dependencies on the target machine or use `backend = \"cpu\"`."
        ),
        call. = FALSE
      )
    }
    warning(
      paste0(
        "CUDA backend requested but unavailable; falling back to CPU. ",
        "Set `cuda_fallback = FALSE` to error instead."
      ),
      call. = FALSE
    )
    backend <- "cpu"
    cuda_requested <- FALSE
  }

  if (!is.null(srule)) {
    if (!is.numeric(srule) || length(srule) != 1L || is.na(srule) || srule < 0.5 || srule >= 1) {
      stop("Early stop matching rule must be in the [0.5, 1) range.", call. = FALSE)
    }
  }

  if (!is.null(rrule)) {
    if (!is.numeric(rrule) || length(rrule) != 1L || is.na(rrule) || rrule < 10) {
      stop("Random sampling early stopping rule must sample at least 10 cases.", call. = FALSE)
    }
    rrule <- as.integer(floor(rrule))
  }

  cont <- normalize_varlist(cont)
  bincat <- normalize_varlist(bincat)

  ccnt <- length(cont)
  bcnt <- length(bincat)
  allvar <- c(cont, bincat)

  if ((ccnt + bcnt) == 0L) {
    stop(
      "No variables provided on which the sample will be drawn. Use `cont` and/or `bincat`.",
      call. = FALSE
    )
  }

  if (anyDuplicated(allvar) > 0L) {
    stop("Variables provided need to be unique.", call. = FALSE)
  }

  if (!all(allvar %in% names(data))) {
    missing_vars <- allvar[!(allvar %in% names(data))]
    stop(sprintf("Unknown variable(s): %s", paste(missing_vars, collapse = ", ")), call. = FALSE)
  }

  mean <- normalize_numlist(mean)
  sd <- normalize_numlist(sd)
  perc <- normalize_numlist(perc)
  cont_dist <- normalize_cont_dist(dist, cont)

  theoretical <- (length(mean) > 0L || length(sd) > 0L || length(perc) > 0L)
  if (!theoretical && !is.null(dist)) {
    warning(
      "`dist` is ignored in population mode (no theoretical targets were provided).",
      call. = FALSE
    )
  }

  if (theoretical && (ccnt != length(mean) || ccnt != length(sd) || bcnt != length(perc))) {
    stop(
      "Details for theoretical distributions do not agree with the variables provided.",
      call. = FALSE
    )
  }

  if (length(perc) > 0L && any(perc <= 0 | perc >= 100)) {
    stop("Percentages must be in the (0, 100) range.", call. = FALSE)
  }

  if (length(sd) > 0L && any(sd <= 0)) {
    stop("Standard deviations must be positive in theoretical sampling.", call. = FALSE)
  }

  if (theoretical && ccnt > 0L) {
    for (i in seq_len(ccnt)) {
      d_i <- cont_dist[i]
      if (d_i == "lognormal") {
        if (mean[i] <= 0 || sd[i] <= 0) {
          stop(
            sprintf(
              "Lognormal target for `%s` requires positive target mean and sd.",
              cont[i]
            ),
            call. = FALSE
          )
        }
      } else if (d_i == "poisson") {
        if (mean[i] <= 0) {
          stop(
            sprintf(
              "Poisson target for `%s` requires a positive target mean (lambda).",
              cont[i]
            ),
            call. = FALSE
          )
        }
        sd_expected <- sqrt(mean[i])
        if (is.finite(sd[i]) && abs(sd[i] - sd_expected) > max(0.05, 0.05 * sd_expected)) {
          warning(
            sprintf(
              "Poisson target for `%s`: provided sd %.4g differs from sqrt(mean)=%.4g; using mean as lambda.",
              cont[i], sd[i], sd_expected
            ),
            call. = FALSE
          )
        }
      }
    }
  }
  all_cont_normal <- (ccnt == 0L || all(cont_dist == "normal"))
  cuda_fastpath <- (cuda_requested && theoretical && ccnt > 0L && bcnt == 0L && !exact && all_cont_normal)
  if (cuda_requested && !cuda_fastpath) {
    warning(
      paste0(
        "CUDA requested, but current model is not supported by the CUDA fast-path ",
        "(requires theoretical sampling, normal continuous variables only, exact = FALSE). ",
        "Falling back to CPU scoring."
      ),
      call. = FALSE
    )
  }

  wght <- normalize_numlist(wght)
  k_total <- ccnt + bcnt

  if (length(wght) > 0L) {
    if (length(wght) != k_total) {
      stop(
        "Number of weights provided must match number of variables; continuous first, then binary/categorical.",
        call. = FALSE
      )
    }
    if (abs(sum(wght) - 100) > 1e-9) {
      stop("Weight scores must add up to 100.", call. = FALSE)
    }
  }

  cont_weights <- if (ccnt > 0L) rep(1, ccnt) else numeric(0)
  bincat_weights <- if (bcnt > 0L) rep(1, bcnt) else numeric(0)
  weights_provided <- length(wght) > 0L

  if (weights_provided) {
    if (ccnt > 0L) {
      cont_weights <- wght[seq_len(ccnt)] / 100
    }
    if (bcnt > 0L) {
      bincat_weights <- wght[(ccnt + 1L):(ccnt + bcnt)] / 100
    }
  }

  for (v in cont) {
    distinct_n <- length(unique(data[[v]][!is.na(data[[v]])]))
    if (distinct_n < 7L) {
      warning(sprintf(
        "Only %d distinct values observed for 'continuous' variable `%s`.",
        distinct_n,
        v
      ), call. = FALSE)
    }
  }

  for (v in bincat) {
    distinct_n <- length(unique(data[[v]][!is.na(data[[v]])]))

    if (theoretical && distinct_n > 2L) {
      stop(sprintf(
        "Theoretical sampling cannot work with categorical variable `%s`; use binary variables.",
        v
      ), call. = FALSE)
    }

    if (distinct_n > 7L) {
      warning(sprintf(
        "%d distinct values observed for 'categorical' variable `%s`; consider fewer categories.",
        distinct_n,
        v
      ), call. = FALSE)
    }
  }

  if (!is.null(retain)) {
    if (!is.character(retain) || length(retain) != 1L) {
      stop("`retain` must be a single column name.", call. = FALSE)
    }

    if (!retain %in% names(data)) {
      stop(sprintf("Retain variable `%s` not found.", retain), call. = FALSE)
    }

    retain_vals <- data[[retain]]
    retain_non_na <- retain_vals[!is.na(retain_vals)]
    retain_unique <- unique(retain_non_na)

    if (length(retain_unique) > 2L || any(!(retain_unique %in% c(0, 1)))) {
      stop(
        "Variable used to continue sampling must be a 0-1 binary with selected cases coded 1.",
        call. = FALSE
      )
    }

    if (identical(retain, "repsample") && force) {
      stop(
        "You cannot force-drop `repsample` and update it through retain simultaneously.",
        call. = FALSE
      )
    }
  }

  if ("repsample" %in% names(data) && !(force || identical(retain, "repsample"))) {
    stop(
      "Column `repsample` already exists. Use `force = TRUE` or `retain = \"repsample\"`.",
      call. = FALSE
    )
  }

  n_total <- nrow(data)

  subset_mask <- rep(TRUE, n_total)
  if (!missing(subset) && !is.null(subset)) {
    subset_expr <- substitute(subset)
    subset_mask <- eval(subset_expr, data, parent.frame())

    if (!is.logical(subset_mask) || length(subset_mask) != n_total) {
      stop("`subset` must evaluate to a logical vector with one value per row.", call. = FALSE)
    }

    subset_mask[is.na(subset_mask)] <- FALSE
  }

  in_mask <- rep(TRUE, n_total)
  if (!is.null(in_rows)) {
    if (!is.numeric(in_rows) || anyNA(in_rows)) {
      stop("`in_rows` must be integer row positions.", call. = FALSE)
    }
    in_rows <- unique(as.integer(in_rows))
    in_rows <- in_rows[in_rows >= 1L & in_rows <= n_total]
    in_mask <- rep(FALSE, n_total)
    in_mask[in_rows] <- TRUE
  }

  selected_scope <- subset_mask & in_mask
  work_idx <- which(selected_scope)

  if (length(work_idx) == 0L) {
    stop("No rows available after applying subset/in_rows filters.", call. = FALSE)
  }

  work <- data[work_idx, , drop = FALSE]

  cc_mask <- stats::complete.cases(work[, allvar, drop = FALSE])
  work <- work[cc_mask, , drop = FALSE]
  work_idx <- work_idx[cc_mask]

  if (nrow(work) == 0L) {
    stop("No rows available after removing missing values in sampling variables.", call. = FALSE)
  }

  if (theoretical && bcnt > 0L) {
    for (v in bincat) {
      vals <- work[[v]]
      if (any(!(vals %in% c(0, 1)))) {
        stop(sprintf(
          "Theoretical sampling requires binary 0/1 variables; `%s` has other values.",
          v
        ), call. = FALSE)
      }
    }
  }

  oldsample <- rep(0, nrow(work))
  if (!is.null(retain)) {
    oldsample <- work[[retain]]
  }

  prevsel <- sum(oldsample == 1, na.rm = TRUE)
  needed <- size - prevsel

  if (prevsel > size) {
    stop("Continuation variable holds more cases than requested sample size.", call. = FALSE)
  }
  if (prevsel == size) {
    stop("Continuation variable holds as many cases as the requested sample size.", call. = FALSE)
  }

  available <- sum(oldsample == 0, na.rm = TRUE)
  if (available < needed) {
    stop("Not enough cases to provide the requested sample.", call. = FALSE)
  }
  if (available == needed) {
    stop("Requested sample matches the available number of cases.", call. = FALSE)
  }

  eligible_new <- if (is.null(retain)) {
    rep(TRUE, nrow(work))
  } else {
    (!is.na(oldsample) & oldsample == 0)
  }

  selected <- (!is.na(oldsample) & oldsample == 1)

  with_seed <- function(seed, code) {
    has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_seed) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }

    on.exit({
      if (has_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)

    set.seed(seed)
    force(code)
  }

  random_selection <- function() {
    samplenum <- as.integer(needed)
    randsel <- as.integer(stata_round(samplenum * randomperc / 100))

    pool <- which(eligible_new & !selected)
    if (randsel > length(pool)) {
      randsel <- length(pool)
    }

    if (randsel > 0L) {
      picks <- sample(pool, size = randsel, replace = FALSE)
      selected[picks] <<- TRUE
    }

    randsel
  }

  randsel <- with_seed(seednum, random_selection())

  samplenum <- as.integer(needed)
  scount <- as.integer(randsel)

  cont_values <- if (ccnt > 0L) lapply(cont, function(v) as.numeric(work[[v]])) else list()
  cont_pop_sorted <- if (ccnt > 0L && !theoretical) lapply(cont_values, sort.int, method = "quick") else list()
  cont_pop_raw <- if (ccnt > 0L && !theoretical) cont_values else list()
  cont_cdf_values <- vector("list", ccnt)

  if (ccnt > 0L && theoretical) {
    cdf_eps <- 1e-12
    for (i in seq_len(ccnt)) {
      vals <- cont_values[[i]]
      d_i <- cont_dist[i]

      if (d_i == "normal") {
        next
      }

      if (d_i == "lognormal") {
        sigma2 <- log1p((sd[i] * sd[i]) / (mean[i] * mean[i]))
        sdlog <- sqrt(sigma2)
        meanlog <- log(mean[i]) - 0.5 * sigma2
        cdf_vals <- stats::plnorm(vals, meanlog = meanlog, sdlog = sdlog)
        cdf_vals <- pmin(1 - cdf_eps, pmax(cdf_eps, as.numeric(cdf_vals)))
        cont_cdf_values[[i]] <- cdf_vals
        next
      }

      if (d_i == "poisson") {
        lambda <- mean[i]
        nonint <- mean(abs(vals - round(vals)) > 1e-8)
        if (is.finite(nonint) && nonint > 0.05) {
          warning(
            sprintf(
              "Poisson target for `%s`: %.1f%% of values are non-integer; using floor(x) for CDF matching.",
              cont[i], 100 * nonint
            ),
            call. = FALSE
          )
        }
        cdf_vals <- stats::ppois(floor(vals), lambda = lambda)
        cdf_vals <- pmin(1 - cdf_eps, pmax(cdf_eps, as.numeric(cdf_vals)))
        cont_cdf_values[[i]] <- cdf_vals
        next
      }
    }
  }

  cont_selected_sorted <- if (ccnt > 0L) {
    lapply(cont_values, function(v) sort.int(v[selected], method = "quick"))
  } else {
    list()
  }
  cont_selected_cdf_sorted <- if (ccnt > 0L && theoretical) {
    lapply(seq_len(ccnt), function(i) {
      if (cont_dist[i] == "normal") {
        numeric(0)
      } else {
        sort.int(cont_cdf_values[[i]][selected], method = "quick")
      }
    })
  } else {
    list()
  }

  bincat_pop_counts <- list()
  bincat_sel_counts <- list()
  bincat_enc <- list()
  bincat_values <- list()
  bincat_sel_success <- numeric(bcnt)

  if (bcnt > 0L) {
    if (!theoretical) {
      for (i in seq_len(bcnt)) {
        vals <- work[[bincat[i]]]
        lvls <- sort(unique(vals))
        enc <- match(vals, lvls)
        pop_counts <- tabulate(enc, nbins = length(lvls))
        sel_counts <- tabulate(enc[selected], nbins = length(lvls))

        bincat_enc[[i]] <- enc
        bincat_pop_counts[[i]] <- pop_counts
        bincat_sel_counts[[i]] <- sel_counts
      }
    } else {
      for (i in seq_len(bcnt)) {
        vals <- as.numeric(work[[bincat[i]]])
        bincat_values[[i]] <- vals
        bincat_sel_success[i] <- sum(vals[selected] == 1)
      }
    }
  }

  current_selected_n <- sum(selected)
  log_floor <- .Machine$double.xmin

  evaluate_candidate <- function(candidate_idx) {
    n1 <- current_selected_n + 1L

    total_score <- 0
    cont_stats <- vector("list", ccnt)
    bin_stats <- vector("list", bcnt)

    if (ccnt > 0L) {
      for (i in seq_len(ccnt)) {
        cand_value <- cont_values[[i]][candidate_idx]

        if (!theoretical) {
          if (exact) {
            x_trial <- c(cont_selected_sorted[[i]], cand_value)
            ks <- suppressWarnings(stats::ks.test(x_trial, cont_pop_raw[[i]], exact = TRUE))
            p <- as.numeric(ks$p.value)
            d <- as.numeric(ks$statistic[[1]])
          } else {
            d <- ks2_d_with_candidate_cpp(
              cont_selected_sorted[[i]],
              cand_value,
              cont_pop_sorted[[i]]
            )
            n_eff <- (n1 * length(cont_pop_sorted[[i]])) / (n1 + length(cont_pop_sorted[[i]]))
            p <- ks_p_asymptotic(d, n_eff)
          }
        } else {
          if (cont_dist[i] == "normal") {
            d <- ks1_norm_d_with_candidate_cpp(
              cont_selected_sorted[[i]],
              cand_value,
              mean[i],
              sd[i]
            )
          } else {
            d <- ks1_cdf_d_with_candidate_cpp(
              cont_selected_sorted[[i]],
              cont_selected_cdf_sorted[[i]],
              cand_value,
              cont_cdf_values[[i]][candidate_idx]
            )
          }
          p <- ks_p_asymptotic(d, n1)
        }

        p_safe <- max(log_floor, p)
        total_score <- total_score - cont_weights[i] * 2 * log(p_safe)
        cont_stats[[i]] <- list(D = d, p = p)
      }
    }

    if (bcnt > 0L) {
      for (i in seq_len(bcnt)) {
        if (!theoretical) {
          trial_counts <- bincat_sel_counts[[i]]
          enc_value <- bincat_enc[[i]][candidate_idx]
          trial_counts[enc_value] <- trial_counts[enc_value] + 1L

          if (exact) {
            if (length(trial_counts) <= 1L) {
              p <- 1
              chi2 <- NA_real_
            } else {
              mat <- cbind(trial_counts, bincat_pop_counts[[i]])
              p <- as.numeric(stats::fisher.test(mat)$p.value)
              chi2 <- NA_real_
            }
          } else {
            ch <- chisq_from_counts(trial_counts, bincat_pop_counts[[i]])
            p <- ch$p
            chi2 <- ch$chi2
          }

          p_safe <- max(log_floor, p)
          total_score <- total_score - bincat_weights[i] * 2 * log(p_safe)

          bin_stats[[i]] <- list(p = p, chi2 = chi2, z = NA_real_)
        } else {
          succ <- bincat_sel_success[i] + as.integer(bincat_values[[i]][candidate_idx] == 1)
          p0 <- perc[i] / 100

          if (exact) {
            p <- as.numeric(stats::binom.test(succ, n1, p = p0)$p.value)
            z <- NA_real_
          } else {
            phat <- succ / n1
            z <- (phat - p0) / sqrt(p0 * (1 - p0) / n1)
            p <- 2 * stats::pnorm(-abs(z))
          }

          p_safe <- max(log_floor, p)
          total_score <- total_score - bincat_weights[i] * 2 * log(p_safe)

          bin_stats[[i]] <- list(p = p, chi2 = NA_real_, z = z)
        }
      }
    }

    if (weights_provided) {
      total_score <- total_score * k_total
    }

    list(score = total_score, cont = cont_stats, bin = bin_stats)
  }

  score_candidate_with_state <- function(candidate_idx,
                                         state_current_selected_n,
                                         state_cont_selected_sorted,
                                         state_cont_selected_cdf_sorted,
                                         state_bincat_sel_counts,
                                         state_bincat_sel_success) {
    n1 <- state_current_selected_n + 1L
    total_score <- 0

    if (ccnt > 0L) {
      for (i in seq_len(ccnt)) {
        cand_value <- cont_values[[i]][candidate_idx]

        if (!theoretical) {
          if (exact) {
            x_trial <- c(state_cont_selected_sorted[[i]], cand_value)
            ks <- suppressWarnings(stats::ks.test(x_trial, cont_pop_raw[[i]], exact = TRUE))
            p <- as.numeric(ks$p.value)
          } else {
            d <- ks2_d_with_candidate_cpp(
              state_cont_selected_sorted[[i]],
              cand_value,
              cont_pop_sorted[[i]]
            )
            n_eff <- (n1 * length(cont_pop_sorted[[i]])) / (n1 + length(cont_pop_sorted[[i]]))
            p <- ks_p_asymptotic(d, n_eff)
          }
        } else {
          if (cont_dist[i] == "normal") {
            d <- ks1_norm_d_with_candidate_cpp(
              state_cont_selected_sorted[[i]],
              cand_value,
              mean[i],
              sd[i]
            )
          } else {
            d <- ks1_cdf_d_with_candidate_cpp(
              state_cont_selected_sorted[[i]],
              state_cont_selected_cdf_sorted[[i]],
              cand_value,
              cont_cdf_values[[i]][candidate_idx]
            )
          }
          p <- ks_p_asymptotic(d, n1)
        }

        p_safe <- max(log_floor, p)
        total_score <- total_score - cont_weights[i] * 2 * log(p_safe)
      }
    }

    if (bcnt > 0L) {
      for (i in seq_len(bcnt)) {
        if (!theoretical) {
          trial_counts <- state_bincat_sel_counts[[i]]
          enc_value <- bincat_enc[[i]][candidate_idx]
          trial_counts[enc_value] <- trial_counts[enc_value] + 1L

          if (exact) {
            if (length(trial_counts) <= 1L) {
              p <- 1
            } else {
              mat <- cbind(trial_counts, bincat_pop_counts[[i]])
              p <- as.numeric(stats::fisher.test(mat)$p.value)
            }
          } else {
            p <- chisq_from_counts(trial_counts, bincat_pop_counts[[i]])$p
          }
        } else {
          succ <- state_bincat_sel_success[i] + as.integer(bincat_values[[i]][candidate_idx] == 1)
          p0 <- perc[i] / 100

          if (exact) {
            p <- as.numeric(stats::binom.test(succ, n1, p = p0)$p.value)
          } else {
            phat <- succ / n1
            z <- (phat - p0) / sqrt(p0 * (1 - p0) / n1)
            p <- 2 * stats::pnorm(-abs(z))
          }
        }

        p_safe <- max(log_floor, p)
        total_score <- total_score - bincat_weights[i] * 2 * log(p_safe)
      }
    }

    if (weights_provided) {
      total_score <- total_score * k_total
    }

    total_score
  }

  psock_stateful_fastpath <- (parallel_mode == "psock" && theoretical && bcnt == 0L && !exact && all_cont_normal)
  psock_cluster <- NULL
  if (parallel_mode == "psock") {
    max_candidates_est <- if (!is.null(rrule)) {
      as.integer(rrule)
    } else {
      as.integer(sum(eligible_new & !selected))
    }
    if (max_candidates_est < psock_parallel_threshold) {
      parallel_mode <- "none"
      parallel_enabled <- FALSE
      psock_stateful_fastpath <- FALSE
    }
  }
  if (parallel_mode == "psock") {
    psock_cluster <- tryCatch(
      parallel::makeCluster(n_cores),
      error = function(e) e
    )

    if (inherits(psock_cluster, "error")) {
      warning(
        sprintf(
          "Failed to start PSOCK cluster (%s). Falling back to serial.",
          conditionMessage(psock_cluster)
        ),
        call. = FALSE
      )
      parallel_mode <- "none"
      parallel_enabled <- FALSE
      psock_cluster <- NULL
    } else {
      on.exit(parallel::stopCluster(psock_cluster), add = TRUE)
      libs <- .libPaths()
      parallel::clusterCall(psock_cluster, function(x) .libPaths(x), libs)

      load_ok <- tryCatch({
        parallel::clusterEvalQ(
          psock_cluster,
          library(RepsampleR, quietly = TRUE, warn.conflicts = FALSE)
        )
        parallel::clusterExport(psock_cluster, "score_candidate_with_state", envir = environment())
        parallel::clusterEvalQ(
          psock_cluster,
          {
            .repsample_psock_chunk <- function(candidate_chunk) {
              vapply(candidate_chunk, function(cand) {
                score_candidate_with_state(
                  cand,
                  psock_state$current_selected_n,
                  psock_state$cont_selected_sorted,
                  psock_state$cont_selected_cdf_sorted,
                  psock_state$bincat_sel_counts,
                  psock_state$bincat_sel_success
                )
              }, numeric(1))
            }
            TRUE
          }
        )

        if (psock_stateful_fastpath) {
          parallel::clusterExport(
            psock_cluster,
            varlist = c(
              "ccnt", "cont_values", "mean", "sd", "cont_weights",
              "weights_provided", "k_total", "log_floor"
            ),
            envir = environment()
          )
          parallel::clusterCall(
            psock_cluster,
            function(n_sel, sorted_list) {
              assign(
                ".reps_state",
                list(
                  current_selected_n = as.integer(n_sel),
                  cont_selected_sorted = sorted_list
                ),
                envir = .GlobalEnv
              )
              TRUE
            },
            current_selected_n,
            cont_selected_sorted
          )
          parallel::clusterEvalQ(
            psock_cluster,
            {
              .repsample_psock_score_chunk_fast <- function(candidate_chunk) {
                st <- get(".reps_state", envir = .GlobalEnv)
                n1 <- st$current_selected_n + 1L
                out <- numeric(length(candidate_chunk))

                for (idx in seq_along(candidate_chunk)) {
                  cand <- candidate_chunk[idx]
                  total <- 0

                  for (i in seq_len(ccnt)) {
                    cand_value <- cont_values[[i]][cand]
                    d <- RepsampleR:::ks1_norm_d_with_candidate_cpp(
                      st$cont_selected_sorted[[i]],
                      cand_value,
                      mean[i],
                      sd[i]
                    )
                    p <- RepsampleR:::ks_p_asymptotic(d, n1)
                    total <- total - cont_weights[i] * 2 * log(max(log_floor, p))
                  }

                  if (weights_provided) {
                    total <- total * k_total
                  }
                  out[idx] <- total
                }
                out
              }

              .repsample_psock_apply_pick <- function(chosen_idx) {
                st <- get(".reps_state", envir = .GlobalEnv)
                for (i in seq_len(ccnt)) {
                  st$cont_selected_sorted[[i]] <- RepsampleR:::insert_sorted_cpp(
                    st$cont_selected_sorted[[i]],
                    cont_values[[i]][chosen_idx]
                  )
                }
                st$current_selected_n <- st$current_selected_n + 1L
                assign(".reps_state", st, envir = .GlobalEnv)
                TRUE
              }
              TRUE
            }
          )
        }
        TRUE
      }, error = function(e) FALSE)

      if (!load_ok) {
        warning(
          "Failed to load RepsampleR on PSOCK workers. Falling back to serial.",
          call. = FALSE
        )
        parallel::stopCluster(psock_cluster)
        psock_cluster <- NULL
        parallel_mode <- "none"
        parallel_enabled <- FALSE
      }
    }
  }

  score_candidates_serial <- function(candidates) {
    scores <- numeric(length(candidates))
    eval_n <- 0L

    for (j in seq_along(candidates)) {
      cand <- candidates[j]
      info <- evaluate_candidate(cand)

      scores[j] <- info$score
      eval_n <- j

      if (!is.null(srule)) {
        p_now <- stats::pchisq(info$score, df = 2 * k_total, lower.tail = FALSE)
        if (p_now >= srule) {
          break
        }
      }
    }

    list(
      evaluated_candidates = candidates[seq_len(eval_n)],
      evaluated_scores = scores[seq_len(eval_n)]
    )
  }

  last_overall_p <- NA_real_
  last_overall_chi2 <- NA_real_
  last_details <- NULL
  
  evaluate_candidates_ordered <- function(candidates) {
    n_candidates <- length(candidates)
    parallel_candidate_threshold <- if (parallel_mode == "psock") {
      psock_parallel_threshold
    } else {
      max(128L, n_cores * 32L)
    }

    if (n_candidates == 0L) {
      return(list(
        evaluated_candidates = integer(0),
        evaluated_scores = numeric(0)
      ))
    }

    if (cuda_fastpath) {
      scores <- numeric(n_candidates)

      for (i in seq_len(ccnt)) {
        candidate_values <- cont_values[[i]][candidates]
        gpu <- cuda_ks1_batch(
          selected_sorted = cont_selected_sorted[[i]],
          candidates = candidate_values,
          mean = mean[i],
          sd = sd[i]
        )

        p_safe <- pmax(log_floor, gpu$p)
        scores <- scores - cont_weights[i] * 2 * log(p_safe)
      }

      if (weights_provided) {
        scores <- scores * k_total
      }

      eval_n <- n_candidates
      if (!is.null(srule)) {
        p_seq <- stats::pchisq(scores, df = 2 * k_total, lower.tail = FALSE)
        hit <- which(p_seq >= srule)
        if (length(hit) > 0L) {
          eval_n <- hit[1L]
        }
      }

      return(list(
        evaluated_candidates = candidates[seq_len(eval_n)],
        evaluated_scores = scores[seq_len(eval_n)]
      ))
    }

    if (!parallel_enabled || n_candidates < parallel_candidate_threshold) {
      return(score_candidates_serial(candidates))
    }

    worker_cores <- min(n_cores, n_candidates)
    scores <- numeric(n_candidates)
    eval_n <- n_candidates

    if (parallel_mode == "fork") {
      if (is.null(srule)) {
        details <- suppressWarnings(parallel::mclapply(
          candidates,
          evaluate_candidate,
          mc.cores = worker_cores,
          mc.preschedule = TRUE,
          mc.set.seed = FALSE,
          mc.cleanup = FALSE,
          mc.allow.recursive = FALSE
        ))
        scores <- vapply(details, function(x) x$score, numeric(1))
      } else {
        batch_size <- max(worker_cores * 4L, 32L)
        eval_n <- 0L

        for (start_i in seq.int(1L, n_candidates, by = batch_size)) {
          stop_i <- min(n_candidates, start_i + batch_size - 1L)
          idx <- start_i:stop_i
          batch_candidates <- candidates[idx]

          batch_details <- suppressWarnings(parallel::mclapply(
            batch_candidates,
            evaluate_candidate,
            mc.cores = min(worker_cores, length(batch_candidates)),
            mc.preschedule = TRUE,
            mc.set.seed = FALSE,
            mc.cleanup = FALSE,
            mc.allow.recursive = FALSE
          ))

          batch_scores <- vapply(batch_details, function(x) x$score, numeric(1))
          scores[idx] <- batch_scores
          batch_p <- stats::pchisq(batch_scores, df = 2 * k_total, lower.tail = FALSE)
          hit <- which(batch_p >= srule)

          if (length(hit) > 0L) {
            eval_n <- idx[hit[1L]]
            break
          }

          eval_n <- stop_i
        }
      }
    } else if (parallel_mode == "psock" && !is.null(psock_cluster)) {
      if (psock_stateful_fastpath) {
        if (is.null(srule)) {
          idx_chunks <- parallel::splitIndices(n_candidates, min(worker_cores, length(psock_cluster)))
          cand_chunks <- lapply(idx_chunks, function(ix) candidates[ix])
          chunk_scores <- parallel::parLapply(
            psock_cluster[seq_along(cand_chunks)],
            cand_chunks,
            function(chunk) .repsample_psock_score_chunk_fast(chunk)
          )
          for (i in seq_along(idx_chunks)) {
            scores[idx_chunks[[i]]] <- chunk_scores[[i]]
          }
        } else {
          batch_size <- max(worker_cores * 16L, 512L)
          eval_n <- 0L

          for (start_i in seq.int(1L, n_candidates, by = batch_size)) {
            stop_i <- min(n_candidates, start_i + batch_size - 1L)
            idx <- start_i:stop_i
            batch_candidates <- candidates[idx]

            idx_chunks <- parallel::splitIndices(length(batch_candidates), min(worker_cores, length(psock_cluster)))
            cand_chunks <- lapply(idx_chunks, function(ix) batch_candidates[ix])
            chunk_scores <- parallel::parLapply(
              psock_cluster[seq_along(cand_chunks)],
              cand_chunks,
              function(chunk) .repsample_psock_score_chunk_fast(chunk)
            )

            batch_scores <- numeric(length(batch_candidates))
            for (k in seq_along(idx_chunks)) {
              batch_scores[idx_chunks[[k]]] <- chunk_scores[[k]]
            }

            scores[idx] <- batch_scores
            batch_p <- stats::pchisq(batch_scores, df = 2 * k_total, lower.tail = FALSE)
            hit <- which(batch_p >= srule)

            if (length(hit) > 0L) {
              eval_n <- idx[hit[1L]]
              break
            }

            eval_n <- stop_i
          }
        }
      } else {
        psock_state <- list(
          current_selected_n = current_selected_n,
          cont_selected_sorted = cont_selected_sorted,
          cont_selected_cdf_sorted = cont_selected_cdf_sorted,
          bincat_sel_counts = bincat_sel_counts,
          bincat_sel_success = bincat_sel_success
        )
        psock_ok <- tryCatch({
          parallel::clusterExport(psock_cluster, "psock_state", envir = environment())
          TRUE
        }, error = function(e) FALSE)

        if (!psock_ok) {
          return(score_candidates_serial(candidates))
        }

        if (is.null(srule)) {
          idx_chunks <- parallel::splitIndices(n_candidates, min(worker_cores, length(psock_cluster)))
          cand_chunks <- lapply(idx_chunks, function(ix) candidates[ix])
          chunk_scores <- parallel::parLapply(
            psock_cluster[seq_along(cand_chunks)],
            cand_chunks,
            function(chunk) .repsample_psock_chunk(chunk)
          )
          for (i in seq_along(idx_chunks)) {
            scores[idx_chunks[[i]]] <- chunk_scores[[i]]
          }
        } else {
          batch_size <- max(worker_cores * 16L, 512L)
          eval_n <- 0L

          for (start_i in seq.int(1L, n_candidates, by = batch_size)) {
            stop_i <- min(n_candidates, start_i + batch_size - 1L)
            idx <- start_i:stop_i
            batch_candidates <- candidates[idx]

            idx_chunks <- parallel::splitIndices(length(batch_candidates), min(worker_cores, length(psock_cluster)))
            cand_chunks <- lapply(idx_chunks, function(ix) batch_candidates[ix])
            chunk_scores <- parallel::parLapply(
              psock_cluster[seq_along(cand_chunks)],
              cand_chunks,
              function(chunk) .repsample_psock_chunk(chunk)
            )

            batch_scores <- numeric(length(batch_candidates))
            for (k in seq_along(idx_chunks)) {
              batch_scores[idx_chunks[[k]]] <- chunk_scores[[k]]
            }

            scores[idx] <- batch_scores
            batch_p <- stats::pchisq(batch_scores, df = 2 * k_total, lower.tail = FALSE)
            hit <- which(batch_p >= srule)

            if (length(hit) > 0L) {
              eval_n <- idx[hit[1L]]
              break
            }

            eval_n <- stop_i
          }
        }
      }
    } else {
      return(score_candidates_serial(candidates))
    }

    list(
      evaluated_candidates = candidates[seq_len(eval_n)],
      evaluated_scores = scores[seq_len(eval_n)]
    )
  }

  deterministic_loop <- function() {
    while (scount < samplenum) {
      pool <- which(eligible_new & !selected)
      if (length(pool) == 0L) {
        stop("Something has gone wrong - no case has been selected.", call. = FALSE)
      }

      if (!is.null(rrule) && length(pool) > rrule) {
        candidates <- sample(pool, size = rrule, replace = FALSE)
      } else {
        candidates <- pool
      }

      eval_res <- evaluate_candidates_ordered(candidates)
      evaluated_candidates <- eval_res$evaluated_candidates
      evaluated_scores <- eval_res$evaluated_scores

      min_score <- min(evaluated_scores)
      tied <- which(evaluated_scores == min_score)

      if (length(tied) == 1L) {
        pick <- tied
      } else {
        pick <- sample(tied, 1L)
      }

      chosen <- evaluated_candidates[pick]
      chosen_info <- evaluate_candidate(chosen)

      selected[chosen] <<- TRUE
      scount <<- scount + 1L
      current_selected_n <<- current_selected_n + 1L

      if (ccnt > 0L) {
        for (i in seq_len(ccnt)) {
          cont_selected_sorted[[i]] <<- insert_sorted_cpp(
            cont_selected_sorted[[i]],
            cont_values[[i]][chosen]
          )
          if (theoretical && cont_dist[i] != "normal") {
            cont_selected_cdf_sorted[[i]] <<- insert_sorted_cpp(
              cont_selected_cdf_sorted[[i]],
              cont_cdf_values[[i]][chosen]
            )
          }
        }
      }

      if (bcnt > 0L) {
        if (!theoretical) {
          for (i in seq_len(bcnt)) {
            enc_value <- bincat_enc[[i]][chosen]
            bincat_sel_counts[[i]][enc_value] <<- bincat_sel_counts[[i]][enc_value] + 1L
          }
        } else {
          for (i in seq_len(bcnt)) {
            bincat_sel_success[i] <<- bincat_sel_success[i] + as.integer(bincat_values[[i]][chosen] == 1)
          }
        }
      }

      if (parallel_mode == "psock" && psock_stateful_fastpath && !is.null(psock_cluster)) {
        sync_ok <- tryCatch({
          parallel::clusterCall(
            psock_cluster,
            function(chosen_idx) .repsample_psock_apply_pick(chosen_idx),
            chosen
          )
          TRUE
        }, error = function(e) FALSE)

        if (!sync_ok) {
          warning(
            "PSOCK state sync failed; disabling PSOCK parallel path for remaining iterations.",
            call. = FALSE
          )
          parallel_enabled <<- FALSE
          parallel_mode <<- "none"
        }
      }

      last_overall_chi2 <<- chosen_info$score
      last_overall_p <<- stats::pchisq(last_overall_chi2, df = 2 * k_total, lower.tail = FALSE)
      last_details <<- chosen_info
    }
  }

  with_seed(seednum, deterministic_loop())

  out <- data
  out$repsample <- 0L
  out$repsample[work_idx[selected]] <- 1L

  r_vals <- numeric(0)

  if (randsel < samplenum) {
    r_list <- list(
      p = last_overall_p,
      chi2 = last_overall_chi2,
      df = 2 * k_total
    )

    if (ccnt > 0L) {
      for (i in seq_len(ccnt)) {
        nm <- cont[i]
        r_list[[paste0(nm, "_D")]] <- last_details$cont[[i]]$D
        r_list[[paste0(nm, "_p")]] <- last_details$cont[[i]]$p
      }
    }

    if (bcnt > 0L) {
      for (i in seq_len(bcnt)) {
        nm <- bincat[i]
        r_list[[paste0(nm, "_p")]] <- last_details$bin[[i]]$p

        if (!theoretical && !exact) {
          r_list[[paste0(nm, "_chi2")]] <- last_details$bin[[i]]$chi2
        }

        if (theoretical && !exact) {
          r_list[[paste0(nm, "_z")]] <- last_details$bin[[i]]$z
        }
      }
    }

    r_vals <- unlist(r_list, use.names = TRUE)
  }

  structure(
    list(
      data = out,
      r = r_vals,
      selected_rows = which(out$repsample == 1L),
      meta = list(
        size = size,
        mode = if (theoretical) "theoretical" else "population",
        seednum = seednum,
        randomperc = randomperc,
        randsel = randsel,
        srule = srule,
        rrule = rrule,
        n_cores = n_cores,
        parallel_enabled = parallel_enabled,
        parallel_mode = parallel_mode,
        psock_stateful_fastpath = psock_stateful_fastpath,
        backend = backend,
        requested_backend = requested_backend,
        cuda_fastpath = cuda_fastpath,
        quality = quality,
        cont_vars = cont,
        bincat_vars = bincat,
        target_mean = if (theoretical) as.numeric(mean) else numeric(0),
        target_sd = if (theoretical) as.numeric(sd) else numeric(0),
        target_perc = if (theoretical) as.numeric(perc) else numeric(0),
        cont_dist = if (ccnt > 0L) stats::setNames(as.character(cont_dist), cont) else character(0),
        exact = exact,
        weighted = weights_provided,
        retain = retain
      )
    ),
    class = "repsample_result"
  )
}

#' Print a `repsample_result`
#'
#' @param x A `repsample_result` object.
#' @param ... Unused.
#'
#' @export
print.repsample_result <- function(x, ...) {
  cat("repsample_result\n")
  cat(sprintf("  Mode: %s\n", x$meta$mode))
  cat(sprintf("  Backend: %s\n", x$meta$backend))
  if (!is.null(x$meta$requested_backend) &&
      is.character(x$meta$requested_backend) &&
      x$meta$requested_backend != x$meta$backend) {
    cat(sprintf("  Requested backend: %s\n", x$meta$requested_backend))
  }
  if (isTRUE(x$meta$parallel_enabled)) {
    cat(sprintf("  Parallel: %s (%d cores)\n", x$meta$parallel_mode, x$meta$n_cores))
  }
  cat(sprintf("  Selected: %d / %d\n", sum(x$data$repsample == 1L), nrow(x$data)))

  if (length(x$r) > 0L) {
    cat(sprintf("  Fisher combined p: %.6g\n", unname(x$r[["p"]])))
    cat(sprintf("  Fisher combined chi2: %.6g\n", unname(x$r[["chi2"]])))
    cat(sprintf("  Degrees of freedom: %.0f\n", unname(x$r[["df"]])))
  } else {
    cat("  No deterministic stage statistics (fully random selection).\n")
  }

  invisible(x)
}
