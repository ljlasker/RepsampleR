test_that("tolerance settings can stop outer search early", {
  set.seed(601)
  dat <- data.frame(x = rnorm(500, 0, 1))

  out <- repsample_search(
    data = dat,
    size = 80,
    cont = "x",
    mean = 1,
    sd = 1,
    n_seeds = 10,
    seed_start = 1,
    n_outer_workers = 1,
    quality = "fast",
    mean_tol = 0.9,
    sd_tol = 0.9,
    ks_tol = 0.4
  )

  expect_s3_class(out, "repsample_search_result")
  expect_true(isTRUE(out$meta$stopped_early_by_tolerance))
  expect_lt(nrow(out$summary), 10)
})

test_that("nearest-neighbor distance controls and caliper options work", {
  set.seed(602)
  dat <- data.frame(
    x = rnorm(600, 0, 1),
    y = rnorm(600, 0, 1)
  )

  out <- repsample_search(
    data = dat,
    size = 60,
    cont = c("x", "y"),
    mean = c(1, -0.5),
    sd = c(1, 0.8),
    dist = c("normal", "normal"),
    method = "nearest",
    nearest_distance = "weighted",
    nearest_feature_weights = c(2, 0.5),
    nearest_caliper = 4,
    n_seeds = 3,
    seed_start = 20,
    n_outer_workers = 1
  )

  expect_s3_class(out, "repsample_search_result")
  expect_equal(out$meta$nearest_distance, "weighted")
  expect_equal(out$best$meta$nearest_distance, "weighted")
  expect_equal(sum(out$best$data$repsample), 60)

  expect_error(
    repsample_easy(
      data = dat,
      size = 40,
      cont = c("x", "y"),
      mean = c(1, -0.5),
      sd = c(1, 0.8),
      dist = c("normal", "normal"),
      method = "nearest",
      nearest_distance = "euclidean",
      nearest_caliper = 1e-8,
      nearest_caliper_strict = TRUE,
      mode = "single",
      seed = 10
    ),
    "No candidate satisfied"
  )
})

test_that("nearest optimal matching and hnsw backend controls behave safely", {
  set.seed(605)
  dat <- data.frame(
    x = rnorm(300, 0, 1),
    y = rnorm(300, 0, 1)
  )

  if (requireNamespace("clue", quietly = TRUE)) {
    out_opt <- repsample_search(
      data = dat,
      size = 50,
      cont = c("x", "y"),
      mean = c(1, -0.5),
      sd = c(1, 0.8),
      dist = c("normal", "normal"),
      method = "nearest",
      nearest_match = "optimal",
      n_seeds = 2,
      n_outer_workers = 1
    )
    expect_s3_class(out_opt, "repsample_search_result")
    expect_equal(out_opt$best$meta$nearest_match, "optimal")
  }

  if (!requireNamespace("RcppHNSW", quietly = TRUE)) {
    expect_warning(
      out_h <- repsample_search(
        data = dat,
        size = 50,
        cont = c("x", "y"),
        mean = c(1, -0.5),
        sd = c(1, 0.8),
        dist = c("normal", "normal"),
        method = "nearest",
        nearest_backend = "hnsw",
        n_seeds = 2,
        n_outer_workers = 1
      ),
      "HNSW backend requested"
    )
    expect_s3_class(out_h, "repsample_search_result")
    expect_equal(out_h$best$meta$nearest_backend, "exact")
  }
})

test_that("repsample_fit unified API runs search and adaptive modes", {
  set.seed(603)
  dat <- data.frame(
    x = rnorm(800, 0, 1),
    y = rnorm(800, 0, 1)
  )

  out_search <- repsample_fit(
    data = dat,
    size = 90,
    cont = c("x", "y"),
    mean = c(1, -0.5),
    sd = c(1, 0.8),
    dist = c("normal", "normal"),
    method = "nearest",
    mode = "search",
    control = list(
      n_seeds = 4,
      n_outer_workers = 1,
      nearest_distance = "mahalanobis"
    )
  )
  expect_s3_class(out_search, "repsample_search_result")
  expect_equal(out_search$meta$nearest_distance, "mahalanobis")

  out_adapt <- repsample_fit(
    data = dat,
    size = 90,
    cont = c("x", "y"),
    mean = c(1, -0.5),
    sd = c(1, 0.8),
    dist = c("normal", "normal"),
    method = "adaptive",
    mode = "search",
    control = list(
      n_seeds = 6,
      adaptive_pilot_n = 2,
      n_outer_workers = 1
    )
  )
  expect_s3_class(out_adapt, "repsample_search_result")
  expect_equal(out_adapt$meta$method, "adaptive")
  expect_true(nrow(out_adapt$meta$pilot) >= 1)
})

test_that("cuda fallback downgrades to CPU when CUDA is unavailable", {
  set.seed(604)
  dat <- data.frame(x = rnorm(300, 0, 1))

  if (!RepsampleR:::cuda_backend_available()) {
    expect_warning(
      out <- repsample(
        data = dat,
        size = 60,
        cont = "x",
        mean = 1,
        sd = 1,
        backend = "cuda",
        cuda_fallback = TRUE,
        quality = "fast"
      ),
      "falling back to CPU"
    )
    expect_equal(out$meta$backend, "cpu")
    expect_equal(out$meta$requested_backend, "cuda")
  } else {
    expect_no_error(
      repsample(
        data = dat,
        size = 60,
        cont = "x",
        mean = 1,
        sd = 1,
        backend = "cuda",
        cuda_fallback = TRUE,
        quality = "fast"
      )
    )
  }
})
