test_that("nearest-neighbor detection behaves as expected", {
  expect_true(RepsampleR:::.use_nearest_neighbor_sampling(
    cont = c("x", "y"),
    bincat = NULL,
    mean = c(1, -0.5),
    sd = c(1, 0.8),
    dist = c("normal", "normal"),
    exact = FALSE
  ))

  expect_false(RepsampleR:::.use_nearest_neighbor_sampling(
    cont = c("x", "y"),
    bincat = "b",
    mean = c(1, -0.5),
    sd = c(1, 0.8),
    dist = c("normal", "normal"),
    exact = FALSE
  ))
})

test_that("repsample_easy single mode supports nearest-neighbor multi-variable matching", {
  set.seed(401)
  dat <- data.frame(
    x = rnorm(1200, 0, 1),
    y = rnorm(1200, 0, 1)
  )

  out <- repsample_easy(
    data = dat,
    size = 120,
    cont = c("x", "y"),
    mean = c(1, -0.5),
    sd = c(1, 0.8),
    dist = c("normal", "normal"),
    mode = "single",
    method = "nearest",
    seed = 10
  )

  expect_s3_class(out, "repsample_result")
  expect_equal(sum(out$data$repsample), 120)
  expect_equal(out$meta$method, "nearest_neighbor")
  expect_equal(out$meta$backend, "nearest_neighbor")
})

test_that("nearest-neighbor single mode is reproducible for fixed seed", {
  set.seed(402)
  dat <- data.frame(
    x = rnorm(800, 0, 1),
    y = rnorm(800, 0, 1)
  )

  out1 <- repsample_easy(
    data = dat,
    size = 80,
    cont = c("x", "y"),
    mean = c(1, -0.5),
    sd = c(1, 0.8),
    dist = c("normal", "normal"),
    mode = "single",
    method = "nearest",
    seed = 99
  )

  out2 <- repsample_easy(
    data = dat,
    size = 80,
    cont = c("x", "y"),
    mean = c(1, -0.5),
    sd = c(1, 0.8),
    dist = c("normal", "normal"),
    mode = "single",
    method = "nearest",
    seed = 99
  )

  expect_equal(out1$selected_rows, out2$selected_rows)
})

test_that("repsample_search supports nearest-neighbor mode", {
  set.seed(403)
  dat <- data.frame(
    x = rnorm(1000, 0, 1),
    y = rnorm(1000, 0, 1)
  )

  out <- repsample_search(
    data = dat,
    size = 100,
    cont = c("x", "y"),
    mean = c(1, -0.5),
    sd = c(1, 0.8),
    dist = c("normal", "normal"),
    method = "nearest",
    n_seeds = 4,
    seed_start = 200,
    n_outer_workers = 1
  )

  expect_s3_class(out, "repsample_search_result")
  expect_equal(nrow(out$summary), 4)
  expect_equal(out$meta$method, "nearest_neighbor")
  expect_equal(out$meta$backend, "nearest_neighbor")
  expect_equal(sum(out$best$data$repsample), 100)
})

test_that("forced nearest mode errors when requirements are not met", {
  set.seed(404)
  dat <- data.frame(
    x = rnorm(500, 0, 1),
    b = rbinom(500, 1, 0.5)
  )

  expect_error(
    repsample_search(
      data = dat,
      size = 80,
      cont = "x",
      bincat = "b",
      mean = 1,
      sd = 1,
      dist = "normal",
      method = "nearest",
      n_seeds = 2
    ),
    "cannot be used"
  )
})
