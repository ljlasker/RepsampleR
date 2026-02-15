test_that("importance-sampling detection behaves as expected", {
  expect_true(RepsampleR:::.use_importance_sampling(
    cont = "x",
    bincat = NULL,
    mean = 1,
    sd = 1,
    dist = "normal",
    exact = FALSE
  ))

  expect_false(RepsampleR:::.use_importance_sampling(
    cont = "x",
    bincat = "b",
    mean = 1,
    sd = 1,
    dist = "normal",
    exact = FALSE
  ))

  expect_false(RepsampleR:::.use_importance_sampling(
    cont = "x",
    bincat = NULL,
    mean = 1,
    sd = 1,
    dist = "normal",
    exact = TRUE
  ))
})

test_that("repsample_easy auto-dispatches to importance sampling in single mode", {
  set.seed(301)
  dat <- data.frame(x = rnorm(3000, 0, 1))

  out <- repsample_easy(
    data = dat,
    size = 300,
    cont = "x",
    mean = 1,
    sd = 1,
    dist = "normal",
    mode = "single",
    method = "auto",
    seed = 13
  )

  expect_s3_class(out, "repsample_result")
  expect_equal(sum(out$data$repsample), 300)
  expect_equal(out$meta$method, "importance_sampling")
  expect_true(is.finite(out$meta$ks))
  expect_true(is.finite(out$meta$ess))
})

test_that("importance-sampling single mode is reproducible for fixed seed", {
  set.seed(302)
  dat <- data.frame(x = rnorm(2500, 0, 1))

  out1 <- repsample_easy(
    data = dat,
    size = 250,
    cont = "x",
    mean = 1,
    sd = 1,
    dist = "normal",
    mode = "single",
    method = "importance",
    seed = 99
  )

  out2 <- repsample_easy(
    data = dat,
    size = 250,
    cont = "x",
    mean = 1,
    sd = 1,
    dist = "normal",
    mode = "single",
    method = "importance",
    seed = 99
  )

  expect_equal(out1$selected_rows, out2$selected_rows)
})

test_that("forced importance mode errors when requirements are not met", {
  set.seed(303)
  dat <- data.frame(
    x = rnorm(300),
    b = rbinom(300, 1, 0.5)
  )

  expect_error(
    repsample_search(
      data = dat,
      size = 60,
      cont = "x",
      bincat = "b",
      mean = 1,
      sd = 1,
      dist = "normal",
      method = "importance",
      n_seeds = 3,
      n_outer_workers = 1
    ),
    "cannot be used"
  )
})

test_that("repsample_search auto-dispatches to importance sampling", {
  set.seed(304)
  dat <- data.frame(x = rnorm(4000, 0, 1))

  out <- repsample_search(
    data = dat,
    size = 400,
    cont = "x",
    mean = 1,
    sd = 1,
    dist = "normal",
    method = "auto",
    n_seeds = 5,
    seed_start = 200,
    n_outer_workers = 1
  )

  expect_s3_class(out, "repsample_search_result")
  expect_equal(out$meta$method, "importance_sampling")
  expect_equal(out$meta$backend, "importance_sampling")
  expect_equal(nrow(out$summary), 5)
  expect_equal(sum(out$best$data$repsample), 400)
})

test_that("auto mode falls back to greedy when overlap is poor", {
  set.seed(305)
  dat <- data.frame(x = rnorm(600, 0, 1))

  out <- NULL
  expect_warning(
    out <- repsample_search(
      data = dat,
      size = 80,
      cont = "x",
      mean = 8,
      sd = 0.1,
      dist = "normal",
      method = "auto",
      n_seeds = 2,
      n_outer_workers = 1,
      randomperc = 20,
      rrule = 120
    ),
    "ESS"
  )

  expect_s3_class(out, "repsample_search_result")
  expect_false(identical(out$meta$backend, "importance_sampling"))
})
