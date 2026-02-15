test_that("population mode returns requested sample and scalars", {
  set.seed(42)
  n <- 140
  dat <- data.frame(
    x = rnorm(n),
    y = rnorm(n, mean = 10, sd = 3),
    b = rbinom(n, 1, 0.4),
    c = sample(1:4, n, replace = TRUE)
  )

  out <- repsample(
    data = dat,
    size = 30,
    cont = c("x", "y"),
    bincat = c("b", "c"),
    seednum = 10
  )

  expect_s3_class(out, "repsample_result")
  expect_equal(sum(out$data$repsample), 30)
  expect_true(all(c("p", "chi2", "df") %in% names(out$r)))
  expect_true(all(c("x_D", "x_p", "y_D", "y_p") %in% names(out$r)))
  expect_true(all(c("b_p", "b_chi2", "c_p", "c_chi2") %in% names(out$r)))
})

test_that("theoretical mode works with binary proportions", {
  set.seed(1)
  n <- 180
  dat <- data.frame(
    x = rnorm(n, 30, 15),
    y = rnorm(n, 7000, 4000),
    b = rbinom(n, 1, 0.3)
  )

  out <- repsample(
    data = dat,
    size = 35,
    cont = c("x", "y"),
    bincat = "b",
    mean = c(30, 7000),
    sd = c(15, 4000),
    perc = 30,
    seednum = 10
  )

  expect_equal(sum(out$data$repsample), 35)
  expect_true(all(c("b_p", "b_z") %in% names(out$r)))
  expect_false("b_chi2" %in% names(out$r))
})

test_that("retain carries cases forward", {
  set.seed(2)
  n <- 120
  dat <- data.frame(
    x = rnorm(n),
    b = rbinom(n, 1, 0.5)
  )

  dat$ret <- NA_real_
  keep_pool <- sample(seq_len(n), 80)
  dat$ret[keep_pool] <- 0
  carried <- keep_pool[1:8]
  dat$ret[carried] <- 1

  out <- repsample(
    data = dat,
    size = 20,
    cont = "x",
    bincat = "b",
    retain = "ret",
    seednum = 17
  )

  expect_equal(sum(out$data$repsample), 20)
  expect_true(all(out$data$repsample[carried] == 1))
})

test_that("force semantics match Stata behavior", {
  set.seed(3)
  dat <- data.frame(
    x = rnorm(90),
    b = rbinom(90, 1, 0.5),
    repsample = 0L
  )

  expect_error(
    repsample(dat, size = 15, cont = "x", bincat = "b"),
    "already exists"
  )

  out <- repsample(dat, size = 15, cont = "x", bincat = "b", force = TRUE)
  expect_equal(sum(out$data$repsample), 15)
})

test_that("weights validation and weighted run", {
  set.seed(4)
  dat <- data.frame(
    x = rnorm(100),
    y = rnorm(100),
    b = rbinom(100, 1, 0.5)
  )

  expect_error(
    repsample(dat, size = 20, cont = c("x", "y"), bincat = "b", wght = c(40, 40)),
    "match number of variables"
  )

  expect_error(
    repsample(dat, size = 20, cont = c("x", "y"), bincat = "b", wght = c(33, 33, 33)),
    "add up to 100"
  )

  out <- repsample(
    dat,
    size = 20,
    cont = c("x", "y"),
    bincat = "b",
    wght = c(50, 30, 20)
  )

  expect_equal(sum(out$data$repsample), 20)
})

test_that("rrule and srule run", {
  set.seed(5)
  dat <- data.frame(
    x = rnorm(220),
    b = rbinom(220, 1, 0.45)
  )

  out <- repsample(
    dat,
    size = 40,
    cont = "x",
    bincat = "b",
    rrule = 25,
    srule = 0.75,
    seednum = 99
  )

  expect_equal(sum(out$data$repsample), 40)
  expect_true(length(out$r) > 0)
})

test_that("randomperc 100 means no deterministic r() stats", {
  set.seed(6)
  dat <- data.frame(
    x = rnorm(100),
    b = rbinom(100, 1, 0.5)
  )

  out <- repsample(dat, size = 20, cont = "x", bincat = "b", randomperc = 100)
  expect_equal(sum(out$data$repsample), 20)
  expect_length(out$r, 0)
})

test_that("exact mode runs and suppresses chi2/z extras as in Stata", {
  set.seed(7)
  dat <- data.frame(
    x = rnorm(60),
    b = rbinom(60, 1, 0.4)
  )

  pop_exact <- repsample(dat, size = 14, cont = "x", bincat = "b", exact = TRUE)
  expect_equal(sum(pop_exact$data$repsample), 14)
  expect_true("b_p" %in% names(pop_exact$r))
  expect_false("b_chi2" %in% names(pop_exact$r))

  th_exact <- repsample(
    dat,
    size = 14,
    cont = "x",
    bincat = "b",
    mean = 0,
    sd = 1,
    perc = 40,
    exact = TRUE
  )

  expect_equal(sum(th_exact$data$repsample), 14)
  expect_true("b_p" %in% names(th_exact$r))
  expect_false("b_z" %in% names(th_exact$r))
})

test_that("n_cores validates and parallel mode matches serial selection", {
  set.seed(8)
  dat <- data.frame(
    x = rnorm(130),
    b = rbinom(130, 1, 0.5)
  )

  expect_error(
    repsample(dat, size = 20, cont = "x", bincat = "b", n_cores = 0),
    "positive integer"
  )

  expect_error(
    repsample(dat, size = 20, cont = "x", bincat = "b", n_cores = 1.5),
    "positive integer"
  )

  out_serial <- repsample(
    dat,
    size = 20,
    cont = "x",
    bincat = "b",
    seednum = 55,
    n_cores = 1
  )

  out_parallel <- suppressWarnings(repsample(
    dat,
    size = 20,
    cont = "x",
    bincat = "b",
    seednum = 55,
    n_cores = 2
  ))

  expect_equal(out_serial$selected_rows, out_parallel$selected_rows)
})

test_that("RStudio macOS route selects PSOCK mode", {
  if (!identical(Sys.info()[["sysname"]], "Darwin")) {
    skip("macOS-only routing test")
  }

  old_rstudio <- Sys.getenv("RSTUDIO", unset = "")
  on.exit(Sys.setenv(RSTUDIO = old_rstudio), add = TRUE)
  Sys.setenv(RSTUDIO = "1")

  set.seed(12)
  dat <- data.frame(x = rnorm(200))

  expect_warning(
    out <- repsample(
      dat,
      size = 30,
      cont = "x",
      mean = 1,
      sd = 1,
      n_cores = 2
    ),
    "Using PSOCK parallel backend"
  )

  expect_true(out$meta$parallel_mode %in% c("psock", "none"))
})

test_that("backend validates and cuda preflight errors without dependencies", {
  set.seed(9)
  dat <- data.frame(
    x = rnorm(80),
    b = rbinom(80, 1, 0.5)
  )

  expect_error(
    repsample(dat, size = 15, cont = "x", bincat = "b", backend = "gpu"),
    "arg"
  )

  out_cpu <- repsample(
    dat,
    size = 15,
    cont = "x",
    bincat = "b",
    backend = "cpu"
  )
  expect_equal(sum(out_cpu$data$repsample), 15)

  if (!(requireNamespace("reticulate", quietly = TRUE) &&
        reticulate::py_module_available("cupy"))) {
    expect_error(
      repsample(dat, size = 15, cont = "x", bincat = "b", backend = "cuda"),
      "CUDA backend requires"
    )
  }
})

test_that("repsample_search finds best fit and returns expected structure", {
  set.seed(21)
  n <- 220
  dat <- data.frame(x = rnorm(n, 0, 1))

  out <- repsample_search(
    data = dat,
    size = 40,
    cont = "x",
    mean = 1,
    sd = 1,
    n_seeds = 4,
    seed_start = 500,
    n_outer_workers = 1,
    backend = "cpu",
    randomperc = 20,
    rrule = 80,
    srule = 0.75
  )

  expect_s3_class(out, "repsample_search_result")
  expect_true(is.list(out$best))
  expect_equal(sum(out$best$data$repsample), 40)
  expect_equal(nrow(out$summary), 4)
  expect_true(all(c("seed", "loss") %in% names(out$summary)))
})

test_that("repsample_search supports outer parallel auto mode", {
  set.seed(22)
  dat <- data.frame(x = rnorm(260))

  out <- suppressWarnings(repsample_search(
    data = dat,
    size = 50,
    cont = "x",
    mean = 1,
    sd = 1,
    n_seeds = 3,
    seed_start = 800,
    n_outer_workers = 2,
    outer_parallel = "auto",
    backend = "cpu",
    randomperc = 25,
    rrule = 90,
    srule = 0.8
  ))

  expect_s3_class(out, "repsample_search_result")
  expect_equal(nrow(out$summary), 3)
  expect_true(out$meta$outer_mode %in% c("serial", "multicore", "psock"))
})

test_that("repsample_search with cuda backend limits outer workers", {
  set.seed(23)
  dat <- data.frame(x = rnorm(180))

  expect_warning(
    tryCatch(
      repsample_search(
        data = dat,
        size = 30,
        cont = "x",
        mean = 1,
        sd = 1,
        method = "greedy",
        n_seeds = 2,
        n_outer_workers = 4,
        backend = "cuda",
        randomperc = 30,
        rrule = 60
      ),
      error = function(e) e
    ),
    "forcing one outer worker"
  )
})

test_that("repsample_search prevents nested parallel oversubscription", {
  set.seed(231)
  dat <- data.frame(x = rnorm(140))

  out <- NULL
  expect_warning(
    out <- repsample_search(
      data = dat,
      size = 24,
      cont = "x",
      mean = 1,
      sd = 1,
      method = "greedy",
      n_seeds = 2,
      seed_start = 900,
      n_outer_workers = 2,
      outer_parallel = "auto",
      backend = "cpu",
      n_cores = 2,
      randomperc = 20,
      rrule = 40
    ),
    "Nested parallelism detected"
  )

  expect_equal(out$best$meta$n_cores, 1)
})

test_that("repsample_search_auto runs staged search in one call", {
  set.seed(24)
  dat <- data.frame(x = rnorm(240))

  out <- repsample_search_auto(
    data = dat,
    size = 40,
    cont = "x",
    mean = 1,
    sd = 1,
    n_stages = 2,
    n_seeds = c(4, 6),
    seed_start = 900,
    top_k = 2,
    refine_radius = 10,
    n_outer_workers = 1,
    backend = "cpu",
    randomperc = 25,
    rrule = 80,
    srule = 0.8
  )

  expect_s3_class(out, "repsample_search_result")
  expect_equal(sum(out$best$data$repsample), 40)
  expect_true(all(c("stage", "seed", "loss") %in% names(out$summary)))
  expect_equal(out$meta$outer_mode, "multistage")
  expect_equal(out$meta$n_stages, 2)
  expect_equal(nrow(out$stage_summary), 2)
})

test_that("repsample_search_auto supports early exit via stop_loss", {
  set.seed(25)
  dat <- data.frame(x = rnorm(160))

  out <- repsample_search_auto(
    data = dat,
    size = 24,
    cont = "x",
    mean = 1,
    sd = 1,
    n_stages = 3,
    n_seeds = c(3, 3, 3),
    n_outer_workers = 1,
    objective = function(x) 0,
    stop_loss = 0,
    randomperc = 30,
    rrule = 50
  )

  expect_s3_class(out, "repsample_search_result")
  expect_true(isTRUE(out$meta$stopped_early))
  expect_equal(out$meta$stopped_stage, 1)
  expect_equal(out$meta$n_stages, 1)
  expect_equal(out$meta$n_stages_requested, 3)
  expect_equal(length(out$stages), 1)
  expect_equal(nrow(out$stage_summary), 1)
})

test_that("quality preset fills tuning defaults when not provided", {
  set.seed(26)
  dat <- data.frame(x = rnorm(200))

  out <- repsample(
    data = dat,
    size = 30,
    cont = "x",
    mean = 1,
    sd = 1,
    quality = "fast"
  )

  expect_equal(sum(out$data$repsample), 30)
  expect_equal(out$meta$quality, "fast")
  expect_equal(out$meta$randomperc, 70)
  expect_equal(out$meta$srule, 0.7)
  expect_true(is.numeric(out$meta$rrule))
  expect_true(out$meta$rrule >= 10)
})

test_that("theoretical lognormal and poisson distribution targets run", {
  set.seed(27)
  n <- 300

  x_lnorm <- rlnorm(n, meanlog = 0.2, sdlog = 0.5)
  dat_lnorm <- data.frame(x = x_lnorm)

  m_l <- mean(x_lnorm)
  s_l <- stats::sd(x_lnorm)
  out_l <- repsample(
    data = dat_lnorm,
    size = 40,
    cont = "x",
    mean = m_l,
    sd = s_l,
    dist = "lognormal",
    seednum = 11,
    quality = "balanced"
  )
  expect_equal(sum(out_l$data$repsample), 40)

  dat_pois <- data.frame(x = rpois(n, lambda = 4))
  out_p <- repsample(
    data = dat_pois,
    size = 40,
    cont = "x",
    mean = 4,
    sd = 2,
    dist = "poisson",
    seednum = 12,
    quality = "balanced"
  )
  expect_equal(sum(out_p$data$repsample), 40)
})

test_that("repsample_easy supports single and search modes", {
  set.seed(28)
  dat <- data.frame(x = rnorm(260))

  one <- repsample_easy(
    data = dat,
    size = 40,
    cont = "x",
    mean = 1,
    sd = 1,
    mode = "single",
    quality = "balanced",
    seed = 10
  )
  expect_s3_class(one, "repsample_result")
  expect_equal(sum(one$data$repsample), 40)

  many <- repsample_easy(
    data = dat,
    size = 40,
    cont = "x",
    mean = 1,
    sd = 1,
    mode = "search",
    quality = "fast",
    n_seeds = 3,
    seed_start = 100
  )
  expect_s3_class(many, "repsample_search_result")
  expect_equal(nrow(many$summary), 3)
})
