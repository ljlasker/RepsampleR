#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (!requireNamespace("RepsampleR", quietly = TRUE)) {
  stop("Please install RepsampleR first (e.g., remotes::install_github(\"ljlasker/RepsampleR\")).", call. = FALSE)
}
suppressPackageStartupMessages(library(RepsampleR))

get_arg <- function(name, default = NULL) {
  key <- paste0("--", name, "=")
  hit <- args[startsWith(args, key)]
  if (length(hit) < 1L) {
    return(default)
  }
  sub(key, "", hit[[1L]], fixed = TRUE)
}

as_int <- function(x, default) {
  if (is.null(x)) {
    return(as.integer(default))
  }
  as.integer(x)
}

as_num <- function(x, default) {
  if (is.null(x)) {
    return(as.numeric(default))
  }
  as.numeric(x)
}

set.seed(as_int(get_arg("seed", "123"), 123))
n <- as_int(get_arg("n", "10000"), 10000)
size <- as_int(get_arg("size", "1000"), 1000)
n_seeds <- as_int(get_arg("n_seeds", "8"), 8)
n_workers <- as_int(get_arg("n_outer_workers", "1"), 1)
outfile <- get_arg("outfile", "benchmark_repsample_methods.csv")

baseline <- data.frame(
  x = rnorm(n, 0, 1),
  y = rnorm(n, 0, 1)
)

target_mean <- c(
  as_num(get_arg("mean_x", "1"), 1),
  as_num(get_arg("mean_y", "-0.5"), -0.5)
)
target_sd <- c(
  as_num(get_arg("sd_x", "1"), 1),
  as_num(get_arg("sd_y", "0.8"), 0.8)
)

methods <- c("importance", "nearest", "greedy")
rows <- list()

for (m in methods) {
  t <- system.time({
    fit <- repsample_fit(
      data = baseline,
      size = size,
      cont = c("x", "y"),
      mean = target_mean,
      sd = target_sd,
      dist = c("normal", "normal"),
      method = m,
      mode = "search",
      control = list(
        n_seeds = n_seeds,
        n_outer_workers = n_workers,
        nearest_backend = "auto",
        nearest_match = "greedy"
      )
    )
  })["elapsed"]

  q <- repsample_quality(
    out = fit$best,
    cont = c("x", "y"),
    mean = target_mean,
    sd = target_sd
  )

  rows[[length(rows) + 1L]] <- data.frame(
    method = m,
    best_seed = as.integer(fit$best_seed),
    best_loss = as.numeric(fit$best_loss),
    quality_loss = as.numeric(q$loss),
    ks = as.numeric(q$ks),
    elapsed_sec = as.numeric(t),
    stringsAsFactors = FALSE
  )
}

out <- do.call(rbind, rows)
out <- out[order(out$quality_loss), , drop = FALSE]
write.csv(out, outfile, row.names = FALSE)
cat("Wrote benchmark results to:", outfile, "\n")
print(out)
