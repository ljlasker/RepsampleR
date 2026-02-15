args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg))
} else {
  normalizePath(getwd())
}
pkg <- normalizePath(file.path(dirname(script_path), ".."))
lib <- .libPaths()[1L]

if ("package:RepsampleR" %in% search()) {
  detach("package:RepsampleR", unload = TRUE, character.only = TRUE)
}

.libPaths(c(lib, .libPaths()))
install.packages(pkg, repos = NULL, type = "source", lib = lib)
library(RepsampleR, lib.loc = lib)

set.seed(123)
baseline <- data.frame(x = rnorm(10000, mean = 0, sd = 1))
target_n <- 1000

out <- repsample_easy(
  data = baseline,
  size = target_n,
  cont = "x",
  mean = 1,
  sd = 1,
  dist = "normal",
  quality = "balanced",
  mode = "single",
  seed = 7
)

picked <- out$data$x[out$data$repsample == 1L]
cat("n =", length(picked), "mean =", mean(picked), "sd =", sd(picked), "var =", var(picked), "\n")
print(out$meta[c("quality", "parallel_mode", "parallel_enabled", "n_cores")])
