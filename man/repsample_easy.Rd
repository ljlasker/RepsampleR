\name{repsample_easy}
\alias{repsample_easy}
\title{Simplified Front-End for Common RepsampleR Runs}
\usage{
repsample_easy(
  data,
  size,
  cont,
  mean,
  sd,
  dist = "normal",
  quality = c("balanced", "fast", "strict", "manual"),
  mode = c("single", "search"),
  seed = 7,
  n_seeds = 16,
  seed_start = 1,
  n_outer_workers = 1,
  n_cores = 1,
  backend = c("cpu", "cuda"),
  method = c("auto", "greedy", "importance", "nearest"),
  nearest_replace = FALSE,
  nearest_distance = c("euclidean", "mahalanobis", "weighted"),
  nearest_backend = c("auto", "exact", "hnsw"),
  nearest_match = c("greedy", "optimal"),
  nearest_feature_weights = NULL,
  nearest_caliper = NULL,
  nearest_caliper_strict = FALSE,
  nearest_optimal_max_size = 2000L,
  nearest_hnsw_k = 64L,
  nearest_hnsw_M = 16L,
  nearest_hnsw_ef_construction = 200L,
  nearest_hnsw_ef_search = 64L,
  mean_tol = NULL,
  sd_tol = NULL,
  ks_tol = NULL,
  perc_tol = NULL,
  cuda_fallback = TRUE,
  ...
)
}
\arguments{
\item{data}{Data frame.}

\item{size}{Desired sample size.}

\item{cont}{Continuous variable names.}

\item{mean}{Target means for \code{cont}.}

\item{sd}{Target standard deviations for \code{cont}.}

\item{dist}{Target distribution form(s): \code{"normal"} (default),
\code{"lognormal"}, \code{"poisson"}.}

\item{quality}{Tuning preset passed to \code{repsample()}: \code{"fast"},
\code{"balanced"} (default), \code{"strict"}, or \code{"manual"}.}

\item{mode}{\code{"single"} runs one seed; \code{"search"} runs multi-seed outer
search with \code{repsample_search()}.}

\item{seed}{Seed for \code{"single"} mode.}

\item{n_seeds}{Number of seeds in \code{"search"} mode.}

\item{seed_start}{Starting seed in \code{"search"} mode.}

\item{n_outer_workers}{Outer workers in \code{"search"} mode.}

\item{n_cores}{Inner per-iteration CPU cores (passed to \code{repsample()}).}

\item{backend}{Compute backend (\code{"cpu"} or \code{"cuda"}).}

\item{method}{Sampling method:
\code{"auto"} (default), \code{"greedy"}, \code{"importance"}, or \code{"nearest"}.}

\item{nearest_replace}{Logical flag used when \code{method = "nearest"}:
\code{FALSE} (default) matches without replacement; \code{TRUE} allows replacement
(duplicate draws are recorded in \code{data$repsample_n}).}

\item{nearest_distance}{Distance metric for \code{method = "nearest"}:
\code{"euclidean"} (default), \code{"mahalanobis"}, or \code{"weighted"}.}

\item{nearest_backend}{Nearest-neighbor search backend:
\code{"auto"} (default), \code{"exact"}, or \code{"hnsw"} (requires optional
package \code{RcppHNSW}).}

\item{nearest_match}{Nearest-neighbor assignment mode:
\code{"greedy"} (default) or \code{"optimal"} (global assignment; requires
optional package \code{clue}).}

\item{nearest_feature_weights}{Optional positive weights for
\code{nearest_distance = "weighted"} (length 1 or length \code{cont}).}

\item{nearest_caliper}{Optional positive caliper(s) on standardized
per-variable absolute distances for \code{method = "nearest"}.}

\item{nearest_caliper_strict}{If \code{TRUE}, error when no candidate satisfies
the caliper for a draw; otherwise fallback to nearest candidate.}

\item{nearest_optimal_max_size}{Maximum allowed \code{size} for
\code{nearest_match = "optimal"} before automatically falling back to greedy.}

\item{nearest_hnsw_k}{Candidate count queried per draw when using HNSW.}

\item{nearest_hnsw_M}{HNSW graph degree parameter.}

\item{nearest_hnsw_ef_construction}{HNSW build-time search width.}

\item{nearest_hnsw_ef_search}{HNSW query-time search width.}

\item{mean_tol}{Optional tolerance on absolute mean error for early stopping
in \code{mode = "search"}.}

\item{sd_tol}{Optional tolerance on absolute SD error for early stopping in
\code{mode = "search"}.}

\item{ks_tol}{Optional tolerance on KS statistic for early stopping in
\code{mode = "search"}.}

\item{perc_tol}{Optional tolerance on absolute percentage-point error for
binary targets in \code{mode = "search"}.}

\item{cuda_fallback}{If \code{TRUE} (default), automatically fall back to CPU
when CUDA is unavailable or a CUDA fast-path cannot be used.}

\item{...}{Additional arguments forwarded to \code{repsample()} (for \code{"single"})
or to \code{repsample_search()} (for \code{"search"}), such as \code{rrule},
\code{srule}, \code{randomperc}, \code{objective}, etc.}
}
\value{
A \code{repsample_result} (\code{mode = "single"}) or
\code{repsample_search_result} (\code{mode = "search"}).
}
\description{
A convenience wrapper for quickly targeting mean/SD (and optional
distribution form) with simple quality presets, without manually setting
\code{randomperc}, \code{rrule}, and \code{srule}.
}
