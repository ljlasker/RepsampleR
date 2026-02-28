\name{repsample_search}
\alias{repsample_search}
\title{Outer Seed Search Preset for \code{repsample()}}
\usage{
repsample_search(
  data,
  size,
  cont = NULL,
  bincat = NULL,
  mean = NULL,
  sd = NULL,
  perc = NULL,
  seeds = NULL,
  n_seeds = 16,
  seed_start = 1,
  n_outer_workers = 1,
  outer_parallel = c("auto", "serial", "multicore", "psock"),
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
  objective = NULL,
  keep_all = FALSE,
  ...
)
}
\arguments{
\item{data}{Data frame.}

\item{size}{Desired sample size.}

\item{cont}{Continuous variable names.}

\item{bincat}{Binary/categorical variable names.}

\item{mean}{Means for continuous variables (theoretical sampling only).}

\item{sd}{Standard deviations for continuous variables (theoretical sampling only).}

\item{perc}{Percentages for binary variables (theoretical sampling only).}

\item{seeds}{Optional explicit seed vector.}

\item{n_seeds}{Number of seeds to search when \code{seeds} is not provided.}

\item{seed_start}{Starting seed when generating seeds.}

\item{n_outer_workers}{Number of workers for outer seed search.}

\item{outer_parallel}{Outer parallel backend: \code{"auto"}, \code{"serial"},
\code{"multicore"}, or \code{"psock"}.}

\item{method}{Sampling method:
\code{"auto"} (default), \code{"greedy"}, \code{"importance"}, or \code{"nearest"}.}

\item{nearest_replace}{Logical flag used when \code{method = "nearest"}:
\code{FALSE} (default) matches without replacement; \code{TRUE} allows replacement
(duplicate draws are recorded in \code{best$data$repsample_n}).}

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

\item{mean_tol}{Optional tolerance on absolute mean error for early stop.}

\item{sd_tol}{Optional tolerance on absolute SD error for early stop.}

\item{ks_tol}{Optional tolerance on KS statistic for early stop.}

\item{perc_tol}{Optional tolerance on absolute percentage-point error for
binary targets in theoretical mode.}

\item{objective}{Optional custom objective function that takes a
\code{repsample_result} and returns a scalar loss (smaller is better).}

\item{keep_all}{If \code{TRUE}, return all per-seed fit objects.}

\item{...}{Additional arguments passed to \code{repsample()}.}
}
\value{
Object of class \code{repsample_search_result} with:
\itemize{
\item \code{best}: best \code{repsample_result}
\item \code{best_seed}: seed with smallest loss
\item \code{best_loss}: smallest loss
\item \code{summary}: data frame of all searched seeds and losses
\item \code{all}: optional full per-seed results when \code{keep_all = TRUE}
\item \code{meta}: search metadata
}
}
\description{
Runs multiple \code{repsample()} calls over different seeds and returns the
best fit under a scalar loss function. This provides a concise preset for
outer parallelization.
}
