\name{repsample_search_auto}
\alias{repsample_search_auto}
\title{One-Command Multi-Stage Search Wrapper for \code{repsample()}}
\usage{
repsample_search_auto(
  data,
  size,
  cont = NULL,
  bincat = NULL,
  mean = NULL,
  sd = NULL,
  perc = NULL,
  n_stages = 2,
  n_seeds = 24,
  seed_start = 1,
  seeds_stage1 = NULL,
  top_k = 4,
  refine_radius = NULL,
  radius_shrink = 0.5,
  allow_repeats = FALSE,
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
  objective = NULL,
  stop_loss = NULL,
  mean_tol = NULL,
  sd_tol = NULL,
  ks_tol = NULL,
  perc_tol = NULL,
  checkpoint_file = NULL,
  resume = FALSE,
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

\item{n_stages}{Number of search stages.}

\item{n_seeds}{Seed budget per stage. Supply one integer to reuse across
all stages, or one value per stage.}

\item{seed_start}{Starting seed for stage 1 when \code{seeds_stage1} is not provided.}

\item{seeds_stage1}{Optional explicit seed vector for stage 1.}

\item{top_k}{Number of best seeds from each stage used to build the next
stage's refinement neighborhood.}

\item{refine_radius}{Refinement radius around each retained seed for
later stages. Provide one value to decay by \code{radius_shrink}, or one
value per refinement stage.}

\item{radius_shrink}{Multiplicative decay for scalar \code{refine_radius}.}

\item{allow_repeats}{If \code{TRUE}, allow the same seed to be re-evaluated in
multiple stages.}

\item{n_outer_workers}{Number of workers for each stage's outer search.}

\item{outer_parallel}{Outer parallel backend for each stage:
\code{"auto"}, \code{"serial"}, \code{"multicore"}, or \code{"psock"}.}

\item{method}{Sampling method:
\code{"auto"} (default), \code{"greedy"}, \code{"importance"}, or \code{"nearest"}.}

\item{nearest_replace}{Logical flag used when \code{method = "nearest"}:
\code{FALSE} (default) matches without replacement; \code{TRUE} allows replacement.}

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

\item{objective}{Optional custom objective function passed to
\code{repsample_search()}.}

\item{stop_loss}{Optional scalar threshold for early exit. If the best
stage loss is less than or equal to \code{stop_loss}, the remaining stages
are skipped.}

\item{mean_tol}{Optional tolerance on absolute mean error for early stop.}

\item{sd_tol}{Optional tolerance on absolute SD error for early stop.}

\item{ks_tol}{Optional tolerance on KS statistic for early stop.}

\item{perc_tol}{Optional tolerance on absolute percentage-point error.}

\item{checkpoint_file}{Optional path to an RDS checkpoint file. If
provided, stage state is saved after each completed stage.}

\item{resume}{If \code{TRUE}, resume from \code{checkpoint_file} when it exists.}

\item{keep_all}{If \code{TRUE}, retain all per-seed fit objects across all stages.}

\item{...}{Additional arguments passed to \code{repsample()}.}
}
\value{
Object of class \code{repsample_search_result}. It follows the same
structure as \code{repsample_search()}, with \code{summary} including a \code{stage}
column and \code{meta$n_stages} describing the multi-stage run. Early-exit
details are returned in \code{meta$stopped_early} and \code{meta$stopped_stage}.
}
\description{
Runs a coarse-to-fine outer seed search in one call by chaining
\code{repsample_search()} stages. Stage 1 performs a broad scan; each later
stage refines around the best seeds from the previous stage.
}
