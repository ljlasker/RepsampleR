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
  objective = NULL,
  stop_loss = NULL,
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

\item{objective}{Optional custom objective function passed to
\code{repsample_search()}.}

\item{stop_loss}{Optional scalar threshold for early exit. If the best
stage loss is less than or equal to \code{stop_loss}, the remaining stages
are skipped.}

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
