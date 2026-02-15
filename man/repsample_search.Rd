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
