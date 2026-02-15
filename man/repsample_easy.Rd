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
