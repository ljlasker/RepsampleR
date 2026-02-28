\name{repsample_fit}
\alias{repsample_fit}
\title{One-Command RepsampleR API with Simplified Defaults}
\usage{
repsample_fit(
  data,
  size,
  cont = NULL,
  bincat = NULL,
  mean = NULL,
  sd = NULL,
  perc = NULL,
  dist = "normal",
  quality = c("balanced", "fast", "strict", "manual"),
  method = c("auto", "adaptive", "greedy", "importance", "nearest"),
  mode = c("auto", "single", "search", "search_auto"),
  control = NULL,
  ...
)
}
\arguments{
\item{data}{Data frame.}

\item{size}{Desired sample size.}

\item{cont}{Continuous variable names.}

\item{bincat}{Binary/categorical variable names.}

\item{mean}{Target means for \code{cont}.}

\item{sd}{Target SDs for \code{cont}.}

\item{perc}{Target percentages for \code{bincat}.}

\item{dist}{Target distribution form(s): \code{"normal"}, \code{"lognormal"}, \code{"poisson"}.}

\item{quality}{Tuning preset for greedy mode: \code{"fast"}, \code{"balanced"},
\code{"strict"}, or \code{"manual"}.}

\item{method}{Sampling method: \code{"auto"}, \code{"adaptive"}, \code{"greedy"},
\code{"importance"}, or \code{"nearest"}.}

\item{mode}{Execution mode: \code{"auto"}, \code{"single"}, \code{"search"}, or \code{"search_auto"}.}

\item{control}{Named list of advanced controls (seeds, parallelism,
nearest options, tolerances, checkpointing, etc.).}

\item{...}{Additional arguments forwarded to lower-level engines.}
}
\value{
A \code{repsample_result} (\code{single}) or \code{repsample_search_result}
(\code{search}, \code{search_auto}, \code{adaptive}).
}
\description{
\code{repsample_fit()} is a high-level entry point that wraps
\code{repsample_easy()}, \code{repsample_search()}, and \code{repsample_search_auto()}
behind a single interface and a \code{control} list.
}
