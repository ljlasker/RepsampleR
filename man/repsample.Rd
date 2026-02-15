\name{repsample}
\alias{repsample}
\title{Greedy Representative Sampling}
\usage{
repsample(
  data,
  size,
  cont = NULL,
  bincat = NULL,
  mean = NULL,
  sd = NULL,
  dist = NULL,
  perc = NULL,
  seednum = 7,
  randomperc = 10,
  srule = NULL,
  rrule = NULL,
  quality = c("manual", "fast", "balanced", "strict"),
  wght = NULL,
  retain = NULL,
  exact = FALSE,
  force = FALSE,
  subset = NULL,
  in_rows = NULL,
  n_cores = 1,
  backend = c("cpu", "cuda")
)
}
\arguments{
\item{data}{Data frame.}

\item{size}{Desired sample size.}

\item{cont}{Continuous variable names.}

\item{bincat}{Binary/categorical variable names.}

\item{mean}{Means for continuous variables (theoretical sampling only).}

\item{sd}{Standard deviations for continuous variables (theoretical sampling only).}

\item{dist}{Target distribution form(s) for \code{cont} variables in theoretical
mode. Supported: \code{"normal"} (default), \code{"lognormal"}, \code{"poisson"}.
You can pass one value, one per \code{cont}, or a named vector/list.}

\item{perc}{Percentages for binary variables (theoretical sampling only), in \verb{(0, 100)}.}

\item{seednum}{Random seed. Default is \code{7}.}

\item{randomperc}{Percentage of random picks at start. Default is \code{10}.}

\item{srule}{Early stopping threshold on Fisher combined p-value in \verb{[0.5, 1)}.}

\item{rrule}{Early stopping candidate subsample size per iteration (minimum \code{10}).}

\item{quality}{Optional tuning preset for faster setup:
\code{"manual"} (default), \code{"fast"}, \code{"balanced"}, \code{"strict"}.
Presets only fill in \code{randomperc} / \code{rrule} / \code{srule} when those are not
explicitly provided.}

\item{wght}{Variable weights that must sum to \code{100}; continuous first, then binary/categorical.}

\item{retain}{Optional existing binary 0/1 column to carry forward selected cases.}

\item{exact}{Use exact tests where Stata does (\code{TRUE/FALSE}).}

\item{force}{Force overwrite of existing \code{repsample} column.}

\item{subset}{Optional logical expression evaluated in \code{data} (Stata \code{if} analog).}

\item{in_rows}{Optional integer row indices to keep (Stata \verb{in} analog).}

\item{n_cores}{Number of CPU cores for candidate scoring in each greedy
iteration. \code{1} (default) is serial. Values greater than \code{1} use
fork-based parallelism in terminal R on Unix-like systems and PSOCK
clusters in environments where forking is unstable (for example,
RStudio on macOS).}

\item{backend}{Compute backend. \code{"cpu"} uses the standard C++/CPU path.
\code{"cuda"} enables an additional CUDA fast-path for theoretical sampling
with continuous variables and \code{exact = FALSE} (uses Python \code{cupy}
through \code{reticulate}).}
}
\value{
An object of class \code{repsample_result} with:
\itemize{
\item \code{data}: input data plus \code{repsample} (\code{0/1})
\item \code{r}: named numeric vector equivalent to Stata \code{r()} scalars (empty if fully random)
\item \code{selected_rows}: row indices selected in the returned data
\item \code{meta}: run metadata
}
}
\description{
\code{repsample()} performs greedy representative sampling with support for
population and theoretical targets, retain/force behavior, exact/asymptotic
testing, variable weights, and early stopping rules.
}
