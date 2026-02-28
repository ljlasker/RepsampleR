\name{repsample_quality}
\alias{repsample_quality}
\title{Compute Match Quality Metrics for a RepsampleR Fit}
\usage{
repsample_quality(
  out,
  cont = NULL,
  bincat = NULL,
  mean = NULL,
  sd = NULL,
  perc = NULL,
  ks_weight = 1
)
}
\arguments{
\item{out}{A \code{repsample_result} object.}

\item{cont}{Continuous variable names used for matching.}

\item{bincat}{Binary variable names used for matching.}

\item{mean}{Target means for \code{cont}.}

\item{sd}{Target standard deviations for \code{cont}.}

\item{perc}{Target percentages for \code{bincat}.}

\item{ks_weight}{Weight on the KS component in the composite loss.}
}
\value{
A named list with components \code{loss}, \code{ks}, \code{cont}, and \code{bincat}.
}
\description{
Returns a compact set of fit-quality diagnostics and a scalar composite loss
that is comparable across methods (\code{greedy}, \code{importance}, \code{nearest}).
}
