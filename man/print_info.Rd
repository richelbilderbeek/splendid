% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main_utils.R
\name{print_info}
\alias{print_info}
\title{Print information about simulated data}
\usage{
print_info(brts, n_0, cond, verbose)
}
\arguments{
\item{brts}{branchin times}

\item{n_0}{starting number of lineages}

\item{cond}{type of conditioning:
\itemize{
  \item cond = 0 no conditiong;
  \item cond = 1 conditions on the survival of crown descendents;
  \item cond = 2 not available;
  \item cond = 3 conditions on the survival of subclade and
  on the other crown descendents in the main clade;
  \item cond = 4 conditions on the survival of the subclade and
  both crown descendents in the main clade;
}}

\item{verbose}{set it to TRUE if you want to see the outputs on screen}
}
\value{
nothing
}
\description{
Print information about simulated data
}
\author{
Giovanni Laudanno
}
