% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_utils.R
\name{sim_use_event}
\alias{sim_use_event}
\title{Update data and time given the event}
\usage{
sim_use_event(data, clade, l_2, event, deltas)
}
\arguments{
\item{data}{contains all the information about the simulated process}

\item{clade}{the id of the clade}

\item{l_2}{the matrix containing the information about how the subclades are
nested into the main clade. See sls_sim.get_standard_l_2() for more info.}

\item{event}{the event occurring in the simulated process at a given time}

\item{deltas}{in the Doob-Gillespie algorithm,
the collection of delta_n and delta_t, which are, respectively,
the change in number of species and
the waiting time for the next event to occur}
}
\value{
data
}
\description{
sim module
}
\author{
Giovanni Laudanno
}
