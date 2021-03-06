% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_convergence.r
\name{plot_convergence}
\alias{plot_convergence}
\title{Plot AutoSpill convergence}
\usage{
plot_convergence(convergence.log, popnegpop.slope, asp)
}
\arguments{
\item{convergence.log}{Dataframe with convergence data of AutoSpill.}

\item{popnegpop.slope}{Matrix with slope errors resulting from calculating
spillover with positive and negative populations. Optional parameter, it
can be \code{NULL}.}

\item{asp}{List with AutoSpill parameters.}
}
\value{
\code{NULL}.
}
\description{
Plots convergence of iterative refinement of the spillover matrix.
}
\references{
Roca \emph{et al}:
AutoSpill: A method for calculating spillover coefficients to compensate
or unmix high-parameter flow cytometry data.
\emph{bioRxiv} 2020.06.29.177196;
\href{https://doi.org/10.1101/2020.06.29.177196}{doi:10.1101/2020.06.29.177196}
(2020).
}
\seealso{
\code{\link{refine.spillover}}, \code{\link{process.posnegpop}},
and \code{\link{get.autospill.param}}.
}
