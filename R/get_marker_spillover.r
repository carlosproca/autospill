# get_marker_spillover.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


#' Get spillover coefficients
#'
#' Calculates spillover coefficients with robust linear models.
#'
#' @param scale.untransformed Logical indicating if calculation has to be run
#'     in untransformed or transformed (bi-exponential) scale.
#' @param flow.gate List of vectors with ids of gated events per sample.
#' @param flow.control List with data and metadata of a set of controls.
#' @param asp List with AutoSpill parameters.
#'
#' @return List of two matrices, with regressions intercepts and coefficients.
#'
#' @references Roca \emph{et al}:
#'     AutoSpill: A method for calculating spillover coefficients to compensate
#'     or unmix high-parameter flow cytometry data.
#'     \emph{bioRxiv} 2020.06.29.177196;
#'     \href{https://doi.org/10.1101/2020.06.29.177196}{doi:10.1101/2020.06.29.177196}
#'     (2020).
#'
#' @seealso \code{\link{gate.flow.data}}, \code{\link{read.flow.control}}, and
#'     \code{\link{get.autospill.param}}.
#'
#' @export

get.marker.spillover <- function( scale.untransformed, flow.gate, flow.control,
    asp )
{
    if ( scale.untransformed )
        expr.data <- flow.control$expr.data.untr
    else
        expr.data <- flow.control$expr.data.tran

    marker.spillover.zero <- rep( 0, flow.control$marker.n )
    names( marker.spillover.zero ) <- flow.control$marker

    marker.spillover <- mclapply( flow.control$sample, function( samp )
    {
        marker.proper <- samp

        marker.proper.expr <- expr.data[
            which( flow.control$event.sample == samp )[ flow.gate[[ samp ]] ],
            marker.proper ]

        marker.expr.n <- length( marker.proper.expr )

        # identify extreme expression values

        expr.trim.n <- round( marker.expr.n * asp$rlm.trim.factor )

        marker.proper.expr.low <- sort( marker.proper.expr )[ expr.trim.n ]
        marker.proper.expr.high <- sort( marker.proper.expr,
            decreasing = TRUE )[ expr.trim.n ]

        marker.spillover.inte <- marker.spillover.zero
        marker.spillover.coef <- marker.spillover.zero

        for ( marker in flow.control$marker )
            if ( marker == marker.proper )
                marker.spillover.coef[ marker ] <- 1.0
            else
            {
                marker.expr <- expr.data[
                    which( flow.control$event.sample == samp )[
                        flow.gate[[ samp ]] ],
                    marker ]

                marker.expr.low <- sort( marker.expr )[ expr.trim.n ]
                marker.expr.high <- sort( marker.expr, decreasing = TRUE )[
                    expr.trim.n ]

                # trim expression values for both markers

                expr.trim.idx <- which (
                    marker.proper.expr > marker.proper.expr.low &
                        marker.proper.expr < marker.proper.expr.high &
                        marker.expr > marker.expr.low &
                        marker.expr < marker.expr.high
                )

                marker.proper.expr.trim <- marker.proper.expr[ expr.trim.idx ]
                marker.expr.trim <- marker.expr[ expr.trim.idx ]

                # fit robust linear model

                spillover.model.result <- fit.robust.linear.model(
                    marker.proper.expr.trim, marker.expr.trim,
                    marker.proper, marker, asp )

                marker.spillover.inte[ marker ] <-
                    spillover.model.result[ 1, 1 ]

                marker.spillover.coef[ marker ] <-
                    spillover.model.result[ 2, 1 ]
            }

        c( marker.spillover.inte, marker.spillover.coef )
    },
    mc.cores = get.worker.process( asp$worker.process.n ) )

    marker.spillover <- do.call( rbind, marker.spillover )
    rownames( marker.spillover ) <- flow.control$marker

    list(
        inte = marker.spillover[ , 1 : flow.control$marker.n ],
        coef = marker.spillover[ , 1 : flow.control$marker.n +
                flow.control$marker.n ]
    )
}

