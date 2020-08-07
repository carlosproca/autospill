# get_marker_spillover_posnegpop.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


#' Get spillover coefficients using positive and negative populations
#'
#' Calculates spillover coefficients using positive and negative populations.
#'
#' @param flow.gate List of vectors with ids of gated events per sample.
#' @param flow.control List with data and metadata of a set of controls.
#' @param asp List with AutoSpill parameters.
#'
#' @return Matrix with spillover coefficients.
#'
#' @references Roca \emph{et al}:
#'     AutoSpill: A method for calculating spillover coefficients to compensate
#'     or unmix high-parameter flow cytometry data.
#'     \emph{bioRxiv} 2020.06.29.177196;
#'     \href{https://doi.org/10.1101/2020.06.29.177196}{doi:10.1101/2020.06.29.177196}
#'     (2020).
#'
#' @seealso \code{\link{get.marker.spillover}}, \code{\link{gate.flow.data}},
#'     \code{\link{read.flow.control}}, and \code{\link{get.autospill.param}}.
#'
#' @export

get.marker.spillover.posnegpop <- function( flow.gate, flow.control, asp )
{
    marker.spillpopu.zero <- rep( 0, flow.control$marker.n )
    names( marker.spillpopu.zero ) <- flow.control$marker

    marker.spillpopu <- mclapply( flow.control$sample, function( samp )
    {
        marker.proper <- samp

        marker.proper.expr <- flow.control$expr.data.tran[
            which( flow.control$event.sample == samp )[ flow.gate[[ samp ]] ],
            marker.proper ]

        marker.expr.n <- length( marker.proper.expr )

        # identify extreme expression values

        expr.trim.n <- round( marker.expr.n * asp$rlm.trim.factor )

        marker.proper.expr.low <- sort( marker.proper.expr )[ expr.trim.n ]
        marker.proper.expr.high <- sort( marker.proper.expr,
            decreasing = TRUE )[ expr.trim.n ]

        marker.spillpopu.coef <- marker.spillpopu.zero

        for ( marker in flow.control$marker )
            if ( marker == marker.proper )
                marker.spillpopu.coef[ marker ] <- 1.0
            else
            {
                marker.expr <- flow.control$expr.data.tran[
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

                # identify distribution modes with 2-means clustering

                marker.proper.expr.trim.cluster <-
                    kmeans( marker.proper.expr.trim, 2 )

                cluster.neg.idx <- ifelse(
                    marker.proper.expr.trim.cluster$centers[ 1 ] <
                        marker.proper.expr.trim.cluster$centers[ 2 ], 1, 2 )
                cluster.pos.idx <- ifelse( cluster.neg.idx == 1, 2, 1 )

                cluster.neg.x <- marker.proper.expr.trim.cluster$centers[
                    cluster.neg.idx ]
                cluster.pos.x <- marker.proper.expr.trim.cluster$centers[
                    cluster.pos.idx ]

                # identify minimum of distribution between modes

                marker.proper.expr.trim.density <-
                    density( marker.proper.expr.trim )

                marker.proper.expr.trim.density.fun <- splinefun(
                    marker.proper.expr.trim.density$x,
                    marker.proper.expr.trim.density$y )

                marker.proper.expr.trim.density.min <- optimize(
                    marker.proper.expr.trim.density.fun,
                    c( cluster.neg.x, cluster.pos.x ) )$minimum

                # identify extreme values at untrimmed sides of modes

                marker.proper.expr.neg.high.bool <- marker.proper.expr.trim <
                    marker.proper.expr.trim.density.min

                if ( sum( marker.proper.expr.neg.high.bool ) > expr.trim.n )
                    marker.proper.expr.neg.high <- sort(
                        marker.proper.expr.trim[
                            marker.proper.expr.neg.high.bool ],
                        decreasing = TRUE )[ expr.trim.n ]
                else
                    marker.proper.expr.neg.high <- sort(
                        marker.proper.expr.trim )[ 2 ]

                marker.proper.expr.pos.low.bool <- marker.proper.expr.trim >
                    marker.proper.expr.trim.density.min

                if ( sum( marker.proper.expr.pos.low.bool ) > expr.trim.n )
                    marker.proper.expr.pos.low <- sort(
                        marker.proper.expr.trim[
                            marker.proper.expr.pos.low.bool ] )[ expr.trim.n ]
                else
                    marker.proper.expr.pos.low <- sort(
                        marker.proper.expr.trim, decreasing = TRUE )[ 2 ]

                # identify indices of negative and positive events

                expr.neg.trim.idx <- expr.trim.idx[
                    marker.proper.expr.trim < marker.proper.expr.neg.high ]

                expr.pos.trim.idx <- expr.trim.idx[
                    marker.proper.expr.trim > marker.proper.expr.pos.low ]

                # get spillover coefficient

                marker.expr.median.diff <-
                    flow.control$transform.inv[[ marker ]](
                        median( marker.expr[ expr.pos.trim.idx ] ) ) -
                    flow.control$transform.inv[[ marker ]](
                        median( marker.expr[ expr.neg.trim.idx ] ) )

                marker.proper.expr.median.diff <-
                    flow.control$transform.inv[[ marker.proper ]](
                        median( marker.proper.expr[ expr.pos.trim.idx ] ) ) -
                    flow.control$transform.inv[[ marker.proper ]](
                        median( marker.proper.expr[ expr.neg.trim.idx ] ) )

                marker.spillpopu.coef[ marker ] <-
                    marker.expr.median.diff / marker.proper.expr.median.diff
            }

        marker.spillpopu.coef
    },
    mc.cores = get.worker.process( asp$worker.process.n ) )

    marker.spillpopu <- do.call( rbind, marker.spillpopu )
    rownames( marker.spillpopu ) <- flow.control$marker

    marker.spillpopu
}

