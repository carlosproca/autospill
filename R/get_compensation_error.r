# get_compensation_error.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


# Returns a list of matrices discribing compensation error, with intercepts,
# coefficients, slopes, and skewness.

get.compensation.error <- function( expr.data.unco, expr.data.comp,
    marker.spillover.unco, scale.untransformed, plot.figure, figure.label,
    flow.gate, flow.control, asp )
{
    marker.spillover.zero <- rep( 0, flow.control$marker.n )
    names( marker.spillover.zero ) <- flow.control$marker

    marker.spillover.comp <- mclapply( flow.control$sample, function( samp )
    {
        marker.proper <- samp

        # get expression data for primary marker

        marker.proper.expr.unco <- expr.data.unco[
            which( flow.control$event.sample == samp )[ flow.gate[[ samp ]] ],
            marker.proper ]

        marker.proper.expr.comp <- expr.data.comp[
            which( flow.control$event.sample == samp )[ flow.gate[[ samp ]] ],
            marker.proper ]

        marker.expr.n <- length( marker.proper.expr.unco )

        # identify extreme expression values

        expr.trim.n <- round( marker.expr.n * asp$rlm.trim.factor )

        marker.proper.expr.low.unco <- sort( marker.proper.expr.unco )[
            expr.trim.n ]
        marker.proper.expr.high.unco <- sort( marker.proper.expr.unco,
            decreasing = TRUE )[ expr.trim.n ]

        marker.proper.expr.low.comp <- sort( marker.proper.expr.comp )[
            expr.trim.n ]
        marker.proper.expr.high.comp <- sort( marker.proper.expr.comp,
            decreasing = TRUE )[ expr.trim.n ]

        marker.spillover.comp.inte <- marker.spillover.zero
        marker.spillover.comp.coef <- marker.spillover.zero
        marker.spillover.comp.slop <- marker.spillover.zero
        marker.spillover.comp.skew <- marker.spillover.zero

        for ( marker in flow.control$marker )
        {
            # get expression data for secondary marker

            marker.expr.unco <- expr.data.unco[
                which( flow.control$event.sample == samp )[
                    flow.gate[[ samp ]] ],
                marker ]

            marker.expr.comp <- expr.data.comp[
                which( flow.control$event.sample == samp )[
                    flow.gate[[ samp ]] ],
                marker ]

            marker.expr.low.unco <- sort( marker.expr.unco )[ expr.trim.n ]
            marker.expr.high.unco <- sort( marker.expr.unco,
                decreasing = TRUE )[ expr.trim.n ]

            marker.expr.low.comp <- sort( marker.expr.comp )[ expr.trim.n ]
            marker.expr.high.comp <- sort( marker.expr.comp,
                decreasing = TRUE )[ expr.trim.n ]

            # trim expression values for both markers

            expr.trim.idx.unco <- which (
                marker.proper.expr.unco > marker.proper.expr.low.unco &
                    marker.proper.expr.unco < marker.proper.expr.high.unco &
                    marker.expr.unco > marker.expr.low.unco &
                    marker.expr.unco < marker.expr.high.unco
            )

            expr.trim.idx.comp <- which (
                marker.proper.expr.comp > marker.proper.expr.low.comp &
                    marker.proper.expr.comp < marker.proper.expr.high.comp &
                    marker.expr.comp > marker.expr.low.comp &
                    marker.expr.comp < marker.expr.high.comp
            )

            marker.proper.expr.trim.unco <- marker.proper.expr.unco[
                expr.trim.idx.unco ]
            marker.expr.trim.unco <- marker.expr.unco[ expr.trim.idx.unco ]

            marker.proper.expr.trim.comp <- marker.proper.expr.comp[
                expr.trim.idx.comp ]
            marker.expr.trim.comp <- marker.expr.comp[ expr.trim.idx.comp ]

            # fit robust linear model with compensated data

            if ( marker == marker.proper )
            {
                marker.spillover.comp.coef[ marker ] <- 1.0
                marker.spillover.comp.slop[ marker ] <- 1.0
            }
            else
            {
                spillover.model.result <- fit.robust.linear.model(
                    marker.proper.expr.trim.comp, marker.expr.trim.comp,
                    marker.proper, marker, asp )

                marker.spillover.comp.inte[ marker ] <-
                    spillover.model.result[ 1, 1 ]

                marker.spillover.comp.coef[ marker ] <-
                    spillover.model.result[ 2, 1 ]

                # get slope and skewness in untransformed scale

                if ( scale.untransformed )
                {
                    marker.spillover.comp.slop[ marker ] <-
                        marker.spillover.comp.coef[ marker ]
                    marker.spillover.comp.skew[ marker ] <-
                        skewness( marker.expr.trim.comp )
                }
                else
                {
                    x.transform.inv <- flow.control$transform.inv[[
                        flow.control$marker.original[
                            match( marker, flow.control$marker ) ] ]]
                    y.transform.inv <- flow.control$transform.inv[[
                        flow.control$marker.original[
                            match( marker.proper, flow.control$marker ) ] ]]

                    y1p <- min( marker.proper.expr.trim.comp )
                    y2p <- max( marker.proper.expr.trim.comp )

                    x1p <- marker.spillover.comp.inte[ marker ] +
                        marker.spillover.comp.coef[ marker ] * y1p
                    x2p <- marker.spillover.comp.inte[ marker ] +
                        marker.spillover.comp.coef[ marker ] * y2p

                    if ( y1p == y2p || x1p == x2p )
                        marker.spillover.comp.slop[ marker ] <- 0
                    else
                    {
                        y1 <- y.transform.inv( y1p )
                        y2 <- y.transform.inv( y2p )

                        x1 <- x.transform.inv( x1p )
                        x2 <- x.transform.inv( x2p )

                        marker.spillover.comp.slop[ marker ] <-
                            marker.spillover.comp.coef[ marker ] *
                            ( x2 - x1 ) * ( y2p - y1p ) /
                            ( ( x2p - x1p ) * ( y2 - y1 ) )
                    }

                    marker.spillover.comp.skew[ marker ] <-
                        skewness( x.transform.inv( marker.expr.trim.comp ) )
                }
            }

            if ( plot.figure && ! is.null( flow.control$figure.scatter.dir ) )
                plot.scatter(
                    expr.data.unco[
                        flow.control$event.sample == samp, marker ],
                    expr.data.unco[
                        flow.control$event.sample == samp, marker.proper ],
                    expr.data.comp[
                        flow.control$event.sample == samp, marker ],
                    expr.data.comp[
                        flow.control$event.sample == samp,marker.proper ],
                    marker.spillover.unco$inte[ marker.proper, marker ],
                    marker.spillover.unco$coef[ marker.proper, marker ],
                    marker.spillover.comp.inte[ marker ],
                    marker.spillover.comp.coef[ marker ],
                    marker.spillover.comp.slop[ marker ],
                    range( c( marker.expr.trim.unco, marker.expr.trim.comp ) ),
                    range( c( marker.proper.expr.trim.unco,
                        marker.proper.expr.trim.comp ) ),
                    samp, marker, marker.proper,
                    scale.untransformed, figure.label,
                    flow.gate, flow.control, asp
                )
        } # marker

        c( marker.spillover.comp.inte, marker.spillover.comp.coef,
            marker.spillover.comp.slop, marker.spillover.comp.skew )
    },
    mc.cores = get.worker.process( asp$worker.process.n ) ) # samp

    marker.spillover.comp <- do.call( rbind, marker.spillover.comp )
    rownames( marker.spillover.comp ) <- flow.control$marker

    list(
        inte = marker.spillover.comp[ , 1 : flow.control$marker.n ],
        coef = marker.spillover.comp[ , 1 : flow.control$marker.n +
                flow.control$marker.n ],
        slop = marker.spillover.comp[ , 1 : flow.control$marker.n +
                2 * flow.control$marker.n ],
        skew = marker.spillover.comp[ , 1 : flow.control$marker.n +
                3 * flow.control$marker.n ]
    )
}

