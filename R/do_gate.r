# do_gate.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


# Returns a vector with the indexes of events inside the initial gate on
# scatter parameters.
#
# Proceeds in three steps:
#    - Defines bound by data trimming
#    - Defines region around the target maximum found in bound
#    - Defines gate around target maximum, only in region
#
# Uses numerical search of maxima over estimated densities.
# Uses Voronoi tesselations to improve density estimation around maxima.

do.gate <- function( gate.data, gate.param, samp, flow.control, asp )
{
    gate.marker <- colnames( gate.data )

    gate.bound <- NULL
    gate.region <- NULL
    gate.boundary <- NULL

    # trim data

    gate.data.x.min <- max( flow.control$expr.data.min,
        min( gate.data[ , 1 ] ) )
    gate.data.x.max <- min( flow.control$expr.data.max,
        max( gate.data[ , 1 ] ) )

    gate.data.y.min <- max( flow.control$expr.data.min,
        min( gate.data[ , 2 ] ) )
    gate.data.y.max <- min( flow.control$expr.data.max,
        max( gate.data[ , 2 ] ) )

    gate.data.trim.x.min <-
        ( 1 - asp$gate.data.trim.factor.x.min ) * gate.data.x.min +
        asp$gate.data.trim.factor.x.min * gate.data.x.max
    gate.data.trim.x.max <-
        ( 1 - asp$gate.data.trim.factor.x.max ) * gate.data.x.min +
        asp$gate.data.trim.factor.x.max * gate.data.x.max

    gate.data.trim.y.min <-
        ( 1 - asp$gate.data.trim.factor.y.min ) * gate.data.y.min +
        asp$gate.data.trim.factor.y.min * gate.data.y.max
    gate.data.trim.y.max <-
        ( 1 - asp$gate.data.trim.factor.y.max ) * gate.data.y.min +
        asp$gate.data.trim.factor.y.max * gate.data.y.max

    # set bound from trimmed data

    gate.bound.x.low <- gate.data.trim.x.min
    gate.bound.x.high <- gate.data.trim.x.max

    gate.bound.y.low <- gate.data.trim.y.min
    gate.bound.y.high <- gate.data.trim.y.max

    gate.bound.data.idx <- which(
        gate.data[ , 1 ] > gate.bound.x.low &
        gate.data[ , 1 ] < gate.bound.x.high &
        gate.data[ , 2 ] > gate.bound.y.low &
        gate.data[ , 2 ] < gate.bound.y.high )

    gate.bound.density <- kde2d(
        gate.data[ gate.bound.data.idx, 1 ],
        gate.data[ gate.bound.data.idx, 2 ],
        asp$gate.bound.density.bw.factor *
            apply( gate.data[ gate.bound.data.idx, ], 2, bandwidth.nrd ),
        n = asp$gate.bound.density.grid.n )

    # get density maxima in bound

    gate.bound.neighbor.idx <- list(
        x = - asp$gate.bound.density.neigh.size :
            asp$gate.bound.density.neigh.size,
        y = - asp$gate.bound.density.neigh.size :
            asp$gate.bound.density.neigh.size )

    gate.bound.density.max.bool <- matrix( FALSE,
        nrow = asp$gate.bound.density.grid.n,
        ncol = asp$gate.bound.density.grid.n )

    for ( x.idx in 1 : asp$gate.bound.density.grid.n )
        for ( y.idx in 1 : asp$gate.bound.density.grid.n )
            gate.bound.density.max.bool[ x.idx, y.idx ] <-
                gate.bound.density$z[ x.idx, y.idx ] >=
                max( gate.bound.density$z[
                    pmax( 0, pmin( asp$gate.bound.density.grid.n,
                        x.idx + gate.bound.neighbor.idx$x ) ),
                    pmax( 0, pmin( asp$gate.bound.density.grid.n,
                        y.idx + gate.bound.neighbor.idx$y ) ) ] )

    gate.bound.density.max.idx <- which( gate.bound.density.max.bool,
        arr.ind = TRUE )

    gate.bound.density.max.n <- nrow( gate.bound.density.max.idx )

    check.critical( gate.bound.density.max.n >= 1,
        paste0( "gate error: no population found in sample bound", samp ) )

    gate.bound.density.max <- data.frame(
        x = gate.bound.density$x[ gate.bound.density.max.idx[ , 1 ] ],
        y = gate.bound.density$y[ gate.bound.density.max.idx[ , 2 ] ],
        z = gate.bound.density$z[ gate.bound.density.max.idx ] )

    gate.bound.density.max <- gate.bound.density.max[
        order( gate.bound.density.max$z, decreasing = TRUE ), ]

    row.names( gate.bound.density.max ) <- NULL
    gate.bound.density.max$num.label <- paste0( " ",
        row.names( gate.bound.density.max ) )

    # locate target maximum in bound, avoiding maxima located near the
    # bottom left corner

    gate.bound.density.max.offset <- 1

    while ( gate.bound.density.max.offset <= gate.bound.density.max.n )
    {
        if ( ( gate.bound.density.max$x[ gate.bound.density.max.offset ] -
                    gate.bound.x.low ) /
                    ( gate.bound.x.high - gate.bound.x.low ) >
                asp$gate.bound.density.max.exclusion.corner ||
            ( gate.bound.density.max$y[ gate.bound.density.max.offset ] -
                    gate.bound.y.low ) /
                    ( gate.bound.y.high - gate.bound.y.low ) >
                asp$gate.bound.density.max.exclusion.corner )
            break

        gate.bound.density.max.offset <- gate.bound.density.max.offset + 1
    }

    check.critical(
        gate.bound.density.max.offset <= gate.bound.density.max.n,
        paste0( "gate error: no good maximum found in sample bound",
            samp ) )

    gate.bound.density.max.target <- asp$gate.bound.density.max.target +
        gate.bound.density.max.offset - 1

    check.critical(
        gate.bound.density.max.target <= gate.bound.density.max.n,
        paste0( "gate error: target maximum not found in sample bound",
            samp ) )

    if ( gate.bound.density.max.n > 1 )
    {
        # get voronoi tesselation for density maxima

        gate.bound.voronoi <- deldir( gate.bound.density.max,
            rw = c( gate.bound.x.low, gate.bound.x.high, gate.bound.y.low,
                gate.bound.y.high ), suppressMsge = TRUE )

        gate.bound.tile <- tile.list( gate.bound.voronoi )

        # get data in the tile of target maximum

        gate.bound.density.max.data.idx <- gate.bound.data.idx[
            sapply( gate.bound.data.idx, function( gbdi )
                which.tile( gate.data[ gbdi, 1 ], gate.data[ gbdi, 2 ],
                    gate.bound.tile ) == gate.bound.density.max.target )
        ]
    }
    else
    {
        gate.bound.voronoi <- NULL
        gate.bound.density.max.data.idx <- gate.bound.data.idx
    }

    if ( gate.param$region.auto )
    {
        # set region from target maximum found in bound

        gate.bound.density.max.x.median <- median( gate.data[
            gate.bound.density.max.data.idx, 1 ] )
        gate.bound.density.max.x.mad <- mad(
            gate.data[ gate.bound.density.max.data.idx, 1 ],
            center = gate.bound.density.max.x.median )

        gate.bound.density.max.y.median <- median( gate.data[
            gate.bound.density.max.data.idx, 2 ] )
        gate.bound.density.max.y.mad <- mad(
            gate.data[ gate.bound.density.max.data.idx, 2 ],
            center = gate.bound.density.max.y.median )

        gate.region.x.low <- max( gate.data.trim.x.min,
            gate.bound.density.max.x.median -
                asp$gate.bound.density.max.mad.factor *
                    gate.bound.density.max.x.mad )
        gate.region.x.high <- min( gate.data.trim.x.max,
            gate.bound.density.max.x.median +
                asp$gate.bound.density.max.mad.factor *
                    gate.bound.density.max.x.mad )

        gate.region.y.low <- max( gate.data.trim.y.min,
            gate.bound.density.max.y.median -
                asp$gate.bound.density.max.mad.factor *
                    gate.bound.density.max.y.mad )
        gate.region.y.high <- min( gate.data.trim.y.max,
            gate.bound.density.max.y.median +
                asp$gate.bound.density.max.mad.factor *
                    gate.bound.density.max.y.mad )
    }
    else
    {
        gate.region.x.low <-
            ( 1 - gate.param$region.factor.x.low ) * gate.bound.x.low +
            gate.param$region.factor.x.low * gate.bound.x.high
        gate.region.x.high <-
            ( 1 - gate.param$region.factor.x.high ) * gate.bound.x.low +
            gate.param$region.factor.x.high * gate.bound.x.high

        gate.region.y.low <-
            ( 1 - gate.param$region.factor.y.low ) * gate.bound.y.low +
            gate.param$region.factor.y.low * gate.bound.y.high
        gate.region.y.high <-
            ( 1 - gate.param$region.factor.y.high ) * gate.bound.y.low +
            gate.param$region.factor.y.high * gate.bound.y.high
    }

    gate.bound <- list(
        density = gate.bound.density,
        density.max = gate.bound.density.max,
        density.max.n = gate.bound.density.max.n,
        density.max.data.idx = gate.bound.density.max.data.idx,
        density.max.target = gate.bound.density.max.target,
        voronoi = gate.bound.voronoi,
        x.low = gate.bound.x.low,
        x.high = gate.bound.x.high,
        y.low = gate.bound.y.low,
        y.high = gate.bound.y.high
    )

    if ( ! is.null( asp$figure.gate.dir ) && asp$gate.plot.stage )
        plot.gate( 1, samp, gate.data, gate.marker, gate.bound,
            gate.region, gate.population, flow.control, asp )

    gate.region.data.idx <- which(
        gate.data[ , 1 ] > gate.region.x.low &
        gate.data[ , 1 ] < gate.region.x.high &
        gate.data[ , 2 ] > gate.region.y.low &
        gate.data[ , 2 ] < gate.region.y.high )

    # get density maxima in region

    gate.region.density <- kde2d(
        gate.data[ gate.region.data.idx, 1 ],
        gate.data[ gate.region.data.idx, 2 ],
        asp$gate.region.density.bw.factor *
            apply( gate.data[ gate.region.data.idx, ], 2, bandwidth.nrd ),
        n = asp$gate.region.density.grid.n )

    gate.region.neighbor.idx <- list(
        x = - asp$gate.region.density.neigh.size :
            asp$gate.region.density.neigh.size,
        y = - asp$gate.region.density.neigh.size :
            asp$gate.region.density.neigh.size )

    gate.region.density.max.bool <- matrix( FALSE,
        nrow = asp$gate.region.density.grid.n,
        ncol = asp$gate.region.density.grid.n )

    for ( x.idx in 1 : asp$gate.region.density.grid.n )
        for ( y.idx in 1 : asp$gate.region.density.grid.n )
            gate.region.density.max.bool[ x.idx, y.idx ] <-
                gate.region.density$z[ x.idx, y.idx ] >=
                max( gate.region.density$z[
                    pmax( 0, pmin( asp$gate.region.density.grid.n,
                        x.idx + gate.region.neighbor.idx$x ) ),
                    pmax( 0, pmin( asp$gate.region.density.grid.n,
                        y.idx + gate.region.neighbor.idx$y ) ) ] )

    gate.region.density.max.idx <- which( gate.region.density.max.bool,
        arr.ind = TRUE )

    gate.region.density.max.n <- nrow( gate.region.density.max.idx )

    check.critical( gate.region.density.max.n >= 1,
        paste0( "gate error: no population found in sample region", samp ) )

    gate.region.density.max <- data.frame(
        x = gate.region.density$x[ gate.region.density.max.idx[ , 1 ] ],
        y = gate.region.density$y[ gate.region.density.max.idx[ , 2 ] ],
        z = gate.region.density$z[ gate.region.density.max.idx ] )

    gate.region.density.max <- gate.region.density.max[
        order( gate.region.density.max$z, decreasing = TRUE ), ]

    row.names( gate.region.density.max ) <- NULL
    gate.region.density.max$num.label <- paste0( " ",
        row.names( gate.region.density.max ) )

    if ( gate.region.density.max.n > 1 )
    {
        # get voronoi tesselation for density maxima

        gate.region.voronoi <- deldir( gate.region.density.max,
            rw = c( gate.region.x.low, gate.region.x.high, gate.region.y.low,
                gate.region.y.high ), suppressMsge = TRUE )

        gate.region.tile <- tile.list( gate.region.voronoi )

        # get data in the tile of largest maximum

        gate.region.density.max.data.idx <- gate.region.data.idx[
            sapply( gate.region.data.idx, function( grdi )
                which.tile( gate.data[ grdi, 1 ], gate.data[ grdi, 2 ],
                    gate.region.tile ) == 1 )
        ]
    }
    else
    {
        gate.region.voronoi <- NULL
        gate.region.density.max.data.idx <- gate.region.data.idx
    }

    gate.region <- list(
        data.idx = gate.region.data.idx,
        density = gate.region.density,
        density.max = gate.region.density.max,
        density.max.n = gate.region.density.max.n,
        density.max.data.idx = gate.region.density.max.data.idx,
        voronoi = gate.region.voronoi,
        x.low = gate.region.x.low,
        x.high = gate.region.x.high,
        y.low = gate.region.y.low,
        y.high = gate.region.y.high
    )

    if ( ! is.null( asp$figure.gate.dir ) && asp$gate.plot.stage ) {
        plot.gate( 2, samp, gate.data, gate.marker, gate.bound, gate.region,
            gate.population, flow.control, asp )
        plot.gate( 3, samp, gate.data, gate.marker, gate.bound, gate.region,
            gate.population, flow.control, asp )
    }

    # thresold data in region around target maximum

    gate.region.max.density <- kde2d(
        gate.data[ gate.region.density.max.data.idx, 1 ],
        gate.data[ gate.region.density.max.data.idx, 2 ],
        asp$gate.region.max.density.bw.factor *
            apply( gate.data[ gate.region.density.max.data.idx, ], 2,
                bandwidth.nrd ),
        n = asp$gate.region.max.density.grid.n )

    gate.region.max.density.interp <- interp.surface( gate.region.max.density,
        gate.data[ gate.region.density.max.data.idx, ] )

    gate.region.max.density.threshold <-
        ( 1 - gate.param$density.threshold ) *
            min( gate.region.max.density.interp ) +
        gate.param$density.threshold * max( gate.region.max.density.interp )

    gate.population.strict.idx <- gate.region.density.max.data.idx[
        gate.region.max.density.interp > gate.region.max.density.threshold ]

    gate.population.strict.idx <- gate.population.strict.idx[
        ! duplicated( data.frame( gate.data[ gate.population.strict.idx, ] ) ) ]

    gate.population.boundary <- convex.hull( tri.mesh(
        gate.data[ gate.population.strict.idx, 1 ],
        gate.data[ gate.population.strict.idx, 2 ] ) )

    gate.population.pip <- point.in.polygon(
        gate.data[ , 1 ], gate.data[ , 2 ],
        gate.population.boundary$x, gate.population.boundary$y )

    gate.population.idx <- which( gate.population.pip != 0 )

    gate.population <- list( boundary = gate.population.boundary )

    if ( ! is.null( asp$figure.gate.dir ) && asp$gate.plot.stage ) {
        plot.gate( 4, samp, gate.data, gate.marker, gate.bound, gate.region,
            gate.population, flow.control, asp )
        plot.gate( 0, samp, gate.data, gate.marker, gate.bound, gate.region,
            gate.population, flow.control, asp )
    }
    else if ( ! is.null( asp$figure.gate.dir ) )
        plot.gate( NULL, samp, gate.data, gate.marker, gate.bound, gate.region,
            gate.population, flow.control, asp )

    gate.population.idx
}

