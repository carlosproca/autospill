# process_posnegpop.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


#' Process data with positive and negative populations
#'
#' Calculates spillover coefficients using positive and negative populations,
#' and the corresponding compensation error.
#'
#' @param marker.spillover.unco.untr List of two matrices, with regressions
#'     intercepts and coefficients, resulting from the initial spillover
#'     calculation in untransformed scale.
#' @param flow.gate List of vectors with ids of gated events per sample.
#' @param flow.control List with data and metadata of a set of controls.
#' @param asp List with AutoSpill parameters.
#'
#' @return List with three elements:
#'     \itemize{
#'         \item{Calculated spillover matrix.}
#'         \item{Calculated compensation matrix.}
#'         \item{Resulting compensation error, a list with four matrices:
#'             intercepts, coefficients, slopes, and skewness.}
#'     }
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

process.posnegpop <- function( marker.spillover.unco.untr, flow.gate,
    flow.control, asp )
{
    # get spillover and compensation matrices

    spillover.posnegpop <- get.marker.spillover.posnegpop( flow.gate,
        flow.control, asp )

    compensation.posnegpop <- solve( spillover.posnegpop )

    spillover.posnegpop.original <- spillover.posnegpop
    rownames( spillover.posnegpop.original ) <- flow.control$marker.original
    colnames( spillover.posnegpop.original ) <- flow.control$marker.original

    compensation.posnegpop.original <- compensation.posnegpop
    rownames( compensation.posnegpop.original ) <- flow.control$marker.original
    colnames( compensation.posnegpop.original ) <- flow.control$marker.original

    # write spillover and compensation matrices

    if ( ! is.null( asp$table.spillover.dir ) )
        write.csv( spillover.posnegpop.original,
            file = file.path( asp$table.spillover.dir,
                sprintf( "%s.csv", asp$spillover.popnegpop.file.name ) ) )

    if ( ! is.null( asp$table.compensation.dir ) )
        write.csv( compensation.posnegpop.original,
            file = file.path( asp$table.compensation.dir,
                sprintf( "%s.csv", asp$compensation.popnegpop.file.name ) ) )

    # plot spillover and compensation matrices, by rows and columns
    # respectively

    if( ! is.null( asp$figure.spillover.dir ) )
    {
        figure.spillover.file.label <- sprintf( "_%s",
            asp$posnegpop.file.label )

        plot.matrix( spillover.posnegpop, TRUE, asp$figure.spillover.dir,
            figure.spillover.file.label, flow.control, asp )
    }

    if( ! is.null( asp$figure.compensation.dir ) )
    {
        figure.compensation.file.label <- sprintf( "_%s",
            asp$posnegpop.file.label )

        plot.matrix( compensation.posnegpop, FALSE,
            asp$figure.compensation.dir, figure.compensation.file.label,
            flow.control, asp )
    }

    # get uncompensated expresion data

    expr.data.unco <- flow.control$expr.data.untr

    # get compensated expression data

    flow.set.comp <- lapply( flow.control$flow.set, compensate,
        compensation( spillover.posnegpop.original ) )

    expr.data.comp <- get.flow.expression.data( flow.set.comp, flow.control )

    check.critical(
        identical( dimnames( expr.data.comp ), dimnames( expr.data.unco ) ),
        "internal error: inconsistent event or dye names in compensated data"
    )

    # get compensation error in compensated data

    compensation.error <- get.compensation.error(
        expr.data.unco, expr.data.comp, marker.spillover.unco.untr,
        TRUE, TRUE, asp$posnegpop.file.label, flow.gate, flow.control, asp
    )

    # get slope error

    slope.error <- compensation.error$slop - diag( flow.control$marker.n )

    if ( asp$verbose )
    {
        posnegpop.delta <- sd( slope.error )
        posnegpop.delta.max <- max( abs( slope.error ) )

        cat( sprintf( "posnegpop, delta %g, delta.max %g\n",
            posnegpop.delta, posnegpop.delta.max ) )
    }

    # save and plot slope error

    if ( ! is.null( asp$table.slope.error.dir ) )
        write.csv( slope.error,
            file = file.path( asp$table.slope.error.dir,
                sprintf( "%s.csv", asp$slope.error.posnegpop.file.name ) ) )

    if ( ! is.null( asp$figure.slope.error.dir ) )
        plot.density.log( slope.error, "compensation error",
            file.path( asp$figure.slope.error.dir,
                sprintf( "%s.png", asp$slope.error.posnegpop.file.name ) ),
            asp )

    # save and plot skewness

    if ( ! is.null( asp$table.skewness.dir ) )
        write.csv( compensation.error$skew,
            file = file.path( asp$table.skewness.dir,
                sprintf( "%s.csv", asp$skewness.posnegpop.file.name ) ) )

    if ( ! is.null( asp$figure.skewness.dir ) )
    {
        if ( ! is.null( flow.control$autof.marker.idx ) )
            spillover.skewness <- compensation.error$skew[
                - flow.control$autof.marker.idx,
                - flow.control$autof.marker.idx ]
        else
            spillover.skewness <- compensation.error$skew

        plot.density.log( spillover.skewness, "spillover skewness",
            file.path( asp$figure.skewness.dir,
                sprintf( "%s.png", asp$skewness.posnegpop.file.name ) ),
            asp )
    }

    list(
        spillover = spillover.posnegpop,
        compensation = compensation.posnegpop,
        error = compensation.error
    )
}

