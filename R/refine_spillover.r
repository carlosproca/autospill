# refine_spillover.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


#' Refine spillover coefficients
#'
#' Refines spillover coefficients iteratively.
#'
#' @param marker.spillover.unco.untr List of two matrices, with regressions
#'     intercepts and coefficients, resulting from the initial spillover
#'     calculation in untransformed scale.
#' @param marker.spillover.unco.tran List of two matrices, with regressions
#'     intercepts and coefficients, resulting from the initial spillover
#'     calculation in transformed scale. Optional parameter used only in
#'     scatter plots, it can be \code{NULL}.
#' @param flow.gate List of vectors with ids of gated events per sample.
#' @param flow.control List with data and metadata of a set of controls.
#' @param asp List with AutoSpill parameters.
#'
#' @return List with four elements:
#'     \itemize{
#'         \item{Spillover matrix at final step.}
#'         \item{Compensation matrix at final step.}
#'         \item{Compensation error at final step, a list with four matrices:
#'             intercepts, coefficients, slopes, and skewness.}
#'         \item{Dataframe with convergence data.}
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

refine.spillover <- function( marker.spillover.unco.untr,
    marker.spillover.unco.tran, flow.gate, flow.control, asp )
{
    # set initial values for iteration variables

    rs.convergence <- FALSE
    rs.exit <- FALSE

    rs.iter <- 0
    rs.iter.last <- FALSE
    rs.iter.width <- floor( log10( asp$rs.iter.max ) ) + 1

    rs.lambda <- asp$rs.lambda.coarse

    rs.delta <- -1.0
    rs.delta.threshold <- asp$rs.delta.threshold.untr

    rs.delta.history <- rep( -1, asp$rs.delta.history.n )

    rs.scale.untransformed <- TRUE

    rs.convergence.log <- data.frame(
        iter = numeric(),
        scale = character(),
        lambda = numeric(),
        delta = numeric(),
        delta.max = numeric(),
        delta.change = numeric(),
        stringsAsFactors = FALSE
    )

    # set initial values for spillover calculation

    spillover.curr <- diag( flow.control$marker.n )
    spillover.update <- marker.spillover.unco.untr$coef -
        diag( flow.control$marker.n )

    while ( ! rs.exit )
    {
        # update spillover matrix and calculate compensation matrix

        spillover.curr <- spillover.curr + spillover.update
        spillover.curr <- sweep( spillover.curr, 1,
            diag( spillover.curr ), "/" )

        compensation.curr <- solve( spillover.curr )

        spillover.curr.original <- spillover.curr
        rownames( spillover.curr.original ) <- flow.control$marker.original
        colnames( spillover.curr.original ) <- flow.control$marker.original

        compensation.curr.original <- compensation.curr
        rownames( compensation.curr.original ) <- flow.control$marker.original
        colnames( compensation.curr.original ) <- flow.control$marker.original

        if ( ( rs.iter == 0 && asp$rs.save.table.initial ) ||
            ( asp$rs.save.table.every > 0 &&
                    rs.iter %% asp$rs.save.table.every == 0 ) ||
            rs.iter.last )
        {
            # save spillover and compensation matrices

            if( ! is.null( asp$table.spillover.dir ) )
            {
                table.spillover.file.name <- ifelse( rs.iter.last,
                    sprintf( "%s.csv", asp$spillover.file.name ),
                    sprintf( "%s_%0*d.csv", asp$spillover.file.name,
                        rs.iter.width, rs.iter ) )

                write.csv( spillover.curr.original,
                    file = file.path( asp$table.spillover.dir,
                        table.spillover.file.name ) )
            }

            if( ! is.null( asp$table.compensation.dir ) )
            {
                table.compensation.file.name <- ifelse( rs.iter.last,
                    sprintf( "%s.csv", asp$compensation.file.name ),
                    sprintf( "%s_%0*d.csv", asp$compensation.file.name,
                        rs.iter.width, rs.iter ) )

                write.csv( compensation.curr.original,
                    file = file.path( asp$table.compensation.dir,
                        table.compensation.file.name ) )
            }

            # plot spillover and compensation matrices, by rows and columns
            # respectively

            if( ! is.null( asp$figure.spillover.dir ) )
            {
                figure.spillover.file.label <- ifelse( rs.iter.last, "",
                    sprintf( "_%0*d", rs.iter.width, rs.iter ) )

                plot.matrix( spillover.curr, TRUE, asp$figure.spillover.dir,
                    figure.spillover.file.label, flow.control, asp )
            }

            if( ! is.null( asp$figure.compensation.dir ) )
            {
                figure.compensation.file.label <- ifelse( rs.iter.last, "",
                    sprintf( "_%0*d", rs.iter.width, rs.iter ) )

                plot.matrix( compensation.curr, FALSE,
                    asp$figure.compensation.dir,
                    figure.compensation.file.label, flow.control, asp )
            }
        }

        # set uncompensated expresion data and spillover

        if ( rs.scale.untransformed ) {
            expr.data.unco <- flow.control$expr.data.untr
            marker.spillover.unco <- marker.spillover.unco.untr
        }
        else {
            expr.data.unco <- flow.control$expr.data.tran
            marker.spillover.unco <- marker.spillover.unco.tran
        }

        # get compensated expression data

        flow.set.comp <- lapply( flow.control$flow.set, compensate,
            compensation( spillover.curr.original ) )

        if ( ! rs.scale.untransformed )
            flow.set.comp <- lapply( flow.set.comp, transform,
                transformList( names( flow.control$transform ),
                    flow.control$transform ) )

        expr.data.comp <- get.flow.expression.data( flow.set.comp,
            flow.control )

        check.critical(
            identical( dimnames( expr.data.comp ),
                dimnames( expr.data.unco ) ),
            "internal error: inconsistent event or dye names in compensated data"
        )

        # get compensation error in compensated data

        if ( ( rs.iter == 0 && asp$rs.plot.figure.initial ) ||
            ( asp$rs.plot.figure.every > 0 &&
                    rs.iter %% asp$rs.plot.figure.every == 0 ) ||
            rs.iter.last )
        {
            plot.scatter.figure <- TRUE
            figure.scatter.file.label <- ifelse( rs.iter.last, "final",
                sprintf( "%0*d", rs.iter.width, rs.iter ) )
        }
        else
        {
            plot.scatter.figure <- FALSE
            figure.scatter.file.label <- NULL
        }

        compensation.error <- get.compensation.error(
            expr.data.unco, expr.data.comp, marker.spillover.unco,
            rs.scale.untransformed, plot.scatter.figure,
            figure.scatter.file.label, flow.gate, flow.control, asp
        )

        # get slope error and update delta variables

        slope.error <- compensation.error$slop - diag( flow.control$marker.n )

        rs.delta.prev <- rs.delta
        rs.delta <- sd( slope.error )

        rs.delta.max <- max( abs( slope.error ) )

        if ( rs.delta.prev >= 0 )
            rs.delta.history[ rs.iter %% asp$rs.delta.history.n + 1 ] <-
            rs.delta - rs.delta.prev
        else
            rs.delta.history[ rs.iter %% asp$rs.delta.history.n + 1 ] <- -1

        rs.delta.change <- mean( rs.delta.history )

        rs.convergence.log[ rs.iter + 1, ] <- list( rs.iter,
            ifelse( rs.scale.untransformed, "linear", "bi-exp" ),
            rs.lambda, rs.delta, rs.delta.max, rs.delta.change )

        if ( asp$verbose )
        {
            cat( sprintf(
                "iter %0*d, %s scale, lambda %.1f, delta %g, delta.max %g, delta.change %g\n",
                rs.iter.width, rs.iter,
                ifelse( rs.scale.untransformed, "linear", "bi-exp" ),
                rs.lambda, rs.delta, rs.delta.max, rs.delta.change ) )
        }

        if ( ( rs.iter == 0 && asp$rs.save.table.initial ) ||
                ( asp$rs.save.table.every > 0 &&
                        rs.iter %% asp$rs.save.table.every == 0 ) ||
                rs.iter.last )
        {
            # save and plot slope error

            if( ! is.null( asp$table.slope.error.dir ) )
            {
                table.slope.error.file.name <- ifelse( rs.iter.last,
                    sprintf( "%s.csv", asp$slope.error.file.name ),
                    sprintf( "%s_%0*d.csv", asp$slope.error.file.name,
                        rs.iter.width, rs.iter ) )

                write.csv( slope.error,
                    file = file.path( asp$table.slope.error.dir,
                        table.slope.error.file.name ) )
            }

            if( ! is.null( asp$figure.slope.error.dir ) )
            {
                figure.slope.error.file.name <- sprintf( "%s%s.png",
                    asp$slope.error.file.name,
                    ifelse( rs.iter.last, "",
                        sprintf( "_%0*d", rs.iter.width, rs.iter ) ) )

                plot.density.log( slope.error, "compensation error",
                    file.path( asp$figure.slope.error.dir,
                        figure.slope.error.file.name ),
                    asp )
            }

            # save and plot skewness

            if( ! is.null( asp$table.skewness.dir ) )
            {
                table.skewness.file.name <- ifelse( rs.iter.last,
                    sprintf( "%s.csv", asp$skewness.file.name ),
                    sprintf( "%s_%0*d.csv", asp$skewness.file.name,
                        rs.iter.width, rs.iter ) )

                write.csv( compensation.error$skew,
                    file = file.path( asp$table.skewness.dir,
                        table.skewness.file.name ) )
            }

            if( ! is.null( asp$figure.skewness.dir ) )
            {
                if ( ! is.null( flow.control$autof.marker.idx ) )
                    spillover.skewness <- compensation.error$skew[
                        - flow.control$autof.marker.idx,
                        - flow.control$autof.marker.idx ]
                else
                    spillover.skewness <- compensation.error$skew

                figure.skewness.file.name <- sprintf( "%s%s.png",
                    asp$skewness.file.name,
                    ifelse( rs.iter.last, "",
                        sprintf( "_%0*d", rs.iter.width, rs.iter ) ) )

                plot.density.log( spillover.skewness, "spillover skewness",
                    file.path( asp$figure.skewness.dir,
                        figure.skewness.file.name ),
                    asp )
            }
        }

        # update iteration variables

        if( rs.scale.untransformed && rs.delta.max < rs.delta.threshold )
        {
            # switch to bi-exponential scale and reset lambda and delta history
            rs.scale.untransformed <- FALSE
            rs.delta.threshold <- asp$rs.delta.threshold.tran
            rs.lambda <- asp$rs.lambda.coarse
            rs.delta <- -1.0
            rs.delta.history <- rep( -1, asp$rs.delta.history.n )
            rs.delta.change <- -1
        }

        if ( rs.delta.change > - asp$rs.delta.threshold.change &&
                rs.lambda == asp$rs.lambda.coarse )
        {
            # reduce lambda and reset delta history
            rs.lambda <- asp$rs.lambda.fine
            rs.delta <- -1.0
            rs.delta.history <- rep( -1, asp$rs.delta.history.n )
            rs.delta.change <- -1
        }

        rs.convergence <- ! rs.scale.untransformed &&
            ( rs.delta.max < rs.delta.threshold ||
                    rs.delta.change > - asp$rs.delta.threshold.change )

        rs.exit <- ( rs.convergence && rs.iter.last ) ||
            ( rs.delta.change > - asp$rs.delta.threshold.change &&
                    rs.scale.untransformed ) ||
            ( ! rs.convergence && rs.iter == asp$rs.iter.max ) ||
            rs.iter > asp$rs.iter.max

        rs.iter.last <- rs.convergence

        rs.iter <- rs.iter + 1

        # update spillover matrix

        spillover.update <- rs.lambda * ( slope.error %*% spillover.curr )
    }

    # save and plot convergence

    if ( ! is.null( asp$table.convergence.dir ) )
        write.csv( rs.convergence.log,
            file = file.path( asp$table.convergence.dir,
                sprintf( "%s.csv", asp$convergence.file.name ) ),
            row.names = FALSE )

    if ( ! is.null( asp$figure.convergence.dir ) )
        plot_convergence( rs.convergence.log, NULL, asp )

    # check convergence

    check.critical( rs.convergence,
        "no convergence in refinement of spillover matrix" )

    list(
        spillover = spillover.curr,
        compensation = compensation.curr,
        error = compensation.error,
        convergence = rs.convergence.log
    )
}

