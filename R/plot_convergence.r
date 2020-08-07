# plot_convergence.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


#' Plot AutoSpill convergence
#'
#' Plots convergence of iterative refinement of the spillover matrix.
#'
#' @param convergence.log Dataframe with convergence data of AutoSpill.
#' @param popnegpop.slope Matrix with slope errors resulting from calculating
#'     spillover with positive and negative populations. Optional parameter, it
#'     can be \code{NULL}.
#' @param asp List with AutoSpill parameters.
#'
#' @return \code{NULL}.
#'
#' @references Roca \emph{et al}:
#'     AutoSpill: A method for calculating spillover coefficients to compensate
#'     or unmix high-parameter flow cytometry data.
#'     \emph{bioRxiv} 2020.06.29.177196;
#'     \href{https://doi.org/10.1101/2020.06.29.177196}{doi:10.1101/2020.06.29.177196}
#'     (2020).
#'
#' @seealso \code{\link{refine.spillover}}, \code{\link{process.posnegpop}},
#'     and \code{\link{get.autospill.param}}.
#'
#' @export

plot_convergence <- function( convergence.log, popnegpop.slope, asp )
{
    convergence.ggdata <- convergence.log
    convergence.ggdata$delta.change <- abs( convergence.ggdata$delta.change )

    if ( ! is.null( popnegpop.slope ) )
    {
        slope.error <- popnegpop.slope - diag( nrow( popnegpop.slope ) )

        posnegpop.delta <- sd( slope.error )
        posnegpop.delta.max <- max( abs( slope.error ) )
    }
    else
    {
        posnegpop.delta <- NULL
        posnegpop.delta.max <- NULL
    }

    plot.xaxis.max  <- ceiling( max( convergence.log$iter ) / 10 ) * 10
    if ( plot.xaxis.max == 10 )
        plot.xaxis.step <- 5
    else if ( plot.xaxis.max <= 50 )
        plot.xaxis.step <- 10
    else
        plot.xaxis.step <- 25

    convergence.ggplot <- ggplot( convergence.ggdata, aes( x = .data$iter,
            shape = .data$scale ) ) +
        labs( x = "Iteration", y = "Absolute compensation error" ) +
        geom_hline( yintercept = asp$rs.delta.threshold.untr,
            size = asp$figure.convergence.line.size,
            linetype = "dashed" ) +
        geom_hline( yintercept = asp$rs.delta.threshold.tran,
            size = asp$figure.convergence.line.size,
            linetype = "dashed" ) +
        geom_hline( yintercept = asp$rs.delta.threshold.change,
            size = asp$figure.convergence.line.size,
            linetype = "dashed" ) +
        geom_point( aes( y = .data$delta ),
            size = asp$figure.convergence.point.size,
            color = asp$convergence.color.delta ) +
        geom_point( aes( y = .data$delta.max ),
            size = asp$figure.convergence.point.size,
            color = asp$convergence.color.delta.max ) +
        geom_point( aes( y = .data$delta.change ),
            size = asp$figure.convergence.point.size,
            color = asp$convergence.color.delta.change ) +
        scale_x_continuous( limits = c( 0, plot.xaxis.max ),
            breaks = seq( 0, plot.xaxis.max, plot.xaxis.step ) ) +
        scale_y_log10() +
        scale_shape_manual( values = c(
            "linear" = asp$convergence.shape.linear,
            "bi-exp" = asp$convergence.shape.biexp,
            "posnegpop" = asp$convergence.shape.posnegpop
        ) ) +
        theme_bw() +
        theme( legend.position = "none",
            plot.margin = margin( asp$figure.margin, asp$figure.margin,
                asp$figure.margin, asp$figure.margin ),
            axis.ticks = element_line( size = asp$figure.panel.line.size ),
            axis.text = element_text( size = asp$figure.axis.text.size ),
            axis.title = element_text( size = asp$figure.axis.title.size ),
            panel.border = element_rect( size = asp$figure.panel.line.size ),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank() )

    if ( ! is.null( posnegpop.delta ) && ! is.null( posnegpop.delta.max ) )
        convergence.ggplot <- convergence.ggplot +
            geom_point( aes( x = 0, y = posnegpop.delta, shape = "posnegpop" ),
                size = asp$figure.convergence.point.size,
                color = asp$convergence.color.delta ) +
            geom_point( aes( x = 0, y = posnegpop.delta.max,
                    shape = "posnegpop" ),
                size = asp$figure.convergence.point.size,
                color = asp$convergence.color.delta.max )

    ggsave( file.path( asp$figure.convergence.dir,
        sprintf( "%s.png", asp$convergence.file.name ) ),
        plot = convergence.ggplot,
        width = asp$figure.width, height = asp$figure.height )
}

