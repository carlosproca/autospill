# plot_scatter.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


# Makes scatter plots of compensated and uncompensated data, with regression
# and reference lines, and a label showing slope error.

plot.scatter <- function(
    expr.data.unco.x, expr.data.unco.y, expr.data.comp.x, expr.data.comp.y,
    spillover.unco.inte, spillover.unco.coef,
    spillover.comp.inte, spillover.comp.coef, spillover.comp.slop,
    marker.limit.range, marker.proper.limit.range, samp, marker, marker.proper,
    scale.untransformed, figure.file.label, flow.gate, flow.control, asp
)
{
    expr.data.ggp <- data.frame(
        x = c( expr.data.unco.x, expr.data.comp.x ),
        y = c( expr.data.unco.y, expr.data.comp.y ),
        z = factor( rep( c( "unco", "comp" ),
                    each = length( expr.data.unco.x ) ),
                levels = c( "unco", "comp" ) ),
        w = rep( 1 : length( expr.data.unco.x ) %in%
                flow.gate[[ samp ]], 2 )
    )

    marker.slope.unco <- 1 / spillover.unco.coef
    marker.intercept.unco <- - spillover.unco.inte / spillover.unco.coef
    marker.x.intercept.unco <- spillover.unco.inte

    marker.slope.comp <- 1 / spillover.comp.coef
    marker.intercept.comp <- - spillover.comp.inte / spillover.comp.coef
    marker.x.intercept.comp <- spillover.comp.inte

    x.transform <- flow.control$transform[[
        flow.control$marker.original[ match( marker, flow.control$marker ) ] ]]
    y.transform <- flow.control$transform[[
        flow.control$marker.original[ match( marker.proper,
            flow.control$marker ) ] ]]

    x.transform.inv <- flow.control$transform.inv[[
        flow.control$marker.original[ match( marker, flow.control$marker ) ] ]]
    y.transform.inv <- flow.control$transform.inv[[
        flow.control$marker.original[ match( marker.proper,
            flow.control$marker ) ] ]]

    marker.limit.adaptor <- diag( 2 ) +
        matrix( c( 0.01 * c( 1, -1 ), 0.05 * c( -1, 1 ) ),
            nrow = 2, byrow = TRUE )
    marker.limit.plot <- as.vector( marker.limit.adaptor %*%
            marker.limit.range )

    marker.proper.limit.adaptor <- diag( 2 ) +
        matrix( c( 0.01 * c( 1, -1 ), 0.05 * c( -1, 1 ) ),
            nrow = 2, byrow = TRUE )
    marker.proper.limit.plot <- as.vector(
        marker.proper.limit.adaptor %*%
            marker.proper.limit.range )

    if ( ! scale.untransformed )
    {
        marker.limit.plot <- x.transform.inv( marker.limit.plot )
        marker.proper.limit.plot <- y.transform.inv(
            marker.proper.limit.plot )
    }

    marker.limit.factor <- log(
        max( abs( marker.proper.limit.plot ) ) /
            max( abs( marker.limit.plot ) ) )

    if ( marker.limit.factor > 1 )
        marker.limit.plot <- as.vector(
            matrix( c(
                ( marker.limit.factor + 1 ) / 2,
                ( 1 - marker.limit.factor ) / 2,
                ( 1 - marker.limit.factor ) / 2,
                ( marker.limit.factor + 1 ) / 2 ),
                nrow = 2, byrow = TRUE ) %*%
                marker.limit.plot )

    if ( ! scale.untransformed )
    {
        marker.limit.plot <- x.transform( marker.limit.plot )
        marker.proper.limit.plot <- y.transform(
            marker.proper.limit.plot )
    }

    for ( scale.same in c( TRUE, FALSE ) )
        if ( scale.same || asp$scatter.plot.scale.other )
        {
            if ( xor( scale.same, ! scale.untransformed ) )
            {
                if ( scale.untransformed ) {
                    gg.scale.x <- scale_x_continuous(
                        limits = marker.limit.plot )
                    gg.scale.y <- scale_y_continuous(
                        limits = marker.proper.limit.plot )
                }
                else {
                    gg.scale.x <- scale_x_continuous(
                        limits = x.transform.inv( marker.limit.plot ) )
                    gg.scale.y <- scale_y_continuous(
                        limits = y.transform.inv( marker.proper.limit.plot ) )
                }
            }
            else
            {
                if ( scale.untransformed ) {
                    marker.limit.plot.untr <- marker.limit.plot
                    marker.limit.plot.tran <- x.transform(
                        marker.limit.plot )
                }
                else {
                    marker.limit.plot.untr <- x.transform.inv(
                        marker.limit.plot )
                    marker.limit.plot.tran <- marker.limit.plot
                }

                marker.breaks.exp <- min( floor( log10( max( abs(
                    marker.limit.plot.untr ) ) /
                        asp$scatter.scale.breaks.coef ) ), 3 )
                marker.breaks.untr <-
                    get.scale.breaks.untransformed(
                        marker.limit.plot.untr,
                        marker.breaks.exp )
                marker.breaks.tran <- x.transform(
                    marker.breaks.untr )

                if ( scale.untransformed ) {
                    marker.proper.limit.plot.untr <-
                        marker.proper.limit.plot
                    marker.proper.limit.plot.tran <- y.transform(
                        marker.proper.limit.plot )
                }
                else {
                    marker.proper.limit.plot.untr <-
                        y.transform.inv( marker.proper.limit.plot )
                    marker.proper.limit.plot.tran <-
                        marker.proper.limit.plot
                }

                marker.proper.breaks.exp <- min( floor( log10( max(
                    abs( marker.proper.limit.plot.untr ) ) /
                        asp$scatter.scale.breaks.coef ) ), 3 )
                marker.proper.breaks.untr <-
                    get.scale.breaks.untransformed(
                        marker.proper.limit.plot.untr,
                        marker.proper.breaks.exp )
                marker.proper.breaks.tran <- y.transform(
                    marker.proper.breaks.untr )

                gg.scale.x <- scale_x_continuous(
                    breaks = marker.breaks.tran,
                    labels = marker.breaks.untr,
                    limits = marker.limit.plot.tran )

                gg.scale.y <- scale_y_continuous(
                    breaks = marker.proper.breaks.tran,
                    labels = marker.proper.breaks.untr,
                    limits = marker.proper.limit.plot.tran )
            }

            if ( marker == marker.proper )
            {
                if ( asp$scatter.ref.line.unco )
                    gg.ref.line.unco <- geom_abline( slope = 1,
                        intercept = 0,
                        color = asp$scatter.ref.line.color,
                        linetype = "dashed",
                        size = asp$scatter.ref.line.size.factor *
                            asp$figure.scatter.line.size )
                else
                    gg.ref.line.unco <- geom_blank()

                gg.ref.line.comp <- geom_abline( slope = 1,
                    intercept = 0,
                    color = asp$scatter.ref.line.color,
                    linetype = "dashed",
                    size = asp$scatter.ref.line.size.factor *
                        asp$figure.scatter.line.size )
            }
            else
            {
                if ( asp$scatter.ref.line.unco )
                {
                    if ( scale.same )
                        marker.x.intercept.unco.plot <- marker.x.intercept.unco
                    else
                        marker.x.intercept.unco.plot <- ifelse(
                            scale.untransformed,
                            x.transform( marker.x.intercept.unco ),
                            x.transform.inv( marker.x.intercept.unco )
                        )

                    gg.ref.line.unco <- geom_vline(
                        xintercept = marker.x.intercept.unco.plot,
                        color = asp$scatter.ref.line.color,
                        linetype = "dashed",
                        size = asp$scatter.ref.line.size.factor *
                            asp$figure.scatter.line.size )
                }
                else
                    gg.ref.line.unco <- geom_blank()

                if ( scale.same )
                    marker.x.intercept.comp.plot <- marker.x.intercept.comp
                else
                    marker.x.intercept.comp.plot <- ifelse(
                        scale.untransformed,
                        x.transform( marker.x.intercept.comp ),
                        x.transform.inv( marker.x.intercept.comp )
                    )

                gg.ref.line.comp <- geom_vline(
                    xintercept = marker.x.intercept.comp.plot,
                    color = asp$scatter.ref.line.color,
                    linetype = "dashed",
                    size = asp$scatter.ref.line.size.factor *
                        asp$figure.scatter.line.size )
            }

            if ( scale.same )
            {
                if ( is.infinite( marker.slope.unco ) )
                    gg.slope.line.unco <- geom_vline(
                        xintercept = marker.x.intercept.unco,
                        color = asp$scatter.expr.color.unco,
                        size = asp$figure.scatter.line.size
                    )
                else
                    gg.slope.line.unco <- geom_abline(
                        slope = marker.slope.unco,
                        intercept = marker.intercept.unco,
                        color = asp$scatter.expr.color.unco,
                        size = asp$figure.scatter.line.size
                    )

                if ( is.infinite( marker.slope.comp ) )
                    gg.slope.line.comp <- geom_vline(
                        xintercept = marker.x.intercept.comp,
                        color = asp$scatter.expr.color.comp,
                        size = asp$figure.scatter.line.size
                    )
                else
                    gg.slope.line.comp <- geom_abline(
                        slope = marker.slope.comp,
                        intercept = marker.intercept.comp,
                        color = asp$scatter.expr.color.comp,
                        size = asp$figure.scatter.line.size )
            }
            else
            {
                gg.slope.line.unco <- geom_blank()
                gg.slope.line.comp <- geom_blank()
            }

            gg.error.label <- annotate(
                "label",
                label = sprintf( "%.4g", spillover.comp.slop ),
                size = asp$figure.scatter.error.label.size,
                x = sum( c( 1 - asp$figure.scatter.error.label.pos.x,
                            asp$figure.scatter.error.label.pos.x ) *
                        gg.scale.x$limits ),
                y = sum( c( 1 - asp$figure.scatter.error.label.pos.y,
                            asp$figure.scatter.error.label.pos.y ) *
                        gg.scale.y$limits ),
                label.size = NA,
                alpha = 0.8
            )

            gg.label <- labs(
                x = flow.control$scatter.and.marker.label[ marker ],
                y = flow.control$scatter.and.marker.label[ marker.proper ]
            )

            gg.scale.color <- scale_color_manual(
                values = c( asp$scatter.expr.color.unco,
                    asp$scatter.expr.color.comp ),
                guide = "none"
            )

            data.ggp <- expr.data.ggp

            if ( ! scale.same )
            {
                if ( scale.untransformed )
                {
                    data.ggp$x <- x.transform( data.ggp$x )
                    data.ggp$y <- y.transform( data.ggp$y )
                }
                else
                {
                    data.ggp$x <- x.transform.inv( data.ggp$x )
                    data.ggp$y <- y.transform.inv( data.ggp$y )
                }
            }

            the.plot <- ggplot( data.ggp, aes( .data$x, .data$y,
                    color = .data$z ) ) +
                geom_point(
                    size = 0.9 * asp$figure.scatter.point.size,
                    stroke = 0.1 * asp$figure.scatter.point.size,
                    alpha = ifelse( data.ggp$w,
                        asp$figure.scatter.alpha.gate.in,
                        asp$figure.scatter.alpha.gate.out ) ) +
                gg.slope.line.unco +
                gg.ref.line.unco +
                gg.slope.line.comp +
                gg.ref.line.comp +
                gg.error.label +
                gg.scale.x +
                gg.scale.y +
                gg.label +
                gg.scale.color +
                theme_bw() +
                theme( plot.margin = margin( asp$figure.margin,
                    asp$figure.margin, asp$figure.margin,
                    asp$figure.margin ),
                    axis.ticks = element_line(
                        size = asp$figure.panel.line.size ),
                    axis.text = element_text(
                        size = asp$figure.scatter.axis.text.size ),
                    axis.title = element_text(
                        size = asp$figure.scatter.axis.title.size ),
                    panel.border = element_rect(
                        size = asp$figure.panel.line.size ),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank() )

            figure.file.name <- sprintf(
                "%s_%s_%s.png", marker, figure.file.label,
                ifelse( xor( scale.same, ! scale.untransformed ),
                    "linear", "bi-exp" )
            )

            ggsave(
                file.path( flow.control$figure.scatter.dir[ samp ],
                    figure.file.name ),
                plot = the.plot, width = asp$figure.width,
                height = asp$figure.height
            )

            if ( asp$make.thumbnail )
            {
                if ( scale.same )
                {
                    if ( marker == marker.proper )
                    {
                        if ( asp$scatter.ref.line.unco )
                            gg.ref.line.unco$aes_params$size <-
                                asp$scatter.ref.line.size.factor *
                                    asp$thumbnail.scatter.line.size

                        gg.ref.line.comp$aes_params$size <-
                            asp$scatter.ref.line.size.factor *
                                asp$thumbnail.scatter.line.size
                    }
                    else
                    {
                        if ( asp$scatter.ref.line.unco )
                            gg.ref.line.unco$aes_params$size <-
                                asp$scatter.ref.line.size.factor *
                                    asp$thumbnail.scatter.line.size

                        gg.ref.line.comp$aes_params$size <-
                            asp$scatter.ref.line.size.factor *
                                asp$thumbnail.scatter.line.size
                    }
                }

                if ( scale.same )
                {
                    if ( is.infinite( marker.slope.unco ) )
                        gg.slope.line.unco$aes_params$size <-
                            asp$thumbnail.scatter.line.size
                    else
                        gg.slope.line.unco$aes_params$size <-
                            asp$thumbnail.scatter.line.size

                    if ( is.infinite( marker.slope.comp ) )
                        gg.slope.line.comp$aes_params$size <-
                            asp$thumbnail.scatter.line.size
                    else
                        gg.slope.line.comp$aes_params$size <-
                            asp$thumbnail.scatter.line.size
                }

                gg.error.label$aes_params$size <-
                    asp$thumbnail.scatter.error.label.size
                gg.error.label$data$x <- sum(
                    c( 1 - asp$thumbnail.scatter.error.label.pos.x,
                            asp$thumbnail.scatter.error.label.pos.x ) *
                    gg.scale.x$limits
                )
                gg.error.label$data$y <- sum(
                    c( 1 - asp$thumbnail.scatter.error.label.pos.y,
                        asp$thumbnail.scatter.error.label.pos.y ) *
                    gg.scale.y$limits
                )

                the.plot <- ggplot( data.ggp, aes( .data$x, .data$y,
                        color = .data$z ) ) +
                    geom_point(
                        size = 0.9 * asp$thumbnail.scatter.point.size,
                        stroke = 0.1 * asp$thumbnail.scatter.point.size,
                        alpha = ifelse( data.ggp$w,
                            asp$thumbnail.scatter.alpha.gate.in,
                            asp$thumbnail.scatter.alpha.gate.out ) ) +
                    gg.slope.line.unco +
                    gg.ref.line.unco +
                    gg.slope.line.comp +
                    gg.ref.line.comp +
                    gg.error.label +
                    gg.scale.x +
                    gg.scale.y +
                    gg.label +
                    gg.scale.color +
                    theme_bw() +
                    theme( plot.margin = margin( asp$thumbnail.margin,
                        asp$thumbnail.margin, asp$thumbnail.margin,
                        asp$thumbnail.margin ),
                        axis.ticks = element_line(
                            size = asp$thumbnail.panel.line.size ),
                        axis.text = element_text(
                            size = asp$thumbnail.scatter.axis.text.size ),
                        axis.title = element_text(
                            size = asp$thumbnail.scatter.axis.title.size ),
                        panel.border = element_rect(
                            size = asp$thumbnail.panel.line.size ),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank() )

                figure.file.name <- sprintf(
                    "%s_%s_%s_thumbnail.png", marker, figure.file.label,
                    ifelse( xor( scale.same, ! scale.untransformed ),
                        "linear", "bi-exp" )
                )

                ggsave(
                    file.path( flow.control$figure.scatter.dir[ samp ],
                        figure.file.name ),
                    plot = the.plot, width = asp$thumbnail.width,
                    height = asp$thumbnail.height
                )
            } # make.thumbnail
        } # scale.same TRUE or FALSE
}

