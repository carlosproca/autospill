# plot_gate.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


# Plots gate, including intermediate steps.

plot.gate <- function( gate.stage, samp, gate.data, gate.marker, gate.bound,
    gate.region, gate.population, flow.control, asp )
{
    if ( is.null( gate.stage ) )
        gate.stage <- -1

    if ( gate.stage <= 1 )
        gate.data.ggp <- data.frame(
            x = gate.data[ , 1 ],
            y = gate.data[ , 2 ],
            z = interp.surface( gate.bound$density, gate.data ) )
    else if ( gate.stage == 2 )
        gate.data.ggp <- data.frame(
            x = gate.data[ gate.bound$density.max.data.idx, 1 ],
            y = gate.data[ gate.bound$density.max.data.idx, 2 ],
            z = interp.surface( gate.bound$density,
                gate.data[ gate.bound$density.max.data.idx, ] ) )
    else if ( gate.stage == 3 )
        gate.data.ggp <- data.frame(
            x = gate.data[ gate.region$data.idx, 1 ],
            y = gate.data[ gate.region$data.idx, 2 ],
            z = interp.surface( gate.region$density,
                gate.data[ gate.region$data.idx, ] ) )
    else    # gate.stage == 4
        gate.data.ggp <- data.frame(
            x = gate.data[ gate.region$density.max.data.idx, 1 ],
            y = gate.data[ gate.region$density.max.data.idx, 2 ],
            z = interp.surface( gate.region$density,
                gate.data[ gate.region$density.max.data.idx, ] ) )

    density.palette <- get.density.palette( gate.data.ggp$z, asp )

    gate.plot <- ggplot( gate.data.ggp, aes( .data$x, .data$y,
            color = .data$z ) ) +
        scale_x_continuous(
            name = flow.control$scatter.and.marker.label[ gate.marker[ 1 ] ],
            breaks = seq( flow.control$expr.data.min,
                flow.control$expr.data.max, asp$data.step ),
            labels = sprintf( "%dK",
                seq( flow.control$expr.data.min, flow.control$expr.data.max,
                    asp$data.step ) / 1e3 ),
            limits = c( flow.control$expr.data.min,
                flow.control$expr.data.max ),
            expand = expansion( asp$figure.gate.scale.expand ) ) +
        scale_y_continuous(
            name = flow.control$scatter.and.marker.label[ gate.marker[ 2 ] ],
            breaks = seq( flow.control$expr.data.min,
                flow.control$expr.data.max, asp$data.step ),
            labels = sprintf( "%dK",
                seq( flow.control$expr.data.min, flow.control$expr.data.max,
                    asp$data.step ) / 1e3 ),
            limits = c( flow.control$expr.data.min,
                flow.control$expr.data.max ),
            expand = expansion( asp$figure.gate.scale.expand ) ) +
        geom_point( size = 0.9 * asp$figure.gate.point.size,
            stroke = 0.1 * asp$figure.gate.point.size, alpha = 0.4 ) +
        scale_color_gradientn( "", labels = NULL, colors = density.palette,
            guide = guide_colorbar( barwidth = asp$figure.gate.bar.width,
                barheight = asp$figure.gate.bar.height ) ) +
        theme_bw() +
        theme( plot.margin = margin( asp$figure.margin, asp$figure.margin,
                asp$figure.margin, asp$figure.margin ),
            legend.margin = margin( asp$figure.gate.bar.margin,
                asp$figure.gate.bar.margin, asp$figure.gate.bar.margin,
                asp$figure.gate.bar.margin ),
            axis.ticks = element_line( size = asp$figure.panel.line.size ),
            axis.text = element_text( size = asp$figure.axis.text.size ),
            axis.title = element_text( size = asp$figure.axis.title.size ),
            panel.border = element_rect( size = asp$figure.panel.line.size ),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank() )

    if ( gate.stage <= 2 )
    {
        gate.bound.ggp <- data.frame(
            x = c(
                gate.bound$x.low,
                gate.bound$x.high,
                gate.bound$x.high,
                gate.bound$x.low,
                gate.bound$x.low
            ),
            y = c(
                gate.bound$y.low,
                gate.bound$y.low,
                gate.bound$y.high,
                gate.bound$y.high,
                gate.bound$y.low
            )
        )

        gate.plot <- gate.plot +
            geom_path( aes( .data$x, .data$y, color = NULL ),
                data = gate.bound.ggp, size = asp$figure.gate.line.size,
                linetype = "dashed" )
    }

    if ( gate.stage != 1 )
    {
        gate.region.ggp <- data.frame(
            x = c(
                gate.region$x.low,
                gate.region$x.high,
                gate.region$x.high,
                gate.region$x.low,
                gate.region$x.low
            ),
            y = c(
                gate.region$y.low,
                gate.region$y.low,
                gate.region$y.high,
                gate.region$y.high,
                gate.region$y.low )
        )

        gate.plot <- gate.plot +
            geom_path( aes( .data$x, .data$y, color = NULL ),
                data = gate.region.ggp, size = asp$figure.gate.line.size )
    }

    if ( gate.stage <= 0 || gate.stage == 4 )
    {
        gate.boundary.ggp <- data.frame(
            x = c( gate.population$boundary$x,
                gate.population$boundary$x[ 1 ] ),
            y = c( gate.population$boundary$y,
                gate.population$boundary$y[ 1 ] )
        )

        gate.plot <- gate.plot +
            geom_path( aes( .data$x, .data$y, color = NULL ),
                data = gate.boundary.ggp, size = asp$figure.gate.line.size )
    }

    if ( gate.stage <= 0 )
    {
        gate.plot <- gate.plot +
            geom_point( data = gate.bound$density.max,
                size = 1.9 * asp$figure.gate.point.size,
                stroke = 0.1 * asp$figure.gate.point.size,
                color = asp$gate.tesselation.color ) +
            geom_text( data = gate.bound$density.max,
                    aes( label = .data$num.label ),
                hjust = 0, vjust = 0, size = asp$figure.axis.text.size / 2.5,
                color = asp$gate.tesselation.color )
    }
    else if ( gate.stage == 1 )
    {
        gate.plot <- gate.plot +
            geom_point( data = gate.bound$density.max,
                size = 1.9 * asp$figure.gate.point.size,
                stroke = 0.1 * asp$figure.gate.point.size,
                color = asp$gate.tesselation.color ) +
            geom_text( data = gate.bound$density.max,
                    aes( label = .data$num.label ),
                hjust = 0, vjust = 0, size = asp$figure.axis.text.size / 2.5,
                color = asp$gate.tesselation.color )

        if ( gate.bound$density.max.n > 1 )
            gate.plot <- gate.plot +
                geom_segment( data = gate.bound$voronoi$dirsgs,
                    mapping = aes( x = .data$x1, y = .data$y1,
                        xend = .data$x2, yend = .data$y2 ),
                    color = asp$gate.tesselation.color,
                    size = asp$figure.gate.line.size )
    }
    else if ( gate.stage == 2 )
    {
        gate.plot <- gate.plot +
            geom_point( data = gate.bound$density.max[
                    gate.bound$density.max.target, ],
                size = 1.9 * asp$figure.gate.point.size,
                stroke = 0.1 * asp$figure.gate.point.size,
                color = asp$gate.tesselation.color ) +
            geom_text( data = gate.bound$density.max[
                    gate.bound$density.max.target, ],
                aes( label = .data$num.label ),
                hjust = 0, vjust = 0, size = asp$figure.axis.text.size / 2.5,
                color = asp$gate.tesselation.color )

        if ( gate.bound$density.max.n > 1 )
            gate.plot <- gate.plot +
                geom_segment( data = gate.bound$voronoi$dirsgs[
                        gate.bound$voronoi$dirsgs$ind1 ==
                            gate.bound$density.max.target |
                        gate.bound$voronoi$dirsgs$ind2 ==
                            gate.bound$density.max.target,
                    ],
                    mapping = aes( x = .data$x1, y = .data$y1, xend = .data$x2,
                        yend = .data$y2 ),
                    color = asp$gate.tesselation.color,
                    size = asp$figure.gate.line.size )
    }
    else if ( gate.stage == 3 )
    {
        gate.plot <- gate.plot +
            geom_point( data = gate.region$density.max,
                size = 1.9 * asp$figure.gate.point.size,
                stroke = 0.1 * asp$figure.gate.point.size,
                color = asp$gate.tesselation.color ) +
            geom_text( data = gate.region$density.max,
                    aes( label = .data$num.label ),
                hjust = 0, vjust = 0, size = asp$figure.axis.text.size / 2.5,
                color = asp$gate.tesselation.color )

        if ( gate.region$density.max.n > 1 )
            gate.plot <- gate.plot +
                geom_segment( data = gate.region$voronoi$dirsgs,
                    mapping = aes( x = .data$x1, y = .data$y1, xend = .data$x2,
                        yend = .data$y2 ),
                    color = asp$gate.tesselation.color,
                    size = asp$figure.gate.line.size )
    }
    else if ( gate.stage == 4 )
    {
        gate.plot <- gate.plot +
            geom_point( data = gate.region$density.max[ 1, ],
                size = 1.9 * asp$figure.gate.point.size,
                stroke = 0.1 * asp$figure.gate.point.size,
                color = asp$gate.tesselation.color ) +
            geom_text( data = gate.region$density.max[ 1, ],
                    aes( label = .data$num.label ),
                hjust = 0, vjust = 0, size = asp$figure.axis.text.size / 2.5,
                color = asp$gate.tesselation.color )

        if ( gate.region$density.max.n > 1 )
            gate.plot <- gate.plot +
                geom_segment( data = gate.region$voronoi$dirsgs[
                        gate.region$voronoi$dirsgs$ind1 == 1 |
                        gate.region$voronoi$dirsgs$ind2 == 1 , ],
                    mapping = aes( x = .data$x1, y = .data$y1, xend = .data$x2,
                        yend = .data$y2 ),
                    color = asp$gate.tesselation.color,
                    size = asp$figure.gate.line.size )
    }

    if ( gate.stage >= 0 )
        ggsave( file.path( asp$figure.gate.dir,
                sprintf( "%s_%d.png", samp, gate.stage ) ),
            plot = gate.plot, width = asp$figure.width,
            height = asp$figure.height )
    else
        ggsave( file.path( asp$figure.gate.dir, sprintf( "%s.png", samp ) ),
            plot = gate.plot, width = asp$figure.width,
            height = asp$figure.height )

    if ( gate.stage <= 0 && asp$make.thumbnail )
    {
        gate.plot <- ggplot( gate.data.ggp, aes( .data$x, .data$y,
                color = .data$z ) ) +
            scale_x_continuous(
                name = flow.control$scatter.and.marker.label[
                    gate.marker[ 1 ] ],
                breaks = seq( flow.control$expr.data.min,
                    flow.control$expr.data.max, asp$data.step ),
                labels = sprintf( "%dK",
                    seq( flow.control$expr.data.min,
                        flow.control$expr.data.max, asp$data.step ) / 1e3 ),
                limits = c( flow.control$expr.data.min,
                    flow.control$expr.data.max ),
                expand = expansion( asp$thumbnail.gate.scale.expand ) ) +
            scale_y_continuous(
                name = flow.control$scatter.and.marker.label[
                    gate.marker[ 2 ] ],
                breaks = seq( flow.control$expr.data.min,
                    flow.control$expr.data.max, asp$data.step ),
                labels = sprintf( "%dK",
                    seq( flow.control$expr.data.min,
                        flow.control$expr.data.max, asp$data.step ) / 1e3 ),
                limits = c( flow.control$expr.data.min,
                    flow.control$expr.data.max ),
                expand = expansion( asp$thumbnail.gate.scale.expand ) ) +
            geom_point( size = 0.9 * asp$thumbnail.gate.point.size,
                stroke = 0.1 * asp$thumbnail.gate.point.size,
                alpha = 0.4 ) +
            scale_color_gradientn( "", labels = NULL, colors = density.palette,
                guide = guide_colorbar(
                    barwidth = asp$thumbnail.gate.bar.width,
                    barheight = asp$thumbnail.gate.bar.height
                ) ) +
            geom_path( aes( .data$x, .data$y, color = NULL ),
                    data = gate.bound.ggp,
                size = asp$thumbnail.gate.line.size, linetype = "dashed" ) +
            geom_path( aes( .data$x, .data$y, color = NULL ),
                data = gate.region.ggp, size = asp$thumbnail.gate.line.size ) +
            geom_path( aes( .data$x, .data$y, color = NULL ),
                data = gate.boundary.ggp,
                size = asp$thumbnail.gate.line.size ) +
            geom_point( data = gate.bound$density.max,
                size = 1.9 * asp$thumbnail.gate.point.size,
                stroke = 0.1 * asp$thumbnail.gate.point.size,
                color = asp$gate.tesselation.color ) +
            geom_text( data = gate.bound$density.max,
                aes( label = .data$num.label ), hjust = 0, vjust = 0,
                size = asp$thumbnail.axis.text.size / 2.5,
                color = asp$gate.tesselation.color ) +
            theme_bw() +
            theme( plot.margin = margin( asp$thumbnail.margin,
                asp$thumbnail.margin, asp$thumbnail.margin,
                asp$thumbnail.margin ),
                legend.margin = margin( asp$thumbnail.gate.bar.margin,
                    asp$thumbnail.gate.bar.margin,
                    asp$thumbnail.gate.bar.margin,
                    asp$thumbnail.gate.bar.margin ),
                axis.ticks = element_line(
                    size = asp$thumbnail.panel.line.size ),
                axis.text = element_text(
                    size = asp$thumbnail.axis.text.size ),
                axis.title = element_text(
                    size = asp$thumbnail.axis.title.size ),
                panel.border = element_rect(
                    size = asp$thumbnail.panel.line.size ),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank() )

        ggsave( file.path( asp$figure.gate.dir,
                sprintf( "%s%s_thumbnail.png",
                    ifelse( gate.stage == 0, "_0", "" ), samp ) ),
            plot = gate.plot, width = asp$thumbnail.width,
            height = asp$thumbnail.height )
    }
}

