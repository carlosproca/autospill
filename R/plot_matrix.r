# plot_matrix.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


# Plots matrix coefficients, by rows or columns.

plot.matrix <- function( the.matrix, by.rows, figure.dir, figure.file.label,
    flow.control, asp )
{
    for ( marker.proper in flow.control$marker )
    {
        if ( is.numeric( flow.control$wavelength ) &&
                all( ! is.na( flow.control$wavelength ) &
                        flow.control$wavelength >= asp$matrix.wavelength.min &
                        flow.control$wavelength <= asp$matrix.wavelength.max ) )
        {
            plot.x.value <- flow.control$wavelength

            plot.x.axis <- scale_x_continuous(
                breaks = plot.x.value,
                labels = flow.control$scatter.and.marker.label[
                    flow.control$marker ],
                sec.axis = sec_axis( ~., breaks = seq(
                    asp$matrix.wavelength.min + 100,
                    asp$matrix.wavelength.max - 100, by = 50 ) )
            )
        }
        else
        {
            plot.x.value <- 1 : flow.control$marker.n

            plot.x.axis <- scale_x_continuous(
                breaks = plot.x.value,
                labels = flow.control$scatter.and.marker.label[
                    flow.control$marker ]
            )
        }

        if ( by.rows )
            plot.y.value <- the.matrix[ marker.proper, ]
        else
            plot.y.value <- the.matrix[ , marker.proper ]

        data.ggp <- data.frame(
            x = plot.x.value,
            y = plot.y.value,
            z = ifelse( flow.control$marker == marker.proper,
                asp$matrix.marker.proper.color, asp$matrix.marker.color )
        )

        matrix.plot <- ggplot( data.ggp, aes( .data$x, .data$y ) ) +
            geom_point( size = 0.9 * asp$figure.matrix.point.size,
                stroke = 0.1 * asp$figure.matrix.point.size,
                color = data.ggp$z ) +
            geom_hline( yintercept = 0, linetype = "dashed",
                size = asp$figure.matrix.line.size ) +
            labs( x = "", y = "" ) +
            plot.x.axis +
            theme_bw() +
            theme( plot.margin = margin( asp$figure.margin,
                asp$figure.margin, asp$figure.margin, asp$figure.margin ),
                legend.position = "none",
                axis.line = element_blank(),
                axis.ticks = element_line(
                    size = asp$figure.panel.line.size ),
                axis.text.x = element_text(
                    size = asp$figure.axis.text.size, angle = 90,
                    hjust = 1, vjust = 0.5 ),
                axis.text.x.top = element_text(
                    size = asp$figure.axis.text.size, hjust = 1,
                    vjust = 0.5 ),
                axis.text.y = element_text(
                    size = asp$figure.axis.text.size ),
                axis.title = element_text( size =
                        asp$figure.axis.title.size ),
                panel.border = element_rect(
                    size = asp$figure.panel.line.size, fill = NA ),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank() )

        figure.file.path <- file.path( figure.dir,
            sprintf( "%s%s.png" , marker.proper, figure.file.label ) )

        ggsave( figure.file.path, plot = matrix.plot, width = asp$figure.width,
            height = asp$figure.height )

        if ( asp$make.thumbnail )
        {
            plot.x.axis$labels <- flow.control$marker.original

            thumbnail.plot <- ggplot( data.ggp, aes( .data$x, .data$y ) ) +
                geom_point( size = 0.9 * asp$thumbnail.matrix.point.size,
                    stroke = 0.1 * asp$thumbnail.matrix.point.size,
                    color = data.ggp$z ) +
                geom_hline( yintercept = 0, linetype = "dashed",
                    size = asp$thumbnail.matrix.line.size ) +
                labs( x = "", y = "" ) +
                plot.x.axis +
                theme_bw() +
                theme( plot.margin = margin( asp$thumbnail.margin,
                    asp$thumbnail.margin, asp$thumbnail.margin,
                    asp$thumbnail.margin ),
                    legend.position = "none",
                    axis.line = element_blank(),
                    axis.ticks = element_line(
                        size = asp$thumbnail.panel.line.size ),
                    axis.text.x = element_text(
                        size = asp$thumbnail.axis.text.size,
                        angle = 90, hjust = 1, vjust = 0.5 ),
                    axis.text.x.top = element_text(
                        size = asp$thumbnail.axis.text.size, hjust = 1,
                        vjust = 0.5 ),
                    axis.text.y = element_text(
                        size = asp$thumbnail.axis.text.size ),
                    panel.border = element_rect(
                        size = asp$thumbnail.panel.line.size,
                        fill = NA ),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank() )

            thumbnail.file.path <- file.path( figure.dir,
                sprintf( "%s%s_thumbnail.png" , marker.proper,
                    figure.file.label ) )

            ggsave( thumbnail.file.path, plot = thumbnail.plot,
                width = asp$thumbnail.width, height = asp$thumbnail.height )
        }
    } # marker.proper
}

