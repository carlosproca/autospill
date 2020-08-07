# plot_density_log.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


# Plots density of log values segreggated by sign, ignoring diagonal
# coefficients and zeroes.

plot.density.log <- function( x.data, x.label, plot.file.path, asp )
{
    # get values avoiding matrix diagonal
    x.data.n <- nrow( x.data )
    x.dat <- x.data[ -( 1 + 0:(x.data.n-1) * (x.data.n+1) ) ]

    # remove other zero values
    x.dat <- x.dat[ x.dat != 0 ]

    x.ggdata <- data.frame( log.abs.x = log10( abs( x.dat ) ),
        sign = factor( sign( x.dat ), levels = c( -1, 0, 1 ) ) )

    scale.color.value <- c( "-1" = asp$density.color.single,
        "1" = asp$density.color.single )

    scale.linetype.value <- c( "-1" = "dashed", "1" = "solid" )

    x.ggplot <- ggplot( x.ggdata, aes( x = .data$log.abs.x, color = .data$sign,
            fill = .data$sign, linetype = .data$sign ) ) +
        geom_density( alpha = 0.1, size = asp$figure.density.line.size ) +
        scale_color_manual( values = scale.color.value ) +
        scale_fill_manual( values = scale.color.value ) +
        scale_linetype_manual( values = scale.linetype.value ) +
        labs( x = sprintf( "log10( Absolute %s )", x.label ), y = "Density" ) +
        theme_bw() +
        theme( plot.margin = margin( asp$figure.margin,
                asp$figure.margin, asp$figure.margin,
                asp$figure.margin ),
            axis.ticks = element_line( size = asp$figure.panel.line.size ),
            axis.text = element_text( size = asp$figure.axis.text.size ),
            axis.title = element_text( size = asp$figure.axis.title.size ),
            panel.border = element_rect( size = asp$figure.panel.line.size ),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none" )

    ggsave( plot.file.path, plot = x.ggplot,
        width = asp$figure.width, height = asp$figure.height )
}


# Plots density of log values, together for initial and final steps of
# autospill, and calculation with positive and negative populations.

plot.density.log.together <- function( x.table, x.label, plot.file.path, asp )
{

    # get values avoiding matrix diagonal
    x.data <- lapply( x.table, function( xt ) {
        xt.n <- nrow( xt )
        xt[ -( 1 + 0:(xt.n-1) * (xt.n+1) ) ]
    } )
    x.data.n <- length( x.data )

    # remove other zero values
    x.data <- lapply( x.data, function( xd )
        xd[ xd != 0 ]
    )

    x.ggdata <- lapply( 1 : x.data.n, function( xd.idx ) {
        data.frame( log.abs.x = log10( abs( x.data[[ xd.idx ]] ) ),
            sign = sign( x.data[[ xd.idx ]] ) * xd.idx )
    } )
    x.ggdata <- do.call( rbind, x.ggdata )
    x.ggdata$sign <- factor( x.ggdata$sign,
        levels = c( 0, rep( 1:x.data.n, each = 2 ) * c( -1, 1 ) ) )

    scale.color.value <- c(
        "-1" = asp$density.color.initial, "1" = asp$density.color.initial,
        "-2" = asp$density.color.final, "2" = asp$density.color.final,
        "-3" = asp$density.color.posnegpop, "3" = asp$density.color.posnegpop
    )

    scale.linetype.value <- c(
        "-1" = "dashed", "1" = "solid",
        "-2" = "dashed", "2" = "solid",
        "-3" = "dashed" , "3"= "solid"
    )

    x.ggplot <- ggplot( x.ggdata, aes( x = .data$log.abs.x, color = .data$sign,
            fill = .data$sign, linetype = .data$sign ) ) +
        geom_density( alpha = 0.1, size = asp$figure.density.line.size ) +
        scale_color_manual( values = scale.color.value ) +
        scale_fill_manual( values = scale.color.value  ) +
        scale_linetype_manual( values = scale.linetype.value ) +
        labs( x = sprintf( "log10( Absolute %s )", x.label ), y = "Density" ) +
        theme_bw() +
        theme( plot.margin = margin( asp$figure.margin,
            asp$figure.margin, asp$figure.margin,
            asp$figure.margin ),
            axis.ticks = element_line( size = asp$figure.panel.line.size ),
            axis.text = element_text( size = asp$figure.axis.text.size ),
            axis.title = element_text( size = asp$figure.axis.title.size ),
            panel.border = element_rect( size = asp$figure.panel.line.size ),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none" )

    ggsave( plot.file.path, plot = x.ggplot,
        width = asp$figure.width, height = asp$figure.height )
}
