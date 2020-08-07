# get_scale_breaks_untransformed.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


# Returns a vector of scale breaks well distributed as powers of 10, inside
# the given range limit.plot.untr and avoiding exponents smaller than mix.exp.

get.scale.breaks.untransformed <- function( limit.plot.untr, min.exp )
{
    if ( limit.plot.untr[ 1 ] > 10^min.exp )
        breaks.exp.min <- ceiling( log10( limit.plot.untr[ 1 ] ) )
    else if ( limit.plot.untr[ 1 ] > 0 )
        breaks.exp.min <- min.exp
    else if ( limit.plot.untr[ 1 ] > -10^min.exp )
        breaks.exp.min <- 0
    else
        breaks.exp.min <- - floor( log10( - limit.plot.untr[ 1 ] ) )

    if ( limit.plot.untr[ 2 ] < -10^min.exp )
        breaks.exp.max <- - ceiling( log10( - limit.plot.untr[ 2 ] ) )
    else if ( limit.plot.untr[ 2 ] < 0 )
        breaks.exp.max <- - min.exp
    else if ( limit.plot.untr[ 2 ] < 10^min.exp )
        breaks.exp.max <- 0
    else
        breaks.exp.max <- floor( log10( limit.plot.untr[ 2 ] ) )

    breaks.exp <- breaks.exp.min : breaks.exp.max

    breaks.exp <- breaks.exp[
        abs( breaks.exp ) == 0 | abs( breaks.exp ) >= min.exp ]

    sapply( breaks.exp, function ( be )
        if ( be > 0 )
            10^be
        else if ( be < 0 )
            -10^-be
        else
            be
    )
}

