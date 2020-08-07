# get_density_palette.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


# Returns a good color palette for a given density distribution.

get.density.palette <- function( dens, asp )
{
    rainbow.palette <- colorRampPalette( asp$density.palette.base.color )(
        asp$density.palette.base.n )

    dens.range <- range( dens, na.rm = TRUE )

    dens.grid <- seq( dens.range[ 1 ], dens.range[ 2 ],
        length.out = asp$density.palette.n )

    density.palette.idx <-
        round( ecdf( dens )( dens.grid ) * asp$density.palette.base.n )

    rainbow.palette[ density.palette.idx ]
}

