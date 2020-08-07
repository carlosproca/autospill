# get_transformation.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


# Returns two lists with transformation parameters per marker, for direct and
# inverse tranformation.

get.transformation <- function( flow.control, asp )
{
    if ( ! is.null( asp$transformation.parameter.file.name ) &&
            file.exists( asp$transformation.parameter.file.name ) )
    {
        transformation.param <- read.csv(
            asp$transformation.parameter.file.name, stringsAsFactors = FALSE )

        check.critical(
            sort( transformation.param$dye ) ==
                sort( flow.control$marker.original ),
            "wrong dye name in tranformation parameters"
        )

        transf <- lapply( flow.control$marker.original, function( mo ) {
            mo.idx <- which( transformation.param$dye == mo )
            flowjo_biexp(
                channelRange = transformation.param$length[ mo.idx ],
                maxValue = transformation.param$max.range[ mo.idx ],
                pos = transformation.param$pos[ mo.idx ],
                neg = transformation.param$neg[ mo.idx ],
                widthBasis = transformation.param$width[ mo.idx ],
                inverse = FALSE
        ) } )

        transf.inv <- lapply( flow.control$marker.original, function( mo ) {
            mo.idx <- which( transformation.param$dye == mo )
            flowjo_biexp(
                channelRange = transformation.param$length[ mo.idx ],
                maxValue = transformation.param$max.range[ mo.idx ],
                pos = transformation.param$pos[ mo.idx ],
                neg = transformation.param$neg[ mo.idx ],
                widthBasis = transformation.param$width[ mo.idx ],
                inverse = TRUE
        ) } )
    }
    else
    {
        transf <- lapply( flow.control$marker.original, function( mo )
            flowjo_biexp(
                channelRange = asp$default.transformation.param$length,
                maxValue = asp$default.transformation.param$max.range,
                pos = asp$default.transformation.param$pos,
                neg = asp$default.transformation.param$neg,
                widthBasis = asp$default.transformation.param$width,
                inverse = FALSE
        ) )

        transf.inv <- lapply( flow.control$marker.original, function( mo )
            flowjo_biexp(
                channelRange = asp$default.transformation.param$length,
                maxValue = asp$default.transformation.param$max.range,
                pos = asp$default.transformation.param$pos,
                neg = asp$default.transformation.param$neg,
                widthBasis = asp$default.transformation.param$width,
                inverse = TRUE
        ) )
    }

    names( transf ) <- flow.control$marker.original
    names( transf.inv ) <- flow.control$marker.original

    list( transf, transf.inv )
}

