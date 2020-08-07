# read_gate_parameter.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


# Returns a list with gate parameters per marker.

read.gate.parameter <- function( flow.control, asp )
{
    if ( ! is.null( asp$gate.parameter.file.name ) &&
            file.exists( asp$gate.parameter.file.name ) )
    {
        gate.param.table <- read.csv( asp$gate.parameter.file.name,
            stringsAsFactors = FALSE )

        check.critical(
            sort( gate.param.table$dye ) ==
                sort( flow.control$marker.original ),
            "wrong dye name in gate parameters" )

        dye.var.idx <- which( names( gate.param.table ) == "dye" )

        gate.param <- lapply( flow.control$marker.original, function( mo ) {
            mo.idx <- which( gate.param.table$dye == mo )
            as.list( gate.param.table[ mo.idx, - dye.var.idx ] )
        } )
    }
    else
        gate.param <- lapply( flow.control$marker.original, function( mo )
            asp$default.gate.param )

    names( gate.param ) <- flow.control$marker.original

    gate.param
}

