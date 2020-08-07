# read_scatter_parameter.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


# Returns a vector with two scatter parameters.

read.scatter.parameter <- function( asp )
{
    if ( ! is.null( asp$scatter.parameter.file.name ) &&
            file.exists( asp$scatter.parameter.file.name ) )
        scatter.parameter <- read.csv( asp$scatter.parameter.file.name,
            stringsAsFactors = FALSE )
    else
        scatter.parameter <- data.frame(
            parameter = asp$default.scatter.parameter,
            stringsAsFactors = FALSE )

    scatter.parameter$parameter
}

