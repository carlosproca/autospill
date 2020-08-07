# create_directory.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


# Creates figure and table directories.
# Returns directories for scatter figures.

create.directory <- function( flow.control, asp )
{
    if ( ! is.null( asp$figure.scatter.dir.base ) ) {
        figure.scatter.dir <- sprintf( "%s_%s", asp$figure.scatter.dir.base,
            flow.control$sample )
        names( figure.scatter.dir ) <- flow.control$sample
    }
    else
        figure.scatter.dir <- NULL

    figure.dir <- c(
        asp$figure.compensation.dir,
        asp$figure.convergence.dir,
        asp$figure.gate.dir,
        asp$figure.skewness.dir,
        asp$figure.slope.error.dir,
        asp$figure.spillover.dir,
        figure.scatter.dir
    )

    table.dir <- c(
        asp$table.compensation.dir,
        asp$table.convergence.dir,
        asp$table.skewness.dir,
        asp$table.slope.error.dir,
        asp$table.spillover.dir
    )

    for ( ftd in c( figure.dir, table.dir ) )
        if ( ! is.null( ftd ) && ! file.exists( ftd ) )
            dir.create( ftd, recursive = TRUE )

    figure.scatter.dir
}

