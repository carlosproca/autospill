# read_marker.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


#' Read markers in a set of flow controls
#'
#' Returns a dataframe and writes a csv file with the common set of markers in
#' a set of single-color controls, together with corrected names in case of
#' presence of forbidden characters.
#'
#' @param control.dir Character string with the directory with the set of
#'     single-color controls.
#' @param control.def.file Character string with the CSV file defining the
#'     names and channels of the single-color controls.
#' @param asp List with AutoSpill parameters.
#'
#' @return Dataframe with the original and corrected names of markers.
#'
#' @references Roca \emph{et al}:
#'     AutoSpill: A method for calculating spillover coefficients to compensate
#'     or unmix high-parameter flow cytometry data.
#'     \emph{bioRxiv} 2020.06.29.177196;
#'     \href{https://doi.org/10.1101/2020.06.29.177196}{doi:10.1101/2020.06.29.177196}
#'     (2020).
#'
#' @seealso \code{\link{get.autospill.param}}.
#'
#' @export

read.marker <- function( control.dir, control.def.file, asp )
{
    # read markers from file if available

    if ( ! is.null( asp$marker.file.name ) &&
            file.exists( asp$marker.file.name ) )
        return( read.table( asp$marker.file.name, sep = ",",
            stringsAsFactors = FALSE ) )

    # read definition of controls

    control <- read.csv( control.def.file, stringsAsFactors = FALSE )

    check.critical( anyDuplicated( control$file.name ) == 0,
        "duplicated filenames in fcs data" )

    # get common set of markers from controls

    flow.set.marker.all <- lapply( control$filename, function( cf )
        colnames( read.FCS( file.path( control.dir, cf ),
            transformation = NULL ) ) )

    flow.set.marker <- flow.set.marker.all[[ 1 ]]

    for( fsm in flow.set.marker.all[ -1 ] )
        flow.set.marker <- intersect( flow.set.marker, fsm )

    # correct marker names

    flow.set.marker.corrected <- flow.set.marker

    for ( fmfc.idx in 1 : nchar( asp$marker.forbidden.char ) )
    {
        fmfc <- substr( asp$marker.forbidden.char, fmfc.idx, fmfc.idx )

        flow.set.marker.corrected <- gsub( fmfc, asp$marker.substitution.char,
            flow.set.marker.corrected, fixed = TRUE )
    }

    # save list of markers

    flow.set.marker.table <- data.frame( flow.set.marker,
        flow.set.marker.corrected, stringsAsFactors = FALSE )

    if ( ! is.null( asp$marker.file.name ) )
        write.table( flow.set.marker.table, file = asp$marker.file.name,
            row.names = FALSE, col.names = FALSE, sep = "," )

    flow.set.marker.table
}

