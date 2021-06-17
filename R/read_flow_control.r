# read_flow_control.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


#' Read set of flow controls
#'
#' Read data and metadata from a set of single-color controls.
#'
#' @param control.dir Character string with the directory with the set of
#'     single-color controls.
#' @param control.def.file Character string with the CSV file defining the
#'     names and channels of the single-color controls.
#' @param asp List with AutoSpill parameters.
#'
#' @return List with the following elements:
#'     \itemize{
#'         \item{antigen}
#'         \item{autof.marker.idx}
#'         \item{event}
#'         \item{event.n}
#'         \item{event.number.width}
#'         \item{event.sample}
#'         \item{expr.data.max}
#'         \item{expr.data.max.ceil}
#'         \item{expr.data.min}
#'         \item{expr.data.tran}
#'         \item{expr.data.untr}
#'         \item{figure.scatter.dir}
#'         \item{flow.set}
#'         \item{gate.parameter}
#'         \item{marker}
#'         \item{marker.n}
#'         \item{marker.original}
#'         \item{sample}
#'         \item{scatter.and.marker}
#'         \item{scatter.and.marker.label}
#'         \item{scatter.and.marker.original}
#'         \item{scatter.parameter}
#'         \item{transform}
#'         \item{transform.inv}
#'         \item{wavelength}
#'     }
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

read.flow.control <- function( control.dir, control.def.file, asp )
{
    # read markers from controls

    flow.set.marker.table <- read.marker( control.dir, control.def.file, asp )

    flow.set.marker <- flow.set.marker.table[[ 1 ]]
    flow.set.marker.corrected <- flow.set.marker.table[[ 2 ]]

    # read definition of controls

    control.table <- read.csv( control.def.file, na.strings = "",
        stringsAsFactors = FALSE )

    check.critical( anyDuplicated( control.table$filename ) == 0,
        "duplicated filenames in fcs data" )

    check.critical( anyDuplicated( control.table$dye ) == 0,
        "same dye assigned to more than one fcs data file" )

    control.table.order <- order( control.table$dye )

    flow.file.name <- control.table$filename[ control.table.order ]
    flow.marker <- control.table$dye[ control.table.order ]
    flow.antigen <- control.table$antigen[ control.table.order ]
    flow.wavelength <- control.table$wavelength[ control.table.order ]

    names( flow.wavelength ) <- flow.marker

    flow.marker.n <- length( flow.marker )

    flow.autof.marker.idx <- which( flow.antigen == asp$antigen.autof )
    if ( length( flow.autof.marker.idx ) != 1 )
        flow.autof.marker.idx <- NULL

    # read scatter parameters

    flow.scatter.parameter <- read.scatter.parameter( asp )

    # set scatter parameters and markers

    flow.scatter.and.marker <- c( flow.scatter.parameter, flow.marker )

    flow.scatter.and.marker.matched.bool <-
        flow.scatter.and.marker %in% flow.set.marker

    if ( ! all( flow.scatter.and.marker.matched.bool ) )
    {
        marker.matched <-
            flow.scatter.and.marker[ flow.scatter.and.marker.matched.bool ]
        flow.scatter.and.marker.unmatched <- paste0(
            sort( setdiff( flow.scatter.and.marker, marker.matched ) ),
            collapse = ", "
        )
        flow.set.unmatched <- paste0(
            sort( setdiff( flow.set.marker, marker.matched ) ),
            collapse = ", "
        )
        error.msg <- sprintf(
            "wrong dye name, not found in fcs data\n\texpected: %s\n\tfound: %s",
            flow.scatter.and.marker.unmatched, flow.set.unmatched
        )
        check.critical( FALSE, error.msg )
    }

    # use corrected names for scatter parameters and markers

    flow.scatter.parameter.original <- flow.scatter.parameter
    flow.scatter.parameter <- flow.set.marker.corrected[
        match( flow.scatter.parameter.original, flow.set.marker ) ]

    flow.marker.original <- flow.marker
    flow.marker <- flow.set.marker.corrected[
        match( flow.marker.original, flow.set.marker ) ]

    flow.scatter.and.marker.original <- flow.scatter.and.marker
    flow.scatter.and.marker <- flow.set.marker.corrected[
        match( flow.scatter.and.marker.original, flow.set.marker ) ]

    check.critical( ! anyDuplicated( flow.scatter.and.marker ),
        "internal error: corrected names for dyes collide" )

    # set labels for scatter parameters and markers

    flow.scatter.and.marker.label <- c( flow.scatter.parameter.original,
        ifelse( ! is.na( flow.antigen ),
            paste0( flow.marker.original, " - ", flow.antigen ),
            flow.marker.original ) )
    names( flow.scatter.and.marker.label ) <- flow.scatter.and.marker

    # set samples

    flow.sample <- flow.marker
    flow.sample.n <- length( flow.sample )

    # read fcs files

    flow.set <- lapply( flow.file.name, function( ff )
        read.FCS( file.path( control.dir, ff ), transformation = NULL ) )

    # get range of fcs data

    flow.set.resolution.read <- sapply( flow.set, function( fs )
        as.numeric( description( fs )[["$P1R" ]] ) )

    check.critical( all( sapply( flow.set.resolution.read, length ) == 1 ),
        "keyword $P1R not found in fcs data" )

    flow.set.resolution <- max( flow.set.resolution.read )

    flow.expr.data.min <- 0
    flow.expr.data.max <- flow.set.resolution
    flow.expr.data.max.ceil <- ceiling( flow.expr.data.max / asp$data.step ) *
        asp$data.step

    # get maximum number of events per sample to adjust event numbering

    flow.sample.event.number.max <- 0

    for ( fs.idx in 1 : flow.sample.n )
    {
        flow.sample.data <- flow.set[[ fs.idx ]]

        flow.sample.event.number <- nrow( exprs( flow.sample.data ) )

        if ( flow.sample.event.number > flow.sample.event.number.max )
            flow.sample.event.number.max <- flow.sample.event.number
    }

    flow.event.number.width <-
        floor( log10( flow.sample.event.number.max ) ) + 1
    flow.event.regexp <- sprintf( "\\.[0-9]{%d}$", flow.event.number.width )

    # make preliminary control info

    flow.control <- list(
        antigen = flow.antigen,
        autof.marker.idx = flow.autof.marker.idx,
        event.number.width = flow.event.number.width,
        expr.data.max = flow.expr.data.max,
        expr.data.max.ceil = flow.expr.data.max.ceil,
        expr.data.min = flow.expr.data.min,
        flow.set = flow.set,
        marker = flow.marker,
        marker.n = flow.marker.n,
        marker.original = flow.marker.original,
        sample = flow.sample,
        scatter.and.marker = flow.scatter.and.marker,
        scatter.and.marker.label = flow.scatter.and.marker.label,
        scatter.and.marker.original = flow.scatter.and.marker.original,
        scatter.parameter = flow.scatter.parameter,
        wavelength = flow.wavelength
    )

    # get expression data

    flow.expr.data.untr <- get.flow.expression.data( flow.set, flow.control )

    # set events

    flow.event <- rownames( flow.expr.data.untr )

    flow.event.n <- length( flow.event )

    flow.event.sample <- sub( flow.event.regexp, "", flow.event )
    flow.event.sample <- factor( flow.event.sample, levels = flow.sample )

    # get transformations

    flow.transform.both <- get.transformation( flow.control, asp )

    flow.transform <- flow.transform.both[[ 1 ]]
    flow.transform.inv <- flow.transform.both[[ 2 ]]

    # get transformed expression data

    flow.set.tran <- lapply( flow.set, transform,
        transformList( names( flow.transform ), flow.transform ) )

    flow.expr.data.tran <- get.flow.expression.data( flow.set.tran,
        flow.control )

    # read gate parameters

    flow.gate.parameter <- read.gate.parameter( flow.control, asp )

    # create figure and table directories

    flow.figure.scatter.dir <- create.directory( flow.control, asp )

    # update control info

    flow.control$event <- flow.event
    flow.control$event.n <- flow.event.n
    flow.control$event.sample <- flow.event.sample
    flow.control$expr.data.tran <- flow.expr.data.tran
    flow.control$expr.data.untr <- flow.expr.data.untr
    flow.control$gate.parameter <- flow.gate.parameter
    flow.control$figure.scatter.dir <- flow.figure.scatter.dir
    flow.control$transform <- flow.transform
    flow.control$transform.inv <- flow.transform.inv

    flow.control
}

