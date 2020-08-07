# get_autospill_param.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


#' Get AutoSpill parameters
#'
#' Returns parameters for running a calculation of compensation with AutoSpill.
#'
#' @param param.set Character string with the name of a parameter set. Posible
#'     values are:
#'     \describe{
#'         \item{\code{"minimal"}}{No figures or tables.}
#'         \item{\code{"final.step"}}{With figures and tables at final step.}
#'         \item{\code{"paper"}}{With all figures and tables used in AutoSpill
#'             paper.}
#'         \item{\code{"website"}}{With all figures and tables used in AutoSpill
#'             website.}
#'     }
#'
#' @return List of AutoSpill parameters. As reference, see source file
#'     \file{get_autospill_param_minimal.r}
#'
#' @references Roca \emph{et al}:
#'     AutoSpill: A method for calculating spillover coefficients to compensate
#'     or unmix high-parameter flow cytometry data.
#'     \emph{bioRxiv} 2020.06.29.177196;
#'     \href{https://doi.org/10.1101/2020.06.29.177196}{doi:10.1101/2020.06.29.177196}
#'     (2020).
#'
#' @export

get.autospill.param <- function( param.set = "minimal" )
{
    get.param.function <- get0( sprintf( "get.autospill.param.%s", param.set ) )

    check.critical( ! is.null( get.param.function ), "bad param set" )

    get.param.function()
}

