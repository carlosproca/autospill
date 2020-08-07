# get_autospill_param_website.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


# Returns parameters for running a calculation of compensation with autospill,
# creating all figures and tables used in autospill website.

get.autospill.param.website <- function()
{
    autosp.param <- get.autospill.param.final.step()

    # filename parameters

    autosp.param$marker.file.name <- "fcs_marker.csv"

    # thumbnail parameters

    autosp.param$make.thumbnail <- TRUE

    autosp.param
}

