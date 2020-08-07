# get_autospill_param_final_step.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


# Returns parameters for running a calculation of compensation with autospill,
# creating figures and tables at final step.

get.autospill.param.final.step <- function()
{
    autosp.param <- get.autospill.param.minimal()

    # directory parameters

    autosp.param$figure.scatter.dir.base <- "figure_scatter"

    autosp.param$figure.gate.dir <- "figure_gate"

    autosp.param$figure.compensation.dir <- "figure_compensation"
    autosp.param$figure.convergence.dir <- "figure_convergence"
    autosp.param$figure.spillover.dir <- "figure_spillover"
    autosp.param$figure.slope.error.dir <- "figure_slope_error"
    autosp.param$figure.skewness.dir <- "figure_skewness"

    autosp.param$table.compensation.dir <- "table_compensation"
    autosp.param$table.convergence.dir <- "table_convergence"
    autosp.param$table.spillover.dir <- "table_spillover"
    autosp.param$table.slope.error.dir <- "table_slope_error"
    autosp.param$table.skewness.dir <- "table_skewness"

    autosp.param
}

