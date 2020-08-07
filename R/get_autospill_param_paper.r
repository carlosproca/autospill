# get_autospill_param_paper.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


# Returns parameters for running a calculation of compensation with autospill,
# creating all figures and tables used in autospill paper.

get.autospill.param.paper <- function()
{
    autosp.param <- get.autospill.param.final.step()

    # table parameters

    autosp.param$rs.save.table.initial <- TRUE

    # graphics parameters

    autosp.param$rs.plot.figure.initial <- TRUE
    autosp.param$gate.plot.stage <- TRUE
    autosp.param$scatter.plot.scale.other <- TRUE

    # parameters for main figures

    autosp.param$figure.width <- 3.25
    autosp.param$figure.height <- 2.60
    autosp.param$figure.margin <- 1.0

    autosp.param$figure.panel.line.size <- 0.5

    autosp.param$figure.axis.text.size <- 7.0
    autosp.param$figure.axis.title.size <- 7.0

    autosp.param$figure.convergence.point.size <- 1.0
    autosp.param$figure.convergence.line.size <- 0.4

    autosp.param$figure.density.line.size <- 0.4

    autosp.param$figure.gate.scale.expand <- 0.01
    autosp.param$figure.gate.point.size <- 0.4
    autosp.param$figure.gate.line.size <- 0.3
    autosp.param$figure.gate.bar.width <- 0.3
    autosp.param$figure.gate.bar.height <- 10.0
    autosp.param$figure.gate.bar.margin <- 0.0

    autosp.param$figure.matrix.point.size <- 1.2
    autosp.param$figure.matrix.line.size <- 0.4

    autosp.param$figure.scatter.alpha.gate.in <- 0.8
    autosp.param$figure.scatter.alpha.gate.out <- 0.2
    autosp.param$figure.scatter.point.size <- 0.6
    autosp.param$figure.scatter.line.size <- 0.5
    autosp.param$figure.scatter.error.label.size <- 3.0
    autosp.param$figure.scatter.error.label.pos.x <- 0.85
    autosp.param$figure.scatter.error.label.pos.y <- 0.05
    autosp.param$figure.scatter.axis.text.size <- 9.0
    autosp.param$figure.scatter.axis.title.size <- 9.0

    autosp.param
}

