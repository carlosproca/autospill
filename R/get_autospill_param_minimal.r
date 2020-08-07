# get_autospill_param_minimal.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


# Returns parameters for running a calculation of compensation with autospill,
# without creating any figures or tables.

get.autospill.param.minimal <- function()
{
    color.pal <- brewer.pal( 9, "Set1" )

    list(
        # general parameters

        verbose = TRUE,

        worker.process.n = 0,

        default.scatter.parameter = c( "FSC-A", "SSC-A" ),

        antigen.autof = "AutoF",

        marker.forbidden.char = " !\"#$%&'()*,/:;?@[\\]^{|}~",
        marker.substitution.char = "-",

        # algorithm parameters

        default.transformation.param = list(
            length = 256,
            max.range = 262144,
            pos = 4.418539922,
            neg = 0,
            width = -100
        ),

        default.gate.param = list(
            density.threshold = 0.33,
            region.auto = TRUE,
            region.factor.x.low = 0.05,
            region.factor.x.high = 0.80,
            region.factor.y.low = 0.05,
            region.factor.y.high = 0.80
        ),

        gate.data.trim.factor.x.min = 0.01,
        gate.data.trim.factor.x.max = 0.99,
        gate.data.trim.factor.y.min = 0.01,
        gate.data.trim.factor.y.max = 0.99,

        gate.bound.density.bw.factor = 3.0,
        gate.bound.density.grid.n = 100,
        gate.bound.density.neigh.size = 3,

        gate.bound.density.max.target = 1,
        gate.bound.density.max.exclusion.corner = 0.05,
        gate.bound.density.max.mad.factor = 3.0,

        gate.region.density.bw.factor = 2.0,
        gate.region.density.grid.n = 100,
        gate.region.density.neigh.size = 2,

        gate.region.max.density.bw.factor = 1.0,
        gate.region.max.density.grid.n = 100,

        rlm.iter.max = 100,
        rlm.trim.factor = 0.01,

        rs.iter.max = 100,

        rs.lambda.coarse = 1.0,
        rs.lambda.fine = 0.1,

        rs.delta.history.n = 10,

        rs.delta.threshold.untr = 1e-2,
        rs.delta.threshold.tran = 1e-4,
        rs.delta.threshold.change = 1e-6,

        # directory parameters

        figure.scatter.dir.base = NULL,

        figure.gate.dir = NULL,

        figure.compensation.dir = NULL,
        figure.convergence.dir = NULL,
        figure.spillover.dir = NULL,
        figure.slope.error.dir = NULL,
        figure.skewness.dir = NULL,

        table.compensation.dir = NULL,
        table.convergence.dir = NULL,
        table.spillover.dir = NULL,
        table.slope.error.dir = NULL,
        table.skewness.dir = NULL,

        # filename parameters

        marker.file.name = NULL,

        gate.parameter.file.name = "fcs_gate_parameter.csv",
        scatter.parameter.file.name = "fcs_scatter_parameter.csv",
        transformation.parameter.file.name =
            "fcs_transformation_parameter.csv",

        convergence.file.name = "autospill_convergence",

        compensation.file.name = "autospill_compensation",
        spillover.file.name = "autospill_spillover",
        slope.error.file.name = "autospill_slope_error",
        skewness.file.name = "autospill_skewness",

        compensation.popnegpop.file.name = "posnegpop_compensation",
        spillover.popnegpop.file.name = "posnegpop_spillover",
        slope.error.posnegpop.file.name = "posnegpop_slope_error",
        skewness.posnegpop.file.name = "posnegpop_skewness",

        posnegpop.file.label = "posnegpop",

        # table parameters

        rs.save.table.initial = FALSE,
        rs.save.table.every = 0,

        # graphics parameters

        rs.plot.figure.initial = FALSE,
        rs.plot.figure.every = 0,

        data.step = 50e3,

        convergence.color.delta = color.pal[ 7 ],    # brown
        convergence.color.delta.max = color.pal[ 5 ],    # orange
        convergence.color.delta.change = color.pal[ 8 ],    # pink
        convergence.shape.linear = "triangle",
        convergence.shape.biexp = "circle",
        convergence.shape.posnegpop = "triangle open",

        density.color.single = "blue3",
        density.color.initial = color.pal[ 3 ],    # green
        density.color.final = color.pal[ 2 ],    # blue
        density.color.posnegpop = color.pal[ 1 ],    # red
        density.palette.n = 1000,
        density.palette.base.n = 1000000,
        density.palette.base.color = c( "blue", "cyan", "green", "yellow",
            "red" ),

        gate.plot.stage = FALSE,
        gate.tesselation.color = "blue3",

        matrix.marker.color = "blue3",
        matrix.marker.proper.color = "red3",
        matrix.wavelength.min = 300,
        matrix.wavelength.max = 900,

        scatter.plot.scale.other = FALSE,
        scatter.expr.color.unco = "blue",
        scatter.expr.color.comp = "black",
        scatter.scale.breaks.coef = 2.0,
        scatter.ref.line.unco = FALSE,
        scatter.ref.line.color = "green3",
        scatter.ref.line.size.factor = 1.8,

        # parameters for main figures

        figure.width = 8.0,
        figure.height = 6.0,
        figure.margin = 4.0,

        figure.panel.line.size = 0.5,

        figure.axis.text.size = 12.0,
        figure.axis.title.size = 12.0,

        figure.convergence.point.size = 2.0,
        figure.convergence.line.size = 0.8,

        figure.density.line.size = 0.2,

        figure.gate.scale.expand = 0.01,
        figure.gate.point.size = 0.8,
        figure.gate.line.size = 0.5,
        figure.gate.bar.width = 1.0,
        figure.gate.bar.height = 25.0,
        figure.gate.bar.margin = 2.0,

        figure.matrix.point.size = 2.5,
        figure.matrix.line.size = 0.8,

        figure.scatter.alpha.gate.in = 0.8,
        figure.scatter.alpha.gate.out = 0.1,
        figure.scatter.point.size = 0.8,
        figure.scatter.line.size = 0.6,
        figure.scatter.error.label.size = 4.0,
        figure.scatter.error.label.pos.x = 0.90,
        figure.scatter.error.label.pos.y = 0.05,
        figure.scatter.axis.text.size = 12.0,
        figure.scatter.axis.title.size = 12.0,

        # thumbnail parameters

        make.thumbnail = FALSE,

        thumbnail.width = 2.0,
        thumbnail.height = 1.5,
        thumbnail.margin = 2.0,

        thumbnail.panel.line.size = 0.2,

        thumbnail.axis.text.size = 5.0,
        thumbnail.axis.title.size = 5.0,

        thumbnail.convergence.point.size = 0.8,
        thumbnail.convergence.line.size = 0.3,

        thumbnail.density.line.size = 0.4,

        thumbnail.gate.scale.expand = 0.01,
        thumbnail.gate.point.size = 0.3,
        thumbnail.gate.line.size = 0.2,
        thumbnail.gate.bar.width = 0.4,
        thumbnail.gate.bar.height = 6.0,
        thumbnail.gate.bar.margin = 0.0,

        thumbnail.matrix.point.size = 1.0,
        thumbnail.matrix.line.size = 0.3,

        thumbnail.scatter.alpha.gate.in = 0.8,
        thumbnail.scatter.alpha.gate.out = 0.1,
        thumbnail.scatter.point.size = 0.3,
        thumbnail.scatter.line.size = 0.2,
        thumbnail.scatter.error.label.size = 1.8,
        thumbnail.scatter.error.label.pos.x = 0.85,
        thumbnail.scatter.error.label.pos.y = 0.05,
        thumbnail.scatter.axis.text.size = 5.0,
        thumbnail.scatter.axis.title.size = 5.0
    )
}

