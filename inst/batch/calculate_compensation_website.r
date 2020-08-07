# calculate_compensation_website.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


# Runs a calculation of compensation with autospill, creating all figures and
# tables used in autospill website.
#
# Requires being called as a batch script, and assumes fixed values for the
# following two variables (see below):
#     control.dir    directory with the set of single-color controls
#     control.def.file    csv file defining the names and channels of the
#         single-color controls


library( autospill )


# set parameters

asp <- get.autospill.param( "website" )

asp


# read flow controls

control.dir <- "./samples/"
control.def.file <- "./fcs_control.csv"

flow.control <- read.flow.control( control.dir, control.def.file, asp )

names( flow.control )

flow.control$antigen
flow.control$autof.marker.idx
str( flow.control$event )
flow.control$event.n
flow.control$event.number.width
table( flow.control$event.sample )
flow.control$expr.data.max
flow.control$expr.data.max.ceil
flow.control$expr.data.min
str( flow.control$expr.data.tran )
str( flow.control$expr.data.untr )
flow.control$figure.scatter.dir
str( flow.control$flow.set, max.level = 1 )
str( flow.control$gate.parameter )
flow.control$marker
flow.control$marker.n
flow.control$marker.original
flow.control$sample
flow.control$scatter.and.marker
flow.control$scatter.and.marker.label
flow.control$scatter.and.marker.original
flow.control$scatter.parameter
str( flow.control$transform )
str( flow.control$transform.inv )
flow.control$wavelength


# gate events before calculating spillover

flow.gate <- gate.flow.data( flow.control, asp )

str( flow.gate )


# get initial spillover matrices from untransformed and transformed data

marker.spillover.unco.untr <- get.marker.spillover( TRUE, flow.gate,
    flow.control, asp )
marker.spillover.unco.tran <- get.marker.spillover( FALSE, flow.gate,
    flow.control, asp )

str( marker.spillover.unco.untr$inte )
str( marker.spillover.unco.untr$coef )

str( marker.spillover.unco.tran$inte )
str( marker.spillover.unco.tran$coef )


# refine spillover matrix iteratively

refine.spillover.result <- refine.spillover( marker.spillover.unco.untr,
    marker.spillover.unco.tran, flow.gate, flow.control, asp )

str( refine.spillover.result )


# output session info

sessionInfo()

