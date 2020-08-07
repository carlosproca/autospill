# calculate_compensation_final_step.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


# Runs a calculation of compensation with autospill, creating figures and
# tables at the final step.
#
# Requires assigning proper values to the variables:
#     control.dir    directory with the set of single-color controls
#     control.def.file    csv file defining the names and channels of the
#         single-color controls


library( autospill )


# set parameters

asp <- get.autospill.param( "final.step" )


# read flow controls

control.dir <- "../fcs_control_data/"
control.def.file <- "../fcs_control_data/fcs_control.csv"

flow.control <- read.flow.control( control.dir, control.def.file, asp )


# gate events before calculating spillover

flow.gate <- gate.flow.data( flow.control, asp )


# get initial spillover matrices from untransformed and transformed data

marker.spillover.unco.untr <- get.marker.spillover( TRUE, flow.gate,
    flow.control, asp )
marker.spillover.unco.tran <- get.marker.spillover( FALSE, flow.gate,
    flow.control, asp )


# refine spillover matrix iteratively

refine.spillover.result <- refine.spillover( marker.spillover.unco.untr,
    marker.spillover.unco.tran, flow.gate, flow.control, asp )

