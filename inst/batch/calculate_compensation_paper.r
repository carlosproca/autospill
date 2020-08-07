# calculate_compensation_paper.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


# Runs a calculation of compensation with autospill, creating all figures and
# tables used in autospill paper.
#
# Requires being called as a batch script with the following two arguments:
#     control.dir    directory with the set of single-color controls
#     control.def.file    csv file defining the names and channels of the
#         single-color controls


library( autospill )


# get directory and csv file with definition of control dataset

args <- commandArgs( TRUE )

if ( length( args ) != 2 ) {
    cat( "ERROR: no arguments with directory and csv file with definition of control dataset",
        file = stderr() )
    stop()
}

control.dir <- args[[ 1 ]]
control.def.file <- args[[ 2 ]]


# set parameters

asp <- get.autospill.param( "paper" )


# read flow controls

flow.control <- read.flow.control( control.dir, control.def.file, asp )


# gate events before calculating spillover

flow.gate <- gate.flow.data( flow.control, asp )


# get initial spillover matrices from untransformed and transformed data

marker.spillover.unco.untr <- get.marker.spillover( TRUE, flow.gate,
    flow.control, asp )
marker.spillover.unco.tran <- get.marker.spillover( FALSE, flow.gate,
    flow.control, asp )


# get spillover and compensation matrices with positive and negative
# populations

spillover.error.posnegpop <- process.posnegpop( marker.spillover.unco.untr,
    flow.gate, flow.control, asp )


# refine spillover matrix iteratively

refine.spillover.result <- refine.spillover( marker.spillover.unco.untr,
    marker.spillover.unco.tran, flow.gate, flow.control, asp )


# plot results together for slope error and skewness

plot_result.together( flow.control, asp )


# replot convergence adding spillover error from the calculation with positive
# and negative populations

plot_convergence( refine.spillover.result$convergence,
    spillover.error.posnegpop$error$slop, asp )


# output session info

sessionInfo()

