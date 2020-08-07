# get_marker.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


# Writes a csv file with the common set of markers in a set of single-color
# controls, together with corrected names in case of presence of forbidden
# characters.
#
# Requires being called as a batch script, and assumes fixed values for the
# following two variables (see below):
#     control.dir    directory with the set of single-color controls
#     control.def.file    csv file defining the names and channels of the
#         single-color controls


library( autospill )


# set parameters

asp <- get.autospill.param()

asp$marker.file.name


# read markers from controls

control.dir <- "./samples/"
control.def.file <- "./fcs_control.csv"

flow.set.marker.table <- read.marker( control.dir, control.def.file, asp )

flow.set.marker.table


# output session info

sessionInfo()

