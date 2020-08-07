#!/bin/sh

# calculate_compensation_paper.sh
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
# Requires the environment variable AUTOSPILL_PATH with the location of the 
# autospill package


# check definition of AUTOSPILL_PATH

if [ -z "$AUTOSPILL_PATH" ]
then echo ' error: environment variable AUTOSPILL_PATH not defined' 1>&2
     exit 2
fi


# check argument with directory and csv file with definition of control dataset

if [ $# -ne 2 ]
then echo ' usage: calculate_compensation_paper.sh  control_dir' \
         ' control_def_file' 1>&2
     exit 2
fi


# calculate compensation

R CMD BATCH --quiet --vanilla "--args $1 $2" \
    "$AUTOSPILL_PATH/batch/calculate_compensation_paper.r" \
    ./calculate_compensation_paper.out

