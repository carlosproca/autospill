#!/bin/sh

# calculate_compensation_website.sh
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


# Runs a calculation of compensation with autospill, creating all figures and
# tables used in autospill website. Converts the final spillover matrix to
# flowjo format.
#
# Requires the environment variable AUTOSPILL_PATH with the location of the 
# autospill package


asp_spillover_dir=./table_spillover
asp_spillover_filebase=autospill_spillover


check_exit_code()
{
    if [ $? -eq 0 ]
    then echo $1 OK
    else echo $1 FAILED
         exit 1
    fi
}


# check definition of AUTOSPILL_PATH

if [ -z "$AUTOSPILL_PATH" ]
then echo ' error: environment variable AUTOSPILL_PATH not defined' 1>&2
     exit 2
fi

asp_batch_dir="$AUTOSPILL_PATH/batch"


# get markers from control files

R CMD BATCH --quiet --vanilla "$asp_batch_dir/get_marker.r" get_marker.out

check_exit_code get_marker


# calculate compensation

R CMD BATCH --quiet --vanilla \
    "$asp_batch_dir/calculate_compensation_website.r" \
    ./calculate_compensation_website.out

check_exit_code calculate_compensation_website


# convert spillover matrix to flowjo mtx format

"$asp_batch_dir/convert_spillover_to_flowjo.py" autospill-spillover \
    $asp_spillover_dir/$asp_spillover_filebase.csv \
    $asp_spillover_dir/$asp_spillover_filebase.mtx \
    > generate_flowjo_spillover.out 2>&1

check_exit_code generate_flowjo_spillover

