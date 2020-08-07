# get_worker_process.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


# Returns the number of worker processes, using the number of available logical
# cores as default value.

get.worker.process <- function( worker.process.n )
{
    ifelse( worker.process.n != 0, worker.process.n, detectCores() )
}

