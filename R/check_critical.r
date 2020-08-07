# check_critical.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


# Checks condition, if not true prints error message and stops execution.

check.critical <- function( condition, error.msg )
{
    if ( ! all( condition ) )
    {
        cat( sprintf( "ERROR: %s\n", error.msg ), file = stderr() )
        stop()
    }
}

