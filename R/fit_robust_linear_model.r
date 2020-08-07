# fit_robust_linear_model.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


# Returns a matrix by rows, with the intercept and p-value, and coefficient and
# p-value, of a robust linear model fitted to the input data.
#
# Reverts to a standard linear model in case of no convergence.

fit.robust.linear.model <- function( x.data, y.data, x.name, y.name, asp )
{
    xy.data <- data.frame( x = x.data, y = y.data )

    xy.model <- rlm( y ~ x, xy.data, maxit = asp$rlm.iter.max )

    if ( xy.model$converged )
    {
        xy.coef <- xy.model$coefficients
        xy.t <- summary( xy.model )$coefficients[ , 3 ]
        xy.df <- summary( xy.model )$df[ 2 ]
        xy.pval <- 2*( pt( abs( xy.t ), xy.df, lower.tail = FALSE ) )
    }
    else
    {
        cat( sprintf( "WARNING: rlm of %s ~ %s did not converge - then using ols\n",
            y.name, x.name ), file = stderr() )

        xy.model <- lm( y ~ x, xy.data )

        xy.coef <- xy.model$coefficients
        xy.pval <- summary( xy.model )$coefficients[ , 4 ]
    }

    res <- cbind( xy.coef, xy.pval )
    dimnames( res ) <- NULL
    res
}

