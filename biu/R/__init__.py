"""
This module provides many functions available in R. This is to make it easier to convert code between R and python. These functions 

At some point, I hope it will contain versions of all the following methods:
http://www.sr.bham.ac.uk/~ajrs/R/r-function_list.html
"""

from ..ops.array import order, pmin, pmax, cummin, cummax

from ..stats.p_adjust import p_adjust

from ..stats.regression import lowess
