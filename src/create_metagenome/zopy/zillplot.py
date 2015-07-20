#! /usr/bin/env python

"""
Module for making "ZILL" ( Zero-Inflated Log-log ) plots
---------------------------------------
Eric Franzosa ( eric.franzosa@gmail.com )
"""

import os, sys, re, glob, argparse
from math import log10
import numpy as np
import matplotlib.pyplot as plt
import zopy.mplutils as mu
from zopy.utils import reader

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_iDashSize  = 75
c_iLogMargin = 4

# ---------------------------------------------------------------
# set up the plotting area
# ---------------------------------------------------------------

def funcZillPlotArea ( ax, aX, aY, logmin=None, logmax=None ):
    """ set up the plotting area """
    # decide how to scale axes
    data_logmin, data_logmax = mu.funcLogMinMax( aX + aY )
    logmin = data_logmin if logmin is None else logmin
    logmax = data_logmax if logmax is None else logmax
    # set up logmin margin ( extra space near x/y axes )   
    logmin_margin = logmin / 10 * c_iLogMargin
    # math to get margin midpoint ( tricky because we're in log space ); this is where the "zeroes" go            
    zeroval = 10 ** ( np.mean( [np.log10( logmin_margin ), np.log10( logmin )] ) )
    # configure axes
    mu.funcSetTickParams( ax )
    ax.set_xscale( "log" )
    ax.set_yscale( "log" )
    ax.set_xlim( logmin_margin, logmax )
    ax.set_ylim( logmin_margin, logmax )
    # adjust minor ticks ( don't want anything outside logmin/logmax )
    aMinorTickLocations = [ logmin + logmin * i for i in range( 1, 10 ) ]
    while aMinorTickLocations[-1] < logmax:
        aTemp = [ 10 * k for k in aMinorTickLocations[-9:] ]
        aMinorTickLocations += aTemp
    ax.set_xticks( aMinorTickLocations, minor=True )
    ax.set_yticks( aMinorTickLocations, minor=True )
    # add a default grid
    mu.funcGrid( ax, xy=True )
    # the point plotter needs to know this
    return zeroval

# ---------------------------------------------------------------
# plot the points
# ---------------------------------------------------------------

def funcIndex ( aX, aY ):
    """ provides the index positions of ( 0,y ) ( x,0 ) and ( 0,0 ) points """
    aIndexZeroX = []
    aIndexZeroY = []
    aIndexJoint = []
    discard = 0
    for i, ( x, y ) in enumerate( zip( aX, aY ) ):
        if x > 0 and y > 0:
            aIndexJoint.append( i )
        elif x > 0 and y == 0:
            aIndexZeroY.append( i )
        elif y > 0 and x == 0:
            aIndexZeroX.append( i )
        else:
            discard += 1
    print >>sys.stderr, "ATTENTION ( ZillPlotIndex ):", "Will be discarding", discard, "( 0,0 ) points"
    return aIndexZeroX, aIndexZeroY, aIndexJoint

def funcSubVector ( aItems, aIndex ):
    """ perform a slice of a vector based on index positions """
    return [aItems[i] for i in aIndex]

def funcZillPlotPoints( ax, aX, aY, zero=None, **kwargs ):
    """ add the point; needs to inherit zero location from Area function """
    aIndexZeroX, aIndexZeroY, aIndexJoint = funcIndex( aX, aY )
    # plot the joint points
    kwargs2 = { k:funcSubVector( v, aIndexJoint ) if isinstance( v, list ) else v for k, v in kwargs.items() }
    mu.funcScatter( ax, funcSubVector( aX, aIndexJoint ), funcSubVector( aY, aIndexJoint ), **kwargs2 )
    # plot the zero-x points ( note kwarg overrides )
    kwargs2 = { k:funcSubVector( v, aIndexZeroX ) if isinstance( v, list ) else v for k, v in kwargs.items() }
    kwargs2["marker"] = "_"
    kwargs2["edgecolor"] = kwargs2["color"]
    kwargs2["s"] = c_iDashSize if "s" not in kwargs2 else kwargs2["s"]
    mu.funcScatter( ax, [zero for k in aIndexZeroX], funcSubVector( aY, aIndexZeroX ), **kwargs2 )
    # plot the zero-y points ( note kwarg overrides )
    kwargs2 = { k:funcSubVector( v, aIndexZeroY ) if isinstance( v, list ) else v for k, v in kwargs.items() }
    kwargs2["marker"] = "|"
    kwargs2["edgecolor"] = kwargs2["color"]
    kwargs2["s"] = c_iDashSize if "s" not in kwargs2 else kwargs2["s"]
    mu.funcScatter( ax, funcSubVector( aX, aIndexZeroY ), [zero for k in aIndexZeroY], **kwargs2 )

# ---------------------------------------------------------------
# main function
# ---------------------------------------------------------------

def zillplot( ax, aX, aY, logmin=None, logmax=None, color="blue", edgecolor="none", alpha=0.5, **kwargs ):
    """ the main function; call this if importing into other scripts """
    # set up plot area and retrieve zero position
    zero = funcZillPlotArea( ax, aX, aY, logmin, logmax )
    # plot points
    funcZillPlotPoints( ax, aX, aY, zero=zero, color=color, edgecolor=edgecolor, alpha=alpha, **kwargs )
