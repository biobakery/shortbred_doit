#! /usr/bin/env python

"""
---------------------------------------
Eric Franzosa ( eric.franzosa@gmail.com )
"""

import os, sys, re, glob, argparse
from math import log10
import numpy as np
import matplotlib.pyplot as plt
import zopy.mplutils as mu

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_fStepWidth = 0.4
c_fHalfStepWidth = c_fStepWidth / 2.0
c_fStepBuffer = 0.1

# ---------------------------------------------------------------
# functions for making step plots
# ---------------------------------------------------------------

def stepplot ( ax, aaData, colors=None, thickness=3.0, alpha=0.5, connect_color="gray", median=False ):
    """
    The rows of aaData are samples
    The columns are a type of step (~group in a barplot)
    """
    iGroups = len( aaData[0] )
    aGroupIndex = mu.funcGroupIndex( ax, iGroups )
    if colors is None:
        colors = ["cornflowerblue" for k in aGroupIndex]
    elif type( colors ) is str:
        colors = [colors for k in aGroupIndex]
    # draw hashes
    for aSteps in aaData:
        for i in range( len( aSteps ) ):
            if aSteps[i] is not None:
                aX = [ aGroupIndex[i] - c_fHalfStepWidth, 
                       aGroupIndex[i] + c_fHalfStepWidth ]
                aY = [ aSteps[i], aSteps[i] ]
                ax.add_line( plt.Line2D( aX, aY, alpha=alpha, lw=thickness, color=colors[i] ) )   
    # draw connections
    for aSteps in aaData:
        for i in range( len( aSteps ) ):
            if i > 0:
                if aSteps[i] is not None and aSteps[i-1] is not None:
                    aX = [ aGroupIndex[i-1] + c_fHalfStepWidth + c_fStepBuffer, 
                           aGroupIndex[i]   - c_fHalfStepWidth - c_fStepBuffer ]
                    aY = [ aSteps[i-1],      aSteps[i]      ]
                    ax.add_line( plt.Line2D( aX, aY, alpha=alpha, lw=thickness / 2.0, color=connect_color ) )
    # typical ax settings
    ax.xaxis.set_ticks( aGroupIndex )
    # add median?
    if median:
        mdata = [[] for k in aaData[0]]
        for i in range( len( aaData[0] ) ):
            for sample in aaData:
                mdata[i].append( sample[i] )
        mdata = [map( np.median, mdata )]
        stepplot( ax, mdata, colors="black", thickness=thickness*2, alpha=1.0, connect_color="black" )
