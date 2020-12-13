
# -*- coding: utf-8 -*-
"""
toolCase.py
"""

from __future__ import print_function

import os, sys, math #, shutil
import numpy as np
import matplotlib.ticker as ticker

from itertools  import product 

sys.dont_write_bytecode = True

"###################################"

def chk_dir(dir):
    check = True
    if os.path.isdir(dir) == False:
        check = False
    return check

"###################################################################################" 

class OOMFormatter(ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_orderOfMagnitude(self, nothing):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin, vmax):
        self.format = self.fformat
        if self._useMathText:
            self.format = '$%s$' % ticker._mathdefault(self.format)

"###################################################################################" 

class OOMFormatter_eng(ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.0f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_orderOfMagnitude(self, nothing):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin, vmax):
        self.format = self.fformat
        if self._useMathText:
            self.format = '$%s$' % ticker._mathdefault(self.format)

"#################################################"

def OOMUp(x):

    if abs(x) <= 1.e-16: return -16

    x = math.log10(x) if x > 0. else math.log10(-x)
    x = math.ceil (x) if x > 0. else math.floor( x)

    return x 

"#################################################"

def OOMDown(x):

    if abs(x) <= 1.e-16: return -16

    x = math.log10(x) if x > 0. else math.log10(-x)
    x = math.floor(x) if x > 0. else math.ceil ( x)

    return x 

"#################################################"

def OOMRoundUp(x):

    x = OOMUp(x)
    x = math.pow(10.0,x) if x > -16 else 0

    return x 

"#################################################"

def OOMRoundDown(x):

    x = OOMDown(x)
    x = math.pow(10.0,x) if x > -16 else 0

    return x 

"##############################################################" "check it out here: ./legacy/modUtils.py"

"Element in nd array `a` closest to the scalar value `a0`"

def find_nearest(a, a0, ax):

    idx = np.argmin(abs(a - a0), axis=ax)

    return idx

"#########################################"

def ProbeAtLocation(z, x, y, x0, y0):

    idx = find_nearest(x, x0, 1)
    idy = find_nearest(y, y0, 0)

    probe_dim = 5 ; probe_stencil = product(range(idy[0], idy[0] + probe_dim), range(idx[0], idx[0] + probe_dim))

    z_probe = np.average(np.array([z[k] for k in probe_stencil])) # = z[idy[0],idx[0]]

    return z_probe

"#####################################################################################" 

#Label line with line2D label data
def labelLine(line,x,label=None,align=True,**kwargs):

    ax = line.axes
    xdata = line.get_xdata()
    ydata = line.get_ydata()

    if (x < xdata[0]) or (x > xdata[-1]):
        print('x label location is outside data range!')
        return

    #Find corresponding y co-ordinate and angle of the line
    ip = 1
    for i in range(len(xdata)):
        if x < xdata[i]:
            ip = i
            break

    y = ydata[ip-1] + (ydata[ip]-ydata[ip-1])*(x-xdata[ip-1])/(xdata[ip]-xdata[ip-1])

    if not label:
        label = line.get_label()

    if align:
        #Compute the slope
        dx = xdata[ip] - xdata[ip-1]
        dy = ydata[ip] - ydata[ip-1]
        ang = degrees(atan2(dy,dx))

        #Transform to screen co-ordinates
        pt = np.array([x,y]).reshape((1,2))
        trans_angle = ax.transData.transform_angles(np.array((ang,)),pt)[0]

    else:
        trans_angle = 0

    #Set a bunch of keyword arguments
    if 'color' not in kwargs:
        kwargs['color'] = line.get_color()

    if ('horizontalalignment' not in kwargs) and ('ha' not in kwargs):
        kwargs['ha'] = 'center'

    if ('verticalalignment' not in kwargs) and ('va' not in kwargs):
        kwargs['va'] = 'center'

    if 'backgroundcolor' not in kwargs:
#        kwargs['backgroundcolor'] = ax.get_facecolor()
        kwargs['backgroundcolor'] = ax.get_axis_bgcolor()

    if 'clip_on' not in kwargs:
        kwargs['clip_on'] = True

    if 'zorder' not in kwargs:
        kwargs['zorder'] = 2.5

    ax.text(x,y,label,rotation=trans_angle,**kwargs)

"#####################################################################################" 

def labelLines(lines,align=True,xvals=None,**kwargs):

    ax = lines[0].axes
    labLines = []
    labels = []

    #Take only the lines which have labels other than the default ones
    for line in lines:
        label = line.get_label()
        if "_line" not in label:
            labLines.append(line)
            labels.append(label)

    if xvals is None:
        xmin,xmax = ax.get_xlim()
        xvals = np.linspace(xmin,xmax,len(labLines)+2)[1:-1]

    for line,x,label in zip(labLines,xvals,labels):
        labelLine(line,x,label,align,**kwargs)

