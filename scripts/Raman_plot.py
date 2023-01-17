#!/usr/bin/env python3
#############################################
# This script is to plot basic figure for Raman data
#############################################

from load import *

def plot_line(X,Y,xlabel='Wavenumber (c$m^{-1}$)',\
        ylabel='Intensity (.u.)',
        ax=None,legend=False):
    # plot spectra as line 
    if ax==None:
        ax0 = None
        fig, ax = plt.subplots(1,1,figsize=(5,5))
    n,d = Y.shape
    labels = Y.index
    for i in range(n):
        ax.plot(X,Y.iloc[i,:],label=labels[i])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if legend==True:
        ax.legend()
    if ax0 == None:
        return ax, fig
    else:
        return ax

def add_vline(ax,xvalues,ylim,**kwargs):
    # add vertical lines to ax
    ymin, ymax = ylim
    if not isinstance(xvalues,list):
        raise TypeError("The x values are not list")
    for x in xvalues:
        ax.axvline(x,ymin,ymax,\
                **kwargs)
    return ax

