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

def add_title(ax,title,style):
    ax.set_title(title+", n={}".format(Y.shape[1]))

def add_vline(ax,xvalues:list,ylim,**kwargs):
    # add vertical lines to ax
    ymin, ymax = ylim
    if not isinstance(xvalues,list):
        raise TypeError("The x values are not list")
    for x in xvalues:
        ax.axvline(x,ymin,ymax,\
                **kwargs)
    return ax

def draw_between(ax,X,Y,):
    """
    shape of Y : n x d
    fill color between +- standard deviation
    """
    if ax==None:
        ax0 = None
        fig, ax = plt.subplots(1,1,figsize=(5,5))
    mean = Y.mean(axis=0)
    std = Y.std(axis=0)
    ax.plot(X,mean,color='k',label=Y.index.values.tolist())
    ax.fill_between(X,mean-std,mean+std,color='r',alpha=.3)
    ax.legend()
    if ax0 == None:
        return ax, fig
    else:
        return ax

def viz_summary(metalevel,summary,droplabel,**kwargs):
    if "repetitive" in metalevel:
        metalevel.remove("repetitive")
    summary.drop(droplabel,axis=1,inplace=True)
    obj = summary.groupby(by=metalevel)
    mean = obj.mean()
    std = obj.std()
    fig,ax = plt.subplots(1,1,**kwargs)
    mean.plot(kind='bar',yerr=std,ax=ax)
    ax.legend()
    plt.xticks(rotation=45,ha='right')
    ax.set_ylabel("Recovery (%)")
    return fig, ax

def viz_phenotype(Y:pd.Series,figsize=(10,10)):
    # Y is saves the phenotypes in index
    # plot figure
    fignum = Y.shape[0]
    rownum = math.ceil(math.sqrt(fignum))
    colnum = math.ceil(fignum/rownum)
    fig,ax = plt.subplots(rownum,colnum,figsize=figsize)
    figcount = 0
    X = Y['wavenumber']
    idx = Y.index.tolist()
    idx.remove("wavenumber")

    for i in idx:
        axi = ax[figcount//colnum,figcount%colnum]
        Y[i].columns = X
        Y[i].T.plot(ax=axi,legend=None)
        axi.set_title(i)
        axi.set_xlabel(None)
        figcount += 1
    plt.close()
    return fig,ax