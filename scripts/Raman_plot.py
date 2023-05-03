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

def draw_between(X,Y,ax=None,xlim=(300,2000),ylim=(0,10),xstep=200):
    """
    plot mean and +/- 1 std for Raman spectra
    shape of Y : n x d
    """
    # upper fig 1
    fig, ax = plt.subplots(2,1,figsize=(12,2),gridspec_kw={'height_ratios': [1, 3]})
    ax0 = ax[0]
    ax0.set_xticks([])
    ax0.set_yticks([])
    ax0.set_xlim(xlim)
    ax0.set_ylim(ylim)
    for i in ['right','left','top','bottom']:
        ax0.spines[i].set_visible(False)
    def xrange(x,xp, y,name,c):
        ax0.fill_between(x,ylim[0],ylim[1],color=c,alpha=0.5)
        ax0.text(xp,y*0.9,name) 
    for tempx,tempx2,tempy,name,c in [[(600,800),575,np.mean(ylim),"Nucleotide",'r'],
                        [(1000,1200),950,np.mean(ylim),"Saccharides",'b'],
                        [(1200,1800),1350,np.mean(ylim),"Lipids/Proteins",'y'],
                        [(1800,2800),2200,np.mean(ylim),"Silent Zone",'k'],
                        [(2800,3100),2850,np.mean(ylim),"C-H bond",'g']
                                    ]:
        if max(tempx) <= max(xlim):
            xrange(tempx,tempx2,tempy, name,c)

    # plot mean and +/- 1std
    ax1 = ax[1]
    ax1.plot(X,Y.mean(),c='#1D7874',label='mean')
    low, high = Y.mean()-Y.std(), Y.mean()+Y.std()
    ax1.fill_between(X,low,high,color='r',alpha=0.5,label='range')
    ax1.set_xticks(np.arange(min(X),max(X),xstep))
    ax1.set_xlabel("Wavenumber (cm$^{-1}$)",fontsize='large')
    ax1.set_ylabel("Intensity (CCD)",fontsize='large')
    ax1.vlines(1172,min(low),max(high),linestyle='--',color='r')
    ax1.vlines(690,min(low),max(high),linestyle='--',color='r')
    ylim1 = ax1.get_ylim()
    ax1.text(1180,ylim1[1]*0.8,'P–O–P')
    ax1.text(700,ylim1[1]*0.8,'P=O')
    ax1.set_xlim(xlim)
    ax1.legend(title='n = {}'.format(Y.shape[0]))

    return fig, ax

def plot_mol(df):
    """
    plot the violin plot for bonds detected by Raman
    """
    import plotly.graph_objects as go
    temp = df.dropna(axis=1,how='all').fillna(0,)
    temp.reset_index(inplace=True,drop=True)

    # Order the columns by height
    temp_mean = temp.mean().sort_values(ascending=False)
    temp = temp[temp_mean.index]

    fig = go.Figure()
    for col in temp.columns:
        # remove values greater than mean + 3 std
        col_mean = temp[col].mean()
        col_std = temp[col].std()
        temp_col = temp[col][(temp[col] <= col_mean + 3*col_std)]
        fig.add_trace(go.Violin(y=temp_col, name=col.replace("$",""), box_visible=True, meanline_visible=True))
    fig.update_layout(title="Intracellular Biomolecular Distributions",
                       yaxis_title="Intensity (CCD)",
                       xaxis_title="Detected Intracellular Biomolecular",
                       showlegend = False,
                       xaxis_tickangle=45,
                       width=1000, height=600,
                       xaxis={"tickprefix": "$", "ticksuffix": "$"}
                    )
    return fig

def meta_boxplot(metadf,by,fig_dict,method="ttest"):
    """
    metadf : dataframe with metadata; n x d
    by : group the metadf by the keywords
    fig_dict : to assign figure relevant parameters
        starti: which column to start for metadf
        endi: which columns to stop
        ncols: how many columns you want for the generated figure
    method : what statistical method to compare differences between 2 data series
    draw boxplot for each feature in metadata
    """
    def add_asterisks(pval):
        if pval < 0.001:
            return '***'
        elif pval < 0.01:
            return '**'
        elif pval < 0.05:
            return '*'
        else:
            return ''
        
    # figure parameter setup
    colnames = metadf.columns
    # continuous
    if fig_dict["continuous"] == "Y":
        starti = fig_dict["starti"]
        endi = fig_dict["endi"]
        if endi == "end":
            endi = len(colnames)
        selected_colnames = colnames[range(starti,endi,1)] # all selected columns for metadf to be computed
    elif fig_dict["continuous"] == "N":
        selected_colnames = fig_dict["allcol"]
        
    ncols = fig_dict["ncols"]
    nrows = math.ceil(len(selected_colnames)/ncols)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=( 5 * ncols,5 * nrows))
    fig.subplots_adjust(hspace=0.9)

    # prepare grouped data
    if by == None:
        grouped = metadf.reset_index().groupby(by="index")
    else:
        grouped = metadf.groupby(by=by)
        
    for i, colnamei in enumerate(selected_colnames): # plot all features
        ax = axes[i//ncols,i%ncols]
        groupi = grouped[colnamei].apply(lambda x: x.reset_index(drop=True)).unstack().T 
        ax = groupi.plot(kind="box",showfliers=fig_dict["showfliers"],
                         ax=ax,boxprops={"facecolor":"b","alpha":0.4,"linewidth":2},
                                                whiskerprops={"linewidth":2},
                                                capprops={"linewidth":2},
                                                medianprops={"linewidth":4,"color":"k"},
                                                patch_artist=True)
        ax.set_title(colnamei)

        colii=groupi.columns
        ncoli = len(colii)
        stri = ""
        for i in range(ncoli):
            for j in range(ncoli):
                if j > i:
                    tempi = groupi.iloc[:,i].dropna().values
                    tempj = groupi.iloc[:,j].dropna().values
                    if method=="ttest":
                        if len(tempi) == len(tempj):
                            t_stat, p_val = scipy.stats.ttest_rel(tempi,tempj)
                        else:
                            t_stat, p_val = scipy.stats.ttest_ind(tempi, tempj, equal_var=False)
                    elif method=="tukey":
                        stati = scipy.stats.tukey_hsd(tempi,tempj)
                        t_stat = stati.statistic[0,1]
                        p_val = stati.pvalue[0,1]
                    asterisks = add_asterisks(p_val)
                    stri += f'{colii[i],colii[j]}: statistic of {method} = {t_stat:.3f}, p-value = {p_val:.3f} {asterisks}\n'
        ax.text(0,-0.7,stri, transform=ax.transAxes)
    return fig

def viz_summary(metalevel,summary,plots,labels,**kwargs):
    if "repetitive" in metalevel:
        metalevel.remove("repetitive")
    #summary.drop(droplabel,axis=1,inplace=True)
    
    obj = summary.groupby(by=metalevel)
    fig,ax = plt.subplots(1,len(plots),**kwargs)
    for i, plot_temp in enumerate(plots):
        mean = obj[plot_temp].mean()
        std = obj[plot_temp].std()
        std[np.isnan(std)] = 0
    
        mean.plot(kind='bar',yerr=std,ax=ax[i])
        ax[i].legend()
        plt.xticks(rotation=45,ha='right')
        ax[i].set_ylabel(labels[i])
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
        Y[i].columns = X # what if no numeric data in?
        if Y[i].empty:
            print('Warning: For %s, the dataFrame is empty'%i)
        else:
            Y[i].T.plot(ax=axi,legend=None)
            axi.set_title(i+", shape=%s"%(str(Y[i].shape))
                        )
            axi.set_xlabel(None)
            figcount += 1
    plt.close()
    return fig,ax