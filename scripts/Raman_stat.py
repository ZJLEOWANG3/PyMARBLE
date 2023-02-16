#!/usr/bin/env python3
#############################################
# This script is to perform statistical analysis for Raman 
# different discussions
#############################################
from load import *

#############################################
# This script is to perform Tukey analysis for ML model data
# boxplot_tukey(acc,(0.96,1.03),title='classification',label='accuracy') as an example
# acc is a dataframe (m x n) : m denotes samples of each n
#############################################
def get_stat(totaln, phenotype,phenotypename):
    """
    Get basic statistics like cell recovery %, cell number, PAO%, GAO%, PHAAO%
    """
    celln,_ = phenotype.shape # 
    phenotypestat = phenotype.mean(axis=0) * 100
    phenotypestat["cell"] = celln
    recovery = celln/totaln * 100
    v = [recovery] + phenotypestat.values.tolist()
    basic_stat = pd.Series(v,index = ['cell recovery']+phenotypename)
    return basic_stat

def boxplot(df,title,ylabel,figsize=(5,5)):
        fig = plt.figure(figsize=figsize)
        def rearrange_df(df):
            #large to small

            id = df.max(axis=0).sort_values(ascending=False).index.to_list()
            df = df[id]
            return df
        df = rearrange_df(df)
        df.boxplot(rot=45)
        plt.title(title,fontsize='x-large')
        plt.ylabel(ylabel,fontsize='x-large')
        loc,_ = plt.xticks()
        labels = df.columns.to_list()
        dict_pos = {}
        for i in range(len(labels)):
            #key: polymer, xloc, ymaxloc
            dict_pos[labels[i]] = (loc[i],df[labels[i]].max())
        return fig,dict_pos

def Tukey(df,xlabel='Raman CCD intensity (a.u.)',figsize=(5,5)):
        def get_tukey_df(df):
            name = df.columns
            df_value, df_treatment = pd.DataFrame(), pd.DataFrame()
            for i in range(df.shape[1]):
                df_value = pd.concat([df_value, df.iloc[:, i]])
                temp_treatment = np.repeat(name[i], repeats=df.shape[0])
                df_treatment = pd.concat([df_treatment, pd.DataFrame(temp_treatment)])
            return df_value, df_treatment
        df_value, df_treatment = get_tukey_df(df)
        data = MultiComparison(df_value, df_treatment)
        results = data.tukeyhsd()
        results_df = table_to_df(results)
        fig, ax = plt.subplots(1,1,figsize=figsize)
        fig = results.plot_simultaneous(ax=ax,figsize=figsize)

        #this is supposed to revise the code with the aim of getting data from fig,ax
        # ax = fig.get_axes()[0]
        # print(ax.get_children()[1].get_xdata())#ax.get_children()
        # print(ax.get_children()[1].get_ydata())
        # print(ax.get_children())
        # print(dir(ax.get_children()[2]))
        # print(ax.get_children()[2].get_ec())
        # # print(dir(ax.get_children()[2]))
        # # print(ax.get_children()[2].get_bounds())
        #
        # # print(ax.get_children())

        plt.xlabel(xlabel, fontsize='x-large')
        return fig[0],results
def table_to_df(results):
        summary = results.summary().data
        df = pd.DataFrame()
        for data in summary[1:]:
            df = pd.concat([df,pd.Series(data)],axis=1)
        df = df.T
        df.columns = summary[0]
        table_df = df
        return table_df

def set_significance(fig,dict_pos,results,ylim):
        table_df = table_to_df(results)
        table_df2 = table_df[table_df['reject']==True]#find the significantly different pair
        g1,g2,pv = table_df2['group1'],table_df2['group2'],table_df2['p-adj']

        for i in range(table_df2.shape[0]):
            polymer1,polymer2 = g1.iloc[i],g2.iloc[i]
            print(polymer1,polymer2)
            x1,y1 = dict_pos[polymer1][0],dict_pos[polymer1][1]
            x2, y2 = dict_pos[polymer2][0],dict_pos[polymer2][1]
            gcf = plt.gcf()
            xmax = max(x1,x2)
            xmin = min(x1,x2)
            xmedian = (xmin+xmax)/2
            ymax = max(y1,y2)
            ax = fig.add_subplot(1,1,1)
            h = ymax*0.02#*(x1+x2)
            ax.plot([xmin,xmin,xmax,xmax],[ymax+0.2*h,ymax+h,ymax+h,ymax+0.2*h])
            if pv.iloc[i]<=0.05:
                ax.text(xmedian,ymax+1.1*h,'*')
            if pv.iloc[i]<=0.005:
                ax.text(xmedian,ymax+1.1*h,'**')
            ax.set_ylim(ylim)
        return fig

def boxplot_tukey(df2,ylim,xtick_label=None,title='boxplot',label='Raman CCD intensity (a.u.)',method='df'):
        #df2 = Raman.Raman_find_polymer.find_multiple_polymer(df.iloc[:,0],df.iloc[:,1:])
        #fig1,fig2 = Raman.Raman_analysis.boxplot_tukey(df2)
        if method=='df':#columns name : polymers (multi-poly) or treatments (single polymer);
            ## each col contains various intensity of polymers
            fig0, dict_pos = boxplot(df2, title,label)
            fig1,results = Tukey(df2,xlabel=label)
            fig2 = set_significance(fig0, dict_pos, results,ylim)
        return fig1,fig2


########################
## This section is to achieve ANOSIM
########################
def Raman_anosim(feature,label):
    # feautre : n x d
    # label : groupings of the single sample
    dm = pdist(feature)
    dm = squareform(dm)
    IDs = np.arange(dm.shape[0])
    dm2 = DistanceMatrix(dm,IDs)
    grouping = label
    # make output deterministic; not necessary for normal use
    np.random.seed(0)
    verbose = anosim(dm2, grouping, permutations=99)
    return verbose

########################
## This section is to compute Jaccard and consistency
########################
def jaccard(a,b):
    """
    compute the jaccard similarity for two clustering results
    J = intersect(A,B)/union(A,B)
    a : 1D array
    b : 1D array
    """
    a_mat = a.reshape(1, -1) == a.reshape(-1, 1)
    b_mat = b.reshape(1, -1) == b.reshape(-1, 1)
    jaccard = (a_mat & b_mat).sum() / (a_mat | b_mat).sum()
    return jaccard

def consistency(a,b):
        return np.equal(a,b).sum()/a.shape[0]
