#!/usr/bin/env python3
##########################################
# This script is to plot clustering heatmap for Raman datasets based on the log-fold change
# with multiple filteration criterion
##########################################
from load import *
import Raman_preprocess

def log_fold(dict_df,order:list,filter,cus_map,rowfontcolor='k',
        topn=6,cb_pos = [0.9,0.72,0.01,0.1],
        fontsize=['small','small'],genus=None,
        base_kw = 'Exp_', T='NT', 
        method='clustermap',
        fs=(6,8),baseline_id = 4,shrink=0.5,row_colors=None,**kwargs):
        """
        input
        dict_df : dictionary [phase_taxonomy] = nxd
        base_kw : baseline for log fold
        order : a list for cluster heatmap visualization orders from t_0 to t_n
        T : whether transverse the figure or not
        topn : how many top chemicals to select; if None, we did not select
        baseline_id : period numbers. E.g., 4 time points for each sample = 4
        genus :  which genus to select
        cus_map: custom colormap
        other keywords is to beautify the visualization

        output
        clustering heatmap figure
        """
        df_l = pd.DataFrame()#save processed df: average intensity x polymers
        for key in order:
            df = dict_df[key]#df: cell x feature
            #df = df[df>=0]#mean for bacteria contain
            df = df[df>=0].fillna(0)#mean for all bacteria
            df2 = df.mean(axis=0)#df2 is series
            if key.startswith(base_kw): # e.g., 'Exp_' to remove zero value in baseline and replaced with 1 to show increase
                df2 = df2.replace(0,1)
            df_l = pd.concat([df_l,df2],axis=1)
        #df_l should now be: polymers x treatment groups or time series
        df_l.columns = [str(v) for v in order]
        df_l = Raman_preprocess.normalize(df_l,baseline_id=baseline_id,method='heatmap')
        temp_df = df_l
        df_l = np.log10(df_l).replace(-np.Inf,-1)
        #get the top n of highest absolute value of chemicals (row)

        def get_top_bond(df_l,topn=topn):
            df_l2 = df_l.abs().mean(axis=1).sort_values(ascending=False)
            idx = df_l2.index[:topn]
            df_l_top = df_l.loc[idx,:] #the top n chemicals with highest changes
            return df_l_top


        def plot_hcmap(df_l,cus_map,shrink,method,fs,T,filter,cb_pos,fontsize,genus=None):
            #filter log change
            print(filter,type(filter))
            if isinstance(filter, str):
                if filter == "pos":
                    df_l = df_l[df_l>0]
                    df_l = df_l.fillna(0)
                    #df_l = df_l.loc[:,df_l.mean(axis=0)>0]
                elif filter == "neg":
                    df_l = df_l[df_l<0]
                    df_l = df_l.fillna(0)
                    #df_l = df_l.loc[:,df_l.mean(axis=0)<0]
            else:
                df_l = df_l.loc[:, df_l.abs().mean(axis=0) > filter]
            #only select certain genus
            if genus != None:
                df_l_name = df_l.columns
                df_l_name_part = [i for i in df_l_name if genus in i]
                df_l = df_l.loc[:, df_l_name_part]
            if T == 'T':
                df_l = df_l.T #transverse or not
            elif T =='NT':
                df_l = df_l
            if method == 'heatmap':
                fig,ax = plt.subplots(1,1,figsize=fs)
                ax = sns.heatmap(df_l.T,ax=ax,\
                             cmap=cus_map,\
                             cbar_kws={"shrink": shrink,'label': 'Log10 Fold Change'},\
                             **kwargs)#linewidths=1, linecolor='black'
            if method == 'clustermap':
                try:
                    #cbar_kws:https://matplotlib.org/stable/api/figure_api.html#matplotlib.figure.Figure.colorbar
                    fig = sns.clustermap(df_l.T,\
                                        cmap=cus_map,\
                                        yticklabels=1,
                                         col_linkage=None,\
                                          cbar_kws={"shrink": shrink,'label':'Log10 Fold Change',},\
                                          cbar_pos = cb_pos,#[0.02,0.85,0.05,0.1],
                                          # col_cluster=False,
                                          figsize=fs,
                                          annot_kws={"size": fontsize[0]},
                                          tree_kws={'linewidth':1},#'colors': row_colors.iloc[:,0] tree_kws is a list of colors
                                          # row_colors=row_colors,#add one extra column for genus
                                          **kwargs,)
                    if T == 'T':
                        ticklabels = fig.ax_heatmap.axes.get_xticklabels()
                        # ticklabels2 = fig.ax_heatmap.axes.get_yticklabels()
                        # ticklabels2.set_yticklabels()
                    elif T == 'NT':
                        ticklabels = fig.ax_heatmap.axes.get_yticklabels()
                        # ticklabels2 = fig.ax_heatmap.axes.get_xticklabels()


                    #set_text_same_length
                    tick_labels = dict([(tick_label,len(str(tick_label))) for tick_label in ticklabels])
                    longest_str = max(tick_labels.values())#make other shorter string have same length
                    rc_dict = dict(zip(row_colors.index, row_colors.iloc[:, 0]))# dict[name] = color

                    #set colored tick label
                    changed_ticklabel = []
                    for tick_label in ticklabels:
                        tick_text = tick_label.get_text()#get the original text
                        changed_label =  str(tick_text) #' '*(longest_str - len(str(tick_text))) + str(tick_text)
                        changed_ticklabel.append(changed_label)
                        # tick_label.set_color(rc_dict[tick_text])
                        # tick_label.set_backgroundcolor(rc_dict[tick_text])
                        tick_label.set_fontsize(fontsize[1])
                        tick_label.set_color(rowfontcolor)
                        #print("rc_dict",rc_dict)
                        #print("tick_text",tick_text)
                        tick_label.set_bbox(dict(edgecolor='black',facecolor=rc_dict[tick_text]
                                                 ))
                        #tick_label._bbox_patch.set(width=100,height=1000)
                        #print("tick_label name",type(tick_label).__name__)
                        #print(help(tick_label.set_bbox))
                        tick_label.set_fontstyle('normal')
                        tick_label.set_ha('left')
                        # tick_label.set_text(s=changed_label)  # change the text but not the one in figure

                    if T == 'T':
                        fig.ax_heatmap.axes.set_xticklabels(changed_ticklabel)
                        fig.ax_heatmap.axes.set_yticklabels(fig.ax_heatmap.axes.get_yticklabels(),
                                                            rotation=45, fontsize=fontsize[0])

                    elif T == 'NT':
                        fig.ax_heatmap.axes.set_yticklabels(changed_ticklabel)
                        fig.ax_heatmap.axes.set_xticklabels(fig.ax_heatmap.axes.get_xticklabels(),
                                                            rotation=45,fontsize=fontsize[0],ha='right')

                except ValueError:
                    print('Value Error for clustermap')
            return fig
        if topn!=None:
            df_l_top = get_top_bond(df_l, topn)
            fig = plot_hcmap(df_l_top,cus_map, shrink, method,fs,T=T,filter=filter,
                          cb_pos=cb_pos,
                          fontsize=fontsize,genus=genus)
            return fig, df_l_top
        else:
            fig = plot_hcmap(df_l,cus_map, shrink, method,fs,T=T,filter=filter,
                          cb_pos=cb_pos,
                          fontsize=fontsize,
                          genus=genus)
            return fig, df_l
