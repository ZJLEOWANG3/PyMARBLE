# this is to show how to run the basic workflow
from load import *
import sys
import os
import math
sys.path.append("/Users/zijianleowang/Desktop/GitHub/RamanomeSpec/scripts")
sys.path.append(os.path.join(".", "../scripts/OPU/opu_analysis"))
import pandas as pd
import json
import matplotlib.pyplot as plt
import Raman_read, Raman_preprocess, Raman_find_polymer, Raman_stat, Raman_plot, Raman_color, Raman_workflow
import imp
import numpy as np
import plotly.io as pio
import opu_analysis_lib as oal # need to install git@github.com:lguangyu/scikit-feature.git
from scipy import stats 

############
## Raman process spectra
############
def process_spec(pathin,pathout,dataloop,metalevel,peakwindow,verbose=True):
    speclist = Raman_read.read_txt(pathin,typpe='dir',dataloop=dataloop,metalevel=metalevel,outdir=pathout)
    n, d = speclist.shape
    stat = pd.DataFrame()
    for i in range(n):
        specp = speclist.iloc[i,:]
        
        print("Begin to process %s"%specp["abspath"])
        pathouti = specp['abspathout']
        # save files
        if os.path.exists(pathouti): # remove exisiting files, since we need to save multiple times, to avoid inf saving
            os.remove(pathouti)
        
        X,Y,Y_BG = Raman_read.read_txt(specp['abspath'])
        if X is None and Y is None and Y_BG is None:
            print("warning: No cell retrieved")
            continue
        else:
            Y,peaks,mol,basic_stat,phenotype_spec = Raman_preprocess.std_wf(X,Y,Y_BG,peakwindow)
        
        # plot figure & save
        idx = phenotype_spec.index.tolist()
        X = phenotype_spec['wavenumber']
        idx.remove('wavenumber')
        fig, ax = Raman_plot.viz_phenotype(phenotype_spec)
        Raman_read.Raman_save(pathouti.replace(".xlsx",".png"),fig,filetype='image',verbose=verbose)
        
        # for group member
        Raman_read.Raman_save(pathouti,peaks,"peak",typpe='single',verbose=verbose)
        Raman_read.Raman_save(pathouti,mol,"molecule",typpe='single',verbose=verbose)
        Raman_read.Raman_save(pathouti,basic_stat,"stat",typpe='single',verbose=verbose)
        Raman_read.Raman_save(pathouti,phenotype_spec,typpe='series',verbose=verbose)

        stat = pd.concat([stat,basic_stat],axis=1)

        # for OPU analysis; the output data shape is n x d
        Raman_read.Raman_save(specp['abspathopu'].replace(".xlsx",".txt"),\
                            phenotype_spec['cell'],\
                            typpe='single',filetype='txt')


    stat.columns = speclist.index
    summary = pd.concat([speclist,stat.T],axis=1)
    Raman_read.Raman_save(pathout+'/summary.xlsx',summary,"summary",typpe='single')
    return summary

############
## Raman average plot
############
def plot_average(metalevel,RamanLevel,Ramanpathin,Ramanpathout,phenotypes):
    
    # init
    dirs = os.listdir(Ramanpathin)
    dirs = [dir for dir in dirs if not dir.startswith(".")]
    colorlist = Raman_color.color(len(dirs)).default_c()
    speclist = Raman_read.read_txt(Ramanpathin,typpe='dir',dataloop=len(metalevel),metalevel=metalevel,outdir=Ramanpathin+"_out")
    PATH_Raman = speclist.groupby(by=RamanLevel)["abspathout"]

    all_violin = pd.DataFrame() # to save dataframe to plot all violinplot based on treatment
    
    all_group = [pd.DataFrame()] * len(phenotypes) # layer 1 contains treatment, which contains different phenotypes
    X = pd.DataFrame()
    for group, path in PATH_Raman: # group here means the metalevel you chose, like CA, CB
        
        group_datas = [0]*len(phenotypes) # list to save different phenotype 
        for i, path_i in enumerate(path):
            # plot molecular distribution for group 
            moldata = pd.read_excel(path_i,"molecule",index_col=0) # n x d
            if i == 0: # init 
                group_mol = moldata
            else:
                # combine mol data
                group_mol = pd.concat([group_mol,moldata],axis=0)
            # plot average raman plot for specific phenotype
            for j, phenotype in enumerate(phenotypes):
                specdata = pd.read_excel(path_i,phenotype,index_col=0,header=None) # n x d
                if i == 0: # init
                    group_datas[j] = specdata 
                else:
                    # combine spec data
                    new_X, df1_new, df2_new = Raman_preprocess.binning(group_datas[j].T,specdata.T) # input and ouput is d x n
                    group_datas[j] = pd.concat([new_X, df1_new.iloc[:,1:], df2_new.iloc[:,1:]],axis=1).T
        
        ### spec data
        pdata = os.path.join(Ramanpathout,"data")
        pfig = os.path.join(Ramanpathout,"fig")
        if not os.path.exists(pdata):
            os.makedirs(pdata)
        if not os.path.exists(pfig):
            os.makedirs(pfig)
        # save spec data with different phenotypes
        for i, group_data in enumerate(group_datas):
            group_data.reset_index(drop=True,inplace=True)
            
            # # save mean
            # allspec_save_data = os.path.join(pdata,"Raman.mean.%s.%s.csv"%(phenotypes[i],group))
            # group_data.to_csv(allspec_save_data, index=False,header=False)
            # # plot&save mean of spec data
            # allspec_save_fig = os.path.join(pfig,"Raman.mean.%s.%s.png"%(phenotypes[i],group))
            # fig, ax = Raman_plot.draw_between(group_data.iloc[0,:],group_data.iloc[1:,],xlim=(380,2000))
            # plt.close()
            # fig.savefig(allspec_save_fig,bbox_inches="tight")
        
        # for each phenotype, it is concatenated into a dataframe with label of metalevel such as treatment WA,WB...
            
            labeled_group_datas = pd.concat([ group_data,
                                        pd.Series([group] * group_data.shape[0]) ], # add metalevel treatment label
                                        axis=1) 
            X = pd.concat([X,labeled_group_datas.iloc[0,:]],axis=1)
            all_group[i] = pd.concat([all_group[i],labeled_group_datas.iloc[1:,:]],axis=0) # save all labeled phenotype dataframe
        
        #### mol data
        mol_save = os.path.join(pfig,"Raman.violin.%s.csv"%group)
        group_mol.to_csv(mol_save,index=False)
        # plot&save mol data
        fig2 = Raman_plot.plot_mol(group_mol)
        pio.write_image(fig2, mol_save.replace(".csv",".png"),scale=3)

        # plot partial molecular boxplots
        temp_group_mol = group_mol.reset_index(names=RamanLevel)
        temp_group_mol[RamanLevel] = group
        all_violin = pd.concat([all_violin,temp_group_mol],axis=0)

    # plot&save partial data
    fig_dict = {"starti":1, "endi":"end", "ncols":3,"continuous":"N",
                "allcol":["O-P-O","polyP","glycogen","PHB","PHB-co-PHV"],
                "showfliers":False}
    by = RamanLevel # if none
    fig = Raman_plot.meta_boxplot(all_violin,by,fig_dict,method="tukey")
    fig.savefig(os.path.join(pfig,"Raman.violin.partmol.png"),bbox_inches="tight",dpi=300)

    # feature analysis for all_group list, list are for various phenotypes
    
    pfeature = os.path.join(Ramanpathout,"feature")
    if not os.path.exists(pfeature):
        os.makedirs(pfeature)
    X0 = X.iloc[:,[0]].T
    X0.reset_index(inplace=True,drop=True)
    for name,data in zip(phenotypes,all_group):
        data.reset_index(inplace=True,drop=True)
        data = pd.concat([X0,data],axis=0)
        data.to_csv(os.path.join(pfeature,name+".csv"),index=False)

############
## OPU analysis
############
def prep_dataset_config(metalevel,OPUlevel,OPUpathin,OPUpathout):
    # get the path and metalevel info into a dataframe for all Raman txt
    speclist = Raman_read.read_txt(OPUpathin,typpe='dir',dataloop=len(metalevel),
                                   metalevel=metalevel,
                                   outdir=os.path.split(OPUpathin)[0] # get the parent dir path for OPU
                                   ) 
    dirs = speclist.groupby(OPUlevel)["abspathopu"].apply(lambda x: "$".join(x)).reset_index()
    n, d = dirs.shape # n denotes how many grouped dirs? such as 2 time points x 3 treatments
    colorlist = Raman_color.color(n).default_c()
    
    # run iteration for config file
    dataset_config = []
    for i in range(n):
        if type(OPUlevel) == list:
            pathid = len(OPUlevel)
        elif type(OPUlevel) == str:
            pathid = 1

        file_strs = dirs.iloc[i,pathid]
        file_list = file_strs.split("$")
        file_list = [i.replace(".xlsx",".txt") for i in file_list]

        temp_dict = {}
        temp_dict["name"] = "".join( dirs[OPUlevel].iloc[i,] )
        temp_dict["color"] = colorlist[i]
        temp_dict['file'] = file_list
        dataset_config.append(temp_dict)
        
    with open(os.path.join(OPUpathout,"dataset_config.json"), "w") as f:
        json.dump(dataset_config, f)
    return dataset_config

# now we can create the OPU analyzer object
def run_opu(dataset_config,params):
	opu_analysis = oal.OPUAnalysis.from_config(dataset_config,
		# we also need extra reconcile parameters to tell the program how to
		# preprocess the data, and more importantly, if need binning to force align
		# the wavenumbers in each dataset (as they can be slightly different each
		# run)
		reconcile_param=dict(
			bin_size=params["load"]["bin_size"],  # use bin-size of 5
			wavenum_low=params["load"]["wavenum_low"],  # extract the wavenumber range >= 400
			wavenum_high=params["load"]["wavenum_high"],  # extract the wavenumber range <= 1800
			normalize=params["load"]["normalize"],  # use Euclidean (l2) normalization
		)
	)
	# run HCA 
	opu_analysis.run_hca(metric=params["HCA"]["metric"], cutoff=params["HCA"]["cutoff"],
						max_n_opus=params["HCA"]["max_n_opus"], opu_min_size=params["HCA"]["opu_min_size"])
	# save label
	opu_analysis.save_opu_labels(params["save"]["opulabel"], delimiter="\t")
	# save spectra data
	opu_analysis.save_opu_collections(params["save"]["opucollection"])
	# plot HCA
	opu_analysis.plot_opu_hca(plot_to=params["save"]["plothca"], dpi=params["save"]["dpi"])
	# plot stack barchart
	opu_analysis.plot_opu_abundance_stackbar(plot_to=params["save"]["plotstackbar"], dpi=params["save"]["dpi"])
	# PCA biplot 
	opu_analysis.plot_opu_abundance_biplot(plot_to=params["save"]["biplot"], dpi=params["save"]["dpi"], method="pca")
	# feature ranking compute, save, plot
	opu_analysis.rank_features("fisher_score")  # 'method=' can be omitted
	opu_analysis.save_opu_feature_rank_table(params["save"]["featuretable"], delimiter="\t")
	opu_analysis.plot_opu_feature_rank(plot_to=params["save"]["plotfeaturerank"], dpi=params["save"]["dpi"])



############
## draw the metadata with statistics
############

# Draw metadata
def draw_meta(fig_dict,pathin,by,pathout,sheet_name):
    metadf = pd.read_excel(pathin,sheet_name=sheet_name)
    fig = Raman_plot.meta_boxplot(metadf,by,fig_dict)
    fig.savefig(pathout,bbox_inches="tight",dpi=300)

# Draw metadata by groups
def draw_meta_grouped(pathin,sheet_name,by,pathout):
    metadf = pd.read_excel(pathin,sheet_name=sheet_name)
    # Group the dataframe by "Sample name"
    grouped = metadf.groupby(by=by)
    # Create a figure with subplots for each group
    ncols=2
    nrows=math.ceil(len(grouped)//ncols)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10 * nrows,5*ncols))

    col_iter = 5
    # Loop over the groups and create a boxplot for each one
    for i, (group_name, group_df) in enumerate(grouped):
        ax = group_df.iloc[:,col_iter:].transform(lambda x: (x - x.mean())/x.std()).boxplot(ax=axes[i//ncols, i%ncols],
                                            boxprops={"facecolor":"b","alpha":0.4,"linewidth":2},
                                            whiskerprops={"linewidth":2},
                                            capprops={"linewidth":2},
                                            medianprops={"linewidth":4,"color":"k"},
                                            patch_artist=True)
        ax.set_title(group_name)
        if i %ncols == 0:
            ax.set_ylabel('Standardized value')
        if i//ncols == nrows-1:
            ax.set_xlabel('Soil chemical profiles')
            # set x-label rotation and horizontal alignment
            labels = ax.get_xticklabels()
            plt.setp(labels, rotation=45, ha='right')
        else:
            ax.set_xticklabels([])
    # Adjust the layout and spacing of the subplots
    fig.tight_layout()
    plt.subplots_adjust(hspace=0.2)
    fig.savefig(pathout,bbox_inches="tight",dpi=300)