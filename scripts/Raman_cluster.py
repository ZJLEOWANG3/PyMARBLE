#!/usr/bin/env python3
#############################################
# This script is to automatically pick optimal clusters based on various metrics
# Visualization includes the PCA biplot and selected wavenumber with average+std+min+max spectra
#############################################

from load import *
import Raman_molecule
import Raman_stat
#############################################
# This section is to customize PCA and Kmeans
#############################################
def customized_PCA(Y,verbose, k=2):
    """
    Customized PCA algorithm to obtain original feature with highest singular value results
    Input:
    Y - DataFrame, Feature with shape n x d
    k - Select first k features
    
    Output:
    transformed data: new_pt
    Sigma matrix w/ singular value: Sigma
    explained variance ratio: explained_variance_ratio_
    PC_dict: the sorting of original features according to its importance
    vh: to show the rotation of original matrix
    """
    
    np.random.seed(0)
    # Change to numpy
    Y = Y.to_numpy()
    n_sample, n_feature = Y.shape
    # Calculate row mean
    Ymean = Y.mean(axis=0)
    
    # Decentering
    YDec = Y - Ymean
    
    # SVD for decentered feature matrix
    ## max \phi \hat(X) \hat(X)^T \phi^T    
    ## SVD, \hat{X} = U Sigma V^H, where U \in R^{nxn}, Sigma \in R^{nxd}, VH \in R^{dxd}
    ## But the svd algorithm can only return s with dimension of min{n,d}, we need to recreate the true Sigma
    u, s, vh = np.linalg.svd(YDec,full_matrices=True)#,full_matrices=False) # nxn, n, dxd, the second n is because it is at most n for features due to rank
    
    # Get variance explained by singular values
    explained_variance_ = (s ** 2) / (n_sample - 1)
    total_var = explained_variance_.sum()
    explained_variance_ratio_ = explained_variance_ / total_var
    singular_values_ = s.copy()  # Store the singular values.
    
    # create m x n Sigma matrix
    Sigma = np.zeros((n_sample, n_feature))
    # populate Sigma with n x n diagonal matrix
    Sigma[:s.shape[0], :s.shape[0]] = np.diag(s)
    
    # flip eigenvectors' sign to enforce deterministic output
#     u, vh = svd_flip(u,vh)
    ## Check SVD
    #str_checkSVD = "SVD Test Pass: {} \n".format(np.allclose((u@Sigma)@vh,YDec,equal_nan=True))
    #verbose += str_checkSVD
    verbose["SVD Test Pass"] = np.allclose((u@Sigma)@vh,YDec,equal_nan=True)
    
    ## Transform original matrix X by U
    new_pt = u@Sigma
    ## Check transformed matrix with Python PCA
    pca = sklearn.decomposition.PCA()#n components, if not set, will keep all
    principalComponents = pca.fit_transform(Y)#n x d', d' is reduced dimension
    ## Check PCA
    #str_checkPCA = "PCA Test Pass: {} \n".format( np.allclose(np.absolute(new_pt[:principalComponents.shape[0],:principalComponents.shape[1]])\
    #                                    ,np.absolute(principalComponents)) )
    verbose["PCA Test Pass"] = np.allclose(np.absolute(new_pt[:principalComponents.shape[0],:principalComponents.shape[1]])\
                                        ,np.absolute(principalComponents))
    #verbose += str_checkPCA    
    ## Find the name of selected feature
    
    important_index = [np.absolute(vh[i]).argmax() for i in range(vh.shape[0])]
    feature_name = list(range(Y.shape[1]))
    PC_dict = {"PC{}".format(i): feature_name[important_index[i]] for i in range(vh.shape[0])}
    
    #return the transformed matrix and Sigma and variance ratio
    return new_pt, Sigma, explained_variance_ratio_, PC_dict, vh, verbose

def customized_Kmeans(Y, kmax=10):
    """
    To Find Kmeans using Silhouette metric
    Input
    Y: data nxd
    kmax: At most kmax clusters, default 10
    
    Output
    cluster_n: optimal cluster number
    labels: n labels of the clusters
    """
    #Find the optimal cluster number
    sil = []
    for k in range(1,kmax+1):
        kmeans = KMeans(n_clusters = k).fit(Y)
        labels = kmeans.labels_
        if k==1:
            sil.append(0)
        else:
            sil.append(silhouette_score(Y,labels,metric='euclidean'))
    index = np.array(sil).argmax()
    cluster_n = range(1,kmax+1)[index]#the optimal clustering number
    
    # Calculate kmeans with initialization of k-means++ for rapid convergence
    np.random.seed(0)
    kmeans = KMeans(n_clusters = cluster_n,init='k-means++',\
                   n_init=10,max_iter=300).fit(Y)
    labels = kmeans.labels_
    return cluster_n, labels


def PCA_biplot(ax,PC_dict,vh,X,mol_df,narrow=2,labels=None,dm='pca'):
    """
    plot the arrow of the PCA rotation, aka VH in SVD 
    input:
    ax: fig ax
    coeff: vh
    X: wavenumber 
    narrow : number of arrows for feature direction
    
    output:
    Selected_Wavelength: list with pairs of idfeature and corresponding wavenumber
    """
    Selected_Wavelength = []
    if dm == 'pca' or isinstance(dm,int):
        vh = vh.T
        scaler = 2
        for i in range(narrow):
            idfeature = PC_dict["PC{}".format(i)]#which feature is selcted 
            ax.arrow(0, 0, vh[idfeature,0]* scaler, vh[idfeature,1]* scaler,color = 'k',alpha = 1)
            if labels is None:
                wavenumber = X[idfeature]
                Selected_Wavelength.append(("PC%i"%i,\
                        idfeature,wavenumber,\
                        Raman_molecule.get_mol(mol_df,wavenumber)))
                ax.text(vh[idfeature,0]* scaler, vh[idfeature,1] * scaler,
                        "{}".format(round(X[idfeature],1)) + " $cm^{-1}$", #"{}-th Feature={}".format(idfeature,X[idfeature]),
                        color = 'g', ha = 'center', va = 'center')
            else:
                ax.text(vh[idfeature,0]* scaler, vh[idfeature,1] * scaler, labels[i], color = 'g', ha = 'center', va = 'center')
    elif dm == 'ica':
        for i in range(narrow):
            idfeature = PC_dict["ICA{}".format(i)]#which feature is selcted 
            if labels is None:
                wavenumber = X.loc[int(idfeature)]
                Selected_Wavelength.append(("ICA%i"%i,\
                        idfeature,wavenumber,\
                        Raman_molecule.get_mol(mol_df,wavenumber)))
    return Selected_Wavelength

def plotPK(Y,dm,labels, EVR,PC_dict, vh,X,narrow,clist, fs,title, mol_df,ax=None):
    """
    Plot the customized PCA & Kmeans
    Input
    Y: data nxd
    labels: n labels
    EVR: explained variance ratio for d features
    PC_dict: the sorting of original features according to its importance
    vh: vh of SVD
    X: wavenumber feature
    clist : color list
    narrow : number of feature arrows to show feature direction
    title : figure title
    
    ouput 
    Selected_Wavelength : selected wavenumber by PCA
    data : data for figure
    """
    if EVR is not None:
        xlabel = "{}1 {:.2f}%".format(dm,EVR[0]*100)
        ylabel = "{}2 {:.2f}%".format(dm,EVR[1]*100)
    else:
        xlabel = "{}1 {}".format(dm,PC_dict["ICA0"]) + " c$m^{-1}$"
        ylabel = "{}2 {}".format(dm,PC_dict["ICA1"]) + " c$m^{-1}$"
    xs = Y[:,0]
    ys = Y[:,1]
    scalex = 1.0/(xs.max() - xs.min())
    scaley = 1.0/(ys.max() - ys.min())
    xs *= scalex
    ys *= scaley
    if ax==None:
        fig,ax = plt.subplots(1,1,figsize=(3,3))
    #create cs : color for each scatter with label
    cs = [clist[i] for i in labels]
    data = pd.concat([pd.Series(xs),pd.Series(ys),pd.Series(cs)],axis=1)
    data.columns = ["xaxis","yaxis","color"]
    ax.scatter(xs,ys,c=cs,linewidths=1,edgecolors='k')
    # wrapped text for title
    wrapped_title = "\n".join(textwrap.wrap(title,40))
    ax.set_title(wrapped_title,fontsize='medium')
    # assign label
    for label in set(labels):
        tempxs = xs[labels==label][0]
        tempys = ys[labels==label][0]
        ax.scatter(tempxs,tempys,c=clist[label],linewidths=1,edgecolors='k',label='cluster {}'.format(label))

    ax.set_xlabel(xlabel,fontsize=fs)
    ax.set_ylabel(ylabel,fontsize=fs)
    #legend
    ax.legend()
    #plot biplot
    Selected_Wavelength = PCA_biplot(ax,PC_dict, vh,X,mol_df,narrow,dm=dm)
         
    if ax == None:
        return Selected_Wavelength, data, fig, ax
    else:
        return Selected_Wavelength, data

#############################################
# This section is to get values of various metrics
#############################################

def get_pred(k, X, clr,aff='nearest_neighbors'):
    if clr=='kmeans':
        km = KMeans(n_clusters=k, init='k-means++',\
                   n_init=10,max_iter=300,random_state=37)
        km.fit(X)
        y_pred = km.predict(X)
    elif clr == 'Spectral':
        sc = SpectralClustering(n_clusters=k,random_state=37,affinity=aff)
        print(X)
        y_pred = sc.fit_predict(X)
    return y_pred

def get_bic_aic(k, X):
    gmm = GaussianMixture(n_components=k, init_params='kmeans')
    gmm.fit(X)
    return gmm.bic(X), gmm.aic(X)

def get_score(k, X, clr = 'kmeans'):#, y):
    """
    Input 
    cluster number k
    feature matrix X, nxd
    clustering method clr
    """
    y_pred = get_pred(k, X, clr)
    bic, aic = get_bic_aic(k, X)
    sil = silhouette_score(X, y_pred)
    db = davies_bouldin_score(X, y_pred)
#     hom = homogeneity_score(y, y_pred)
#     com = completeness_score(y, y_pred)
#     vms = v_measure_score(y, y_pred)
    cal = calinski_harabasz_score(X, y_pred)
    return k, bic, aic, db, sil, cal#, hom, com, vms, 


def get_opt(data,clr,verbose,aff='nearest_neighbors'):
    """
    data : nxd
    clr : classifier method name
    verbose : to save the results
    """
    df = pd.DataFrame([get_score(k, data,clr) for k in range(2, 11)],
                  columns=['k', 'BIC', 'AIC', 
                           'davies','silhouette','calinski'])
    opt_n_dict = optimal_metric(df)
    verbose["Methods"] = "PCA + {}".format(clr)
    verbose["Metric results"] = opt_n_dict
    #verbose += "Methods: PCA + {} \n".format(clr)
    #verbose += "Metric results: {} \n".format(opt_n_dict)
    # pick the clusters with largest numbers of metric
    opt_n = max(opt_n_dict,key = lambda k:len(opt_n_dict[k]))
    verbose["Optimal cluster number"] = opt_n
    #verbose += "Optimal cluster number: {} \n".format(opt_n)
    labels = get_pred(opt_n, data, clr,aff)
    return opt_n_dict, opt_n, labels, verbose

# plot the EVR
def plot_EVR(df,title):
    fig, ax = plt.subplots(1,1,figsize=(3,3))
    ax.plot(df["cumEVR"],df["jaccard"],label="Jaccard")
    ax.plot(df["cumEVR"],df["consistency"],label="Consistency")
    ax.vlines(0.8,0,1.1,linestyles='dashed')
    ax.text(0.65,1.05,"80%")
    ax.set_ylabel("Similarity",fontsize='large')
    ax.set_xlabel("EVR (%)",fontsize='large')
    ax.set_title(title,fontsize='medium')

    ax2 = ax.twinx()
    ax2.plot(df["cumEVR"],df["PCAcomps"],'-ko',markerfacecolor='white',markeredgecolor='k')
    ax2.set_ylabel("PCA Components")
    ax.legend(loc='center left')
    return fig

# Get optimal cluster number based on optimal metric
def optimal_metric(df):
    """
    Input dataframe with columns of cluster, bic, aic, silhouette, davies, calinski
    the Higher Silhouette and Calinski, The better
    others are lower the better
    
    Output: Optimal dictionary key=k, value=(metric,value)
    """
    def pick_hl(s,col_k,hl):
        """
        Input dataframe series, Series of cluster k, hl: maximum or minimum
        
        Output optimal k and maximum metric value
        """
        if hl=='h':
            idm = s.idxmax()
            extrema = s.max()
        elif hl=='l':
            idm = s.idxmin()
            extrema = s.min()
        try:
            optk = col_k[idm]
            return optk, extrema
        except UnboundLocalError:
            print("UnboundLocalError")
        
    
    dict_df = {}
    for metric in df.columns[1:]:
        if metric == 'silhouette' or metric== 'calinski':
            hl='h'
        else:
            hl='l'
        optk, extrema = pick_hl(df[metric],df['k'],hl)
        optk = int(optk)
        if optk not in dict_df.keys():
            dict_df[optk] = (metric,extrema)
        else:
            dict_df[optk] += (metric,extrema)
    return dict_df


######################
# For the label switch
######################
#hard vote
"""
def switch_cluster(a,b):
    # a is reference
    # b is comparison
    # for binary
    b1 = np.where((b==0)|(b==1), b^1, b)
    b_score = (a==b).sum()
    b1_score = (a==b1).sum()
    if b1_score > b_score:
        return a, b1
    else:
        return a, b
"""

def switch(a,b):
    # switch
    # since the ref is fixed, so we cannot use get_ref here
    # a as ref, b as comp
    # final goal is to ouput the b with changed labels
    
    # get the reference and comparison
#     def get_ref(a,b):
#         na = np.unique(a)
#         nb = np.unique(b)

#         # define ref & comp
#         if na.shape[0] > nb.shape[0]:
#             ref = np.copy(a)
#             comp = np.copy(b)
#         else:
#             ref = np.copy(b)
#             comp = np.copy(a)
#         return ref, comp

    # get consistency
    #def compute_consistency(a,b):
    #    return round(np.equal(a,b).sum()/a.shape[0],5)
    
#     a,b = get_ref(a,b) # a as ref, b as comp
    uniqueb = np.unique(b)
    # permute labels of b
    n = uniqueb.shape[0]
    perm = set(permutations(np.arange(n)))
    
    # get index for each unique in b
    idb = []
    for bi in uniqueb: # bi - i-th unique label in b
        idb.append([b == bi])

    # compute consistency for all perm
    # initialize new b
    b_new = []
    b_temp = np.ones(b.shape,dtype=np.int8)
    consistency = []
    
    for  labels in perm:
        for (i, idx) in enumerate(idb):
            # generate the complete temp b
            b_temp[tuple(idx)[0]] = labels[i]
        #b_temp[idx] = labels[i]
        temp_consistency = round(Raman_stat.consistency(a,b_temp),5)
        b_new.append(b_temp.copy())
        consistency.append(temp_consistency)
        
    consistency = np.asarray(consistency)
    b_new = np.asarray(b_new)
    max_con = np.max(consistency)
    idx = np.argwhere(consistency==max_con)
    b_new = b_new[idx]  
    return b_new,max_con
#############################################
# This section is to visulize 3 figures and generate related data
#############################################
# modified viz_clr
def viz_clr(X,Y,title,mol_df, clr='kmeans',cm=None, pca='pca',aff='nearest_neighbors',fs='x-large'):
    """
    visualize the clusters with defined optimal n after PCA
    Input 
    X: Wavenumbers for Raman
    Y: feature matrix (n x d)
    title : figure title for plotPK
    clr: clustering methods
    cm : cluster match; if None, do not match
    pca: use PCA for clustering, pca:use PCA or F-False, if int I then pick the fist I PCA component; T-ica: use ICA to compute
    aff: affinity for Spectral clustering. rbf by default, or nearest_neighbors
    
    Output:
    verbose : documents of script verbose, previously is a string but currently revised to be a dictionary
    opt_n_dict: optimal cluster number and related metric value
    Selected_Wavelength: selected wavelength by PCA
    """
    # plot spectra
    def plot_spectra(ax,X,Y,clist=['k','r','b'],YON='Y'):
        """
        input
        ax : plot in this ax
        X : Raman wavenumber
        Y : Raman features a.u. (n x d)
        clist : color list for plot
        YON : Y to set x,y label
        
        ouput
        data : dataframe for this plotting
        """
        fs = 'x-large'
        Ymean = Y.mean(axis=0)
        Ystd = Y.std(axis=0)
        Ymax = Y.loc[Y.apply(np.linalg.norm, axis=1).idxmax(),:]#the spectra with max norm
        Ymin = Y.loc[Y.apply(np.linalg.norm, axis=1).idxmin(),:]
        Ylow = Ymean-Ystd
        Yhigh = Ymean+Ystd
        data = pd.concat([X,Ymean,Ystd,Ymax,Ymin,Ylow,Yhigh],axis=1,ignore_index=True)
        data.columns = ["X","Ymean","Ystd","Ymax","Ymin","Ylow1sig","Yhigh1sig"]
        ax.plot(X,Ymean,c=clist[0],lw=0.5)
        ax.fill_between(X,Ylow,Yhigh,facecolor=clist[1],alpha=0.3)
        ax.plot(X,Ymin,c=clist[2],lw=0.1)
        ax.plot(X,Ymax,c=clist[2],lw=0.1)

        if YON=='Y':
            ax.set_xlabel("Wavenumber ($cm^{-1}$)",fontsize=fs)
            ax.set_ylabel("Intensity (a.u.)",fontsize=fs)
        if YON=='N':
            ax.axes.yaxis.set_ticklabels([])
        return data
    # plot selected spectra based on labels
    def plot_sep(X,Y,labels,clist,Yrange):
        """
        input 
        X : wavenumber·
        Y : Raman feature (n x d)
        labels : clustering labels
        clist : color list for labels
        Yrange : consistency of Ylim
        
        output
        fig : figure
        data1 : PCA cluster figure data, only save the X, Ymean, Ystd, Ymax, Ymin, Ylow1sig, Yhigh1sig; \
                read_csv(index_col=0)
        data2 : total spectra data, same format as data1
        data_dict : data dictionary of the figure, [label] = (X,Y)
        labels : the labels for clusters
        """
        label_uniq = set(labels)
        n = len(label_uniq) # number of clusters
        fig,ax = plt.subplots(1,n,figsize=(n*3,3))
        i = 0 # count number, for labels
        data_dict = {}
        ratio = {}
        for label in label_uniq:
            idx = np.argwhere(labels==label).reshape(-1)
            c = clist[int(label)]
            if i==0:
                data = plot_spectra(ax[int(label.tolist())],X,Y.iloc[idx,:],['k',c,c],'N1')
            elif i >=1:
                data = plot_spectra(ax[int(label)],X,Y.iloc[idx,:],['k',c,c],'N')
            #data_dict[label] = [data,Y.iloc[idx,:].T]
            data_dict[label] =pd.concat([data, Y.iloc[idx,:].T],axis=1)
            ax[int(label)].set_ylim(Yrange) # same Y range
            ax[int(label)].set_title('cluster {}, n = {}'.format(label,len(idx)))
            ratio[str(label)] = data_dict[label].shape[1]
            i += 1
        ax[0].set_ylabel( "Intensity (a.u.)",fontsize=fs)
        fig.text(0.5,-0.02,"Wavenumber ($cm^{-1}$)",fontsize=fs,ha='center')
        return fig, data_dict, ratio#,mol_df,war

    def sort_labels(Y,labels):
        # to make labels in order of max Y values from small to larger 
        Y = pd.DataFrame(Y)
        label_uniq = set(labels)
        data = []
        for label in label_uniq:
            idx = np.argwhere(labels==label).reshape(-1)
            temp_Y = Y.iloc[idx,:].mean().max()
            data.append([label,temp_Y,idx])
        data = pd.DataFrame(data)
        data = data.sort_values(by=1)
        data.iloc[:,0] = np.arange(len(label_uniq))
        labels2 = np.zeros(labels.shape,dtype=np.int32)
        for row in range(data.shape[0]):
            labels2[data.iloc[row,2]] = data.iloc[row,0]
        return labels2

    verbose = {}
     

    # pick PCA components
    if isinstance(pca,int):
        new_pt, Sigma, EVR, PC_dict,vh,verbose = customized_PCA(Y,verbose)# nxd,
        data = new_pt[:,:int]
    elif pca=='F':
        data = Y
    elif pca=='pca':
        new_pt, Sigma, EVR, PC_dict,vh,verbose = customized_PCA(Y,verbose)# nxd,
        data = new_pt
    elif pca=="ica":
        Y.columns = Y.columns.astype(str)
        clf = FastICA(random_state=38,whiten='unit-variance')
        data = clf.fit_transform(Y)
        name = clf.get_feature_names_out()
        Sigma, vh, EVR = None, None, None
        PC_dict = {}
        for i in range(data.shape[1]):
            PC_dict["ICA{}".format(i)] = name[i].split("fastica")[1]
        
        
    # pick optimal cluster number
    opt_n_dict, opt_n, labels, verbose = get_opt(data,clr,verbose)
    # sort labels by Y
    labels = sort_labels(data,labels)
    """
    # original switch for binary
    if isinstance(cm,np.ndarray) and np.unique(labels).shape[0]==2:
        # we switch the labels
        _, labels = switch_cluster(cm, labels)
    """
    def plot_main(cm, labels, X, Y, pca, verbose):
        # switch labels
        if isinstance(cm,np.ndarray):
            labels, cons = switch(cm,labels)
            if len(labels.shape)==3:
                labels = labels[0][0]
        clist = Raman_color.color(10).default_c()# consistency of plotting color
        
        # spectra plot
        fig1, ax1 = plt.subplots(1,1,figsize=(3,3))
        fig2, ax2 = plt.subplots(1,1,figsize=(3,3))
        data2 = plot_spectra(ax2,X,Y)# data for fig2
        ax2.set_title("total spectra, n={}".format(Y.shape[0]))
        Yrange = [data2.iloc[:,1:].min().min(),data2.iloc[:,1:].max().max()]# Range for figure Y axis
        fig3, data_dict,ratio = plot_sep(X,Y,labels,clist,Yrange)#data for fig3
        verbose["ratio"] = ratio

        if pca != 'F':
            narrow = 2
            # plot clusters after PCA/ICA
            Selected_Wavelength, data1 = plotPK(data,pca,labels,EVR,PC_dict,vh,X,narrow,clist, fs,title,mol_df, ax1)
            if isinstance(pca,int):
                name = 'PCA'
            else:
                name = pca.upper() # PCA or ICA
            verbose["Selected wavenumber by {}".format(name)] = Selected_Wavelength
            # add features to spectra plot
            i = 0 # for the alignment of y axis of text
            d_text = 0 # for the distance of two text if there is a wrap
            for (PCA_i,idf,v,mol) in Selected_Wavelength:
                # idf : the i-th feature in Raman
                # v : the selected wavenumber of Raman
                # PCA_i : the i-th PCA components
                if mol==None:
                    mol = "unknown"
                ax2.vlines(v,Yrange[0],Yrange[1],linestyles='dashed',lw=1,color='b')
                # plot which PCA : which wavenumber and molecular name
                wrap_thre = 20
                if len(mol)>=wrap_thre:
                    wrap_mol = "\n".join(textwrap.wrap(mol,wrap_thre)) 
                    d_text += 0.1
                else:
                    wrap_mol = mol
                ax2.text(0.03,0.9-i*0.05-d_text,\
                        "{}:{}".format(PCA_i, str(round(v,1))) + "$\mathregular{cm^{-1}}$,"+wrap_mol,\
                        transform=ax2.transAxes,fontsize='x-small')
                i += 1

        # Plot predicted UMAP & T-SNE
        fig4,ax4 = plt.subplots(1,2,figsize=(6,3))#,gridspec_kw={'width_ratios': [1.5, 2.8]})
        plt.subplots_adjust(wspace=0.35)
        for (method,i) in zip(["umap",'tsne'],np.arange(2)):
            pred_clf(Y,method,data1.iloc[:,2].to_numpy(),ax4[i],clist) 
        fig4.suptitle("Predicted Dimension Reduction Clustering \n {}".format(title),fontsize='medium')
        return data1, data2, data_dict, verbose, fig1, fig2, fig3, fig4
    # plot
    data1, data2, data_dict, verbose, fig1, fig2, fig3, fig4 = plot_main(cm, labels, X, Y, pca, verbose)
    return opt_n_dict, verbose, [fig1,fig2,fig3,fig4], [data1,data2,data_dict, labels]
"""
        if pca=='T':
            narrow = 2
            # plot clusters
            Selected_Wavelength, data1 = plotPK(new_pt,labels,EVR,PC_dict,vh,X,narrow,clist, fs,title,mol_df, ax1)
            verbose["Selected wavenumber by PCA"] = Selected_Wavelength

            # Plot predicted UMAP & T-SNE
            fig4,ax4 = plt.subplots(1,2,figsize=(6,3))#,gridspec_kw={'width_ratios': [1.5, 2.8]})
            plt.subplots_adjust(wspace=0.35)
            for (method,i) in zip(["umap",'tsne'],np.arange(2)):
                pred_clf(Y,method,data1.iloc[:,2].to_numpy(),ax4[i],clist) 
            fig4.suptitle("Predicted Dimension Reduction Clustering \n {}".format(title),fontsize='medium')

            # add features to spectra plot
            i = 0 # for the alignment of y axis of text
            d_text = 0 # for the distance of two text if there is a wrap
            for (PCA_i,idf,v,mol) in Selected_Wavelength:
                # idf : the i-th feature in Raman
                # v : the selected wavenumber of Raman
                # PCA_i : the i-th PCA components
                if mol==None:
                    mol = "unknown"
                ax2.vlines(v,Yrange[0],Yrange[1],linestyles='dashed',lw=1,color='b')
                # plot which PCA : which wavenumber and molecular name
                wrap_thre = 20
                if len(mol)>=wrap_thre:
                    wrap_mol = "\n".join(textwrap.wrap(mol,wrap_thre)) 
                    d_text += 0.1
                else:
                    wrap_mol = mol
                ax2.text(0.03,0.9-i*0.05-d_text,\
                        "{}:{}".format(PCA_i, str(round(v,1))) + "$\mathregular{cm^{-1}}$,"+wrap_mol,\
                        transform=ax2.transAxes,fontsize='x-small')
                i += 1
"""
        
    

'''
def viz_clr(X,Y,title,mol_df, clr='kmeans',cm=None, pca='T',aff='nearest_neighbors',fs='x-large'):
    """
    visualize the clusters with defined optimal n after PCA
    Input 
    X: Wavenumbers for Raman
    Y: feature matrix (n x d)
    title : figure title for plotPK
    clr: clustering methods
    cm : cluster match; if None, do not match
    pca: use PCA for clustering, T-True or F-False, if int I then pick the fist I PCA component
    aff: affinity for Spectral clustering. rbf by default, or nearest_neighbors
    
    Output:
    verbose : documents of script verbose, previously is a string but currently revised to be a dictionary
    opt_n_dict: optimal cluster number and related metric value
    Selected_Wavelength: selected wavelength by PCA
    """
    # plot spectra
    def plot_spectra(ax,X,Y,clist=['k','r','b'],YON='Y'):
        """
        input
        ax : plot in this ax
        X : Raman wavenumber
        Y : Raman features a.u. (n x d)
        clist : color list for plot
        YON : Y to set x,y label
        
        ouput
        data : dataframe for this plotting
        """
        fs = 'x-large'
        Ymean = Y.mean(axis=0)
        Ystd = Y.std(axis=0)
        Ymax = Y.loc[Y.apply(np.linalg.norm, axis=1).idxmax(),:]#the spectra with max norm
        Ymin = Y.loc[Y.apply(np.linalg.norm, axis=1).idxmin(),:]
        Ylow = Ymean-Ystd
        Yhigh = Ymean+Ystd
        data = pd.concat([X,Ymean,Ystd,Ymax,Ymin,Ylow,Yhigh],axis=1,ignore_index=True)
        data.columns = ["X","Ymean","Ystd","Ymax","Ymin","Ylow1sig","Yhigh1sig"]
        ax.plot(X,Ymean,c=clist[0],lw=0.5)
        ax.fill_between(X,Ylow,Yhigh,facecolor=clist[1],alpha=0.3)
        ax.plot(X,Ymin,c=clist[2],lw=0.1)
        ax.plot(X,Ymax,c=clist[2],lw=0.1)

        if YON=='Y':
            ax.set_xlabel("Wavenumber ($cm^{-1}$)",fontsize=fs)
            ax.set_ylabel("Intensity (a.u.)",fontsize=fs)
        if YON=='N':
            ax.axes.yaxis.set_ticklabels([])
        return data
    # plot selected spectra based on labels
    def plot_sep(X,Y,labels,clist,Yrange):
        """
        input 
        X : wavenumber·
        Y : Raman feature (n x d)
        labels : clustering labels
        clist : color list for labels
        Yrange : consistency of Ylim
        
        output
        fig : figure
        data1 : PCA cluster figure data, only save the X, Ymean, Ystd, Ymax, Ymin, Ylow1sig, Yhigh1sig; \
                read_csv(index_col=0)
        data2 : total spectra data, same format as data1
        data_dict : data dictionary of the figure, [label] = (X,Y)
        labels : the labels for clusters
        """
        label_uniq = set(labels)
        n = len(label_uniq) # number of clusters
        fig,ax = plt.subplots(1,n,figsize=(n*3,3))
        i = 0 # count number, for labels
        data_dict = {}
    
        for label in label_uniq:
            idx = np.argwhere(labels==label).reshape(-1)
            c = clist[int(label)]
            if i==0:
                data = plot_spectra(ax[int(label.tolist())],X,Y.iloc[idx,:],['k',c,c],'N1')
            elif i >=1:
                data = plot_spectra(ax[int(label)],X,Y.iloc[idx,:],['k',c,c],'N')
            #data_dict[label] = [data,Y.iloc[idx,:].T]
            data_dict[label] =pd.concat([data, Y.iloc[idx,:].T],axis=1)
            ax[int(label)].set_ylim(Yrange) # same Y range
            ax[int(label)].set_title('cluster {}, n = {}'.format(label,len(idx)))
            i += 1
        ax[0].set_ylabel( "Intensity (a.u.)",fontsize=fs)
        fig.text(0.5,-0.02,"Wavenumber ($cm^{-1}$)",fontsize=fs,ha='center')
        return fig, data_dict#,mol_df,war

    def sort_labels(Y,labels):
        # to make labels in order of max Y values from small to larger 
        Y = pd.DataFrame(Y)
        label_uniq = set(labels)
        data = []
        for label in label_uniq:
            idx = np.argwhere(labels==label).reshape(-1)
            temp_Y = Y.iloc[idx,:].mean().max()
            data.append([label,temp_Y,idx])
        data = pd.DataFrame(data)
        data = data.sort_values(by=1)
        data.iloc[:,0] = np.arange(len(label_uniq))
        labels2 = np.zeros(labels.shape,dtype=np.int32)
        for row in range(data.shape[0]):
            labels2[data.iloc[row,2]] = data.iloc[row,0]
        return labels2

    verbose = {}
    new_pt, Sigma, EVR, PC_dict,vh,verbose = customized_PCA(Y,verbose)# nxd, 

    # pick PCA components
    if isinstance(pca,int):
        data = new_pt[:,:int]
    elif pca=='F':
        data = Y
    elif pca=='T':
        data = new_pt
    # pick optimal cluster number
    opt_n_dict, opt_n, labels, verbose = get_opt(data,clr,verbose)
    # sort labels by Y
    labels = sort_labels(data,labels)
    """
    # original switch for binary
    if isinstance(cm,np.ndarray) and np.unique(labels).shape[0]==2:
        # we switch the labels
        _, labels = switch_cluster(cm, labels)
    """
    if isinstance(cm,np.ndarray):
        labels, cons = switch(cm,labels)
        if len(labels.shape)==3:
            labels = labels[0][0]
    clist = Raman_color.color(10).default_c()# consistency of plotting color

    # spectra plot
    fig1, ax1 = plt.subplots(1,1,figsize=(3,3))
    fig2, ax2 = plt.subplots(1,1,figsize=(3,3))
    data2 = plot_spectra(ax2,X,Y)# data for fig2
    ax2.set_title("total spectra, n={}".format(Y.shape[0]))
    Yrange = [data2.iloc[:,1:].min().min(),data2.iloc[:,1:].max().max()]# Range for figure Y axis
    fig3, data_dict = plot_sep(X,Y,labels,clist,Yrange)#data for fig3
    if pca=='T':
        narrow = 2
        # plot clusters
        Selected_Wavelength, data1 = plotPK(new_pt,labels,EVR,PC_dict,vh,X,narrow,clist, fs,title,mol_df, ax1)
        verbose["Selected wavenumber by PCA"] = Selected_Wavelength

        # Plot predicted UMAP & T-SNE
        fig4,ax4 = plt.subplots(1,2,figsize=(6,3))#,gridspec_kw={'width_ratios': [1.5, 2.8]})
        plt.subplots_adjust(wspace=0.35)
        for (method,i) in zip(["umap",'tsne'],np.arange(2)):
            pred_clf(Y,method,data1.iloc[:,2].to_numpy(),ax4[i],clist) 
        fig4.suptitle("Predicted Dimension Reduction Clustering \n {}".format(title),fontsize='medium')

        #verbose += "Selected wavenumber by PCA: {} \n".format(Selected_Wavelength)

        # add features to spectra plot
        i = 0 # for the alignment of y axis of text
        d_text = 0 # for the distance of two text if there is a wrap
        for (PCA_i,idf,v,mol) in Selected_Wavelength:
            # idf : the i-th feature in Raman
            # v : the selected wavenumber of Raman
            # PCA_i : the i-th PCA components
            if mol==None:
                mol = "unknown"
            ax2.vlines(v,Yrange[0],Yrange[1],linestyles='dashed',lw=1,color='b')
            # plot which PCA : which wavenumber and molecular name
            wrap_thre = 20
            if len(mol)>=wrap_thre:
                wrap_mol = "\n".join(textwrap.wrap(mol,wrap_thre)) 
                d_text += 0.1
            else:
                wrap_mol = mol
            ax2.text(0.03,0.9-i*0.05-d_text,\
                    "{}:{}".format(PCA_i, str(round(v,1))) + "$\mathregular{cm^{-1}}$,"+wrap_mol,\
                    transform=ax2.transAxes,fontsize='x-small')
            i += 1
#             ax2.arrow(v,data.max().max(),0,-3,width=5)
        return opt_n_dict, Selected_Wavelength, verbose, [fig1,fig2,fig3,fig4], [data1,data2,data_dict, labels]
    elif pca=='F':
        return opt_n_dict, verbose, [fig1,fig2,fig3,fig4], [data1,data2,data_dict, labels]
'''

########################################
# This section is to generate figures for predicted UMAP
########################################
def pred_clf(Y,method,colors,ax,clist):
    ####
    # plot the predicted UMAP or T-SNE with labels from the classifier

    ####
    if method=='umap':
        clf = UMAP(random_state=38)
        Y2 = clf.fit_transform(Y)
    elif method=="tsne":
        clf = TSNE(random_state=38,perplexity = min(30,Y.shape[0]-1))
        Y2 = clf.fit_transform(Y)
    ax.scatter(Y2[:,0],Y2[:,1],c=colors,edgecolors='k')
    
    # set legend
    for c in np.unique(colors):
        idx = np.argwhere(colors==c)[0]
        ax.scatter(Y2[idx,0],Y2[idx,1],c=c, edgecolors='k', label = "cluster {}".format(clist.index(c)))
    ax.legend()
    ax.set_xlabel("{} 1".format(method),fontsize='x-large')
    ax.set_ylabel("{} 2".format(method),fontsize='x-large')

