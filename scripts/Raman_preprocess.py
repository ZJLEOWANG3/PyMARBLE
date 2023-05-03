#!/usr/bin/env python3
##############################################
# This library is to preprocess Raman spectroscopy
# It includes background subtraction, smoothing, 
# baseline correction, binning, normalization
# required packages include pandas, scipy, 
##############################################
from load import *
sys.path.append("/Users/zijianleowang/Desktop/GitHub/RamanomeSpec/scripts")
import Raman_stat

def std_wf(X,Y,Y_BG,peakwindow):
    """
    standard preprocessing workflow
    X : wavenumber
    Y : n x d, intensity
    Y_BG : n x d, intensity for background spectra
    """
    totaln, d = Y.shape
    Y, Xid = remove_burnt(X,Y.T) # Xid :: removed particle index
    
    with open(os.path.join(".",'../data/molecule_dict.json'),'r') as f:
        molecule_dict = json.load(f)
    with open(os.path.join(".",'../data/molecule_win.json'),'r') as f:
        molecule_wn = json.load(f)
    if Y.shape[1] != 0: # this means there is non-burnt cell
        # clean data
        Y = smooth(Y)
        Y = bg_subtraction(Y,Y_BG.T)
        Y = baseline(Y) #d x n
        peaks,mol,phenotype,phenotypename = Raman_find_polymer.get_all_peak(\
            X,Y.T,mol_dict=molecule_dict,window=molecule_wn,wid=peakwindow)
        # phenotype shape is : n x phenotypename
        
        # get the stat and save it
        phenotypebool = phenotype.loc[phenotype['cell']==1,:] # get only cells
        cellspec = Y.loc[:,phenotypebool.index] # d x n
        kl = Raman_stat.comp_KL(cellspec) # chaos by KL divergence
        wd = Raman_stat.comp_WD(cellspec) # chaos by Wasserstein distance
        
        # celln, percentage of cell recovery, PAO, GAO, PHAAO
        basic_stat = Raman_stat.get_stat(totaln, phenotypebool, phenotypename) 
        basic_stat = pd.concat([basic_stat,pd.Series([kl,wd],index=["KL chaos","Wasserstain distance"])])
        
        # save them
        phenotype_spec = [X]
        for i in phenotypename:
            phenotype_speci = cellspec.loc[:,phenotypebool[i]==1].T
            phenotype_spec.append(phenotype_speci)# save n x d
        phenotype_spec = pd.Series(phenotype_spec,index=["wavenumber"]+phenotypename)
    else:
        phenotype_spec = None
        print("All cell burnt")

    return Y,peaks,mol, basic_stat, phenotype_spec

def remove_burnt(X,Y0,threshold=4000):
    """
    find the burnt cells in Raman with raw data X and Y0
    threshold :: the highest intensity value accepted
    Y0 :: d x n
    burnt cell marker peak: 1496-"1589"-1698; 1285-"1361"-1430
    """
    Xid1, _ = Raman_find_polymer.find_peak(X,Y0.T,(1560,1610))
    Xid2, _ = Raman_find_polymer.find_peak(X,Y0.T,(1330,1390))
    Xid = np.intersect1d(Xid1,Xid2) # contains both of the peaks
    Y2 = Y0.copy()
    Y2.columns = np.arange(Y2.shape[1])
    Y2.drop(Xid,axis=1,inplace=True)

    Y2 = Y2[Y2<threshold]    #remove super high intensity
    Y2.dropna(axis=1,how='any',inplace=True)
    #Y2.columns = np.arange(Y2.shape[1])
    return Y2,Xid # dxn


def bg_subtraction(Y,Y_BG,reset=True):
        """
        d x n
        Method to subtract background spectra from sample spectra.
        Y is sample spectra, Y_BG is background spectra
        shape should be features x samples
        Input Type should be Pandas DataFrame
        Output Type is Pandas DataFrame
        """
        if (not isinstance(Y,pd.DataFrame) and \
                            not isinstance(Y_BG,pd.DataFrame)) and \
           (not isinstance(Y,pd.Series) and not isinstance(Y_BG,pd.Series)):
            raise TypeError('Input of bg subtraction is not pandas dataframe')
        if isinstance(Y_BG,pd.Series):
            Y_BG_ave = Y_BG
        else: 
            Y_BG_ave = Y_BG.mean(axis=1)
        if reset==True:
            Y.index = range(Y.shape[0])
            Y_BG_ave.index = range(Y.shape[0])
        else:
            Y_BG_ave.index = Y.index
        Y2 = Y.subtract(Y_BG_ave,axis=0)
        return Y2

def smooth(Y,wl=5,po=3):
        """
        d x n
        Method to smooth the spectra.
        Y is sample spectra.
        Shape should be features x samples
        Input Type should be Pandas DataFrame
        Output Type is Pandas DataFrame
        """
        
        if len(Y.shape)==1:
            Y = pd.Series(Y)
        else:
            Y = pd.DataFrame(Y)
        if (not isinstance(Y,pd.DataFrame)) and (not isinstance(Y,pd.Series)):
            raise TypeError('The input is not pandas dataframe')
        if isinstance(Y,pd.Series):
            Y_smoothed  = pd.Series(sf(Y,window_length=wl,polyorder=po))
        elif isinstance(Y,pd.DataFrame):
            col = Y.columns
            Y_smoothed = pd.DataFrame(sf(Y,wl,po,axis=1))
            Y_smoothed.columns = col
        return Y_smoothed

def baseline(Y,method='IModPoly',degree=20):
        """
        d x n
        Method for baseline corrections. Here, we provide\
        three methods, including Modified polynomial, Improved modified\
        polynomial, ZhangFit, etc.
        method name: 'ModPoly', 'IModPoly','ZhangFit'
        Y is sample spectra.
        Shape should be features x samples
        Input Type should be Pandas DataFrame
        Output Type is Pandas DataFrame
        """
        if len(Y.shape)==1:
            Y = pd.Series(Y)
        else:
            Y = pd.DataFrame(Y)
        if (not isinstance(Y,pd.DataFrame)) and (not isinstance(Y,pd.Series)):
            raise TypeError('The input is not pandas dataframe')
                    
        if isinstance(Y,pd.Series):
            base_obj = BaselineRemoval(Y)
            if method=='ModPoly':
               Y_baselined = pd.Series(base_obj.ModPoly(degree))
            if method=='IModPoly':
               Y_baselined = pd.Series(base_obj.IModPoly(degree))
            if method=='ZhangFit':
               Y_baselined = pd.Series(base_obj.ZhangFit())
            
        elif isinstance(Y,pd.DataFrame):
            Y_baselined = pd.DataFrame()
            for i in range(Y.shape[1]):
                Y_smoothed = Y.iloc[:,i]
                base_obj = BaselineRemoval(Y_smoothed) 
                if method=='ModPoly':
                    Modpoly = base_obj.ModPoly(degree)
                    Y2 = Modpoly
                if method=='IModPoly':
                    IModpoly = base_obj.IModPoly(degree)
                    Y2 = IModpoly
                if method=='ZhangFit':
                    Zhangfit = base_obj.ZhangFit()
                    Y2 = Zhangfit
                Y_baselined = pd.concat([Y_baselined,\
                        pd.Series(Y2,name=Y.iloc[:,i].name)],axis=1)
        
        return Y_baselined

def augment(df,Y,number,typpe='gaussian',s=1):
    """
    df : d,n, np.ndarray
    Y : label
    number: samples generated
    typpe : noise type
    s : stepsize of shift
    how many samples to be generated via three types of augmentation
    add Gaussian noise, Lorentzian noise, or simply shift left or right for one stepsize (~2 cm-1)
    """
    import torch
    if isinstance(df,np.ndarray):
        df = df.copy()
    elif isinstance(df,torch.Tensor):
        df = df.clone()
    Y = Y.copy()
    d,n = df.shape
    np.random.seed(38)
    shift = np.random.choice([-s,s],number)
    sampleid = np.random.choice(np.arange(n),number)
    if typpe=='gaussian':
        noise = np.random.normal(0,1,(d,number))# gaussian
        df[:,sampleid] = np.add(df[:,sampleid],noise)
    elif typpe=='lorenz':
        noise = np.random.standard_cauchy((d,number))
        df[:,sampleid] = np.add(df[:,sampleid],noise)
    elif typpe=='shift':
        df[s:,sampleid] = df[:-s,sampleid]
    X = df[:,sampleid]
    Y = Y[sampleid]

    """
    for i in sampleid:
        if typpe!='shift':
            if typpe=='gaussian':
                noise = np.random.normal(0,1,d)# gaussian
            if typpe=='lorenz':
                rv = cauchy()
                x = np.arange(d)
                noise = rv.pdf(x)
            df[:,i] += noise
        if typpe=='shift':
            # shift left or right randomly
            df[s:,i] = df[:-s,i]
        df[:,i] = Raman_preprocess.smooth(df[:,i])
        df[:,i] = Raman_preprocess.baseline(df[:,i]) # baseline correction
    X = df[:,sampleid]
    Y = Y[sampleid]
    """
    return X, Y

def augmentcombo(X:np.ndarray,Y:np.ndarray,ratio:np.ndarray,v:int,s=1):
    """
    X : feature
    Y : label
    the data augmentation combo is to 
    v : the least sample amounts for all the label
    add more data with a ratio of different strategies ratios = [] to a critical value v
    shift, gaussian, lorenz
    """
    n,d = X.shape
    Xnew0 = np.array([]).reshape(d,-1)
    Ynew0 = pd.array([])
    
    labels,counts = np.unique(Y,return_counts=True)
    for labeli,counti in zip(labels[counts<v],counts[counts<v]):
        idn = Y==labeli
        Xi = X[idn,:]
        Yi = Y[idn]
        newn = v-counti # number of new spectra to be generated
        n = np.ceil(ratio*newn).astype(int) # number of generated spectra for different typpe
        Xnew = np.array([]).reshape(d,-1)
        Ynew = pd.array([])
        for i,typpe in enumerate(['gaussian','lorenz','shift']):
            xn,yn = augment(Xi.T,Yi,n[i],typpe,s) 
            Xnew = np.concatenate([Xnew,xn],axis=1) #dxn
            Ynew = np.concatenate([Ynew,yn])
        Xnew0 = np.concatenate([Xnew0,Xnew],axis=1) #dxn
        Ynew0 = np.concatenate([Ynew0,Ynew])
        # add this to total new X
    return Xnew0.T,Ynew0.T

def binning(df1,df2):
        """
        input shape df1 and df2 : d x n
        This is used to align different spectra with various step size of wavelength from shared start and stop
        key idea is to: use df1 as the reference, df2 will be binned into that size as best;
        df1 - sample spectra; df2 - background spectra
        :return: pruned
        """
        
        def bin(df, new_X):
            """
            X is the series to be modified; new_X is series as the reference for modification
            data is the DF shape of feature x sample
            """
            X = df.iloc[:,0]
            data = df.iloc[:,1:]
            bins = pd.concat([pd.Series(2*new_X[0]-new_X[1]), new_X], ignore_index=True)  # left + stepsize
            digitized = np.digitize(X, bins)
            new_data = pd.DataFrame()
            for i in set(digitized):
                bin_means = data.iloc[digitized == i, :].mean(axis=0)
                new_data = pd.concat([new_data, bin_means], axis=1, ignore_index=True)
            #shape the new_data by new_X
            new_data = new_data.T.iloc[:new_X.shape[0],:]
            new = pd.concat([new_X, new_data], axis=1, ignore_index=True).fillna(0)
            # print(new_X.shape,new_data.shape,new)
            return new
        X1 = df1.iloc[:,0]
        X2 = df2.iloc[:,0]
        X1_feature = (X1.min(),X1.max(),X1.diff().max())#X min, X max, max of X step size
        X2_feature = (X2.min(), X2.max(), X2.diff().max())  # X min, X max, max of X step size
        X_feature = (
                    max(X1_feature[0],X2_feature[0]),
                    min(X1_feature[1],X2_feature[1]),
                    max(X1_feature[2],X2_feature[2])
                     )
        new_X = pd.Series(np.arange(X_feature[0],X_feature[1]+X_feature[2],X_feature[2]))
        t0 = time.time()
        df1_new = bin(df1,new_X)
        t1 = time.time()
        df2_new = bin(df2,new_X)
        t2 = time.time()
        return new_X,df1_new,df2_new

def lowerbound(df,LB=1):
        # LB: lower bound of the dataframe
        # df: dataframe, n x d
        minv = df.min(axis=0)#minimum value
        df = df.add((-1)*minv+LB,axis=1)
        return df

def normalize(df,baseline_id=None,method='spectra'):
        """
        :param baseline_id: the columns between two baseline id are the one normalized in terms of smaller id.
        For example, [0,4,7] as baseline_id; the columns [1,2,3] normalized in terms of 0, [5,6] by 4, [8,...] by 7
        :param method:
        :return:
        """
        if method == 'L1':
            # normalize based on L1
            df = df.div(df.sum(axis=1),axis=0) # n x d
            return df
        if method == "L2":
            df2 = df.applymap(lambda x: x**2)
            df = df.div(df2.sum(axis=1)**(0.5),axis=0)
            return df
        if method=='spectra':
            # map all data to 0 and 1
            def get_max_min(df):
                Temp_max = df.max().max()
                Temp_min = df.min().min()
                return Temp_max,Temp_min

            def normalize(df,max,min,mode='maxmin'):
                df2 = (df-min)/(max-min)
                return df2
            max,min = get_max_min(df)
            df2 = normalize(df,max,min)
        if method=='heatmap':

            if isinstance(baseline_id,list):
                if len(baseline_id)==1:
                    base_df = df.iloc[:,baseline_id[0]]
                    df2 = df
                    for i in range(df.shape[1]):
                        df2.iloc[:,i] = df.iloc[:,i].divide(base_df)
                elif len(baseline_id)>1:
                    df2 = pd.DataFrame()
                    for i in range(len(baseline_id)):
                        i0 = baseline_id[i]
                        i1 = baseline_id[i+1]
                        base_df = df.iloc[:, i0]
                        for col in range(i0+1,i1,1):
                            series = df.iloc[:, col].divide(base_df)
                            df2 = pd.concat([df2,series],axis=1)
            elif isinstance(baseline_id,int):
                temp_baseline_id = [i for i in range(0,df.shape[1],baseline_id)]
                #copy from the len >1
                df2 = pd.DataFrame()
                for i in range(len(temp_baseline_id)):
                    i0 = temp_baseline_id[i]
                    try:
                        i1 = temp_baseline_id[i + 1]
                    except:#the last item
                        i1 = df.shape[1]
                    base_df = df.iloc[:, i0]
                    for col in range(i0 + 1, i1, 1):
                        series = df.iloc[:, col].divide(base_df)
                        series.name = df.columns[col]
                        df2 = pd.concat([df2, series], axis=1)
        return df2

def lowerbound(df,LB=1):
    # LB: lower bound of the dataframe
    # df: dataframe, n x d
    minv = df.min(axis=0)#minimum value
    df = df.add((-1)*minv+LB,axis=1)
    return df

def normalizeeach(df,L=1):
    # normalize each sample
    # L: 1 means L1 normalization
    # df: dataframe, n x d
    if L==1:
        df = df.div(df.sum(axis=1),axis=0)
    return df

def reject_outliers(data, m=1.5):
    return data[abs(data - np.mean(data)) < m * np.std(data)]
