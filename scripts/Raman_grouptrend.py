#!/usr/bin/env python3
##########################################
# This script is to generate group trend with different groups of same samples
# with multiple filteration criterion
##########################################
from load import *
class grouptrend:
    """
    to detect the trend within the dataframe group based on keywords (aka level)
    """
    def __init__(self,db:pd.DataFrame,kws:list):
        
        """
        db format is n x (d, kws level)
        """
        self.db = db
        self.kws = kws # ['phase','tax']
        
    def lowerbound(self,df,LB=1):
        # LB: lower bound of the dataframe
        # df: dataframe, n x d
        minv = df.min(axis=0)#minimum value
        df = df.add((-1)*minv+LB,axis=1)
        return df
    
    def normalize(self,df,L=1):
        # L: 1 means L1 normalization
        # df: dataframe, n x d
        if L==1:
            df = df.div(df.sum(axis=1),axis=0)
        return df
    
    def compute_chaos(self,df,method='kl',cluster=True):
        """
        compute the system chaos within the dataframe df, shape of n x d
        """
        n, d = df.shape
        klist = []# kl list
        if method=='kl':
            func = scipy.special.kl_div
        if method=='ws':
            func = wasserstein_distance
            
        # all value above 0
        df = self.lowerbound(df,1)
        # L1 normalization : summation = 1
        df = self.normalize(df,1)
        # compute all pairs of them - too slow - do not recommend
        if cluster == "All":
            for i in range(n):
                for j in range(i+1,n):
                    xi = df.iloc[i,:]
                    yi = df.iloc[j,:]
                    klv = func(xi,yi) # KL value
                    klv = klv[np.isfinite(klv)]
                    klv = np.sum(klv)
                    klist.append(klv)
        if cluster == 0:
            xi = df.iloc[0,:]
            for j in range(1,n):
                yi = df.iloc[j,:]
                klv = func(xi,yi) # KL value
                klv = klv[np.isfinite(klv)]
                klv = np.sum(klv)
                klist.append(klv)
        # cluster first, then compute
        if cluster == True:
            kmeans = KMeans(n_clusters=2, random_state=0).fit(df.to_numpy()) # cluster = 2, based on previous results by default
            label = kmeans.labels_
            clusters = []
            for labeli in np.unique(label):
                dfi = df.iloc[label==labeli,:]
                clusters.append(dfi)
            df1 = clusters[0]
            df2 = clusters[1]
            n1 = df1.shape[0]
            n2 = df2.shape[0]
            if n1 > n2:
                nlarge = n1
                nsmall = n2
                dflarge = df1
                dfsmall = df2
            else:
                nlarge = n2
                nsmall = n1
                dflarge = df2
                dfsmall = df1
            for i in range(nsmall):
                for j in range(i+1,nlarge):
                    klv = func(dfsmall.iloc[i,:],dflarge.iloc[j,:])
                    klv = klv[np.isfinite(klv)]
                    klv = np.sum(klv)
                    klist.append(klv)
        klist = np.array(klist)
        ave = np.mean(klist)
        std = np.std(klist)
        return klist # ave, std
    
    
    def transform_db(self,method='kl',cluster=True,level=1):
        """
        transform the original database into a groupby format with ascending order of index based on i-th level
        compute mean and chaos for each group within the groupby
        """
        db = self.db
        kws = self.kws
        groupby = db.groupby(kws)
        mean = groupby.mean()
        chaos = groupby.apply(self.compute_chaos,method=method,cluster=cluster) # tuple of mean and std 
        
        idx = mean.index.tolist()
        idx.sort(key=lambda x:x[level])
        self.mean = mean
        self.chaos = chaos
        return mean, chaos
    
    def get_fold(self,base='Exp',groupn=4):
        """
        df: multi-indexed dataframe
        base: who is the baseline in the multi-indexed dataframe
        groupn: every group contains n samples, here n = 4 by default
        """
        df = self.mean
        n, d = df.shape
        head = np.arange(0,n,groupn)
        fold = df.div(df.xs(base),axis=0).fillna(0) # the fold change
        fold.replace([np.inf,-np.inf],0,inplace=True) # remove inf
        self.fold = fold
        return fold

    def process_fold(self,kw,fillin=0):
        """
        compute the log-10 fold change with sign to indicate increase/decrease of certain signal
        fold is multi-indexed dataframe, with baseline of all 1
        kw is the level key of the dataframe fold
        fillin is the value to fill for NAN, np.inf, -np.inf; by default = 0
        return foldabs
        """
        fold = self.fold
        dfi = fold.xs(kw) # 'S1'
        sign = np.sign(dfi)
        foldabs = np.log10(dfi.abs())
        foldabs.replace([np.inf,-np.inf],0,inplace=True) # absolute of fold-change
        foldabs = sign.multiply(foldabs)# positive and negative
        return foldabs
    
    def unstack_chaos(self):
        chaos = self.chaos
        mean = chaos.unstack().T.apply(lambda x: [y[0] for y in x])
        std = chaos.unstack().T.apply(lambda x: [y[1] for y in x])
        self.chaosm = mean # mean of chaos
        self.chaoss = std # std of chaos
        return mean, std
