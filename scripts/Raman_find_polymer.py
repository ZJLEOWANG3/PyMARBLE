#!/usr/bin/env python3
#############################################
# This script is to find single or multiple polymers
#############################################

from load import *

def find_peak(X,Y,w=(1150,1200),wid=np.arange(1,10)):
    """
    find single peak within w tuple
    n x d
    X :: pd.Series, wavenumber
    Y :: pd.DataFrame, n x d, samples x feature intensity
    w :: tuple, (x inf lim,x sup lim)
    find peak 
    """
    if not isinstance(Y,pd.DataFrame):
        raise TypeError("This spectra is not pd.DataFrame, please transform Series to DataFrame")
    n, d = Y.shape
    Xid = [] # selected i-th samples 
    Xpeak = [] # selected peak wavenumber
    Ypeak = [] # selected peak intensity
    for i in range(n):
        Yi = Y.iloc[i,:] # i-th sample, all intensity Y
        peakind = scipy.signal.find_peaks_cwt(Yi, widths=wid)# peak id
        Xi = X[peakind] # selected peak X
        Yi = Yi[peakind] # selected intensity Y
        id2 = np.logical_and(Xi>w[0],Xi<w[1]) # selected peak id for the targeted molecules
        if id2.any():
            peakid = peakind[id2]
            ypeak = Yi[peakid].tolist()[0]
            if ypeak>0.5:
                Xpeak.append(X[peakid].tolist()[0])
                Ypeak.append(ypeak)
                Xid.append(i)
    Y2 = Y.iloc[Xid,:]
    Y2.columns = X.to_numpy().tolist()
    return Xid, Y2.T, peakind# Y2 :: n x d


def find_single_polymer(X, Y, w, size=5):
        """
        input
        X: wavelength
        Y: intensity dxn
        w: polymer peak
        size: window size
        
        output
        df_multi_poly: pd
        """
        #to generate intensity for certain polymers
        #output Y2 is a series of intensity of a certain polymer
        def find_poly(X, Y, w, size):
            idx1 = X[np.logical_and(X > w - size, X < w + size)].index
            int1 = Y.loc[idx1, :].sum(0) / len(idx1)
            return int1

        Y2 = find_poly(X, Y, w, size)
        return Y2

def find_multiple_polymer(X,Y,
                          mol_dict=None,
                          size=5,plot='Y',title='Boxplot'):
        """
        mol_dict : the dictionary with Raman wavenumber ~ polymers
        if None, find a location for the json file
        if str, load the file
        if dict, it is the customized dictionary
        """
        #to generate multiple intensity for certain polymers
        #input Y is pandas dataframe, shape = feature x sample
        #output df is pd.DataFrame, shape = intensities of sample x polymers
        if mol_dict==None:
            mol_file = '../../DATA/molecule_dict.json'
        elif isinstance(mol_dict,str):
            mol_file = mol_dict
        elif isinstance(mol_dict,dict):
            molecule_dict = mol_dict
        try:
            with open(mol_file,'r') as mol_f:
                molecule_dict = json.load(mol_f)
        except FileNotFoundError:
            print("File Not Found, Using user-customized dictioanary")
        df_multi_poly = pd.DataFrame()
        for polymer,w in molecule_dict.items():
            Y2 = find_single_polymer(X,Y,w,size=size)
            Y2.rename(polymer,inplace=True)
            df_multi_poly = pd.concat([df_multi_poly,Y2],axis=1)
        if plot=='Y':
            ax = df_multi_poly.boxplot(rot=90, return_type='axes')
            ax.set_title(title)
            return df_multi_poly,ax
        else:
            return df_multi_poly

