#!/usr/bin/env python3
#############################################
# This script is to find single or multiple polymers
#############################################

from load import *

def find_peak(X,Y,w=(1150,1200),wid=np.arange(1,30)):
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
    Xid = [] # selected i-th samples containing peaks within window w
    Xpeak = [] # selected peak wavenumber
    Ypeak = [] # selected peak intensity
    for i in range(n):
        Yi = Y.iloc[i,:] # i-th sample, all intensity Y
        peakind = scipy.signal.find_peaks_cwt(Yi, widths=wid)# peak id
        Xi = X[peakind] # found peak wavenumbers
        Yi = Yi[peakind] # selected intensity Y
        ii = Yi>0
        if len(ii)==0:
            continue
        Xi = Xi[ii]
        Yi = Yi[ii]
        id2 = np.logical_and(Xi>w[0],Xi<w[1]) # selected peak id for the targeted molecules
        if id2.any(): # found molecule
            peakid = peakind[id2]
            ypeak = Yi[peakid].tolist()[0]
            if ypeak>0.5:
                Xpeak.append(X[peakid].tolist()[0])
                Ypeak.append(ypeak)
                Xid.append(i)
    Y2 = Y.iloc[Xid,:]
    Y2.columns = X.to_numpy().tolist()
    # return Xid : the row number for samples who contains the peak within w
    # return Y2 : bacteria contains the peak
    # return Y3 : bacteria not contains the peak
    return Xid, Y2.T # Y2 :: n x d

def get_all_peak(X,Y,window=5,
                phenotypename = ['cell','PCO','EBPRPAO','GAO1','GAO2','PHBAO','PHBVAO'],
                wid=np.arange(1,30),mol_dict=None):
    """
    n x d
    get all the peaks and identified molecules
    phenotype :: boolean matrix to mention who is the phenotype with col name in phenotypename
    """
    if not isinstance(Y,pd.DataFrame):
        raise TypeError("This spectra is not pd.DataFrame, please transform Series to DataFrame")
    
    # load dict
    if mol_dict==None:
            mol_file = os.path.join(",",'../data/molecule_dict.json')
    elif isinstance(mol_dict,str):
        mol_file = mol_dict
    elif isinstance(mol_dict,dict):
        molecule_dict = mol_dict
    # get peak
    n, d = Y.shape
    peak,mol = [], []
    
    nphenotype = len(phenotypename) # count n types of phenotypes
    phenotype = pd.DataFrame(np.zeros([n,nphenotype]),columns=phenotypename) # dataframe to save samples x phenotype
    for i in range(n): # for each single cell within this dataset such as drop
        Yi = Y.iloc[i,:]
        peakind = scipy.signal.find_peaks_cwt(Yi, widths=wid)# peak id
        Xi = X[peakind] # found peak wavenumbers
        Yi = Yi[peakind] # selected intensity Y
        ii = Yi>0 # intensity above 0 is required
        if len(ii)==0:
            continue
        Xi = Xi[ii]
        Yi = Yi[ii]
        #################
        # get molecule name
        moli = {}
        for polymer,wn in molecule_dict.items():
            if isinstance(window,int):
                tempbool = Xi.between(wn - window,wn + window) # whether it contains a certain molecule
            elif isinstance(window,dict): # customize window for each polymer
                tempbool = Xi.between(wn - window[polymer],wn + window[polymer]) # whether it contains a certain molecule
            if tempbool.any():
                moli[polymer] = Yi[tempbool.tolist()].values[0]
        mol.append(moli)

        ###### count the statistics
        keysi = moli.keys()

        # only if it contains DNA, it is denoted as cell
        if "DNA/RNA,adenine" in keysi:
            phenotype['cell'][i] += 1
        
        # GAO
        if 'glycogen' in keysi:
            phenotype['GAO1'][i] += 1
        if 'glycogen' in keysi and not ('polyP' in keysi and 'O-P-O' in keysi):
            phenotype['GAO2'][i] += 1
        
        # PAO
        if 'polyP' in keysi and 'O-P-O' in keysi:
            phenotype['PCO'][i] += 1
        if 'polyP' in keysi and 'O-P-O' in keysi and ('PHB-co-PHV' in keysi or 'PHB' in keysi):
            phenotype['EBPRPAO'][i] += 1

        if 'PHB-co-PHV' in keysi:
            phenotype['PHBVAO'][i] += 1
        if 'PHB' in keysi:
            phenotype['PHBAO'][i] += 1
        ##############

        peaki = pd.concat([Xi,Yi],axis=1)
        peaki.columns = ['peakwn','peakint'] # peak wavenumber, peak intensity
        peak.append(peaki)

    peak, mol = pd.DataFrame(peak),pd.DataFrame(mol)
    peak = peak.fillna(0)
    mol = mol.fillna(0)

    peak.index = Y.index
    mol.index = Y.index
    phenotype.index = Y.index

    return peak, mol, phenotype, phenotypename


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

