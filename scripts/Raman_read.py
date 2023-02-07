#!/usr/bin/env python3
#############################################
# This script is to read Raman txt files
#############################################

from load import *

def read_BG(path):
        """
        read separate background txt files
        """
        df = pd.read_csv(path,header=None,sep='\t')
        X_BG = df.iloc[0,1:]
        Y_BG = df.iloc[1:,1:]
        return X_BG,Y_BG

def read_txt(path,typpe='point',BG=10,dataloop=3,metalevel=None,outdir=None):
        """
        Method to read txt and generate pandas Dataframe in the shape of\
        features x samples
        Default parameters are Raman type and number of background spectra
        'point' means point mapping, 'shape' means shape mapping.
        Return is X,Y for shape; X,Y,Y_BG for point
        dataloop and metalevel is for typpe==dir only, denoting dire layers and each layer's meaning
        such as 3 layers, ['treatments','time points','repetitive']
        outdir :: the output directory name
        n x d
        """
        
        if typpe == 'shape':
            df = pd.read_csv(path,sep='\t',header=None)
            X = df.iloc[0,2:]
            Y = df.iloc[1:,2:]
            return X,Y
        if typpe == 'point': 
            df = pd.read_csv(path,sep='\t',header=None)
            n, d = df.shape 
            if n > BG:
                X = df.iloc[0,1:]
                Y = df.iloc[1:(n-BG),1:]
                Y_BG = df.iloc[(n-BG):,1:]
                return X,Y,Y_BG
            else:
                return None, None, None
        if typpe == "dir":
            # walk through the directory to obtain the spectra list
            pathall = [] # save all targeted txt path
            pathallout = [] # generate all output txt path
            metacol = [] # save all meta columns info for the txt such as [treatment conditions, time points, drop]
            for root, dirs, files in os.walk(path):
                if root[len(path):].count(os.sep) < dataloop:
                    for f in files:
                        if not f.startswith("."):
                            temppath = os.path.join(root,f) # the abs path for the files in the given dataloop
                            pathall.append(temppath)
                            
                            pathallout.append(temppath.replace(path.split('/')[-1],outdir).replace(".txt",".xlsx"))
                            tempmetacol = temppath.split("/")[-dataloop:]
                            metacol.append(tempmetacol)

            metacol = pd.DataFrame(metacol)
            metacol.columns = metalevel

            pathall = pd.DataFrame(pathall)
            pathall.columns = ['abspath']

            pathallout = pd.DataFrame(pathallout)
            pathallout.columns = ['abspathout']
            # generate output path
            
            
            df = pd.concat([pathall,pathallout,metacol],axis=1) #
            return df

def combine_sig_bg(ps:str,pb:str,pout:str)->pd.DataFrame:
    """
    combine the signal and background if your signal and background spectra are in separate files
    ps :: path for spectra
    pb :: path for background spectra
    pout :: path to save your combined dataset
    
    """
    sig = pd.read_csv(ps,sep='\t',header=None,index_col=0)
    bg = pd.read_csv(pb,sep='\t',header=None,index_col=0).iloc[1:,:]
    new = pd.concat([sig,bg],axis=0)
    idx = pd.DataFrame(np.arange(new.shape[0]))
    new.index = idx.index
    new = pd.concat([idx,new],axis=1)
    new.index = np.arange(new.shape[0])
    new.to_csv(pout,sep='\t',header=None,index=False)
    return new

def Raman_save(path,data,tabname,typpe='list'):
    """
    path is os.path.join(dir,filename)
    """

    dir = "/".join(path.split('/')[:-1])
    if not os.path.exists(dir):
        os.makedirs(dir)

    if os.path.exists(path):
        mode = 'a'
    else:
        mode = 'w'

    with pd.ExcelWriter(path, mode=mode, engine="openpyxl") as writer:
        
        if typpe=='list':
            # to save a list of data
            for df,tabnamei in zip(data,tabname):
                df.to_excel(writer,sheet_name=tabnamei)
        if typpe == 'single':
            data.to_excel(writer,sheet_name=tabname)