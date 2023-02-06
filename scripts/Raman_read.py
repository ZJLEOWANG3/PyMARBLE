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

def read_txt(path,type='point',BG=10):
        """
        Method to read txt and generate pandas Dataframe in the shape of\
        features x samples
        Default parameters are Raman type and number of background spectra
        'point' means point mapping, 'shape' means shape mapping.
        Return is X,Y for shape; X,Y,Y_BG for point
        n x d
        """
        if type == 'shape':
            df = pd.read_csv(path,sep='\t',header=None)
            X = df.iloc[0,2:]
            Y = df.iloc[1:,2:]
            return X,Y
        if type == 'point': 
            df = pd.read_csv(path,sep='\t',header=None)
            X = df.iloc[0,1:]
            Y = df.iloc[1:(df.shape[0]-BG),1:]
            Y_BG = df.iloc[(df.shape[0]-BG):,1:]
            return X,Y,Y_BG
        if type == "dir":
            # read all the scRaman data in the directory
            X = pd.DataFrame()
            Y = pd.DataFrame()
            for sc in os.scandir(path):
                if ".DS" not in sc.name:
                    datai = pd.read_csv(sc.path,sep='\t',header=None)
                    Xi = datai.iloc[:,0]
                    Xi.name = sc.name
                    Yi = datai.iloc[:,1]
                    Yi.name = sc.name

                    X = pd.concat([X,Xi],axis=1)
                    Y = pd.concat([Y,Yi],axis=1)
            return X,Y
