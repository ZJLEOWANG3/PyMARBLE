#!/usr/bin/env python3
#############################################
# This script is to load molecular dictionary, 
# Given wavenumber and return molecular name
# To be developed : a section to update dictionary
#############################################

from load import *
#############################################
# This section is to load molecular dictionary
#############################################
def load(mol_path="../DATA/molecule_dict.json"):
    """
    load molecular dictionary and transform it to datafrane
    """
    with open(mol_path,'r') as f:
        mol_dict = json.load(f)
    mol_df = pd.DataFrame.from_dict([mol_dict]).T
    mol_df.columns = ["wv"]
    return mol_df

#############################################
# This section is to get molecular name based on given wavenumber
#############################################
def get_mol(mol_df, wavenumber,size=5):
    """
    get the selected molecular list
    mol_df : molecular dataframe; index : mol name, columns "wv" : wavenumber to the assigned mol name
             typically the results from Raman_molecule.load
    wavenumber : wavenumber used to match the molecules
    size : window size for the matching

    return mol : list
    """
    col_name = mol_df.columns.values[0]
    def closest(lst, K):
        # find the closest number to K in list lst
        return lst[min(range(len(lst)), key = lambda i: abs(lst[i]-K))]
    mol = mol_df[np.logical_and(mol_df[col_name]>=wavenumber-size,\
            mol_df["wv"]<=wavenumber+size)].index.values
    if len(mol)==0:
        mol = None
    elif len(mol)==1:
        mol = mol.tolist()[0]
    elif len(mol)>=2:
        # get the molecule with closest wavenumber
        candidate = mol_df.loc[mol]
        candidate_wv = candidate["wv"].values
        closest_wv = closest(candidate_wv,wavenumber)
        mol = mol_df[mol_df[col_name]==closest_wv].index.values.tolist()[0]
    return mol


