#!/bin/python3

# to conduct feature selection for SCRS
# Feature selection to only remain the most informative and non-redundant information in a dataset
# https://scikit-learn.org/stable/modules/classes.html#module-sklearn.feature_selection
# a simple tutorial: https://machinelearningmastery.com/feature-selection-with-real-and-categorical-data/
# compute the 
import Raman_preprocess
import pandas as pd
import numpy as np
import os
import pickle
from scipy.special import kl_div, rel_entr
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.feature_selection import mutual_info_classif as MIC
from sklearn.feature_selection import SelectFromModel
from sklearn.linear_model import Lasso
from sklearn.preprocessing import label_binarize
from sklearn.ensemble import AdaBoostClassifier as ABC
import math
from sklearn.feature_selection import GenericUnivariateSelect, chi2
import skfeature
from skfeature.function.similarity_based import lap_score,trace_ratio, reliefF, fisher_score
from skfeature.utility.construct_W import construct_W
import matplotlib.ticker as ticker
from sklearn.preprocessing import LabelEncoder

def get_label_pair(labels:list):
    # label can have n categories, we need to have index for all pairs of categories for downstream feature selection
    unique_labels = set(labels)

    label_indices = {label: [i for i, l in enumerate(labels) if l == label] for label in unique_labels}

    pairs = itertools.combinations(unique_labels, 2)
    name = list(pairs)
    pairs = itertools.combinations(unique_labels, 2)
    result = [(label_indices[pair[0]], label_indices[pair[1]]) for pair in pairs]
    
    return result, name

def Featuremodels(X,label,model="Adaboost",norm=None):
    """
    tempY is Raman feature: n x d
    templabel is a pd.Series/1D np.ndarray
    norm: L1, L2
    """
    if norm != None:
        X = Raman_preprocess.normalize(X,method=norm).to_numpy()
    if type(X) != np.ndarray:
        X = X.values

    le = LabelEncoder()
    label = le.fit_transform(label)

    if model == "MI":
        metric = MIC(X,label)
    elif model == 'Lasso':
        metric = SelectFromModel(estimator=Lasso(alpha=1)).fit(X, label).estimator_.coef_
    elif model == 'Adaboost':
        metric = SelectFromModel(estimator=ABC(n_estimators=200,random_state=38)).fit(X, label).\
        estimator_.feature_importances_
    elif model == "lap":
        metric = lap_score.lap_score(X,W = construct_W(X))
    elif model == "rfe":
        metric = reliefF.reliefF(X,label)
        #metric = trace_ratio.trace_ratio(X,templabel,n_selected_features=X.shape[0])
    elif model == "fisher":
        metric = fisher_score.fisher_score(X,label)
    return metric

def Featureplot(df,X,cbarname='Relief Score',figsize=(10,10)):
    # df is dataframe for feature ranking scores
    # X is wavenumber for SCRS
    fig2, ax2 = plt.subplots(1,1,figsize=figsize)
    sns.heatmap(df,ax=ax2,
                cmap='YlGnBu',
                cbar_kws={'label':cbarname})
    ax2.figure.axes[-1].yaxis.label.set_size(20)
    ax2.set_yticklabels(labels=ax2.get_yticklabels(),fontsize=15)
    start, end = ax2.get_xlim()
    ax2.xaxis.set_ticks(np.linspace(start, end, 10))
    ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
    labels = []
    for i,v in zip(ax2.get_xticklabels(),np.linspace(start, end+1, 10)):
        new = round((v-start)/(end-start)*(1600-400)+400,0)
        i.set_text(new)
        labels.append(i.get_text())
    ax2.set_xticklabels(labels = labels,fontsize=15,rotation=45,ha='right')
    ax2.hlines(np.arange(ax2.get_ylim()[1],ax2.get_ylim()[0],1), *ax2.get_xlim(),\
               linewidth=.5,color='k'
              )
    ax2.set_xlabel("Wavenumber c$m^{-1}$",fontsize=20)
    ax2.set_ylabel("Sample",fontsize=20)

    # get the wavenumber with highest score
    idx = df.T.idxmax(axis=0)
    maxWave = X[idx].astype(float).round(1) # the wavenumber of the max score
    maxWave.index = df.index
    # add double y axis to show the highest peak
    ax3 = ax2.twinx()
    # ax3.set_aspect("equal")
    ax3.set_ylim([0,ax2.get_ylim()[1]])
    ax3.set_yticks(ax2.get_yticks()-.5)
    ax3.set_yticklabels(maxWave.iloc[::-1], fontsize=10)
    ax3.tick_params(top=False)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
    ax3.spines['left'].set_visible(False)
    ax2.text(0.98,1.04,"  Feature c$m^{-1}$ \nw/ Highest Score",transform=ax2.transAxes,fontsize='large')
    
    return fig2