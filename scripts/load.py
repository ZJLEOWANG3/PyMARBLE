#!/usr/bin/env python3
##########################################
# Load all the packages for other packages
##########################################
import pandas as pd
import torch
import os
import re
from scipy.signal import savgol_filter as sf
from BaselineRemoval import BaselineRemoval
import time
from lmfit.models import GaussianModel, VoigtModel, LinearModel, ConstantModel
from lmfit import Minimizer
import lmfit.models as lmfitM
import math
import lmfit
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import f_oneway, pearsonr, entropy, kstest, spearmanr, kendalltau, cauchy, wasserstein_distance
from scipy.stats import multivariate_normal as mvn
from statsmodels.stats.multicomp import pairwise_tukeyhsd, MultiComparison
#from statsmodels.stats.multicomp import MultiComparison
import seaborn as sns
from sklearn.manifold import TSNE
#import Raman
from bioinfokit.visuz import cluster
# from cycler import cycler
import random
import imp
from tqdm import tqdm
from matplotlib.colors import ListedColormap
from cycler import cycler
import matplotlib
import ete3
import sklearn
import scipy
from umap.umap_ import UMAP
from itertools import permutations, combinations, cycle
from collections import defaultdict
import itertools
import Bio
import skbio
import time
from Bio import Phylo
import json
import plotly.graph_objects as go
import textwrap
from matplotlib import colors, rc
from sklearn.utils.extmath import fast_logdet, randomized_svd, svd_flip
from sklearn.metrics import silhouette_score, davies_bouldin_score, r2_score, roc_auc_score, accuracy_score, roc_curve, auc
from sklearn.metrics import homogeneity_score, completeness_score, v_measure_score, calinski_harabasz_score
from sklearn.mixture import GaussianMixture
from sklearn.cluster import SpectralClustering, AgglomerativeClustering, KMeans
from sklearn.ensemble import RandomForestClassifier as RF
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize, LabelBinarizer
from sklearn.multiclass import OneVsRestClassifier
from sklearn.decomposition import FastICA
import warnings
import pickle
from PIL import Image, ImageDraw
import operator
from scipy.cluster.hierarchy import dendrogram, to_tree
from skbio.stats.distance import anosim
from scipy.spatial.distance import pdist, squareform
from skbio import DistanceMatrix
import sys
sys.path.append("/Users/zijianleowang/Desktop/Projects_in_Cornell/Raman Library/RamanSpec/CODE/scripts")
import Raman_find_polymer, Raman_color, Raman_chmap, Raman_preprocess,Raman_cluster, Raman_figmerge, Raman_stat, Raman_tree
import Raman_read
import Raman_molecule 

