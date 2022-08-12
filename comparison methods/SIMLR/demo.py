# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 09:15:57 2020

@author: Administrator
"""

import os
import sys
sys.path.insert(0,os.path.abspath('..'))
import time
import numpy as np
import pandas as pd
import SIMLR
from SIMLR import helper
from sklearn import metrics
from sklearn.metrics.cluster import adjusted_rand_score as ari
from sklearn.metrics.cluster import normalized_mutual_info_score as nmi

X = pd.read_csv('yan.csv',header=None)
X = np.array(X)
X = X.transpose()

label = pd.read_csv('yan_label.csv')
y=np.array(label)
label = y.ravel() 

c = label.max() # number of clusters
### if the number of genes are more than 500, we recommend to perform pca first!
start_main = time.time()
if X.shape[1]>500:
    X = helper.fast_pca(X,500)
else:
    X = X.todense()
start_main = time.time()
simlr = SIMLR.SIMLR_LARGE(c, 30, 0); ###This is how we initialize an object for SIMLR. the first input is number of rank (clusters) and the second input is number of neighbors. The third one is an binary indicator whether to use memory-saving mode. you can turn it on when the number of cells are extremely large to save some memory but with the cost of efficiency.
S, F,val, ind = simlr.fit(X)
julei = simlr.fast_minibatch_kmeans(F,c)
print('NMI value is %f \n' % nmi(julei.flatten(),label.flatten()))
print('ARI value is %f \n' % ari(julei.flatten(),label.flatten()))
print('HOM value is %f \n' % metrics.homogeneity_score(julei,label))
print("AMI: %0.3f"% metrics.adjusted_mutual_info_score(label, julei))