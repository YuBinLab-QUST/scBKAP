# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 08:49:43 2020

@author: Administrator
"""

import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from sklearn import metrics
from sklearn.metrics.cluster import adjusted_rand_score as ari
from sklearn.metrics.cluster import normalized_mutual_info_score as nmi

X = pd.read_csv('yan/yan.csv',header=None)
X = np.array(X)
X = X.transpose()

label = pd.read_csv('sim_seu.csv')
y=np.array(label)
label = y.ravel() 

tsne = TSNE()
A = tsne.fit_transform(X) #进行数据降维,并返回结果

c = label.max()
kk = KMeans(n_clusters=c)
julei = kk.fit(A)
julei = julei.labels_

print('NMI value is %f \n' % nmi(julei.flatten(),label.flatten()))
print('ARI value is %f \n' % ari(julei.flatten(),label.flatten()))
print('HOM value is %f \n' % metrics.homogeneity_score(julei,label))
print('AMI value is %f \n' % metrics.adjusted_mutual_info_score(label, julei))


