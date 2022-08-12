# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 09:29:19 2020

@author: Administrator
"""

import numpy as np
import pandas as pd
from lsnmf import Lsnmf
from sklearn import metrics
from sklearn.metrics.cluster import normalized_mutual_info_score as nmi
from sklearn.metrics.cluster import adjusted_rand_score as ari

X = pd.read_csv(('yan.csv'),header=None)
X = np.array(X)

label = pd.read_csv('yan_label.csv')
y=np.array(label)
label = y.ravel() 

c = label.max()
lsnmf = Lsnmf(X, max_iter=50, rank=c)
lsnmf_fit = lsnmf()
sm = lsnmf_fit.summary()

a = sm['predict_samples']
julei = a[0]
y=np.array(julei)
julei =y.ravel()

print('NMI value is %f \n' % nmi(julei.flatten(),label.flatten()))
print('ARI value is %f \n' % ari(julei.flatten(),label.flatten()))
print('HOM value is %f \n' % metrics.homogeneity_score(julei,label))
print('AMI value is %f \n' % metrics.adjusted_mutual_info_score(label, julei))
