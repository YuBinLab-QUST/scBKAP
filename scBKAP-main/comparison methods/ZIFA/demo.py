import numpy as np
from ZIFA import ZIFA
from ZIFA import block_ZIFA
import pandas as pd
from sklearn.cluster import KMeans
from sklearn import metrics
from sklearn.metrics.cluster import adjusted_rand_score as ari
from sklearn.metrics.cluster import normalized_mutual_info_score as nmi

# This gives an example for how to read in a real data called input.table. 
# genes are columns, samples are rows, each number is separated by a space. 
# If you do not want to install pandas, you can also use np.loadtxt: https://docs.scipy.org/doc/numpy/reference/generated/numpy.loadtxt.html
X = pd.read_csv('yan/yan.csv',header=None)
X = np.array(X)
X = X.transpose()

label = pd.read_csv('yan/yan_label.csv')
y=np.array(label)
label = y.ravel() 

Z, model_params = block_ZIFA.fitModel(X, 5)

c = label.max()
kk = KMeans(n_clusters=c)
julei = kk.fit(Z)
julei = julei.labels_

print('NMI value is %f \n' % nmi(julei.flatten(),label.flatten()))
print('ARI value is %f \n' % ari(julei.flatten(),label.flatten()))
print('HOM value is %f \n' % metrics.homogeneity_score(julei,label))
print('AMI value is %f \n' % metrics.adjusted_mutual_info_score(label, julei))

