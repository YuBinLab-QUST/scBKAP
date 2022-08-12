# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 16:47:33 2019

@author: Administrator
"""

from numpy import *
import matplotlib.pyplot as plt
 
# 加载本地数据
def loadDataSet(fileName):
    dataMat = []
    fr = open(fileName)
    for line in fr.readlines():
        curLine = line.strip().split('\t')
        fltLine = list(map(float,curLine)) # 映射所有数据为浮点数
        dataMat.append(fltLine)
    return dataMat
 
# 欧式距离计算
def distEclud(vecA, vecB):
    return sqrt(sum(power(vecA - vecB, 2))) # 格式相同的两个向量做运算
 
# 中心点生成 随机生成最小到最大值之间的值
def randCent(dataSet, k):
    n = shape(dataSet)[1]
    centroids = mat(zeros((k,n))) # 创建中心点，由于需要与数据向量做运算，所以每个中心点与数据得格式应该一致（特征列）
    for j in range(n):  # 循环所有特征列，获得每个中心点该列的随机值
        minJ = min(dataSet[:,j])
        rangeJ = float(max(dataSet[:,j]) - minJ)
        centroids[:,j] = mat(minJ + rangeJ * random.rand(k,1)) # 获得每列的随机值 一列一列生成
    return centroids
 
# 返回 中心点矩阵和聚类信息
def kMeans(dataSet, k, distMeas=distEclud, createCent=randCent):
    m = shape(dataSet)[0]
    clusterAssment = mat(zeros((m,2))) # 创建一个矩阵用于记录该样本 （所属中心点 与该点距离）
    centroids = createCent(dataSet, k)
    clusterChanged = True
    while clusterChanged:
        clusterChanged = False  # 如果没有点更新则为退出
        for i in range(m):
            minDist = inf; minIndex = -1
            for j in range(k):  # 每个样本点需要与 所有 的中心点作比较
                distJI = distMeas(centroids[j,:],dataSet[i,:])  # 距离计算
                if distJI < minDist:
                    minDist = distJI; minIndex = j
            if clusterAssment[i,0] != minIndex: # 若记录矩阵的i样本的所属中心点更新，则为True，while下次继续循环更新
                clusterChanged = True
            clusterAssment[i,:] = minIndex,minDist**2   # 记录该点的两个信息
        # print(centroids)
        for cent in range(k): # 重新计算中心点
            # print(dataSet[nonzero(clusterAssment[:,0] == cent)[0]]) # nonzero返回True样本的下标
            ptsInClust = dataSet[nonzero(clusterAssment[:,0].A==cent)[0]] # 得到属于该中心点的所有样本数据
            centroids[cent,:] = mean(ptsInClust, axis=0) # 求每列的均值替换原来的中心点
    return centroids, clusterAssment
 
def biKmeans(dataSet, k, distMeas=distEclud):
    m = shape(dataSet)[0]
    clusterAssment = mat(zeros((m,2)))  # 保存数据点的信息（所属类、误差）
    centroid0 = mean(dataSet, axis=0).tolist()[0]   # 根据数据集均值获得第一个簇中心点
    centList =[centroid0] # 创建一个带有质心的 [列表]，因为后面还会添加至k个质心
    for j in range(m):
        clusterAssment[j,1] = distMeas(mat(centroid0), dataSet[j,:])**2     # 求得dataSet点与质心点的SSE
    while (len(centList) < k):
        lowestSSE = inf
        for i in range(len(centList)):
            ptsInCurrCluster = dataSet[nonzero(clusterAssment[:,0].A==i)[0],:] # 与上面kmeans一样获得属于该质心点的所有样本数据
            # 二分类
            centroidMat, splitClustAss = kMeans(ptsInCurrCluster, 2, distMeas)  # 返回中心点信息、该数据集聚类信息
            sseSplit = sum(splitClustAss[:,1])  # 这是划分数据的SSE    加上未划分的 作为本次划分的总误差
            sseNotSplit = sum(clusterAssment[nonzero(clusterAssment[:,0].A!=i)[0],1])   # 这是未划分的数据集的SSE
            if (sseSplit + sseNotSplit) < lowestSSE:    # 将划分与未划分的SSE求和与最小SSE相比较 确定是否划分
                bestCentToSplit = i         # 得出当前最适合做划分的中心点
                bestNewCents = centroidMat  # 划分后的两个新中心点
                bestClustAss = splitClustAss.copy() # 划分点的聚类信息
                lowestSSE = sseSplit + sseNotSplit
        bestClustAss[nonzero(bestClustAss[:,0].A == 1)[0],0] = len(centList)   # 由于是二分，所有只有0，1两个簇编号，将属于1的所属信息转为下一个中心点
        bestClustAss[nonzero(bestClustAss[:,0].A == 0)[0],0] = bestCentToSplit # 将属于0的所属信息替换用来聚类的中心点
        centList[bestCentToSplit] = bestNewCents[0,:].tolist()[0] # 与上面两条替换信息相类似，这里是替换中心点信息，上面是替换数据点所属信息
        centList.append(bestNewCents[1,:].tolist()[0])
        clusterAssment[nonzero(clusterAssment[:,0].A == bestCentToSplit)[0],:] = bestClustAss # 替换部分用来聚类的数据的所属中心点与误差平方和 为新的数据
    return mat(centList), clusterAssment
 
#centList,clusterAssment = biKmeans(data_phate,7)
#
#julei = mat(loadDataSet('label.txt'))
#y=np.array(julei)
#julei = y.ravel() 
#print("NMI: %0.3f"% metrics.normalized_mutual_info_score(label, julei))
#print("ARI: %0.3f"% metrics.adjusted_mutual_info_score(label, julei))