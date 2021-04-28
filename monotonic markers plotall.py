

# this code plots the medians of the 11 markers as a function of age, showing the median trend
# and how the CR group appears younger


import math
from scipy.stats import hypergeom
import statsmodels
import statsmodels.stats.multitest
import statistics
from scipy.stats import mannwhitneyu
from scipy.stats import ttest_ind
from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sns
from scipy.cluster import hierarchy
import scipy.spatial.distance as ssd


datasetF = pd.read_csv('mediansmergemonotonicOYAL', header = None)
values0 = datasetF.values
lstval0 = values0.tolist()
print(lstval0[0]) # [0.038070301, 0.0, 0.414794819, 0.037548572999999995] - yey
print(len(lstval0)) # 11 yey

# plotting for the first marker (34) only
mrkr34 = [lstval0[0][1], lstval0[0][0], lstval0[0][2], lstval0[0][3]]
mrkrnu2 = [lstval0[1][1], lstval0[1][0], lstval0[1][2], lstval0[1][3]]
mrkrnu3 = [lstval0[2][1], lstval0[2][0], lstval0[2][2], lstval0[2][3]]
mrkrnu4 = [lstval0[3][1], lstval0[3][0], lstval0[3][2], lstval0[3][3]]
x = [8, 20, 24, 24]
colors0 = ['blue', 'blue', 'blue', 'red']
for ii in range (0, len(mrkr34)):
    plt.scatter(x[ii],mrkr34[ii], color = colors0[ii])
plt.show() #

shape0 = ['o', 'v', 'p', '*']
# 4 markers at the same time
fig01 = plt.figure()
for ii in range(0, 4): # # 1 is bigger sale, so lates ingnore for now
    # if ii != 1 and ii != 0:
        tmp = [lstval0[ii][1], lstval0[ii][0], lstval0[ii][2], lstval0[ii][3]]
        for jj in range(0, len(tmp)):
            plt.scatter(x[jj], tmp[jj], color=colors0[jj], marker = shape0[ii])
plt.show()

fig, axs = plt.subplots(2,2, figsize=(15, 6) )# , facecolor='w', edgecolor='k')
# fig.subplots_adjust(hspace = .5, wspace=.001)
axs[0, 0].scatter(x,mrkr34)


axs[0, 1].scatter(x,mrkrnu2)
axs[1, 0].scatter(x,mrkrnu3)
axs[1, 1].scatter(x,mrkrnu4)

fig.tight_layout()
plt.show() # yey
# https://jakevdp.github.io/PythonDataScienceHandbook/04.08-multiple-subplots.html
# now in a loop, making the CR in red
fig = plt.figure()
fig.subplots_adjust(hspace=0.4, wspace=0.4)
for ii in range(1,12):
    plt.subplot(3, 4, ii)
    tmp = [lstval0[ii-1][1], lstval0[ii-1][0], lstval0[ii-1][2], lstval0[ii-1][3]]
    for jj in range(0, len(tmp)):
        plt.scatter(x[jj], tmp[jj], color=colors0[jj])
plt.show() # yey

