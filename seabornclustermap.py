

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sns
from scipy.cluster import hierarchy
import scipy.spatial.distance as ssd




df1 = pd.read_csv('mergedscaleddistances', header = None)
distMatrix = df1.iloc[:, :] # [:45, :45] # also for 45 animals (O/Y/AL) [:, :].values
distArray = ssd.squareform(distMatrix)
# define linkage object
distLinkage = hierarchy.linkage(distArray)
gp1 = [1] * 15
gp2 = [2] * 15
gp3 = [3] * 15
gp4 = [4] * 14
gpnames0 = gp1 + gp2 + gp3 + gp4
gpnames1 = pd.Series(gpnames0)

lut = dict(zip(gpnames1.unique(), 'rbgk'))

row_colors = gpnames1.map(lut)

sns.clustermap(distMatrix, row_linkage=distLinkage, col_linkage=distLinkage, cmap = 'YlGnBu', xticklabels = 1, yticklabels = 1, row_colors = row_colors, col_colors = row_colors) #, cmap = 'mako') # mako heatmap  xticklabels=True,

# plt.savefig('clustermtrx.pdf')
plt.show()

