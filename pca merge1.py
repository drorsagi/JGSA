
# generating the PCA of the merged data


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
import matplotlib.colors
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sns
from scipy.cluster import hierarchy
import scipy.spatial.distance as ssd



df1 = pd.read_csv('mergedpolarlipid', header = None)
arr0 = df1.iloc[:, :].values
arr0t = arr0.transpose()

scaler = StandardScaler() # scaling the data (metabolites now in columns!) due to orders of magnitude differences in values between metabolites
scaled_data1 = scaler.fit_transform(arr0t)
print(scaled_data1.mean(axis = 0)) # 4.44089210e-16 - yey (practically mean = 0)
scaled_data = scaled_data1.transpose() # returning to the original matrix, but the metabolites are scaled
scaledlst = scaled_data.tolist()
arr0lst = arr0.tolist()
gp1 = [1] * 15
gp2 = [2] * 15
gp3 = [3] * 15
gp4 = [4] * 14
gpnames0 = gp1 + gp2 + gp3 + gp4
mergeandgrps1 = arr0lst + [gpnames0]
mergeandgrpsscaled = scaledlst + [gpnames0]
mergeandgrps = np.array(mergeandgrps1)
mergedforpca = mergeandgrps.transpose()
mergedscaledarr = np.array(mergeandgrpsscaled)
mergedscaledpca = mergedscaledarr.transpose()


xdf2 = mergedscaledpca[:,0:-1]
ydf2 = mergedscaledpca[:,-1]


colors11 = {1: 'red', 2: 'blue', 3: 'green', 4: 'black'}
colorskey = np.array(['red', 'blue', 'green', 'black'])

marker={1:'o', 2:'v', 3: 'p', 4: '*'} # marker styles - https://matplotlib.org/3.1.1/api/markers_api.html
#



pca2 = PCA(n_components=2)
pca2.fit(xdf2)
print(pca2.explained_variance_ratio_) #
X_pca2 = pca2.transform(xdf2)

p31 = [x[0] for x in X_pca2]
p32 = [x[1] for x in X_pca2]


# all 59 animals with full metabolomics profiling, including legend, colors - blue, red ,green, black for young, old, AL, CR/TF
# p31, p32 are the calculated PCA coordinates

fig101 = plt.figure()

oldlay = plt.scatter(p31[0:15], p32[0:15], marker='o', color=colorskey[0])
younglay = plt.scatter(p31[15:30], p32[15:30], marker='v', color=colorskey[1])
oldestlay = plt.scatter(p31[30:45], p32[30:45], marker='p', color=colorskey[2])
tflay = plt.scatter(p31[45:59], p32[45:59], marker='*', color=colorskey[3])
#
plt.legend((younglay, oldlay, oldestlay, tflay),
           ('8 months', '20 month', '24 months', '24 months TF'),
           scatterpoints=1,
           loc='upper left',
           ncol=2,
           fontsize=6)
# fig101.savefig('pca all 59.pdf') #
plt.show()



datasetF = pd.read_csv('markersOYALmergemonotonicML', header = None) # yey

array = datasetF.values
X00 = array[:, 0:-1] # the columns are the markers, rows are animals, so can do straight up scaling
Y00 = array[:, -1]

scaler = StandardScaler() # scaling the data (metabolites now in columns!) due to orders of magnitude differences in values between metabolites
scaledX00 = scaler.fit_transform(X00)


pca3 = PCA(n_components=2)
pca3.fit(scaledX00)
print(pca3.explained_variance_ratio_) # for O/Y/AL [0.60457381 0.09442841]
X_pca3 = pca3.transform(scaledX00)

p41 = [x[0] for x in X_pca3]
p42 = [x[1] for x in X_pca3]

selected_animals = np.arange(0, 45)
fig131 = plt.figure()
medyoung0 = statistics.median(p41[15:30])
medold0 = statistics.median(p41[:15])
medAL0 = statistics.median(p41[30:45])
for x ,y ,label in zip(p41, p42, Y00):
    plt.scatter(x, y, color=colors11[label], marker = marker[label], s=25) # s for point size
for ii, txt in enumerate(selected_animals):
    if ii == 30:
        plt.annotate('*', (p41[ii], p42[ii]))
plt.axvline(x=medyoung0, ymin=0, ymax=1, color = 'blue', linewidth = 1)
plt.axvline(x=medold0, ymin=0, ymax=1, color = 'red', linewidth = 1)
plt.axvline(x=medAL0, ymin=0, ymax=1, color = 'green', linewidth = 1)

# fig131.savefig('pca 45 11 biomarkers.pdf')
plt.show()





