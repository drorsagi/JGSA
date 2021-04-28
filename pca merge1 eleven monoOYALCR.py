
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



colors11 = {1: 'red', 2: 'blue', 3: 'green', 4: 'black'}


marker={1:'o', 2:'v', 3: 'p', 4: '*'} # marker styles - https://matplotlib.org/3.1.1/api/markers_api.html

datasetF = pd.read_csv('markersOYALCRmergemonotonic', header = None) # yey

array = datasetF.values
X100 = array[:, 0:-1] # the columns are the markers, rows are animals, so can do straight up scaling
Y100 = array[:, -1]

scaler = StandardScaler() # scaling the data (metabolites now in columns!) due to orders of magnitude differences in values between metabolites
scaledX100 = scaler.fit_transform(X100)


pca4 = PCA(n_components=2)
pca4.fit(scaledX100)
print(pca4.explained_variance_ratio_) # [0.55488318 0.10822533]
X_pca4 = pca4.transform(scaledX100)
fig23 = plt.figure() # (figsize=(9, 2), dpi=300)
p51 = [x[0] for x in X_pca4]
p52 = [x[1] for x in X_pca4]
medyoung0 = statistics.median(p51[15:30])
medold0 = statistics.median(p51[:15])
medAL0 = statistics.median(p51[30:45])
medTRF0 = statistics.median(p51[45:])
for x ,y ,label in zip(p51, p52, Y100):
    plt.scatter(x, y, color=colors11[label], marker = marker[label], s=25) # s for point size

plt.axvline(x=medyoung0, ymin=0, ymax=1, color = 'blue', linewidth = 1)
plt.axvline(x=medold0, ymin=0, ymax=1, color = 'red', linewidth = 1)
plt.axvline(x=medAL0, ymin=0, ymax=1, color = 'green', linewidth = 1)
plt.axvline(x=medTRF0, ymin=0, ymax=1, color = 'black', linewidth = 1)

plt.show() # yey

