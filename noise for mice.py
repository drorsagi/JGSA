

# calculates noise from files of matabolomics data (intensities) for different ages

def getnoise(std0, mean0): # returns the noise for a list of std/mean. for each
    noise0 = []
    for ii in range(0, len(std0)):
        if mean0[ii] != 0:
            noise0.append(std0[ii]/mean0[ii])
        elif mean0[ii] == 0:
            noise0.append(0)
    return noise0



def calcnoise (samplst0): # list of samples/animals, each element is a list of metabolite values. each element is a feature/metabolite, 'columns' are sample no.

    mean0 = []
    std0 = []
    for ii in samplst0:
        mean0.append(statistics.mean(ii))
        std0.append(statistics.pstdev(ii))
    # print(len(meanold0), len(stdold0), meanold0[0], stdold0[0])  #
    # noise
    noise0 = getnoise(std0, mean0)
    return noise0




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
import scipy.stats as stats
import math


liver0 = pd.read_csv('liver aging metabolomics.csv') # yey liver aging metabolomics.csv, muscle has only 263 features, and the medians are too low - 0.19/0.16 the histogram seems unreliable
liverval = liver0.values

# each line is a metabolite, each column sample, first 8 columns young, last 7 old

livervallst = liverval.tolist()
print(len(livervallst), len(livervallst[0])) # 473 15
meanliver = statistics.mean(livervallst[0][0:8]) # young mice
stdliver = (statistics.pstdev(livervallst[0][0:8]))
noiseliver0 = stdliver/meanliver
print(noiseliver0) # 0.4254147884064086

liveryoung = liverval[:, :8]
liverold = liverval[:, 8:]
liveryounglst = liveryoung.tolist()
liveroldlst = liverold.tolist()

noiseyoung = calcnoise (liveryounglst)
print(statistics.median(noiseyoung)) # 0.28518980816321543
print(len(noiseyoung)) #
# 473
noiseold = calcnoise (liveroldlst)
print(statistics.median(noiseold)) # 0.31140751616615536

# plot the histogram of young vs old
young00 = noiseyoung #
old00 =  noiseold #
medyoung = statistics.median(noiseyoung)
medold = statistics.median(noiseold)
fig, ax = plt.subplots()
bins = np.linspace(0, 1.5, 75)

plt.hist(young00, bins, alpha=0.5, color = 'blue', label='Young')
plt.hist(old00, bins, alpha=0.5, color = 'red', label='Old')

plt.axvline(x=medyoung, ymin=0, ymax=1, color = 'blue', linewidth = 4)
plt.axvline(x=medold, ymin=0, ymax=1, color = 'red', linewidth = 4)

ax.tick_params(axis='x', labelsize=22)
ax.tick_params(axis='y', labelsize=22)

plt.show()



