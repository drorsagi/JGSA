

# clculates noise from files of matabolomics data (intensities) for different ages

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


dataunder25 = pd.read_csv('age25andunder',header = None) # - yey
dataunder251 = dataunder25.values
dataarr = np.array(dataunder251)
# each line is sample, each column metabolite (feature) - to calculate noise need to transpose
dataarrt = dataarr.transpose()
dataunder25lst = dataarrt.tolist()
print(len(dataunder25lst), len(dataunder25lst[0])) # 7211 127 (features/metabolites, samples)
meanunder250 = statistics.mean(dataunder25lst[0])
stdunder250 = (statistics.pstdev(dataunder25lst[0]))
noiseunder250 = stdunder250/meanunder250
print(noiseunder250) # 0.18127580919418668

noiseunder25 = calcnoise (dataunder25lst)
print(statistics.median(noiseunder25)) # 0.33858318707113455

dataover58 = pd.read_csv('age58andover',header = None) # - yey
dataover581 = dataover58.values
dataarr58 = np.array(dataover581)
# each line is sample, each column metabolite (feature) - to calculate noise need to transpose
dataarr58t = dataarr58.transpose()
dataover58lst = dataarr58t.tolist()
print(len(dataover58lst), len(dataover58lst[0]))

noiseover58 = calcnoise (dataover58lst)
print(statistics.median(noiseover58)) # 0.3855169352374945

# ages 39 to 45, but there are too many samples there...
dataage38to45 = pd.read_csv('ages38to45',header = None) # - yey
dataage38to451 = dataage38to45.values
dataarr38to45 = np.array(dataage38to451)
# each line is sample, each column metabolite (feature) - to calculate noise need to transpose
dataarr38to45t = dataarr38to45.transpose()
dataage38to45lst = dataarr38to45t.tolist()
print(len(dataage38to45lst), len(dataage38to45lst[0])) # 7211 873

noiseages38to45 = calcnoise (dataage38to45lst)
print(statistics.median(noiseages38to45)) # 0.35223152892525894


# plot the histogram of under25 vs above 58
under25 = noiseunder25 # noiseunder25cutoof
over58 =  noiseover58 # noiseover58cutoof
medunder25 = statistics.median(noiseunder25)
medover58 = statistics.median(noiseover58)
fig, ax = plt.subplots()
bins = np.linspace(0, 1.5, 200) # np.linspace(0, 3, 50), https://www.educative.io/edpresso/what-is-the-linspace-method-in-numpy?utm_source=Google%20AdWords&aid=5082902844932096&utm_medium=cpc&utm_campaign=kb-dynamic-edpresso&gclid=Cj0KCQiAnb79BRDgARIsAOVbhRoykVPp0cG9AS-gvzzJX7jS2jZO0VbVP0-R3jPueOMvz6gEUR8L3PkaAuHZEALw_wcB

plt.hist(under25, bins, alpha=0.5, color = 'blue', label='Young')
plt.hist(over58, bins, alpha=0.5, color = 'red', label='Old')

plt.axvline(x=medunder25, ymin=0, ymax=1, color = 'blue', linewidth = 4)
plt.axvline(x=medover58, ymin=0, ymax=1, color = 'red', linewidth = 4)

ax.tick_params(axis='x', labelsize=22)
ax.tick_params(axis='y', labelsize=22)

plt.show()



