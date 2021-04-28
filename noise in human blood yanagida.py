


def strytofloat(lst): # converts a nested list where every element is a str into float values
    ans = []
    for ii in lst:
        tmpval = [float(v) for v in ii]
        ans.append(tmpval)
    return ans


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

    noise0 = getnoise(std0, mean0)
    return noise0


def sepyoungold(arr0): # separates an array to young and old lists
    young0 = []
    old0 = []
    for ii in arr0:
        if 'youth' in str(ii[0]):
            young0.append(ii[1:])
        if 'elder' in str(ii[0]):
            old0.append(ii[1:])
    return young0, old0





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




dataset = pd.read_csv('individual variability in human blood python.csv', header = None)
print(dataset.head()) #

datasettmp = np.array(dataset.iloc[:, :])
print(datasettmp[0, 0]) # Person01_plasma_youth_POS - yey
print('youth' in str(datasettmp[0, 0])) # True - yey
datasetyoungold = datasettmp.transpose() # now the beginning of each line is young or old

young0 = []
old0 = []
for ii in datasetyoungold:
    if 'youth' in str(ii[0]):
        young0.append(ii[1:])
    if 'elder' in str(ii[0]):
        old0.append(ii[1:])
print(len(young0), len(young0[0]), young0[0][0], young0[4][0]) # 15 80 2960993.854 4074829.215
print(len(old0), len(old0[0]), old0[0][0], old0[4][0]) # 15 80 8893910.805 1331638.655

young0val = strytofloat(young0)
old0val = strytofloat(old0)

# each elements is all the metabolite under the same person, so need to transpose
youngarr = np.array(young0val)
youngarrt = youngarr.transpose()
youngmetval = youngarrt.tolist()
oldarr = np.array(old0val)
oldarrt = oldarr.transpose()
oldmetval = oldarrt.tolist()

# calculating noise

noiseyoung = calcnoise (youngmetval)
noiseold = calcnoise (oldmetval)

# excluding values above 1.5 and below 0.1 - they are outlayers
noiseyoung01_10 = [ii for ii in noiseyoung if  0.1 < ii < 1.0]
noiseold01_10 = [ii for ii in noiseold if  0.1 < ii < 1.0]
print(statistics.median(noiseyoung01_10), statistics.median(noiseold01_10)) # 0.36542056832692066 0.4515744670405618
print(statistics.mean(noiseyoung01_10), statistics.mean(noiseold01_10)) # 0.411115281409268 0.4717689607958367

fig, ax = plt.subplots()

medyoung01_10 = statistics.median(noiseyoung01_10)
medold01_10 = statistics.median(noiseold01_10)

bins = np.linspace(0, 1.5, 30)

plt.hist(noiseyoung01_10, bins, alpha=0.5, color = 'blue', label='01')
plt.hist(noiseold01_10, bins, alpha=0.5, color = 'red', label='02')

plt.axvline(x=medyoung01_10, ymin=0, ymax=1, color = 'blue', linewidth = 4)
plt.axvline(x=medold01_10, ymin=0, ymax=1, color = 'red', linewidth = 4)

ax.tick_params(axis='x', labelsize=22)
ax.tick_params(axis='y', labelsize=22)

plt.show()



# calculate p-value by random permutation of the labels
datasettmpt = datasettmp.transpose() # now each line is a sample


datasetperm = np.random.permutation(datasettmpt)
# # separating young and old - wrong, just permute and pick 15 of each
# yo00 = sepyoungold(datasetperm)
# print(len(yo00[0]), len(yo00[0][0]), yo00[0][0][0], yo00[0][4][0]) # 15 80 1619368.353 4532454.002 - yey

youngperm = datasetperm[:15, 1:]
oldperm = datasetperm[15:, 1:]

# transpose to make every line a metabolite
youngpermval = strytofloat(youngperm)
oldpermval = strytofloat(oldperm)

# each elements is all the metabolite under the same person, so need to transpose
youngpermarr = np.array(youngpermval)
youngpermarrt = youngpermarr.transpose()
youngmetpermval = youngpermarrt.tolist()
oldpermarr = np.array(oldpermval)
oldpermarrt = oldpermarr.transpose()
oldmetpermval = oldpermarrt.tolist()
print(len(youngmetpermval), len(youngmetpermval[0]), youngmetpermval[0][0], youngmetpermval[4][0], youngmetpermval[0][4]) # 80 15 3480322.414 1500016.991 1619368.353 - yey
print(len(oldmetpermval), len(oldmetpermval[0]), oldmetpermval[0][0], oldmetpermval[4][0], oldmetpermval[0][4]) # 80 15 3540203.498 651416.4129 7862473.221 - yey

# calculate noise, values between 0.1 to 1
noiseyoungperm = calcnoise (youngmetpermval)
noiseoldperm = calcnoise (oldmetpermval)

# excluding values above 1.5 and below 0.1
noiseyoungperm01_10 = [ii for ii in noiseyoungperm if  0.1 < ii < 1.0]
noiseoldperm01_10 = [ii for ii in noiseoldperm if  0.1 < ii < 1.0]
print(statistics.median(noiseyoungperm01_10), statistics.median(noiseoldperm01_10)) # 0.40173270037488773 0.48904883317306297
print(statistics.mean(noiseyoungperm01_10), statistics.mean(noiseoldperm01_10)) # 0.4486176246837722 0.473013270217896
# # run 100 times
# counter0 = 0
# for ii in range (0, 1000):
#     datasetperm = np.random.permutation(datasettmpt)
#
#     youngperm = datasetperm[:15, 1:]
#     oldperm = datasetperm[15:, 1:]
#
#     # transpose to make every line a metabolite
#     youngpermval = strytofloat(youngperm)
#     oldpermval = strytofloat(oldperm)
#
#     # each elements is all the metabolite under the same person, so need to transpose
#     youngpermarr = np.array(youngpermval)
#     youngpermarrt = youngpermarr.transpose()
#     youngmetpermval = youngpermarrt.tolist()
#     oldpermarr = np.array(oldpermval)
#     oldpermarrt = oldpermarr.transpose()
#     oldmetpermval = oldpermarrt.tolist()
#
#     # calculate noise, values between 0.1 to 1
#     noiseyoungperm = calcnoise(youngmetpermval)
#     noiseoldperm = calcnoise(oldmetpermval)
#
#     # excluding values above 1.5 and below 0.1
#     noiseyoungperm01_10 = [ii for ii in noiseyoungperm if 0.1 < ii < 1.0]
#     noiseoldperm01_10 = [ii for ii in noiseoldperm if 0.1 < ii < 1.0]
#     if (statistics.median(noiseoldperm01_10) - statistics.median(noiseyoungperm01_10)) > 0.09: # (medold01_10 - medyoung01_10):
#         counter0 = counter0 + 1
# print(counter0) # 10 permutations - 0, 100 permutations - 9, 100 permutations - 3, 1000 permutations - 26 (pval = 0.026)




