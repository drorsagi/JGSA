

# generates the ditance matrix for clustering


def caldistance(lst1, lst2): # calculates the L1 distance between lst1, lst2. they represent two layers (chcknum1/2) with heir list of metabolomic values
    dist0 = 0
    for ii in range(0, len(lst1)):
        tmp = abs(lst1[ii] - lst2[ii])
        dist0 = dist0 + tmp
    return dist0


def strtofloat0(lst): # returns a list of float, given al ist of numbers in str type
    ans = []
    for ii in lst:
        tmp = [float(kk) for kk in ii]
        ans.append(tmp)
    return ans

def getvalues0(polarval0, start1, len1, start2, len2): # returns lists with the values of the groups that we want to compare
    grp1 = polarval0[:, start1:(start1+len1)]
    grp2 = polarval0[:, start2:(start2+len2)]
    grp1lst = grp1.tolist()
    grp2lst = grp2.tolist()
    return grp1lst, grp2lst

def averageingroupdis(dim0, mtrx0): # calculates the average distance within the same group (array) - sum of all elements, devided by (number of elements minus diagonal)
# dim0 is the number of rows/columns in the array.
# This is a symmetrical matrix with diagonal = 0. element ij is the distance between animal i and j (in the same group)
    sm0 = np.sum(mtrx0)
    numbrofelmnts = ((dim0*dim0) - dim0)
    ans = sm0/numbrofelmnts
    return ans

def averageoutgroupdis(dim0, mtrx0): # calculates the average distance beween two groups (array) - sum of all elements, devided by number of elements
# dim0 is the number of rows/columns in the array, here the diagonal has no meaning - each row is one group and each column is a second group.
# element ij is the distance between animal i and j (in the different groups!)
    sm0 = np.sum(mtrx0)
    numbrofelmnts = ((dim0*dim0))
    ans = sm0/numbrofelmnts
    return ans

def buildidsmatrx(distarr, perm0): # receives the original distance matrix/array and the permutation vector, builds the permutated matrix
    permdist0 = []
    for ii in perm0:
        permtmp = []
        for jj in perm0:
            tmp = distarr[ii, jj] # need to have the indices starting from 0!
            permtmp.append(tmp)
        # print('permtmp', permtmp)
        permdist0.append(permtmp)
    return permdist0

def originaldistmtrx(distarry): # receives the two-group metabolomics data, generates the distance matrix(list)
    distlstot0 = []
    for ii in range(0, len(distarry)):
        rowdist = []
        for jj in range(0, len(distarry)):
            tmpdist = caldistance(distarry[ii], distarry[jj])
            rowdist.append(tmpdist)
        distlstot0.append(rowdist)
    return distlstot0


def generatepairgroup(group01, group02): # generates the distance matrix (array) for the group01-group02 pair
    group01arr = np.array(group01)
    group01arrt = group01arr.transpose()
    print(len(group01arrt), len(group01arrt[0]))  #
    group01lst0 = group01arrt.tolist()
    group02arr = np.array(group02)
    group02arrt = group02arr.transpose()
    print(len(group02arrt), len(group02arrt[0]))  #
    group02lst0 = group02arrt.tolist()
    group0102lst0 = group01lst0 + group02lst0
    print(len(group0102lst0), len(group0102lst0[0])) #
    distlst0 = originaldistmtrx(group0102lst0) # generating the distance matrix (array)
    print(len(distlst0), len(distlst0[0]), distlst0[0][0], distlst0[0][1], distlst0[1][1], distlst0[1][0])

    return distlst0

def ingpdis(gpnum, gpsize, distmtrx): # receives the distance matrix(list), returns the intragroup distance of gpnum
    distmtrxarr = np.array(distmtrx)
    if gpnum == 1: # always size 15
        tmpdistmtrxarr = distmtrxarr[0:gpsize, 0:gpsize]
        sm0 = np.sum(tmpdistmtrxarr)
        numbrofelmnts = ((gpsize * gpsize) - gpsize)
        ans = sm0 / numbrofelmnts
    if gpnum == 2: # should work for size 15 as well as 14
        tmpdistmtrxarr = distmtrxarr[15:, 15:] # starts with No. 15 always
        sm0 = np.sum(tmpdistmtrxarr)
        numbrofelmnts = ((gpsize * gpsize) - gpsize) # goos for size 15 and 14 - this is the matrix size
        ans = sm0 / numbrofelmnts
    return ans

def outgpdis(gset, gpsize, distmtrx): # receives the distance matrix(list), returns the intergroup distance of gset
    distmtrxarr = np.array(distmtrx)
    if gset[1] != 3:
        tmpdistmtrxarr = distmtrxarr[0:gpsize, gpsize:]
        sm0 = np.sum(tmpdistmtrxarr)
        numbrofelmnts = (gpsize * gpsize)
        ans = sm0 / numbrofelmnts
    elif gset[1] == 3:
        tmpdistmtrxarr = distmtrxarr[0:gpsize, gpsize:]
        sm0 = np.sum(tmpdistmtrxarr)
        numbrofelmnts = (gpsize * (gpsize-1))
        ans = sm0 / numbrofelmnts
    return ans



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

# need to do standard scaling (mean 0, std 1) because of the 10 orders of magnitude difference between the values in polar to lipid.


merge0 = pd.read_csv('mergedpolarlipid', header = None) # the metabolomics data for the merged polar & lipid file
print(merge0.head()) # yey, value for first metabolite for animal 58 is 0.002377, corresponds to C80. last animal in the lipid and merged file!!
merge0val0 = merge0.iloc[:, :].values #
merge0val0t = merge0val0.transpose()
# Feature Scaling
scaler = StandardScaler()
merge0valt = scaler.fit_transform(merge0val0t)


print(merge0valt[0,0]) # -1.962701125923673 - scaled, (before scaling - 0.000390325, after scaling the columns (layers...) - -0.17130511615095023)
# print(merge0valt.mean(axis = 0))  # 0 (4.44089210e-16), yey
merge0val = merge0valt.transpose()

valuesoldyoung = getvalues0(merge0val, 0, 15, 15, 15) # returns lists with the values of the groups that we want to compare
oldtmp1 = valuesoldyoung[0]
oldval = strtofloat0(oldtmp1)
youngtmp1 = valuesoldyoung[1]
youngval = strtofloat0(youngtmp1)
valuesALCR = getvalues0(merge0val, 30, 15, 45, 14)
altmp1 = valuesALCR[0]
alval = strtofloat0(altmp1)
crtmp1 = valuesALCR[1]
crval = strtofloat0(crtmp1)
print(len(oldval), len(oldval[0]), oldval[0][0]) # 434 15 -1.962701125923673 (before scaling 0.000390325) - yey

oldvalarr = np.array(oldval)
oldvalarrt = oldvalarr.transpose()
print(len(oldvalarrt), len(oldvalarrt[0])) # 15 434
oldvallst0 = oldvalarrt.tolist()
youngvalarr = np.array(youngval)
youngvalarrt = youngvalarr.transpose()
print(len(youngvalarrt), len(youngvalarrt[0])) # 15 434
youngvallst0 = youngvalarrt.tolist()
alvalarr = np.array(alval)
alvalarrt = alvalarr.transpose()
print(len(alvalarrt), len(alvalarrt[0])) # 15 434
alvallst0 = alvalarrt.tolist()
crvalarr = np.array(crval)
crvalarrt = crvalarr.transpose()
print(len(crvalarrt), len(crvalarrt[0])) # 14 434
crvallst0 = crvalarrt.tolist()
allanimalsvallst0 = oldvallst0 + youngvallst0 + alvallst0 + crvallst0
print(len(allanimalsvallst0), len(allanimalsvallst0[0])) # 59 434 - yey
# print(np.array(allanimalsvallst0).mean(axis = 0)) # 4.44089210e-16 - yey!



# dist1_2 = caldistance([1,2,3], [1,2,4])
# print(dist1_2) # 1 - yey
dist1_1 = caldistance(allanimalsvallst0[0], allanimalsvallst0[0])
dist1_2 = caldistance(allanimalsvallst0[0], allanimalsvallst0[1])
print(dist1_1, dist1_2) # 0.0 360.7584529430399 - looking good!
# going over all possible pairs
distlstot = []
for ii in range(0, len(allanimalsvallst0)):
    rowdist = []
    for jj in range(0, len(allanimalsvallst0)):
        tmpdist = caldistance(allanimalsvallst0[ii], allanimalsvallst0[jj])
        rowdist.append(tmpdist)
    distlstot.append(rowdist)
print(len(distlstot), len(distlstot[0]), distlstot[0][0], distlstot[0][1], distlstot[1][1], distlstot[1][0]) # 59 59 0.0 360.7584529430399 0.0 360.7584529430399
# distlstot is the matrix/list that consist all the distances between animals!
# # write to file the scaled distances matrix
# # # writing to file
# ff1 = open('mergedscaleddistances', 'w')
# for ii in distlstot:
#     for jj in range(0, len(ii)-1):
#         ff1.write(str(ii[jj]) + ',')
#     ff1.write(str(ii[len(ii)-1]) + '\n')
# ff1.close() # yey!!







