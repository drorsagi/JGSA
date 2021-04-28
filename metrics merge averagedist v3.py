



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

def getstderror(in_out, distarr): # calcuates the std (population; pstdev) of distance matrix distarr, in_out 0 means 0s in the diagonal exclude, 1 intergroup count all elements
    distlst0 = distarr.tolist()
    # print('dist', distlst0)
    if in_out == 0: # distarr represents intragroup distances
        elements0 = []

        for ii in distlst0:
            for jj in ii:
                if jj != 0:
                    elements0.append(jj)
    elif in_out == 1: # distarr represents intergroup distances
        elements0 = []
        for ii in distlst0:
            for jj in ii:
                if jj != 0:
                    elements0.append(jj)
    std0 = statistics.pstdev(elements0)
    # print(len(elements0)) # - yey
    return std0

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



print('intergroup distance')
merge0 = pd.read_csv('mergedpolarlipid', header = None) # the metabolomics data for the merged polar & lipid file
print(merge0.head())
merge0val0 = merge0.iloc[:, :].values #
# Feature Scaling
scaler = StandardScaler()
merge0valt = scaler.fit_transform(merge0val0.transpose()) # the scaling is for the columns, so we scale the transpose matrix (col = metabolites)
merge0val = merge0valt.transpose()




print(merge0val[0,0]) # -1.962701125923673, wrongly scaling the animals gave -0.17130511615095023 - scaled, (before scaling - 0.000390325)
# print(merge0val.mean(axis = 1)) # 4.44089210e-16 pratically zero, axis 1 is the lines in merge0val, which are the metabolites - yey!


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
print(len(oldval), len(oldval[0]), oldval[0][0]) # 434 15 -1.962701125923673 (not scaled 0.000390325) - yey

oldvalarr = np.array(oldval)
oldvalarrt = oldvalarr.transpose()
print(len(oldvalarrt), len(oldvalarrt[0])) # 15 434
oldvallst0 = oldvalarrt.tolist()
youngvalarr = np.array(youngval)
youngvalarrt = youngvalarr.transpose()
print(len(youngvalarrt), len(youngvalarrt[0])) # 15 434
youngvallst0 = youngvalarrt.tolist()
oldyoungvallst0 = oldvallst0 + youngvallst0
print(len(oldyoungvallst0), len(oldyoungvallst0[0])) # 30 434 - yey

# dist1_2 = caldistance([1,2,3], [1,2,4])
# print(dist1_2) # 1 - yey
dist1_1 = caldistance(oldyoungvallst0[0], oldyoungvallst0[0])
dist1_2 = caldistance(oldyoungvallst0[0], oldyoungvallst0[1])
print(dist1_1, dist1_2) # 0.0 360.7584529430399 (not scaled 0.0, 239666801.601786) - looking good!
# going over all possible pairs
distlstot = []
for ii in range(0, len(oldyoungvallst0)):
    rowdist = []
    for jj in range(0, len(oldyoungvallst0)):
        tmpdist = caldistance(oldyoungvallst0[ii], oldyoungvallst0[jj])
        rowdist.append(tmpdist)
    distlstot.append(rowdist)
print(len(distlstot), len(distlstot[0]), distlstot[0][0], distlstot[0][1], distlstot[1][1], distlstot[1][0]) # 30 30 0.0 360.7584529430399 0.0 360.7584529430399 (not scaled 30 30 0.0 239666801.601786 0.0 239666801.601786)
# distlstot is the matrix/list that consist all the distances between old/young groups!
# intragroup average distance - find all intragroup pairs, find their distance, average over them
# the first group is 'old', has 15 members

distmpmtrx = [[0,3,3], [3,0,3], [3,3,0]] # average distance 3
averpairdist = averageingroupdis(3, distmpmtrx)
print(averpairdist) # 3.0 - yey
distlstotarr = np.array(distlstot)
olddistlst = distlstotarr[0:15, 0:15]
print(len(olddistlst), len(olddistlst[0])) # 15 15 - yey
oldaverdist = averageingroupdis(15, olddistlst)
print(oldaverdist) # 392.91898409453125
youngdistlst = distlstotarr[15:, 15:]
print(len(youngdistlst), len(youngdistlst[0])) # 15 15 - yey
print(youngdistlst) #
youngaverdist = averageingroupdis(15, youngdistlst)
print(youngaverdist) # 319.49663046450587
# intergroup distance
oldyoungdistlst = distlstotarr[0:15, 15:]
print(len(oldyoungdistlst), len(oldyoungdistlst[0])) # 15 15
oldyoungaverdist = averageoutgroupdis(15, oldyoungdistlst)
print(oldyoungaverdist) # 421.90849721217114
# # now permuting the labels
# permuting the labels, then pick the corresponding distances from the original matrix. For example, if no. 2 is now # 14,
# then all the distances between no.'i' and no. 2 are replaced by the distances between 'i' and 14 (which is the new 2)
range0 = np.arange(30) # the total size of groups old and young
print('range0', range0) # [ 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29] - yey
rangetoldyoung = range0.tolist()
print('range old-young', rangetoldyoung) # [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29]
permuteoldyoung = np.random.permutation(rangetoldyoung)
print(permuteoldyoung)
# [24 14 15 11  2  1  5 27 17 22  4 21 26 18 10  9  0  3 16 25 12 19 29 13 7  8 20 23  6 28] - yey
# generating the permuted distance matrix
# replacing line 1 with ietm 24th distance from the elements in permuteoldyoung
# looking back at the 3X3 from line 122 (distmpmtrx), some made up permutation
tmppermvector = [1,0,2]
permdist00 = buildidsmatrx(np.array(distmpmtrx), tmppermvector) #
print(permdist00)
# [[0, 3, 3], [3, 0, 3], [3, 3, 0]] - yey
# another test
distmpmtrx1 = [[0,3,1], [3,0,3], [1,3,0]]
permdist00 = buildidsmatrx(np.array(distmpmtrx1), tmppermvector) #
print(permdist00)
# [[0, 3, 3], [3, 0, 1], [3, 1, 0]] - yey
# average distance
averpairdist1 = averageingroupdis(3, distmpmtrx1)
print(averpairdist1) # 2.3333333333333335 - yey (14/6)
print(14/6) # 2.3333333333333335
permdist01 = buildidsmatrx(distlstotarr, permuteoldyoung)
print(len(permdist01[0]), permdist01[0]) # 30, looking good,  - verify
print(len(permdist01), permdist01[1]) # 30, 0s at the diagonal, symmetrical, right dimensions
print(permdist01[2])
# [270.8677751129098, 305.8514809467108, 0.0, 348.1166469342564, 258.079996631177, 292.04749239194604, 438.42205445653906, 407.0033386470068, 292.9168883580974, 288.7996323894231, 232.4099903233905, 240.46271273601133, 430.75493738250066, 305.47679639380016, 442.3428277323883, 357.15540929300744, 487.3632904565961, 300.0764278217727, 480.46630910068916, 356.8892532859533, 413.1425710018274, 275.90971380295446, 360.7584529430399, 246.76716306521368, 630.8520351722232, 368.2628080244128, 255.8729917475618, 240.05998770802532, 344.3910595689221, 290.5095701466949]
# now calculating the intra and intergroup distance
permdistarr = np.array(permdist01)
oldpermdist = permdistarr[0:15, 0:15]

oldpermaverdist = averageingroupdis(15, oldpermdist)
print(oldpermaverdist) # 395.5409465220694 (not scaled 245619507.76641896)
youngpermdist = permdistarr[15:, 15:]

youngpermaverdist = averageingroupdis(15, youngpermdist)
print(youngpermaverdist) # 389.157765335123 ( not scaled 254807718.7233826)
# intergroup distance
oldyoungdpermdist = permdistarr[0:15, 15:]
print(len(oldyoungdistlst), len(oldyoungdistlst[0])) # 15 15
oldyoungpermaverdist = averageoutgroupdis(15, oldyoungdistlst)
print(oldyoungpermaverdist) # 421.90849721217114
# now run 1000 simulations
oldoldyoungdiff = abs(oldyoungaverdist - oldaverdist)
youngoldyoungdiff = abs(oldyoungaverdist - youngaverdist)

print('10000 simulations')
ind0 = 0
for ii in range(0, 1): # 10000
    permuteoldyoung = np.random.permutation(rangetoldyoung)
    permdistoy0 = buildidsmatrx(distlstotarr, permuteoldyoung)
    permdistarr = np.array(permdistoy0)

    oldpermdist = permdistarr[0:15, 0:15]
    oldpermaverdist = averageingroupdis(15, oldpermdist)
    # print(oldpermaverdist)  #
    youngpermdist = permdistarr[15:, 15:]
    youngpermaverdist = averageingroupdis(15, youngpermdist)
    # print(youngpermaverdist)  #
    # intergroup distance
    oldyoungdpermdist = permdistarr[0:15, 15:]
    # print(len(oldyoungdistlst), len(oldyoungdistlst[0]))  # 15 15 - yey
    oldyoungpermaverdist = averageoutgroupdis(15, oldyoungdpermdist)
    # print(oldyoungpermaverdist)
    if (oldyoungpermaverdist > oldpermaverdist) and (oldyoungpermaverdist > youngpermaverdist):
        ind0 = ind0 + 1
print('number of permutations where intergroup (old) distnce is the largest', ind0)
# out of 10 permutations - 2, 100 permutations - 9 # without abs - 10 perm 0, 100 perm - 19, for the abs 10 permutations - 1, 100 permutations - 31...

# so the p-value that the the intregroup distance is the largest (as in the original case) is ~ 0.08 - not (borderline) significant....


# checking specifically Y-AL
pairgroup = generatepairgroup(youngval, alval) # generates the distance matrix of the youngval, alval pair
youngdist = ingpdis(1, 15, pairgroup)
aldist = ingpdis(2, 15, pairgroup)
young_aldist = outgpdis([1,2], 15, pairgroup)
print('young AL average distance', [youngdist, aldist, young_aldist]) # [319.49663046450587, 464.5228587656814, 456.42932070823764]


print('QC originaldistmtrx')
youngvalarr = np.array(youngval)
youngvalarrt = youngvalarr.transpose()
print(len(youngvalarrt), len(youngvalarrt[0])) # 15 434
youngvallst0 = youngvalarrt.tolist()
alvalarr = np.array(alval)
alvalarrt = alvalarr.transpose()
print(len(alvalarrt), len(alvalarrt[0])) # 15 434
alvallst0 = alvalarrt.tolist()
youngalvallst0 = youngvallst0 + alvallst0
print(len(youngalvallst0), len(youngalvallst0[0])) # 30 434 - yey
youngaldist0 = originaldistmtrx(youngalvallst0)
print(len(youngaldist0), len(youngaldist0[0]), youngaldist0[0][0], youngaldist0[0][1], youngaldist0[0][2], youngaldist0[1][0], youngaldist0[1][1], youngaldist0[1][2], youngaldist0[2][0], youngaldist0[2][1], youngaldist0[2][2])
# 30 30 0.0 305.6724255319 243.01211653528443 305.6724255319 0.0 356.78533900468716 243.01211653528443 356.78533900468716 0.0 - looking good

# add error bars, statistics.pstdev(data, mu=None)
distmtrxyoungal = np.array(pairgroup)
youngerror = getstderror(0, distmtrxyoungal[:15, :15])
alrror = getstderror(0, distmtrxyoungal[15:, 15:])
youngalerror = getstderror(1, distmtrxyoungal[15:, :15])
print(youngerror, alrror, youngalerror) # The stdev over the complete symmetricl matrix is the same as

# checking specifically Y-old
pairgroup1 = generatepairgroup(youngval, oldval) # generates the distance matrix of the youngval, oldval pair
youngdist = ingpdis(1, 15, pairgroup1)
olddist = ingpdis(2, 15, pairgroup1)
young_olddist = outgpdis([1,2], 15, pairgroup1)
print('young old average distance', [youngdist, olddist, young_olddist]) # [319.49663046450587, 392.91898409453125, 421.90849721217114]
# add error bars, statistics.pstdev(data, mu=None)
distmtrxyounold = np.array(pairgroup1)
youngerror = getstderror(0, distmtrxyounold[:15, :15])
oldrror = getstderror(0, distmtrxyounold[15:, 15:])
youngolderror = getstderror(1, distmtrxyounold[15:, :15])

# checking specifically old-oldest
pairgroup2 = generatepairgroup(oldval, alval) # generates the distance matrix of the youngval, oldval pair
olddist = ingpdis(1, 15, pairgroup2)
aldist = ingpdis(2, 15, pairgroup2)
old_aldist = outgpdis([1,2], 15, pairgroup2)
print('old al average distance', [olddist, aldist, old_aldist]) # [319.49663046450587, 392.91898409453125, 421.90849721217114]
# add error bars, statistics.pstdev(data, mu=None)
distmtrxoldal = np.array(pairgroup2)
olderror = getstderror(0, distmtrxoldal[:15, :15])
alrror = getstderror(0, distmtrxoldal[15:, 15:])
oldalerror = getstderror(1, distmtrxoldal[15:, :15])



groups00 = ['8', '24', '8-24']
x_pos = np.arange(len(groups00))
distances0 = [youngdist, aldist, young_aldist]
error = [youngerror, alrror, youngalerror]
# Build the plot
fig, ax = plt.subplots()
ax.bar(x_pos, distances0, yerr=error, align='center', alpha=0.5, ecolor='black', capsize=10)

ax.set_ylabel('Distances [AU]')
ax.set_xticks(x_pos)
ax.set_xticklabels(groups00)
ax.set_title('distance within group vs distance between groups')
ax.yaxis.grid(True)

# Save the figure and show
# plt.tight_layout()
# plt.savefig('bar_plot_with_error_bars.png')
plt.show() # yey

# plotting bar graph with standard deviations https://pythonforundergradengineers.com/python-matplotlib-error-bars.html
# groups01 = ['Young', 'Old', 'Young-Old']
groups01 = ['8 months', '20 months', '8-20 months']
x_pos = np.arange(len(groups01))
distances01 = [youngdist, olddist, young_olddist]
error = [youngerror, oldrror, youngolderror]
# Build the plot
fig, ax = plt.subplots()
ax.bar(x_pos, distances01, yerr=error, align='center', alpha=0.5, ecolor='black', capsize=10)
# ax.set_ylabel('Coefficient of Thermal Expansion ($\degree C^{-1}$)')
ax.set_ylabel('Distances [AU]')
ax.set_xticks(x_pos)
ax.set_xticklabels(groups01)
ax.set_title('distance within group vs distance between groups')
ax.yaxis.grid(True)

# Save the figure and show
# plt.tight_layout()
# plt.savefig('inter to intra group distance 8-20.pdf')
plt.show() # yey

# plotting bar graph with standard deviations https://pythonforundergradengineers.com/python-matplotlib-error-bars.html
# groups01 = ['Old', 'AL', 'Old-AL']
groups02 = ['20 months', '24 monhts', '20-24 months']
x_pos = np.arange(len(groups02))
distances02 = [olddist, aldist, old_aldist]
error = [olderror, alrror, oldalerror]
# Build the plot
fig, ax = plt.subplots()
ax.bar(x_pos, distances02, yerr=error, align='center', alpha=0.5, ecolor='black', capsize=10)
# ax.set_ylabel('Coefficient of Thermal Expansion ($\degree C^{-1}$)')
ax.set_ylabel('Distances [AU]')
ax.set_xticks(x_pos)
ax.set_xticklabels(groups02)
ax.set_title('distance within group vs distance between groups')
ax.yaxis.grid(True)

# Save the figure and show
# plt.tight_layout()
# plt.savefig('inter to intra group distance 20-24.pdf')
plt.show() # yey


# for the paper, same bar graph without titles
groups00 = ['Young', 'Oldest(AL)', 'Young-Oldest(AL)']
x_pos = np.arange(len(groups00))
distances0 = [youngdist, aldist, young_aldist]
error = [youngerror, alrror, youngalerror]
# Build the plot
fig, ax = plt.subplots()
ax.bar(x_pos, distances0, yerr=error, align='center', alpha=0.5, ecolor='black', capsize=10)
# ax.set_ylabel('Coefficient of Thermal Expansion ($\degree C^{-1}$)')
ax.set_ylabel('Distances [AU]')
ax.set_xticks(x_pos)
ax.set_xticklabels(groups00)
ax.set_title('distance within group vs distance between groups')
ax.yaxis.grid(True)
plt.show()
