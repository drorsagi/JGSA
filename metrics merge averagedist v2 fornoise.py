

# calculating the intragroup average distance, showing that is is significantly increases with age,
# TRF is lower than 24-month-old AL. Make a figure

# note about scaling - the metabolom matrix has the metabolites in rows (animals col),
# and scaling normalizes the columns, so we need to scale the transpose of the metabolom matrix.

# Checking old-young groups. need to expand the code for all 6 pairs!
# Due to orders of magnitude difference in the reads between different metabolites, we need to scale the data. each metabolite is normmalized to have mean 0, std 1.
# this program calculates the average distance within a group and between groups by calculating the the L1 metrics between all possible pairs (of layers)
# across all metabolites. we then have a matrix of distances where the diagonal is 0.
# we can calculate the average distance for the group and between groups. Next we keep the distances but permute the names, and test
# what is the probability to have intra-group average distance higher than the intergroup one
# specifically, for each pair (a,b) we will compare a to ab, b to ab, and the demand is that the distance is >= original distance
# present the data on a bar graph, p-value according to permutations




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

def ingspdistance(distmtrx): # receives the distance matrix(list), returns the intragroup distances all 4 groups; O/Y/AL/CR
    distmtrxarr = np.array(distmtrx)
    permdistances = []
    for ii in range(0, 4):
        if ii != 3: # always size 15
            tmpdistmtrxarr = distmtrxarr[(ii*15):((ii+1)*15), (ii*15):((ii+1)*15)]
            sm0 = np.sum(tmpdistmtrxarr)
            numbrofelmnts = (210) # 15*15 - 15
            permdistances.append(sm0 / numbrofelmnts)
        elif ii == 3:
            tmpdistmtrxarr = distmtrxarr[45:59, 45:59]
            sm0 = np.sum(tmpdistmtrxarr)
            numbrofelmnts = (182)  # 14*14 - 14
            permdistances.append(sm0 / numbrofelmnts)
    return permdistances

def calcsum0(old, young, al): # recieves three group distances, check monotonic trend between young/old/al. if monotonic, returns the sum of abs(difference).
    if (young < old < al):
        ans = abs(old - young) + abs(old - al)
    else:
        ans = 0
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

def make1dlist(array): # recieves a 2-d array returns a 1-d list of all elements
    onedlist = []
    for ii in array.tolist():
        onedlist = onedlist + ii
    return onedlist


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
from scipy.stats import stats

# need to do standard scaling (mean 0, std 1) because of the 10 orders of magnitude difference between the values in polar to lipid.



print('intergroup distance')
merge0 = pd.read_csv('mergedpolarlipid', header = None) # the metabolomics data for the merged polar & lipid file
print(merge0.head())
merge0val0 = merge0.iloc[:, :].values #
# Feature Scaling
scaler = StandardScaler()
merge0valt = scaler.fit_transform(merge0val0.transpose()) # the scaling is for the columns, so we scale the transpose matrix (col = metabolites)
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
print(len(oldval), len(oldval[0]), oldval[0][0]) # 434 15 v (not scaled 0.000390325) - yey

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
print(oldaverdist) # 392.91898409453125 (not scaled 246372927.80372265) - no errors, verify
youngdistlst = distlstotarr[15:, 15:]
print(len(youngdistlst), len(youngdistlst[0])) # 15 15 - yey
print(youngdistlst) # looking good - symmetrical with 0's in the diagonal
youngaverdist = averageingroupdis(15, youngdistlst)
print(youngaverdist) # 319.49663046450587 (not scaled 198695619.0538811)

# # now permuting the labels
# permuting the labels, then pick the corresponding distances from the original matrix. For example, if no. 2 is now # 14,
# then all the distances between no.'i' and no. 2 are replaced by the distances between 'i' and 14 (which is the new 2)



# looking at the ingroup average distances, and getting p-value for the difference
getdistmtrx = originaldistmtrx(merge0val.transpose().tolist()) # calculates the distance matrix for all 59 animals
print(len(getdistmtrx), len(getdistmtrx[-1]), getdistmtrx[0][0], getdistmtrx[0][1], getdistmtrx[0][2], getdistmtrx[1][0], getdistmtrx[1][1], getdistmtrx[1][2], getdistmtrx[2][0], getdistmtrx[2][1], getdistmtrx[2][2]) #
# 59 59 0.0 360.7584529430399 385.0680196051907 360.7584529430399 0.0 270.8677751129098 385.0680196051907 270.8677751129098 0.0 - yey, same as mergedscaleddistances

distances = [392.91898409453125, 319.49663046450587, 464.5228587656814, 366.00724671998097] # O,Y,AL,TRF
difftot = abs(distances[0] - distances[1]) + abs(distances[2] - distances[0]) # difftot = (O-Y) + (AL-O)
# running permutations
range0 = np.arange(59)
perm0 = np.random.permutation(range0)
permdisttmp = buildidsmatrx(np.array(getdistmtrx), perm0)
# calculating the permuted-ingroup distances
qcingpsdist = ingspdistance(getdistmtrx) # QC ingroup distances
print(qcingpsdist) # [392.91898409453125, 319.49663046450587, 464.5228587656814, 366.00724671998097] - yey
permingpdist = ingspdistance(permdisttmp)
print(permingpdist) # for one permutation - [555.9641151107118, 405.52211994358356, 424.97273571968947, 399.20275047630844] - not monotonic, good
diffsum0 = calcsum0(permingpdist[0], permingpdist[1], permingpdist[2])
print(diffsum0) # 0 - yey
# lets run 1000 permutations
# ind0 = 0
# for ii in range(0, 1): # 100, 1000
#     perm0 = np.random.permutation(range0)
#     permdisttmp = buildidsmatrx(np.array(getdistmtrx), perm0)
#     permingpdist = ingspdistance(permdisttmp)
#     diffsum0 = calcsum0(permingpdist[0], permingpdist[1], permingpdist[2])
#     if (diffsum0 > difftot) and (permingpdist[3] < permingpdist[2]):
#         ind0 = ind0 + 1
# print(ind0) # 1 perm - 0, 10 perm - 0, 100 perm - 0, 1000 perm - 5, 1000 perm - 5 (p_value = 0.005) - yey

# barplot

groups00 = ['8 months', '20 months', '24 months AL', '24 months TRF']
x_pos = np.arange(len(groups00))
distances0 = [distances[1], distances[0], distances[2], distances[3]]

# Build the plot
fig, ax = plt.subplots()
ax.bar(x_pos, distances0, align='center', color='blue', capsize=10)
# ax.set_ylabel('Coefficient of Thermal Expansion ($\degree C^{-1}$)')
# ax.set_ylabel('Average Ingroup Distances [AU]')
ax.set_xticks(x_pos)
ax.set_xticklabels(groups00) #, size = 0)
ax.set_title('Average Ingroup Distance With Age') # ax.set_title('Average Ingroup Distance With Age')
# ax.yaxis.grid(True)
# ax.tick_params(axis='x', labelsize=22)
ax.tick_params(axis='y', labelsize=22)
# plt.savefig('intra group distance.pdf')
plt.show()


