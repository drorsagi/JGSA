


# plots the three effect size graphs

# generates a dictionary of metaboline name and average value
def generatedct(names, values):
    dct0 = {}
    for ii in range(0, len(names)):
        dct0[names[ii]] = values[ii]

    return dct0

def findmaxdif(lst):
    tmp = [lst[0][1], abs(int(lst[0][0]) - int(lst[0][2]))]
    for ii in lst:
        diff0 = abs(int(ii[0]) - int(ii[2]))
        if diff0 > tmp[1]:
            tmp = [ii[1], diff0]
    return tmp

def makedct(excep):
    dct = {}
    for ii in excep:
        dct[ii[1]] = abs(int(ii[0]) - int(ii[2]))
    return dct

def getnoise(std0, mean0): # returns the noise for a list of std/mean. for each
    noise0 = []
    for ii in range(0, len(std0)):
        if mean0[ii] != 0:
            noise0.append(std0[ii]/mean0[ii])
        elif mean0[ii] == 0:
            noise0.append(0)
    return noise0

def ttestlst(lst1, lst2): # calculates the t-test from the valuse of the metabolites
    ttest0 = []
    for ii in range(0, len(lst1)):
        temp = ttest_ind(lst1[ii], lst2[ii])
        ttest0.append(temp[1])
    return ttest0

def getindices(lst): # return the indices having t < 0.05
    ind0 = []
    for ii in range(0, len(lst)):
        if lst[ii] < 0.05:
            ind0.append(ii)
    return ind0

def strtofloat0(lst): # returns a list of float, given al ist of numbers in str type
    ans = []
    for ii in lst:
        tmp = [float(kk) for kk in ii]
        ans.append(tmp)
    return ans


def getvalues0(lstval0, start1, len1, start2, len2): # returns lists with the values of the groups that we want to compare
    arrval0 = np.array(lstval0)
    grp1 = arrval0[:, start1:(start1+len1)]
    grp2 = arrval0[:, start2:(start2+len2)]
    grp1lst = grp1.tolist()
    grp2lst = grp2.tolist()
    return grp1lst, grp2lst


#
def getMNpval00(lst1, lst2): # calculates the mann-whitney p-val over the two lists, each element in the lists is a list of values to compare
    pval00 = []
    # print(len(lst1), len(lst2)) # 140 - good
    for ii in range(0, len(lst1)):
        if (sum(lst1[ii]) + sum(lst2[ii]) == 0): # if all values are 0's give pval = 1, for fdr purposes
            print(ii)
            pval00.append(1)
        if (lst1[ii] != lst2[ii]) and ((sum(lst1[ii]) + sum(lst2[ii])) != 0): # if lst1[ii] is not the same length as lst2[ii] need to check directly not all 0's, so look at sum
            anstmp = mannwhitneyu(lst1[ii], lst2[ii], alternative = 'two-sided')
            pval00.append(anstmp[1])
    return pval00

def fdrandmarkers(lst): # receives a list of p-vaslues, returns the list of markers that passes FDR
    fdr0 = statsmodels.stats.multitest.multipletests(lst, alpha=0.1, method='fdr_bh', is_sorted=False, returnsorted=False)
    mrks = []  # the markers passing FDR
    jj = 1
    for ii in fdr0[0]:
        if ii == True:
            mrks.append(jj)
        jj = jj + 1
    return mrks

def getoverlap(lst1,lst2): # returns the numebr of elemnts that apper in the two lists
    ind = 0
    for ii in lst1:
        if ii in lst2:
            ind = ind + 1
    return ind

def commonmrkrs (lst1, lst2): # returns the elemnts that apper in the two lists
    common0 = []
    for ii in lst1:
        if ii in lst2:
            common0.append(ii)
    return common0

def findmedian(lst): # gets a list, every element is a list of values, returns a list of medians
    med0 = []
    for ii in lst:
        tmp = statistics.median(ii)
        med0.append(tmp)
    return med0

def checkdirec0(lst): # receives the list of medians and checks which elements has the same trend from young to old
    indlst = []
    ind = 0
    for ii in lst: # here ii[0] id pld, ii[1] young, ii[2] AL, ii[3] CR
        # if (ii[1] < ii[0] < ii[3]) or (ii[1] > ii[0] > ii[3]): # CR
        if (ii[1] < ii[0] < ii[2]) or (ii[1] > ii[0] > ii[2]): # only looking at AL
            indlst.append(ind)
        ind = ind + 1
    return indlst

def permdata0(lst): # receives a list, where each element is a list of values (metabolite values) and returns permutation of the tags!!!!
    lstarr = np.array(lst)
    lstarrt = lstarr.transpose()
    tmp = np.random.permutation(lstarrt)
    tmpt = tmp.transpose()
    ans = tmpt.tolist()
    return ans



def sortpolarind(indiceslipid, indicespolar, polarlst0): # receives the originall ist of polar markers, returns the list in same order as lipid markers (same order of animals)
    indiceslipid1 = []
    indicespolar1 = []
    for ii in range (0, len(indicespolar)): # standardizing all animals' names/indices to be uppercase
        indiceslipid1.append(indiceslipid[ii].upper())
        indicespolar1.append(indicespolar[ii].upper())

    polarlst0arr = np.array(polarlst0)
    polarlst0arrt = polarlst0arr.transpose()
    polarlst11  = polarlst0arrt.tolist()

    polarsorted = []
    for ii in range (0, len(indicespolar)):
        if indiceslipid1[ii] == indicespolar1[ii]:
            polarsorted.append(polarlst11[ii])
        else:
            tmp0 = indiceslipid1[ii]
            tmp1 = indicespolar1.index(tmp0)
            polarsorted.append(polarlst11[tmp1])
    polarsortedarr = np.array(polarsorted)
    polarsortedT = polarsortedarr.transpose()
    polarsortedfin = polarsortedT.tolist()
    return polarsortedfin

def writetofile0(markers0, groups, data0, filename0): # writes to file the requested markers for ML, with the proper animal groups - O/Y/AL or O/Y/AL-CR, data is mergelst0, the list of all metabolomics data

    mrkrsforML = []
    for ii in range(0, len(data0)):
        if (ii + 1) in markers0:
            mrkrsforML.append(data0[ii])
    if groups == 3:
        mrkrsarr = np.array(mrkrsforML)
        mrkrs45 = mrkrsarr[:, 0:45]
        mrkrsfin = mrkrs45.tolist()
    elif groups == 4:
        mrkrsfin = mrkrsforML[:]

    # write to file as transpose with group (O,Y,AL,CR = 1,2,3,4) name
    # adding groups' names
    gp1 = [1] * 15
    gp2 = [2] * 15
    gp3 = [3] * 15
    gp4 = [4] * 14
    if groups == 3:
        gpnames0 = gp1 + gp2 + gp3
    elif groups == 4:
        gpnames0 = gp1 + gp2 + gp3 + gp4

    mrkrsforML1 = mrkrsfin + [gpnames0]
    # print(len(mrkrsforML1), len(mrkrsforML1[0]))  #
    mrkrsforMLarr = np.array(mrkrsforML1)
    mrkrsforMLT = mrkrsforMLarr.transpose()
    # print(len(mrkrsforMLT), len(mrkrsforMLT[0]))  #
    mrkrsforMLlst0 = mrkrsforMLT.tolist()
    # writing to file
    ff1 = open(filename0, 'w')
    for ii in mrkrsforMLlst0:
        for jj in range(0, len(ii)-1):
            ff1.write(str(ii[jj]) + ',')
        ff1.write(str(ii[len(ii)-1]) + '\n')
    ff1.close() # yey!!


def getuval(agegp1, agegp2, monotonicmrkrs): # calculates the U values and th MW pvalues between the nested lists
    u_val = []
    p_val = []
    for ii in range(0, len(agegp1)):
        if (ii + 1) in monotonicmrkrs:
            tmp0 = mannwhitneyu(agegp1[ii], agegp2[ii], alternative='two-sided')
            u_val.append(tmp0[0])
            p_val.append(tmp0[1])
    return u_val, p_val

def getuvaldirect(monotonicmrkrs, direct0, minuvalOAL):
    minvaldirect = []
    for ii in range(0, len(monotonicmrkrs)):
        if direct0[ii] == 'up':
            minvaldirect.append(minuvalOAL[ii])
        elif direct0[ii] == 'down':
            minvaldirect.append(-minuvalOAL[ii])
    return minvaldirect






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



merge0 = pd.read_csv('mergedpolarlipid', header = None) # the metabolomics data for the merged polar & lipid file
print(merge0.head()) # yey, value for first metabolite for animal 58 is 0.002377, corresponds to C80. last animal in the lipid and merged file!!
merge0val0 = merge0.iloc[:, :].values #

print(merge0val0[0,0]) # 0.000390325 - yey


valuesoldyoung = getvalues0(merge0val0, 0, 15, 15, 15) # returns lists with the values of the groups that we want to compare
oldtmp1 = valuesoldyoung[0]
oldval = strtofloat0(oldtmp1)
youngtmp1 = valuesoldyoung[1]
youngval = strtofloat0(youngtmp1)
valuesALCR = getvalues0(merge0val0, 30, 15, 45, 14)
altmp1 = valuesALCR[0]
alval = strtofloat0(altmp1)
crtmp1 = valuesALCR[1]
crval = strtofloat0(crtmp1)

anstmp11 = mannwhitneyu(oldtmp1[1], youngtmp1[1], alternative = 'two-sided') # now two tailed MW
print(anstmp11, anstmp11[0]) # two tailed - MannwhitneyuResult(statistic=62.0, pvalue=0.03808828420629213), 62.0 - yey
monotonicmrkrs = [34, 43, 70, 77, 106, 152, 166, 218, 232, 242, 355]
# old vs young
uval = []
pval0 = []
for ii in range(0, len(oldtmp1)):
    if (ii + 1) in monotonicmrkrs:
        tmp0 = mannwhitneyu(oldtmp1[ii], youngtmp1[ii], alternative = 'two-sided')
        uval.append(tmp0[0])
        pval0.append(tmp0[1])
print(uval) # [195.0, 161.0, 57.0, 174.0, 209.0, 171.0, 40.0, 43.0, 178.0, 171.0, 176.0] - verify
print(pval0) # [0.0001702401160966836, 0.04648668639528603, 0.02253107159485839, 0.01140098291155517, 6.836812536001056e-05, 0.01614027989360373, 0.0028226386759533724, 0.004209946028667896, 0.007016199234239618, 0.01614027989360373, 0.00897202481281498]
# need to verify (https://www.socscistatistics.com/tests/mannwhitney/default2.aspx)
# numner 34 - The U-value is 30 (195 is the complementary to 225; n1*n2), The p-value is .00068 ? maybe due to 0s
# number 70 - The U-value is 57, The p-value is .0226 - yey
# number 43 - The U-value is 64 (161 is the complementary to 225; n1*n2), The p-value is .0466 - yey
# number 77 - The U-value is 51 (174 is the complementary to 225), The p-value is .0114 - yey
# number 218 (78 on the lipid list) - The U-value is 43, The p-value is .00424 - close (mine is 0.004209)
# number 355 (215 on the lipid list) - The U-value is 49 (176 is the complementary to 225), The p-value is .00906 - a bit different than the 0.00897...
# change uval to minimum U (min(uval, 225-uval))
minuval = []
for ii in uval:
    minuval.append(min(ii, 225-ii))
print(minuval) # [30.0, 64.0, 57.0, 51.0, 16.0, 54.0, 40.0, 43.0, 47.0, 54.0, 49.0]
direct0 = ['up', 'up', 'down', 'up', 'up', 'up', 'down', 'down', 'up', 'up', 'up'] # the directionality of the marker with age, 'up', 'down'
minvaldirect = []
for ii in range(0, len(monotonicmrkrs)):
    if direct0[ii] == 'up':
        minvaldirect.append(minuval[ii])
    elif direct0[ii] == 'down':
        minvaldirect.append(-minuval[ii])
print(minvaldirect) # [30.0, 64.0, -57.0, 51.0, 16.0, 54.0, -40.0, -43.0, 47.0, 54.0, 49.0] - yey
effectsize0 = [ii/225 for ii in minvaldirect]



# bar graph version 1.0
groups00 = monotonicmrkrs[:]
x_pos = np.arange(len(groups00))
auc0 = effectsize0[:]

#https://stackoverflow.com/questions/34001751/python-how-to-increase-reduce-the-fontsize-of-x-and-y-tick-labels
# Build the plot
fig, ax = plt.subplots()
ax.bar(x_pos, auc0, align='center', color = 'red', alpha=0.5, capsize=10)

ax.set_ylabel('Mann Whitney effect size [AU]', size = 22)
ax.set_xlabel('Marker No.', size = 22)
ax.set_xticks(x_pos)
ax.set_xticklabels(groups00)
ax.tick_params(axis='x', labelsize=22)
ax.tick_params(axis='y', labelsize=22)
ax.set_title('Effect size, Old vs Young', size = 22)
ax.yaxis.grid(True)

# fig.savefig('effect size young-old.pdf')
plt.show() # yey
# effect size mann whitney
# https://www.leeds.ac.uk/educol/documents/00002182.htm#:~:text=Effect%20size%20is%20a%20simple,confounding%20this%20with%20sample%20size.
# https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test#Effect_sizes


# young vs AL
uvalYAL = []
pval0YAL = []
for ii in range(0, len(youngtmp1)):
    if (ii + 1) in monotonicmrkrs:
        tmp0 = mannwhitneyu(youngtmp1[ii], altmp1[ii], alternative = 'two-sided')
        uvalYAL.append(tmp0[0])
        pval0YAL.append(tmp0[1])
print(uvalYAL) # [19.0, 25.0, 212.0, 36.0, 1.0, 4.0, 218.0, 215.0, 12.0, 7.0, 23.0]
print(pval0YAL) # [2.824728717671631e-05, 0.0003078634611799441, 4.0199730927471705e-05, 0.001619713575230349, 4.143220097537536e-06, 7.477207640048691e-06, 1.3294722762606463e-05, 2.3290004782067234e-05, 3.356755084154065e-05, 1.3294722762606463e-05, 0.0002228932910846752]
# looking good, the smaller the min(U, 225) the better (smaller) p-value
minuvalYAL = []
for ii in uvalYAL:
    minuvalYAL.append(min(ii, 225-ii))
print(minuvalYAL) # [19.0, 25.0, 13.0, 36.0, 1.0, 4.0, 7.0, 10.0, 12.0, 7.0, 23.0]
# direct0 = ['up', 'up', 'down', 'up', 'up', 'up', 'down', 'down', 'up', 'up', 'up'] # the directionality of the marker with age, 'up', 'down'
minvaldirectYAL = []
for ii in range(0, len(monotonicmrkrs)):
    if direct0[ii] == 'up':
        minvaldirectYAL.append(minuvalYAL[ii])
    elif direct0[ii] == 'down':
        minvaldirectYAL.append(-minuvalYAL[ii])
print(minvaldirectYAL) # [19.0, 25.0, -13.0, 36.0, 1.0, 4.0, -7.0, -10.0, 12.0, 7.0, 23.0]
effectsize0YAL = [ii/225 for ii in minvaldirectYAL]

# bar graph version 1.0
groups00 = monotonicmrkrs[:]
x_pos = np.arange(len(groups00))
auc0YAL = effectsize0YAL[:]

# Build the plot
fig, ax = plt.subplots()
ax.bar(x_pos, auc0YAL, align='center', color = 'red', alpha=0.5, capsize=10)

ax.set_ylabel('Mann Whitney effect size [AU]', size = 22)
ax.set_xlabel('Marker No.', size = 22)
ax.set_xticks(x_pos)
ax.set_xticklabels(groups00)
ax.tick_params(axis='x', labelsize=22)
ax.tick_params(axis='y', labelsize=22)
ax.set_title('effect size, Young vs Oldest', size = 22)
ax.yaxis.grid(True)
# fig.savefig('effect size young-oldest.pdf')
plt.show()

# old vs AL
uvalOAL = getuval(oldtmp1, altmp1, monotonicmrkrs)
print(uvalOAL[0]) # [57.0, 46.0, 190.0, 61.0, 2.0, 40.0, 203.0, 206.0, 27.0, 27.0, 62.0]
print(uvalOAL[1]) # [0.022231691998559965, 0.006189824426537995, 0.0014040790204992788, 0.03439744889250225, 5.052703540469845e-06, 0.0028226386759533724, 0.00018919298065443355, 0.00011457127866286101, 0.00042246757669391977, 0.00042246757669391977, 0.03808828420629213]

minuvalOAL = []
for ii in uvalOAL[0]:
    minuvalOAL.append(min(ii, 225-ii))
print(minuvalOAL) # [57.0, 46.0, 35.0, 61.0, 2.0, 40.0, 22.0, 19.0, 27.0, 27.0, 62.0]
minuvaldirectOAL = getuvaldirect(monotonicmrkrs, direct0, minuvalOAL)

# minvaldirectOAL = []
# for ii in range(0, len(monotonicmrkrs)):
#     if direct0[ii] == 'up':
#         minvaldirectOAL.append(minuvalOAL[ii])
#     elif direct0[ii] == 'down':
#         minvaldirectOAL.append(-minuvalOAL[ii])
print(minuvaldirectOAL) # [57.0, 46.0, -35.0, 61.0, 2.0, 40.0, -22.0, -19.0, 27.0, 27.0, 62.0] - yey
# markers 232 (92 in lipid list), 242 - verified, both have U = 27
effectsize0OAL = [ii/225 for ii in minuvaldirectOAL]

# bar graph version 1.0
groups00 = monotonicmrkrs[:]
x_pos = np.arange(len(groups00))
auc0OAL = effectsize0OAL[:]

# Build the plot
fig, ax = plt.subplots()
ax.bar(x_pos, auc0OAL, align='center', color = 'red', alpha=0.5, capsize=10)

ax.set_ylabel('Mann Whitney effect size [AU]', size = 22)
ax.set_xlabel('Marker No.', size = 22)
ax.set_xticks(x_pos)
ax.set_xticklabels(groups00)
ax.tick_params(axis='x', labelsize=22)
ax.tick_params(axis='y', labelsize=22)
ax.set_title('effect size, Old vs Oldest', size = 22)
ax.yaxis.grid(True)
# fig.savefig('effect size old-oldest.pdf')
plt.show()







