
# This program calculates the p-value for the noise to increase with age (and TRF being smaller than AL).
# can read directly the merged file and generate the O/Y/AL/CR groups




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


    # perm = []
    # for ii in lst:
    #     tmp = np.random.permutation(ii)
    #     perm.append(tmp)
    # return perm

def sortpolarind(indiceslipid, indicespolar, polarlst0): # receives the originall ist of polar markers, returns the list in same order as lipid markers (same order of animals)
    indiceslipid1 = []
    indicespolar1 = []
    for ii in range (0, len(indicespolar)): # standardizing all animals' names/indices to be uppercase
        indiceslipid1.append(indiceslipid[ii].upper())
        indicespolar1.append(indicespolar[ii].upper())
    # print(indiceslipid1[30], indicespolar1[30]) # L1 L65 - yey
    # transposing the values, to be able to exchange columns
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

def transform0(lst0, nn): # recieves a list of metabolites, each metabolites with the corresponding 59 values, returns the list nn animals, each with it metabolites' values
    lst0arr = np.array(lst0)
    lst0arrt = lst0arr.transpose()
    lst0arrtnn = lst0arrt[:nn, :]
    lst0nn = lst0arrtnn.tolist()
    return lst0nn

def calcnoise (animallst0): # list of animals, each with their metabolites. need to transose to list of metabolites with the corresponding animal-values
    tmp = np.array(animallst0)
    tmpt = tmp.transpose()
    metvals = tmpt.tolist()
    mean0 = []
    std0 = []
    for ii in metvals:
        mean0.append(statistics.mean(ii))
        std0.append(statistics.pstdev(ii))
    # print(len(meanold0), len(stdold0), meanold0[0], stdold0[0])  #
    # noise
    noise0 = getnoise(std0, mean0)
    return noise0

def calcsum0(old, young, al): # recieves three median values, check monotonic trend between young/old/al. if monotonic, returns the sum of abs(difference).
    if (young < old < al):
        ans = abs(old - young) + abs(old - al)
    else:
        ans = 0
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
import scipy.stats as stats
import math

print('test') # yey
lst = [1,2,3]
arr = np.array(lst)
print(arr) # yey

# generating the merged file
polar = pd.read_csv('chicken59  polar final.csv')
print(polar.head()) # yey
polarval = polar.iloc[1:, 1:].values # first two lines are the classes (old, young, AL, CR = 1,2,3,4), note - the first line is the header!. first  column is the metabolite names
polarnames = polar.iloc[1:, 0] # array of metabolites names
print(polarval[0,0]) # 0.000390325 - yey
print(type(polarval[0,0])) # <class 'str'> - need to change to float to work with the numbers!
print(type(polarnames)) # <class 'pandas.core.series.Series'>
print(polarnames[1]) # 1,3-Dimethyluracil -yey, starts from 1 not 0...
polarlst0 = polarval.tolist()

lipid = pd.read_csv('chicken59 lipidomics.csv')
print(lipid.head()) # yey
lipidval = lipid.iloc[1:295, 2:].values # first two lines are the classes (old, young, AL, CR = 1,2,3,4), note - the first line is the header!. first column is the metabolite names, second fattyAcids
lipidnames = lipid.iloc[1:295, 0] # array of metabolites names, there is 'e- ether' after the list ends
print(lipidval[0,0], lipidval[-1,0]) # 2317.147449, 72051.23054 - yey
print(type(lipidval[0,0])) # <class 'str'> - need to change to float to work with the numbers!
print(type(lipidnames)) # <class 'pandas.core.series.Series'>
print(lipidnames[1]) # AcCa(16:0) - yey, starts from 1 not 0...
lipidlst0 = lipidval.tolist()

# before merging need to sort the columns so they match. For example, in the original files the first AL animas is L1 lipidomics, L65 polar
indicespolartmp = polar.iloc[0, 1:]
indicespolar = indicespolartmp.tolist()
indiceslipidtmp = lipid.iloc[0, 2:]
indiceslipid = indiceslipidtmp.tolist()
print(len(indicespolar), len(indiceslipid), indicespolar[0], indicespolar[-1], indiceslipid[0], indiceslipid[-1]) # 59 59 O1 C16 O1 C80 - yey

# we want everything to be sorted like the lipid column, which is ascending order in animal value
sortedlpolar0 = sortpolarind(indiceslipid, indicespolar, polarlst0)
print(lipidlst0[0][15], lipidlst0[0][30], lipidlst0[0][45], lipidlst0[0][-1]) # 4503.415405 4840.658791 0 361.7154458 (Y106, L1, C9, C80)
print(sortedlpolar0[0][15], sortedlpolar0[0][30], sortedlpolar0[0][45], sortedlpolar0[0][-1]) # 0.001339296 0.000973921 0.002258703 0.002377372 (y106, L1, C9, C80) - yey!


mergelst0 = sortedlpolar0 + lipidlst0
print('QC merged metabolomics list')
print(len(mergelst0), len(mergelst0[0]), len(mergelst0[-1])) # 434 59 59 (434 = 140 + 294) - yey
print(mergelst0[0][0], mergelst0[0][-1], mergelst0[-1][0], mergelst0[140][0], mergelst0[140][-1]) # 0.000390325 0.002377372 72051.23054 2317.147449 361.7154458 (O1[0], C80[0], O1[0] - lipid, C80[0] - lipid)
print('QC merged metabolomics list - done') # yey

# write the merged metabolomics data into file
# # writing to file
# ff1 = open('mergedpolarlipid', 'w')
# for ii in mergelst0:
#     for jj in range(0, len(ii)-1):
#         ff1.write(str(ii[jj]) + ',')
#     ff1.write(str(ii[len(ii)-1]) + '\n')
# ff1.close() # yey!!

# generating the four groups (Old, Young, AL, CR)
print('old vs young')
lstofval0 = getvalues0(mergelst0, 0, 15, 15, 15)
oldtmp00 = lstofval0[0]
oldtmp1 = strtofloat0(oldtmp00)
youngtmp00 = lstofval0[1]
youngtmp1 = strtofloat0(youngtmp00)
print(len(oldtmp1[0]), len(oldtmp1), oldtmp1[0][1]) # 15 434 0.001642452 - yey
print(len(youngtmp1[0]), len(youngtmp1), youngtmp1[0][1]) # 15 434 0.002295563 - yey


# now AL vs CR
lstofval0 = getvalues0(mergelst0, 30, 15, 45, 14)
altmp00 = lstofval0[0]
altmp1 = strtofloat0(altmp00)
crtmp00 = lstofval0[1]
crtmp1 = strtofloat0(crtmp00)
print('AL vs CR')
print(len(altmp1[0]), len(altmp1), altmp1[0][1]) # 15 434 0.003053628 (L2) - yey
print(len(crtmp1[0]), len(crtmp1), crtmp1[0][1]) # 14 434 0.004467691 (C10) - yey


print('now noise')
print(len(oldtmp1[0]), len(oldtmp1), oldtmp1[0][1]) # 15 434 0.001642452 - yey
# noise of oldtmp1[0]
mean11 = statistics.mean(oldtmp1[0])
std11 = statistics.pstdev(oldtmp1[0]) # pstdev like numpy std (np.std) uses population standard deviation
print(mean11, std11) # 0.002699757133333333 0.0009042609951173658 - yey
print(oldtmp1[0]) # yey
# [0.000390325, 0.001642452, 0.002573045, 0.003296381, 0.002684807, 0.003102332, 0.00335887, 0.00239444, 0.003543983, 0.002630396, 0.004403292, 0.00211244, 0.002844785, 0.002143229, 0.00337558]
# calculating mean, std for oldtmp1 - the metabolomics of old (20 months) layers
meanold0 = []
stdold0 = []
for ii in oldtmp1:
    meanold0.append(statistics.mean(ii))
    stdold0.append(statistics.pstdev(ii))
print(len(meanold0), len(stdold0), meanold0[0], stdold0[0]) # 434 434 0.002699757133333333 0.0009042609951173658 - yey
# noise
noiseold0 = getnoise(stdold0, meanold0)
print(len(noiseold0), noiseold0[0], noiseold0[-1]) # 434 0.3349416078774811 0.2712327744537025
#

print(max(noiseold0)) # 2.616972086836892
maxold0 = max(noiseold0)
print(noiseold0.index(maxold0)) # 344 - metabolite # 345, lipid #205 - PE(32:1)
print(meanold0[344], stdold0[344], stdold0[344]/meanold0[344]) # 12406.739446 32468.090818840203 2.616972086836892 (it has many 0's, check directly via calculator.net - yey)
# print(oldtmp1[40]) # calculating mean, std (via calculator.net) - 0.00264958088886 0.0064650557210965 - yey
# # [0.002240716, 0.02221791, 0.015284045, 0.0, 0.0, 9.60713e-07, 0.0, 0.0, 0.0, 0.0, 8.16199e-08, 0.0, 0.0, 0.0, 0.0]
# last time serotonin (metabolite 121) was noisy, now -
print('serotonin', noiseold0[120], oldtmp1[120][0]) # serotonin 1.406580269491114 0.021594036 - yey
# # sort the metabolites by noise, nameslst is the list of names
polarnameslst = polarnames.tolist()
lipidnameslst = lipidnames.tolist()
nameslst = polarnameslst + lipidnameslst
print('metabolites\' names', len(nameslst), nameslst[0], nameslst[139], nameslst[140], nameslst[-1]) # 434 1,3-Dimethyluracil Xylitol AcCa(16:0) TG(20:1_18:1_18:2) - yey
# old (20 months animals)
print(nameslst[344]) # PE(32:1) - yey
dctnoiseold = generatedct(nameslst, noiseold0)
print(dctnoiseold['PE(32:1)']) # 2.616972086836892 - yey
# sorting
srtdctnoiseold = sorted(dctnoiseold.items(), key=lambda x: x[1], reverse=True)
print(srtdctnoiseold[0]) # () - yey
print(srtdctnoiseold[0:5]) #
# [('PE(32:1)', 2.616972086836892), ('Chlorogenic acid', 2.4400295715742764), ('O-Phosphoethanolamine', 2.090977909717578), ('PE(34:2e)', 1.6305181820394963), ('Cer(t17:1_23:0)', 1.5792544543873577)]
# comparing polar to lipid noise
medianpolarold = statistics.median(noiseold0[:140])
medianlipidold = statistics.median(noiseold0[140:])
medianmergeold = statistics.median(noiseold0)
print('medians old noise, polar vs. lipid', medianpolarold, medianlipidold) # 0.3604356973962525 0.3200333654988953
print('medians old noise, merged', medianmergeold) # 0.3285230351687166

# calculating mean, std for youngtmp1 - the metabolomics of young (8 months) layers
print('young ', len(youngtmp1[0]), len(youngtmp1), youngtmp1[0][1]) # 15 434 0.002295563 - yey
meanyoung0 = []
stdyoung0 = []
for ii in youngtmp1:
    meanyoung0.append(statistics.mean(ii))
    stdyoung0.append(statistics.pstdev(ii))
print(len(meanyoung0), len(stdyoung0), meanyoung0[0], stdyoung0[0]) # 434 434 0.0025791786 0.0007090162115257826
# # verify
# print(youngtmp1[0])
# # [0.001339296, 0.002295563, 0.002089409, 0.003943424, 0.001694564, 0.002791911, 0.002133811, 0.003432982, 0.00314774, 0.003139989, 0.003359987, 0.001970572, 0.002633137, 0.001940667, 0.002774627]
# # mean, (population) std according to calculator.net - 0.0025791786, 0.00070901621152578 - yey!
# # noise
noiseyoung0 = getnoise(stdyoung0, meanyoung0)
dctnoiseyoung = generatedct(nameslst, noiseyoung0)
srtdctnoiseyoung = sorted(dctnoiseyoung.items(), key=lambda x: x[1], reverse=True)
print(srtdctnoiseyoung[0]) # ('Ascorbic acid', 3.741657386773941)
print(srtdctnoiseyoung[0:5]) # PE(32:1) - #3, Chlorogenic acid - #4,
# [('Ascorbic acid', 3.741657386773941), ("Uridine 5'-monophosphate (5'-UMP)", 3.741657386773941), ('PE(32:1)', 2.573072728631829), ('Chlorogenic acid', 2.154764637849187), ('PE(16:0p_22:5)', 2.0793479295802437)]
# comparing polar to lipid noise for youngs
medianpolaryoung = statistics.median(noiseyoung0[:140])
medianlipidyoung = statistics.median(noiseyoung0[140:])
medianmergeyoung = statistics.median(noiseyoung0)
print('medians young noise, polar vs. lipid', medianpolaryoung, medianlipidyoung) # 0.31544633590536497 0.25547051251099906
print('medians young noise, merged', medianmergeyoung) # 0.27458448113092015


# noise for 24-month-old group
print('AL ', len(altmp1[0]), len(altmp1), altmp1[0][1]) # 15 434 0.003053628
meanal0 = []
stdal0 = []
for ii in altmp1:
    meanal0.append(statistics.mean(ii))
    stdal0.append(statistics.pstdev(ii))
print(len(meanal0), len(stdal0), meanal0[0], stdal0[0]) # 434 434 0.0027756481333333334 0.001440753181040892
#
# # noise
noiseal0 = getnoise(stdal0, meanal0)
dctnoiseal = generatedct(nameslst, noiseal0)
srtdctnoiseal = sorted(dctnoiseal.items(), key=lambda x: x[1], reverse=True)
print(srtdctnoiseal[0]) # ('Shikimic acid', 3.7416573867739418)
print(srtdctnoiseal[0:5]) #
# [('Shikimic acid', 3.7416573867739418), ('Quinolinate', 3.225309094188683), ('Cer(t18:1_26:0)', 2.9857612294991003), ('PE(34:2e)', 2.7340710087326188), ('PE(32:1)', 2.5751117291565166)]
# comparing polar to lipid noise for AL
medianpolaral = statistics.median(noiseal0[:140])
medianlipidal = statistics.median(noiseal0[140:])
medianmergeal = statistics.median(noiseal0)
print('medians AL noise, polar vs. lipid', medianpolaral, medianlipidal) # 0.37728186732191554 0.37545131606213256
print('medians AL noise, merged', medianmergeal) # 0.375682973867254

# noise for TRF group
print('CR ', len(crtmp1[0]), len(crtmp1), crtmp1[0][1]) # 14 434 0.004467691
meancr0 = []
stdcr0 = []
for ii in crtmp1:
    meancr0.append(statistics.mean(ii))
    stdcr0.append(statistics.pstdev(ii))
print(len(meancr0), len(stdcr0), meancr0[0], stdcr0[0]) # 434 434 0.0029796691428571427 0.0015460467837488366
#
# # noise
noisecr0 = getnoise(stdcr0, meancr0)
dctnoisecr = generatedct(nameslst, noisecr0)
srtdctnoisecr = sorted(dctnoisecr.items(), key=lambda x: x[1], reverse=True)
print(srtdctnoisecr[0]) # ('Cer(d18:0_18:0)', 3.6055512754639896)
print(srtdctnoisecr[0:5]) #
# [('Cer(d18:0_18:0)', 3.6055512754639896), ('Cer(t18:1_26:0)', 3.6038375540435927), ('Shikimic acid', 2.4776036064743443), ('Chlorogenic acid', 2.2687588665755682), ('PE(32:1)', 2.0180666803139107)]
# comparing polar to lipid noise for CR
medianpolarcr = statistics.median(noisecr0[:140])
medianlipidcr = statistics.median(noisecr0[140:])
medianmergecr = statistics.median(noisecr0)
print('medians CR noise, polar vs. lipid', medianpolarcr, medianlipidcr) # 0.2642105513892977 0.29127257020496633
print('medians CR noise, merged', medianmergecr) # 0.28561794525565576 - lower than AL!!

# the noise increases monotonically with age at the AL regime. TRF is lower than the AL counterpart,
# meaning TRF animals are physiologically younger!
# noise summary: [Old, Young, AL, CR] = [0.3285230351687166, 0.27458448113092015, 0.375682973867254, 0.28561794525565576]

# need a p-value to this, via permuting the labels, then check sum of differences >= to non permuted.
diff00 = abs(medianmergeold - medianmergeyoung) + abs(medianmergeold - medianmergeal)
print(diff00) # 0.10109849273633387 - yey
# converting to floats
mergelst0val = strtofloat0(mergelst0)
mergelst045 = transform0(mergelst0val, 45) # returns a list of 45 animals, each with its 434 metabolite values
print(len(mergelst045), len(mergelst045[0]), mergelst045[0][0], mergelst045[0][-1]) # 45 434 0.000390325 (O1 first polar) 72051.23054 (O1 last lipid) - yey
perm00 = np.random.permutation(mergelst045)
noisepermold = calcnoise(perm00[0:15])
medianpermnoiseold = statistics.median(noisepermold)
noisepermyoung = calcnoise(perm00[15:30])
medianpermnoiseyoung = statistics.median(noisepermyoung)
noisepermal = calcnoise(perm00[30:45])
medianpermnoiseal = statistics.median(noisepermal)
print(medianpermnoiseold, medianpermnoiseyoung, medianpermnoiseal) # 0.3283630083249569 0.3859786513345388 0.3617745837771553

# QC
oldtmp1QC = np.array(oldtmp1).transpose().tolist()
noiseoldQC = calcnoise(oldtmp1QC)
medianoldnoiseQC = statistics.median(noiseoldQC)
print('QC calcnoise old', medianoldnoiseQC) # 0.3285230351687166 - yey


mergelst059 = transform0(mergelst0val, 59) # returns a list of 59 animals, each with its 434 metabolite values
print(len(mergelst059), len(mergelst059[0]), mergelst059[0][0], mergelst059[0][-1], mergelst059[-1][-1]) # 59 434 0.000390325 72051.23054 85901.02017 - yey
# ind = 0
# for ii in range (0, 10): # 1000
#     perm00 = np.random.permutation(mergelst059)
#     noisepermold = calcnoise(perm00[0:15])
#     medianpermnoiseold = statistics.median(noisepermold)
#     noisepermyoung = calcnoise(perm00[15:30])
#     medianpermnoiseyoung = statistics.median(noisepermyoung)
#     noisepermal = calcnoise(perm00[30:45])
#     medianpermnoiseal = statistics.median(noisepermal)
#     noisepermcr = calcnoise(perm00[45:])
#     medianpermnoisecr = statistics.median(noisepermcr)
#     diffsumperm = calcsum0(medianpermnoiseold, medianpermnoiseyoung, medianpermnoiseal)
#     if (diffsumperm >= diff00) and (medianpermnoisecr < medianpermnoiseal):
#         ind = ind + 1




old00 = noiseold0
young00 = noiseyoung0
al00 = noiseal0
cr00 = noisecr0
fig, ax = plt.subplots()
bins = np.linspace(0, 1.5, 75)

plt.hist(young00, bins, alpha=0.5, color = 'blue', label='Young')
plt.hist(old00, bins, alpha=0.5, color = 'red', label='Old')
plt.hist(al00, bins, alpha=0.5, color = 'green', label='Oldest')
plt.hist(cr00, bins, alpha=0.5, color = 'orange', label='TF')
plt.axvline(x=medianmergeold, ymin=0, ymax=1, color = 'red', linewidth = 4)
plt.axvline(x=medianmergeyoung, ymin=0, ymax=1, color = 'blue', linewidth = 4)
plt.axvline(x=medianmergeal, ymin=0, ymax=1, color = 'green', linewidth = 4)
plt.axvline(x=medianmergecr, ymin=0, ymax=1, color = 'orange', linewidth = 4)
ax.tick_params(axis='x', labelsize=22)
ax.tick_params(axis='y', labelsize=22)


# plt.savefig('noise YOALTRF.pdf') # YOALTRF ALTRF
plt.show()



# plotting with error bars
# reminder - noise is already defined
# old00 = noiseold0
# young00 = noiseyoung0
# al00 = noiseal0
# cr00 = noisecr0
samplesize0 = len(old00)
print(samplesize0) # 434 - yey
stdevold0 = statistics.pstdev(old00)
stdevyoung0 = statistics.pstdev(young00)
stdeval0 = statistics.pstdev(al00)
stdevcr0 = statistics.pstdev(cr00)
errorold0 = stdevold0/math.sqrt(samplesize0)
erroryoung0 = stdevyoung0/math.sqrt(samplesize0)
erroral0 = stdeval0/math.sqrt(samplesize0)
errorcr0 = stdevcr0/math.sqrt(samplesize0)

# plotting - https://stackoverflow.com/questions/22364565/python-pylab-scatter-plot-error-bars-the-error-on-each-point-is-unique
xcoor = [8, 20, 24, 24]
ycoor = [medianmergeyoung, medianmergeold, medianmergeal, medianmergecr]
error00 = [erroryoung0, errorold0, erroral0, errorcr0]
fig, ax = plt.subplots()
# plt.errorbar(xcoor, ycoor, yerr=error00, fmt='ro', capsize = 3) # plotting all data at once, no loop is needed
for ii in range(0, 4):
    if ii < 3:
        plt.errorbar(xcoor[ii], ycoor[ii], yerr=error00[ii], fmt='bo', capsize = 3)
    elif ii == 3:
        plt.errorbar(xcoor[ii], ycoor[ii], yerr=error00[ii], fmt='ro', capsize=3)
ax.tick_params(axis='x', labelsize=22)
ax.tick_params(axis='y', labelsize=22)

# plt.savefig('noise medians.pdf')
plt.show() # looking good




