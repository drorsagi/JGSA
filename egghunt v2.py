



def unify(lst1, names1, lst2, names2): # makes one list of lst1 and lst2, with the same chicken index
    ans = []
    for ii in range (0, len(lst1)):
        ind0 = names2.index(names1[ii])
        tmp = lst1[ii] + lst2[ind0]
        ans.append(tmp)
    return ans



def converttovalues(arr): # arr = array, each element is a list, each element in the list is a string like '31'
    # gp1tCmkr1val = [float(ii) for ii in gp1tCmkr1]
    values0 = []
    for ii in arr:
        tmp = [float(jj) for jj in ii]
        values0.append(tmp)
    return values0



def getnoise(std0, mean0):
    noise0 = []
    for ii in range(0, len(std0)):
        if mean0[ii] != 0:
            noise0.append(std0[ii]/mean0[ii])
        elif mean0[ii] == 0:
            noise0.append(0)
    return noise0

def permutarr0(array0): # permuting the columns of an array, returns a nested list of the column-permutated array
    permcol = []
    for ii in range(0, array0.shape[1]): # array.shape returns (#row, #columns)
        tmp = np.random.permutation(array0[:, ii])
        permcol.append(tmp)
    return permcol

def fraction0(permutationlst, noise):
    pval0 = []
    for jj in range(0, len(permutationlst[0])):
        ind0 = 0
        for ii in permutationlst:
            # if ii[jj] >= noise[jj]: # identifies the most noisiest!
            if ii[jj] <= noise[jj]: # identifies the least noisiest!
                ind0 = ind0 + 1
        pval0.append((ind0 + 1) / 1001)
    return pval0


from sklearn.linear_model import LinearRegression
import statsmodels
import statsmodels.stats.multitest
import statistics
from scipy.stats import ttest_ind
from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


namestmp = []
ff0 = open('egg last three months AL names.csv')
for line in ff0:
    tmp = line.split(',')
    namestmp.append(tmp[0])
ff0.close()
print(len(namestmp)) # 28, this includes the header
# elminating the header
names0 = namestmp[1:]
print(len(names0), names0[0], names0[26]) # 27, 32, 135 - yey

values = []
nameorder = []
ff1 = open('egg summary1.csv')
for line in ff1:
    tmp = line.split(',')
    if tmp[0] in names0:
        nameorder.append(tmp[0])
        tmp1 = tmp[-1].rstrip('\n')
        tmp[-1] = tmp1
        values.append(tmp[1:]) # number of eggs for each of the first year

ff1.close()


tmpstrvalues = np.array(values)
valuesarr = converttovalues(tmpstrvalues) # valuesarr is actually a nested list, but np.mean/std can handle it
chickmean = np.mean(valuesarr, axis = 1)
print(chickmean[0]) # 29.0769230769 - yey!! same as the excel file
timemean = np.mean(valuesarr, axis = 0)
print(timemean[0]) # 12.7407407407 - correct, the first month, april, is only 23 days
print(timemean) # the surviving animals are above average
# [ 12.74074074  30.33333333  29.66666667  30.25925926  29.7037037
#   28.40740741  28.77777778  27.51851852  29.07407407  30.07407407
#   26.48148148  29.25925926  27.62962963]

# look at last three months - the survivers
values1 = []
nameorder1 = []
ii = 0
ff2 = open('last three months eggs.csv')
for line in ff2:
    tmp = line.split(',')
    if ii > 0:
        values1.append(tmp[0:-1]) # number of eggs
        tmp1 = tmp[-1].rstrip('\n')
        nameorder1.append(tmp1) # chicken indices
    ii = ii + 1
ff2.close()

# now we have to unify the lists, where each line is the smae animal
newstrvalues = unify(values, nameorder, values1, nameorder1)
print(len(newstrvalues), newstrvalues[0]) # 27 ['13', '31', '30', '31', '31', '29', '32', '29', '31', '31', '28', '32', '30', '30', '31', '30'] - yey
tmpstrvaluesunify = np.array(newstrvalues)
valuesunifyarr = converttovalues(tmpstrvaluesunify) # returns a nested list
chicktotmean = np.mean(valuesunifyarr, axis = 1)
print(chicktotmean[0]) # 29.3125
timetotmean = np.mean(valuesunifyarr, axis = 0)
print(timetotmean[0]) # 12.740
print(timetotmean)
# [12.74074074 30.33333333 29.66666667 30.25925926 29.7037037  28.40740741
#  28.77777778 27.51851852 29.07407407 30.07407407 26.48148148 29.25925926
#  27.62962963 25.07407407 26.48148148 25.92592593] - looking good
# standard deviation, then noise
timetotstd = np.std(valuesunifyarr, axis = 0)
print(timetotstd)
# [0.58325984 1.15470054 0.90267093 0.92666637 1.46096911 2.00479535
#  2.84583299 2.09709447 1.90371811 0.93989463 1.64137184 2.1010155
#  3.47635556 4.56172422 3.88129771 3.8386526 ] - looking good
noisetimetot = getnoise(timetotstd, timetotmean)
print(noisetimetot)
# [0.0457791155465803, 0.0380670507157995, 0.0304271101297227, 0.030624225222024005, 0.04918474544467834, 0.07057297837144386, 0.09888994961286252, 0.07620666307882674, 0.06547820267972687, 0.031252654040516054, 0.061981873583665296, 0.07180685875282448, 0.12581983946656636, 0.18192991713680548, 0.14656648678387776, 0.14806231451517918]
# plot the noise
age0 = []
for ii in range(0, 17): # age in months, starts from months 0, which is end of apr 2018, age 8 month (sep-apr), May 2019 does not exist in the data
    if ii != 13:
        age0.append(ii + 8)

plt.scatter(age0[:], noisetimetot[:], c='b', marker='x', label='noise with age')

plt.xlabel('age')
plt.ylabel('noise')
plt.legend(loc='upper left')
plt.show()
# how good is linear regression?
corr0 = pearsonr(age0[:], noisetimetot[:]) # pearsonr(age0[1:], noisetimetot[1:])
print('pearsonr=', corr0) # pearsonr= (0.814465079087986, 0.00012233384888687764),

# linear regression
age1 = np.array(age0[:])
noisetimetot1 = np.array(noisetimetot[:])
# print('linear regression')
X = age1.reshape(-1, 1) # values converts it into a numpy array
Y = noisetimetot1.reshape(-1, 1) # -1 means that calculate the dimension of rows, but have 1 column
linear_regressor = LinearRegression()  # create object for the class
linear_regressor.fit(X, Y)  # perform linear regression
Y_pred = linear_regressor.predict(X)  # make predictions
#To retrieve the intercept:
print(linear_regressor.intercept_) # [-0.04066762], no april [1:] - [-0.05283488]
#For retrieving the slope:
print(linear_regressor.coef_) # 0.00766268, no april [1:] - 0.00831026
plt.scatter(X, Y)
plt.plot(X, Y_pred, color='red')
# plt.savefig('noise in egg layed.pdf')
plt.show()
print('regression done')

