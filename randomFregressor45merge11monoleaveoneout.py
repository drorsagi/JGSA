
# calculates the predicted age of the 45 AL hens based on real age data using leave one out cross validation
# and rndom forest regressor/regression



def rmse(score):
    rmse = np.sqrt(-score)
    print(f'rmse= {"{:.2f}".format(rmse)}')

def get_whisker(tmp,dataset):
    whisker = []
    for quantile,data in zip(tmp,dataset):
        data = np.array(data)
        q1 = quantile[0]
        median = quantile[1]
        q3 = quantile[2]
        iqr = q3 - q1
        upper = q3 + 1.5 * iqr
        upper = np.clip(upper,q3,data.max())
        lower = q1 - 1.5 * iqr
        lower = np.clip(lower,data.min(),q1)
        whisker.append((upper,lower))
    return whisker


import matplotlib.pyplot as plt
from sklearn import metrics
import statistics
import pandas
import numpy as np
import numpy.random
from sklearn import model_selection
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.model_selection import cross_val_score
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import cross_validate
from sklearn import datasets
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
from scipy.stats import mannwhitneyu


datasetF = pandas.read_csv('markersOYALmergemonotonicMLrealage', header = None) # real age as age group, good for regression estimates

array = datasetF.values
X = array[:, 0:-1]

Y = array[:, -1]

rf = RandomForestRegressor(n_estimators = 300, random_state = 0) # n_estimators = 1000, random_state = 42


from sklearn.model_selection import LeaveOneOut

# create loocv procedure
cvout = LeaveOneOut()
# # enumerate splits
y_true, y_pred = list(), list()
for train_ix, test_ix in cvout.split(X):
    print(test_ix)
    # split data
    X_train, X_test = X[train_ix, :], X[test_ix, :]
    y_train, y_test = Y[train_ix], Y[test_ix]
    # fit model
    model = RandomForestRegressor(random_state=1)
    model.fit(X_train, y_train)
    # evaluate model
    yhat = model.predict(X_test)
    # store
    y_true.append(y_test[0])
    y_pred.append(yhat[0])
# calculate accuracy
# acc = accuracy_score(y_true, y_pred)
# print('Accuracy: %.3f' % acc)
# # Accuracy: 0.867
print('true, predictions', y_true, y_pred)
# true - [20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0]
# predictions - [17.96, 15.92, 18.56, 17.76, 16.2, 19.4, 20.04, 20.24, 14.4, 19.16, 17.0, 11.2, 21.8, 19.76, 20.08, 8.0, 10.04, 10.04, 8.12, 8.76, 18.84, 8.36, 9.48, 8.6, 19.32, 10.64, 9.32, 8.36, 8.48, 10.92, 15.4, 20.56, 23.8, 23.36, 23.24, 24.0, 23.84, 23.32, 22.96, 24.0, 23.2, 24.0, 23.64, 23.56, 23.96]

# test MW significance between AL24 months to TRF (24 month old)
al24 = [15.4, 20.56, 23.8, 23.36, 23.24, 24.0, 23.84, 23.32, 22.96, 24.0, 23.2, 24.0, 23.64, 23.56, 23.96]
trf24 = [22.02666667, 22.18666667, 20.88, 21.92, 12.37333333, 11.02666667, 11.29333333, 22.85333333, 23.86666667, 21.52, 22.58666667, 21.65333333, 23.86666667, 21.85333333]
mwsig0 = mannwhitneyu(al24, trf24, alternative='two-sided')
print(mwsig0) # MannwhitneyuResult(statistic=170.0, pvalue=0.004851291024865207)

# plotting violin plot of the data
dataset0 = [al24, trf24]
fig, ax = plt.subplots()

pos = [1, 2]
label0 = ['AL', 'TRF']
tmp0 = [np.percentile(data,[25,50,75]) for data in dataset0]
whisker = get_whisker(tmp0,dataset0)

vp = ax.violinplot(dataset0, showmeans = True, positions=[1, 2])

for body in vp['bodies']:
    body.set_facecolor('red')
    body.set_edgecolor('black')
    body.set_alpha(1)
vp['cmaxes'].set_color('black')
vp['cmins'].set_color('black')
vp['cbars'].set_color('black')

ax.scatter(pos,[quantile[1] for quantile in tmp0], marker='o',color='white',s=30,zorder=3)
ax.vlines(pos,[quantile[0] for quantile in tmp0],[quantile[2] for quantile in tmp0],color='black',linestyle='-',lw=5)
ax.vlines(pos,[bound[0] for bound in whisker],[bound[1] for bound in whisker],color='black',linestyle='-',lw=2)

# Add title
ax.set_title('Violin Plot')
ax.set_xticks(pos)
ax.set_xticklabels(label0)


plt.show() # yey

# plot all 59 points as scatter plot, true age vs. predicted age
colorskey = np.array(['blue', 'red', 'green', 'black'])
age8 = [8]*15
age20 = [20]*15
age24 = [24]*29

al20 = [17.96, 15.92, 18.56, 17.76, 16.2, 19.4, 20.04, 20.24, 14.4, 19.16, 17.0, 11.2, 21.8, 19.76, 20.08]
al8 = [8.0, 10.04, 10.04, 8.12, 8.76, 18.84, 8.36, 9.48, 8.6, 19.32, 10.64, 9.32, 8.36, 8.48, 10.92]
med8 = statistics.median(al8)
med20 = statistics.median(al20)
med24 = statistics.median(al24)
medtrf = statistics.median(trf24)
medians0 = [med8, med20, med24, medtrf]
agemeds = [8, 20, 24, 24]
print(medians0) # [9.32, 18.56, 23.56, 21.886666665]
fig = plt.figure()

younglay = plt.scatter(age8, al8, marker='o', color='gold')
oldlay = plt.scatter(age20, al20, marker='o', color='gold')
oldestlay = plt.scatter(age24[0:15], al24, marker='o', color='gold')
trflay = plt.scatter(age24[15:], trf24, marker='p', color='green')
mediansplt = plt.scatter(agemeds[:3], medians0[:3], marker='x', s=100, color='black')
medianstrf = plt.scatter(agemeds[3], medians0[3], marker='x', s=100, color='magenta')


plt.show()





