
# at the bottom we have also leave one out cross validation


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
from sklearn.model_selection import KFold
import math

datasetF = pandas.read_csv('markersOYALmergemonotonicML', header = None)

array = datasetF.values
X = array[:, 0:-1]

Y = array[:, -1]

classifier = RandomForestClassifier(random_state=1) # random_state=seed

# kfold cross validation:
kf =KFold(n_splits=10, shuffle=True, random_state=1)

score1 = cross_val_score(classifier, X, Y, cv= kf, scoring="accuracy")
print('Accuracy: %.3f (%.3f)' % (statistics.mean(score1), statistics.pstdev(score1)))

# Accuracy: 0.865 (0.112)

# leaveoneout

from sklearn.model_selection import LeaveOneOut

cvout = LeaveOneOut()
# # enumerate splits
y_true, y_pred = list(), list()
for train_ix, test_ix in cvout.split(X):
    # print(test_ix)
    # split data
    X_train, X_test = X[train_ix, :], X[test_ix, :]
    y_train, y_test = Y[train_ix], Y[test_ix]
    # fit model
    model = RandomForestClassifier(random_state=1)
    model.fit(X_train, y_train)
    # evaluate model
    yhat = model.predict(X_test)
    # store
    y_true.append(y_test[0])
    y_pred.append(yhat[0])
# calculate accuracy
acc = accuracy_score(y_true, y_pred)
print('Accuracy: %.3f' % acc)
# Accuracy: 0.867
print('true, predictions', y_true, y_pred)
# true - [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
# pred - [1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 3.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.0, 2.0, 2.0, 2.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
