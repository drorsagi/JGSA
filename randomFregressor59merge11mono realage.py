

# this code tries to predict the age group of the CR/TRF group based on the three AL using random forest resression
# groups, based on the 11 markers. as age we use the real age (8,20,24) of the 4 groups




import math
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
from sklearn.metrics import plot_confusion_matrix
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error # https://www.kite.com/python/answers/how-to-take-root-mean-square-error-(rmse)-in-python


datasetF = pandas.read_csv('markersOYALCRmergemonotonicrealage', header = None) # yey

array = datasetF.values
X = array[:, 0:-1]

Y = array[:, -1]

X1_train = X[:45, :]
Y1_train = Y[:45]

X1_validation = X[45:, :]
Y1_validation = Y[30:44] # all TRF ages are 24


# using random forest regression

rf = RandomForestRegressor(n_estimators = 300, random_state = 0)
rf.fit(X1_train, Y1_train)
predictions = rf.predict(X1_validation)
print('regressor', predictions, type(X1_validation), Y1_validation )
# regressor [22.02666667 22.18666667 20.88       21.92       12.37333333 11.02666667
#  11.29333333 22.85333333 23.86666667 21.52       22.58666667 21.65333333
#  23.86666667 21.85333333] <class 'numpy.ndarray'> [24. 24. 24. 24. 24. 24. 24. 24. 24. 24. 24. 24. 24. 24.]
# Always younger! only two animals appear in the same chronological age as the control AL 24 months old - yey
# 3 around 10, 5 more under 22, 6 above 22.
print(statistics.median(predictions)) # 21.886666666666667 - 10% younger.









