



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

datasetF = pandas.read_csv('markersOYALmergemonotonicML', header = None) # yey

array = datasetF.values
X = array[:, 0:-1]
# Y1 = array[:, -1]
# Y = np.random.permutation(Y1)
Y = array[:, -1]
validation_size = 0.2 #
seed = 1 # 44, 12, 1, 14, 15 # 101 # np.random.seed(151) #
X_train, X_validation, Y_train, Y_validation = model_selection.train_test_split(X, Y, test_size=validation_size, random_state=seed)

# print(X_train[:, 0:3])

# summary - consistently wrong: AL = [L1, L66],  O = [O14, O15], Y = [Y111, Y115] - 6/45 wrong
# O14 - old predicted young, O15 - old predicted older (24 months)
# Y111, 115 - young predicted old, L1, L66- older predicted old (once even L1 predicted young)




classifier = RandomForestClassifier(n_estimators=300, random_state=0) # random_state=seed
model21 = classifier.fit(X_train, Y_train)
preds21 = classifier.predict(X_validation)
print('WTF', preds21, type(X_validation), Y_validation )
print(confusion_matrix(Y_validation,preds21))
print(accuracy_score(Y_validation, preds21))



# # let's systematically try 20 seeds
# print('loop')
# accur0 = []
# for ii in range (0, 20):
#     seed = ii
#     X_train, X_validation, Y_train, Y_validation = model_selection.train_test_split(X, Y, test_size=validation_size, random_state=seed)
#
#     # print(X_train[:, 0:3])
#
#     # # Feature Scaling
#     # from sklearn.preprocessing import StandardScaler
#     #
#     # sc = StandardScaler()
#     # X_train = sc.fit_transform(X_train)
#     # X_validation = sc.transform(X_validation)
#
#      # X_train = X[0:4, :]
#     # print(X_train[:, 0:3]) # yey
#     classifier = RandomForestClassifier(n_estimators=300, random_state=0)  # random_state=seed
#     model21 = classifier.fit(X_train, Y_train)
#     preds21 = classifier.predict(X_validation)
#     print('WTF', preds21, type(X_validation), Y_validation)
#     print(confusion_matrix(Y_validation, preds21))
#     print(accuracy_score(Y_validation, preds21))
#     accurtmp = accuracy_score(Y_validation, preds21)
#     accur0.append(accurtmp)
# print(accur0, statistics.median(accur0), statistics.mean(accur0)) # 0.8888888888888888 0.8444444444444444
# # [0.3333333333333333, 0.8888888888888888, 0.8888888888888888, 0.8888888888888888, 0.8888888888888888, 0.8888888888888888, 0.7777777777777778, 0.7777777777777778, 0.8888888888888888, 0.8888888888888888, 0.7777777777777778, 0.8888888888888888, 0.8888888888888888, 1.0, 1.0, 1.0, 0.8888888888888888, 0.8888888888888888, 0.6666666666666666, 0.7777777777777778] 0.8888888888888888 0.8444444444444444








