





from itertools import cycle

# from sklearn import svm, datasets

# from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
# from scipy import interp - use np.interp!!!!!!!!!!! (line 131)
# import scipy.interpolate # https://stackoverflow.com/questions/45634350/missing-interpolate-in-scipy-0-17/45634483
from sklearn.metrics import roc_auc_score




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


datasetF = pandas.read_csv('markersOYALmergemonotonicML', header = None) # yey


array = datasetF.values
X = array[:, 0:-1]
# Y1 = array[:, -1]
# Y = np.random.permutation(Y1)
Y1 = array[:, -1] # run in paralel the standard code
Y = array[:, -1]
# Binarize the output
Y = label_binarize(Y, classes=[1, 2, 3])
n_classes = Y.shape[1]
validation_size = 0.2 #
seed = 1 # 10 #1 # 2 # 13 is roc = 1! # 44, 12, 1, 14, 15 # 101 # np.random.seed(151) #
X_train, X_validation, Y_train, Y_validation = model_selection.train_test_split(X, Y, test_size=validation_size, random_state=seed)
X1_train, X1_validation, Y1_train, Y1_validation = model_selection.train_test_split(X, Y1, test_size=validation_size, random_state=seed)


classifier = RandomForestClassifier(n_estimators=300, random_state=0) # random_state=seed
model21 = classifier.fit(X1_train, Y1_train)
preds21 = classifier.predict(X1_validation)

classifier1 = OneVsRestClassifier(RandomForestClassifier(n_estimators=300, random_state=0))
y_score = classifier1.fit(X_train, Y_train).predict(X_validation) #decision_function(X_validation)

fpr = dict()
tpr = dict()
roc_auc = dict()
for i in range(n_classes):
    fpr[i], tpr[i], _ = roc_curve(Y_validation[:, i], y_score[:, i])
    roc_auc[i] = auc(fpr[i], tpr[i])

# Compute micro-average ROC curve and ROC area
fpr["micro"], tpr["micro"], _ = roc_curve(Y_validation.ravel(), y_score.ravel())
roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])


print('WTF', preds21, type(X1_validation), Y1_validation )
print(confusion_matrix(Y1_validation,preds21))
print(accuracy_score(Y1_validation, preds21))

plot_confusion_matrix(classifier, X1_validation, Y1_validation, cmap=plt.cm.Blues, normalize = 'all') # https://scikit-learn.org/stable/auto_examples/model_selection/plot_confusion_matrix.html

plt.show() #



##############################################################################
# Plot ROC curves for the multilabel problem
# ..........................................
# Compute macro-average ROC curve and ROC area

# First aggregate all false positive rates
all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))

# Then interpolate all ROC curves at this points
mean_tpr = np.zeros_like(all_fpr)
for i in range(n_classes):
    mean_tpr += np.interp(all_fpr, fpr[i], tpr[i]) # interp(all_fpr, fpr[i], tpr[i])

# Finally average it and compute AUC
mean_tpr /= n_classes

fpr["macro"] = all_fpr
tpr["macro"] = mean_tpr
roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

# Plot all ROC curves
# plt.figure()
lw = 2
fig, ax = plt.subplots()
plt.plot(fpr["micro"], tpr["micro"],
         label='micro-average ROC curve (area = {0:0.2f})'
               ''.format(roc_auc["micro"]),
         color='deeppink', linestyle=':', linewidth=4)

plt.plot(fpr["macro"], tpr["macro"],
         label='macro-average ROC curve (area = {0:0.2f})'
               ''.format(roc_auc["macro"]),
         color='navy', linestyle=':', linewidth=4)

colors = cycle(['aqua', 'darkorange', 'cornflowerblue'])
for i, color in zip(range(n_classes), colors):
    plt.plot(fpr[i], tpr[i], color=color, lw=lw,
             label='ROC curve of class {0} (area = {1:0.2f})'
             ''.format(i, roc_auc[i]))

plt.plot([0, 1], [0, 1], 'k--', lw=lw)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate', size = 22)
plt.ylabel('True Positive Rate', size = 22)
ax.tick_params(axis='x', labelsize=22)
ax.tick_params(axis='y', labelsize=22)

plt.legend(loc="lower right")
# fig.savefig('multi-class roc curve.pdf')
plt.show()









