# JGSA
python codes for first paper
Summary for python codes for Journal of Gerontology Series A paper (JGSA)
Ransom forest algorithm:
Regressor to predict real age of hens: 
Violin plot of TRF vs 24-month-old AL, and all 59 ages: All 45 AL values were calculated using leaveoneout model, based on the 11 aging biomarkers.
randomFregressor45merge11monoleaveoneout
Prediction of TRF ages: randomFregressor59merge11monorealage
The regressor uses all 45 AL dataset to generate the predictive model and then predicts the ages of TRF individuals based on the 11 aging biomarkers.
Classifier to predict the correct age group of AL individuals:
randomF45merge11monokfoldcross: uses 10fold cross validation as well as leaveoneout model to predict accuracy of the 11 aging biomarkers. Accuracy was 0.865 and 0.867 respectively. Leaveoneout (L-O-O) also matched for each individual hen its predicted age group. L-O-O predictions were in agreement with my method.
randomF45merge11mono/monoplot: calculates the accuracy of the 11 biomarkers using 25 different seeds, and 9 animals in each validation cohort. 6 hens were mispredicted (11%), same hens as in the L-O-O model. The monoplot version plots a typical confusion matrix and 3 ROC curves for the three age groups. 
Figures:
Fig1: PCA: pca merge1 – PCA of all 59 animals under all metabolites (sup Fig1) and 45 AL animals under the 11 aging biomarkers. 
pca merge1 eleven mono AYALCR: plots the PCA of the 59 animals showing a younger profile for the TRF group under the 11 aging biomarkers.
Fig 2: arrangedata merge59 polarlipid pval noise – calculates (and plots) the noise for all 4 groups, and for the 24-month-old AL vs TRF. Also plots the median noise values for each group. 
metrics merge averagedist v2 for noise – calculates (and plots) the average ingroup distance as a function of age, showing the it increases with age, but TRF are more compact with respect to AL counterparts.
Fig3: noise for human, noise in huma blood yanagida, noise for mice.
Fig 4: egghunt v2.
Supp 1: pca merge1, seabornclustermap (the distance matrix is built in metrics merge averagedistmtrx),  metrics merge averagedist v3.
Supp 2: arrangedata merge59 polarlipid bargraph (effect size for u-values).
Supp 3: corresponding randomF files (mono plot and regressor45...leaveoneout).
Supp 4: monotonic markers plottall (medians of 11 markers).
Supp 5: randomwalk, randomwalk personal.
