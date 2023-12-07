# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 16:44:11 2023

@author: trpen
"""


pip install GEOparse


####################
## Load libraries ##
####################

import GEOparse as gp
import pandas as pd
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, auc, roc_curve, RocCurveDisplay, confusion_matrix
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import ComplementNB
from sklearn.neural_network import MLPClassifier
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

################################################
######            LOAD DATASETS           ######
################################################

###############################################################################
#########                          3 PARTS                             ########
#########   (1) Train on best genes from cortex, validate on spinal    ########
#########   (2) Train on best genes from spinal, validate on cortex    ########
#########   (3) Train on NovaSeq gene expression, validate on HiSeq    ########
#########                                                              ########
###############################################################################


##Part 1) Train on best genes from cortex, validate using expression in the spinal cord

## Load up data into pandas
#####  - Training = Gene expression data in cortex
###        -80% used for training, 20% used for test
###    - want to train on 80% because we want more training data to validate classifier on spinal cord data!
#  - Validation = Gene expression data in spinal cord 

##FIRST PART: spinal cord + cortex phenotype data
cortex = pd.read_csv("C:/Users/trpen/OneDrive/Desktop/ALS Classification/ML/Cortex_Phenotype.csv",header=0, index_col=0) 
spinal = pd.read_csv("C:/Users/trpen/OneDrive/Desktop/ALS Classification/ML/SpinalCord_Phenotype.csv",header=0, index_col=0) 

###############################################################################
#######look at phenotype info (platform, subtype, etc) from each dataset#######
###############################################################################
##Cortex, capitalize 'platform' in the actual csv file 
cortex_info = cortex[['Subject', 'Platform','Subtype','PLSSubtype']] #585 samples

##Spinal
spinal_info = spinal[['Subject', 'Platform','Subtype','PLSSubtype']] #519 samples


#Phenotype data:
#set samples as the index, then have subtypes as the phenotype data (kinda redundant), can skip this step
cortexPhenos = cortex_info[['Subject','Platform','PLSSubtype']]
cortexPhenos = cortexPhenos.set_index('Subject')
spinalPhenos = spinal_info[['Subject','Platform','PLSSubtype']]
spinalPhenos = spinalPhenos.set_index('Subject')

#determine number of ox vs notox subtypes
print(cortexPhenos['PLSSubtype'].value_counts()) #notox = 346, ox = 239
print(spinalPhenos['PLSSubtype'].value_counts()) #notox = 381, ox = 138


###############################################################################
#####Load Gene expression data for spinal & cortex
##       Train on cortex, validate on spinal
##       3 best genes

gexp3c_cortex = pd.read_csv("C:/Users/trpen/OneDrive/Desktop/ALS Classification/ML/HiSeqHoldouot_RPKM_Expression_CortexSamples_3Gene.csv",header=0, index_col=0) 
gexp3c_spinal = pd.read_csv("C:/Users/trpen/OneDrive/Desktop/ALS Classification/ML/HiSeqHoldouot_RPKM_Expression_SpinalSamples_3Gene.csv",header=0, index_col=0) 

##Training data = 3 gene expression in cortex
##Validation data = 3 gene expression in spinal cord

###############################################################################
#####                    TRAINING THE CLASSIFIERS                         #####
###############################################################################

#set data to train classifiers on, (two methods needed, because CNB method needs values [0-1], but since we are not using CNB we technically dont need the 'scaled_mm' part)
gexpTrain_c = gexp3c_cortex.T
gexpTrainc_scld = StandardScaler().fit_transform(gexpTrain_c)
gexpTrainc_scld_mm = MinMaxScaler().fit_transform(gexpTrainc_scld)

## Generate labels that are in the set ['Ox', 'NotOx']
# This will bulid a classifier that can discriminate Ox from the other subtypes
# uses OX label for OX subtype, anything not OX set to 'other'
CortexSubtypes = pd.Series(['OX' if i=='OX' else 'other' for i in cortexPhenos['PLSSubtype']])
print(CortexSubtypes.value_counts()) #should be same as before: 346 other, 239 OX

## Set up classification comparison for 100 fold cross-validation
# Can change kFolds to change number of iterations
# Store out results in a list nested inside a dictionary
kFolds = 100
performanceResults1 = {'KNN':[],'SVM':[],'LDA':[],'RF':[],'MLP':[]}
for i in range(kFolds):
    ## Split up data into train and test sets using random integer seed
    gexpData_train, gexpData_test, disState_train, disState_test = train_test_split(gexpTrainc_scld, CortexSubtypes, test_size=0.2, random_state=i) #20% of data is going into test set (test size)
    gexpData_train_mm, gexpData_test_mm, disState_train_mm, disState_test_mm = train_test_split(gexpTrainc_scld_mm, CortexSubtypes, test_size=0.2, random_state=i)
    ## Train and test the classifiers
    # k-nearest neighbors (KNN)
    knn_clf = KNeighborsClassifier(8)
    knn_clf.fit(gexpData_train, disState_train)
    knn_pred = knn_clf.predict(gexpData_test)
    performanceResults1['KNN'].append(classification_report(disState_test, knn_pred, output_dict=True, zero_division=0))
    # Support vector machines (SVM)
    svm_clf = SVC(kernel='linear', C = 0.025)
    svm_clf.fit(gexpData_train, disState_train)
    svm_pred = svm_clf.predict(gexpData_test)
    performanceResults1['SVM'].append(classification_report(disState_test, svm_pred, output_dict=True, zero_division=0))
    # Linear discriminant analysis (LDA)
    lda_clf = LinearDiscriminantAnalysis()
    lda_clf.fit(gexpData_train, disState_train)
    lda_pred = lda_clf.predict(gexpData_test)
    performanceResults1['LDA'].append(classification_report(disState_test, lda_pred, output_dict=True, zero_division=0))
    # Random Forest ensemble (RF)
    rf_clf = RandomForestClassifier()
    rf_clf.fit(gexpData_train, disState_train)
    rf_pred = rf_clf.predict(gexpData_test)
    performanceResults1['RF'].append(classification_report(disState_test, rf_pred, output_dict=True, zero_division=0))
    # # Complement Naive Bayes (CNB)
    # cnb_clf = ComplementNB()
    # cnb_clf.fit(gexpData_train_mm, disState_train)
    # cnb_pred = cnb_clf.predict(gexpData_test_mm)
    # performanceResults1['CNB'].append(classification_report(disState_test_mm, cnb_pred, output_dict=True, zero_division=0))
    #Multilayer perceptron (MLP)
    mlp_clf = MLPClassifier(hidden_layer_sizes=(100,100,100), alpha=1e-5, learning_rate_init=(0.0001))
    mlp_clf.fit(gexpData_train, disState_train)
    mlp_pred = mlp_clf.predict(gexpData_test)
    performanceResults1['MLP'].append(classification_report(disState_test, mlp_pred, output_dict=True, zero_division=0))


## Build figure to describe classifiers: using 'f1-score'
# f1 score = 2/(1/recall+1/precision) = tp/(tp + 0.5 * (fp + fn))
# tp = true positive, fp = false positive, tn = true negative, fn = false negative
f1_scores = {'KNN':{},'SVM':{},'LDA':{},'RF':{},'MLP':{}}
plotMe = pd.DataFrame(columns=['f1_score','method','class'])
for disState in list(set(CortexSubtypes))+['weighted avg']: # ignoring macro average
    for clf in performanceResults1.keys():
        f1_scores[clf][disState] = sum([performanceResults1[clf][i][disState]['f1-score'] for i in range(kFolds)])/kFolds
        if not disState=='weighted avg':
            plotMe = pd.concat([plotMe,pd.DataFrame({'f1_score':[performanceResults1[clf][i][disState]['f1-score'] for i in range(kFolds)], 'method':[clf for i in range(kFolds)], 'class':[disState for i in range(kFolds)]})])

# Print out the matrix of f1-scores for classifiers
print(pd.DataFrame(f1_scores))

# Dataframe with f1 scores copied and pasted into Excel, then the median value for each set of 100 f1 scores was calculated for each classifier (dataframe has f1 scores from all 100 rounds of cross validation from each classifier)

# Make a plot of the f1-scores to compare across all 100 cross-validated f1-scores
sns.boxplot(x='class', y='f1_score', hue='method', data=plotMe)
plt.legend(loc='center right', bbox_to_anchor =(1.18, 0.5), fontsize = 9)
plt.savefig("f1_pt1_Cortex_Spinal-600dpi.svg", format='svg', dpi = 600, bbox_inches = 'tight')
plt.savefig("f1_pt1_Cortex_Spinal-300dpi.svg", format='svg', dpi = 300,bbox_inches = 'tight')
plt.show()


## Calculate the positive predictive value (PPV) and negative predictive value (NPV)
#####    DID NOT USE THIS INFO     #####
# PPV = precision(OX) = tp / (tp + fp)
# NPV = prcision(other) = tn / (tn + fn)
# tp = true positive, fp = false positive, tn = true negative, fn = false negative
ppvNpv2 = {'KNN':{},'SVM':{},'LDA':{},'RF':{}, 'MLP':{}}
for disState in list(set(CortexSubtypes))+['weighted avg']: # ignoring macro average
    for clf in performanceResults1.keys():
        ppvNpv2[clf][disState] = sum([performanceResults1[clf][i][disState]['precision'] for i in range(kFolds)])/kFolds

# Print out the matrix of postive and negative predictive power for classifiers
print(pd.DataFrame(ppvNpv2))


## Build final classifiers by training on all combined data
# final classifiers are built using 100% of training data (versus splitting it into 80% train, 20% test)

###SVM
svm_clf = SVC(kernel='linear', C = 0.025, probability=True)
svm_clf.fit(gexpTrainc_scld, CortexSubtypes)

##LDA
lda_clf = LinearDiscriminantAnalysis()
lda_clf.fit(gexpTrainc_scld, CortexSubtypes)

#KNN
knn_clf = KNeighborsClassifier(8)
knn_clf.fit(gexpTrainc_scld, CortexSubtypes)

# #CNB
# cnb_clf = ComplementNB()
# cnb_clf.fit(gexpTrainc_scld_mm, CortexSubtypes)

#RF
rf_clf = RandomForestClassifier()
rf_clf.fit(gexpTrainc_scld, CortexSubtypes)

#MLP
mlp_clf = MLPClassifier(hidden_layer_sizes=(100,100,100), alpha=1e-5, learning_rate_init=(0.0001))
mlp_clf.fit(gexpTrainc_scld, CortexSubtypes)

##############
## Validate ##
##############

# Set up data for validation dataset
# In part 1, validating with spinal cord data
gexpValidations = gexp3c_spinal.T
gexpValidations_scld = StandardScaler().fit_transform(gexpValidations)
gexpValidations_scld_mm = MinMaxScaler().fit_transform(gexpValidations_scld)

# Set up disease states for validation dataset
## Generate labels that are in the set ['Ox', 'NotOx']
SpinalSubtypes = pd.Series(['OX' if i=='OX' else 'other' for i in spinalPhenos['PLSSubtype']])
print(SpinalSubtypes.value_counts()) #381 other, 138 OX

# Get results
# use classifiers on validation data
#SVM
svm_pred = svm_clf.predict(gexpValidations_scld)
svm_pred_prob = svm_clf.predict_proba(gexpValidations_scld)
print(pd.DataFrame(classification_report(SpinalSubtypes, svm_pred, output_dict=True, zero_division=0)))
print(confusion_matrix([int(i=='other') for i in SpinalSubtypes], [int(i=='other') for i in svm_pred]))

#LDA
lda_pred = lda_clf.predict(gexpValidations_scld)
lda_pred_prob = lda_clf.predict_proba(gexpValidations_scld)
print(pd.DataFrame(classification_report(SpinalSubtypes, lda_pred, output_dict=True, zero_division=0)))
print(confusion_matrix([int(i=='other') for i in SpinalSubtypes], [int(i=='other') for i in lda_pred]))

#KNN
knn_pred = knn_clf.predict(gexpValidations_scld)
knn_pred_prob = knn_clf.predict_proba(gexpValidations_scld)
print(pd.DataFrame(classification_report(SpinalSubtypes, knn_pred, output_dict=True, zero_division=0)))
print(confusion_matrix([int(i=='other') for i in SpinalSubtypes], [int(i=='other') for i in knn_pred]))

# RF
rf_pred = rf_clf.predict(gexpValidations_scld)
rf_pred_prob = rf_clf.predict_proba(gexpValidations_scld)
print(pd.DataFrame(classification_report(SpinalSubtypes, rf_pred, output_dict=True, zero_division=0)))
print(confusion_matrix([int(i=='other') for i in SpinalSubtypes], [int(i=='other') for i in rf_pred]))

# # CNB
# cnb_pred = cnb_clf.predict(gexpValidationc_scld_mm)
# cnb_pred_prob = cnb_clf.predict_proba(gexpValidationc_scld_mm)
# print(pd.DataFrame(classification_report(SpinalSubtypes, cnb_pred, output_dict=True, zero_division=0)))
# print(confusion_matrix([int(i=='other') for i in SpinalSubtypes], [int(i=='other') for i in cnb_pred]))

#MLP
mlp_pred = mlp_clf.predict(gexpValidations_scld)
mlp_pred_prob = mlp_clf.predict_proba(gexpValidations_scld)
print(pd.DataFrame(classification_report(SpinalSubtypes, mlp_pred, output_dict=True, zero_division=0)))
print(confusion_matrix([int(i=='other') for i in SpinalSubtypes], [int(i=='other') for i in mlp_pred]))


##################################################################
#########          Plot ROC curve for validation         #########
##################################################################

#without CNB since it was not a good classifier:
fig, ax = plt.subplots(figsize=(6, 6))
RocCurveDisplay.from_predictions([int(i=='other') for i in SpinalSubtypes], knn_pred_prob[:,1], name='KNN', ax=ax)
RocCurveDisplay.from_predictions([int(i=='other') for i in SpinalSubtypes], svm_pred_prob[:,1], name='SVM', ax=ax)
RocCurveDisplay.from_predictions([int(i=='other') for i in SpinalSubtypes], lda_pred_prob[:,1], name='LDA', ax=ax)
RocCurveDisplay.from_predictions([int(i=='other') for i in SpinalSubtypes], rf_pred_prob[:,1], name='RF', ax=ax)
RocCurveDisplay.from_predictions([int(i=='other') for i in SpinalSubtypes], mlp_pred_prob[:,1], name='MLP', ax=ax)
plt.savefig("ROC_pt1_Cortex_Spinal-600dpi.svg", format='svg', dpi = 600)
plt.savefig("ROC_pt1_Cortex_Spinal-300dpi.svg", format='svg', dpi = 300)
plt.show()

###  If you want to plot curves on seperate plots (for best classifiers)
###  SVM example -->
fig, ax = plt.subplots(figsize=(6, 6))
RocCurveDisplay.from_predictions([int(i=='other') for i in SpinalSubtypes], svm_pred_prob[:,1], name='SVM', ax=ax)
plt.show()

###############################################################################
###############################################################################


####### Seperate scripts used for figures/ info for parts 2 and 3, but have generally the same format            #######  (parts below are just not the most updated versions)



#################################
#######      Part 2       #######
#################################


#Train on best genes from spinal, validate using expression in the cortex 
#review phenotype data, make sure its the same as before :)
#number of ox vs notox subtypes
print(cortexPhenos['PLSSubtype'].value_counts()) #notox = 346, ox = 239
print(spinalPhenos['PLSSubtype'].value_counts()) #notox = 381, ox = 138

###############################################################################
##### Gene expression data remains same:

gexp3c_cortex = pd.read_csv("C:/Users/trpen/OneDrive/Desktop/ALS Classification/ML/HiSeqHoldouot_RPKM_Expression_CortexSamples_3Gene.csv",header=0, index_col=0) 
gexp3c_spinal = pd.read_csv("C:/Users/trpen/OneDrive/Desktop/ALS Classification/ML/HiSeqHoldouot_RPKM_Expression_SpinalSamples_3Gene.csv",header=0, index_col=0) 

##### Spinal now becomes training data, cortex used for validation
gexpTrain_s = gexp3c_spinal
gexpValidationc = gexp3c_cortex 

###############################################################################
#####                   TRAINING THE CLASSIFIERS RND2                     #####
###############################################################################

#set data to train classifiers on, (two methods needed, because CNB method needs values [0-1])
gexpTrain_s = gexp3c_spinal.T
gexpTrains_scld = StandardScaler().fit_transform(gexpTrain_s)
gexpTrains_scld_mm = MinMaxScaler().fit_transform(gexpTrains_scld)

#Repeat: set spinal subtypes for the training:
## Generate labels that are in the set ['Ox', 'NotOx']
# This will bulid a classifier that can discriminate Ox from the other subtypes
SpinalSubtypes = pd.Series(['OX' if i=='OX' else 'other' for i in spinalPhenos['PLSSubtype']])
print(SpinalSubtypes.value_counts())

## Set up classification comparison for 100 fold cross-validation
# Store out second set of performance results in dictionary
kFolds = 100
performanceResults2 = {'KNN':[],'SVM':[],'LDA':[],'RF':[],'MLP':[]}
for i in range(kFolds):
    ## Split up data into train and test sets using random integer seed
    gexpData_train, gexpData_test, disState_train, disState_test = train_test_split(gexpTrains_scld, SpinalSubtypes, test_size=0.2, random_state=i) #20% of data is going into test set (test size)
    gexpData_train_mm, gexpData_test_mm, disState_train_mm, disState_test_mm = train_test_split(gexpTrains_scld_mm, SpinalSubtypes, test_size=0.2, random_state=i)
    ## Train and test the classifiers
    # k-nearest neighbors (KNN)
    knn_clf = KNeighborsClassifier(8)
    knn_clf.fit(gexpData_train, disState_train)
    knn_pred = knn_clf.predict(gexpData_test)
    performanceResults2['KNN'].append(classification_report(disState_test, knn_pred, output_dict=True, zero_division=0))
    # Support vector machines (SVM)
    svm_clf = SVC(kernel='linear', C = 0.025)
    svm_clf.fit(gexpData_train, disState_train)
    svm_pred = svm_clf.predict(gexpData_test)
    performanceResults2['SVM'].append(classification_report(disState_test, svm_pred, output_dict=True, zero_division=0))
    # Linear discriminant analysis (LDA)
    lda_clf = LinearDiscriminantAnalysis()
    lda_clf.fit(gexpData_train, disState_train)
    lda_pred = lda_clf.predict(gexpData_test)
    performanceResults2['LDA'].append(classification_report(disState_test, lda_pred, output_dict=True, zero_division=0))
    # Random Forest ensemble (RF)
    rf_clf = RandomForestClassifier()
    rf_clf.fit(gexpData_train, disState_train)
    rf_pred = rf_clf.predict(gexpData_test)
    performanceResults2['RF'].append(classification_report(disState_test, rf_pred, output_dict=True, zero_division=0))
    # # Complement Naive Bayes (CNB)
    # cnb_clf = ComplementNB()
    # cnb_clf.fit(gexpData_train_mm, disState_train)
    # cnb_pred = cnb_clf.predict(gexpData_test_mm)
    # performanceResults2['CNB'].append(classification_report(disState_test_mm, cnb_pred, output_dict=True, zero_division=0))
    #Multilayer perceptron (MLP)
    mlp_clf = MLPClassifier(hidden_layer_sizes=(100,100,100), alpha=1e-5, learning_rate_init=(0.0001))
    mlp_clf.fit(gexpData_train, disState_train)
    mlp_pred = mlp_clf.predict(gexpData_test)
    performanceResults2['MLP'].append(classification_report(disState_test, mlp_pred, output_dict=True, zero_division=0))


## Build figure to describe second set of classifiers: using 'f1-score'

f1_scores2 = {'KNN':{},'SVM':{},'LDA':{},'RF':{},'MLP':{}}
plotMe2 = pd.DataFrame(columns=['f1_score','method','class'])
for disState in list(set(SpinalSubtypes))+['weighted avg']: # ignoring macro average
    for clf in performanceResults2.keys():
        f1_scores2[clf][disState] = sum([performanceResults2[clf][i][disState]['f1-score'] for i in range(kFolds)])/kFolds
        if not disState=='weighted avg':
            plotMe2 = pd.concat([plotMe2,pd.DataFrame({'f1_score':[performanceResults2[clf][i][disState]['f1-score'] for i in range(kFolds)], 'method':[clf for i in range(kFolds)], 'class':[disState for i in range(kFolds)]})])

# Print out the matrix of f1-scores for classifiers
print(pd.DataFrame(f1_scores2))

# Make a plot of the f1-scores to compare across all 100 cross-validated f1-scores
sns.boxplot(x='class', y='f1_score', hue='method', data=plotMe2)
plt.savefig("f1_pt2_Spinal_Cortex-600dpi.svg", format='svg', dpi = 600)
plt.savefig("f1_pt2_Spinal_Cortex-300dpi.svg", format='svg', dpi = 300)
plt.show()


## Calculate the positive predictive value (PPV) and negative predictive value (NPV)
ppvNpv3 = {'KNN':{},'SVM':{},'LDA':{},'RF':{}, 'MLP':{}}
for disState in list(set(SpinalSubtypes))+['weighted avg']: # ignoring macro average
    for clf in performanceResults2.keys():
        ppvNpv3[clf][disState] = sum([performanceResults2[clf][i][disState]['precision'] for i in range(kFolds)])/kFolds

# Print out the matrix of postive and negative predictive power for classifiers
print(pd.DataFrame(ppvNpv3))


## Build final classifiers by training on all spinal cord data
###SVM
svm_clf2 = SVC(kernel='linear', C = 0.025, probability=True)
svm_clf2.fit(gexpTrains_scld, SpinalSubtypes)

##LDA
lda_clf2 = LinearDiscriminantAnalysis()
lda_clf2.fit(gexpTrains_scld, SpinalSubtypes)

#KNN
knn_clf2 = KNeighborsClassifier(8)
knn_clf2.fit(gexpTrains_scld, SpinalSubtypes)

# #CNB
# cnb_clf2 = ComplementNB()
# cnb_clf2.fit(gexpTrain2s_scld_mm, SpinalSubtypes)

#RF
rf_clf2 = RandomForestClassifier()
rf_clf2.fit(gexpTrains_scld, SpinalSubtypes)

#MLP
mlp_clf2 = MLPClassifier(hidden_layer_sizes=(100,100,100), alpha=1e-5, learning_rate_init=(0.0001))
mlp_clf2.fit(gexpTrains_scld, SpinalSubtypes)

##############
## Validate ##
##############

# Set up data for validation dataset
gexpValidationc = gexp3c_cortex.T
gexpValidationc_scld = StandardScaler().fit_transform(gexpValidationc)
gexpValidationc_scld_mm = MinMaxScaler().fit_transform(gexpValidationc_scld)

# Set up disease states for validation dataset
CortexSubtypes = pd.Series(['OX' if i=='OX' else 'other' for i in cortexPhenos['PLSSubtype']])
print(CortexSubtypes.value_counts()) #346 other, 239 OX

# Get results
#SVM
svm_pred2 = svm_clf2.predict(gexpValidationc_scld)
svm_pred_prob2 = svm_clf2.predict_proba(gexpValidationc_scld)
print(pd.DataFrame(classification_report(CortexSubtypes, svm_pred2, output_dict=True, zero_division=0)))
print(confusion_matrix([int(i=='other') for i in CortexSubtypes], [int(i=='other') for i in svm_pred2]))

#LDA
lda_pred2 = lda_clf2.predict(gexpValidationc_scld)
lda_pred_prob2 = lda_clf2.predict_proba(gexpValidationc_scld)
print(pd.DataFrame(classification_report(CortexSubtypes, lda_pred2, output_dict=True, zero_division=0)))
print(confusion_matrix([int(i=='other') for i in CortexSubtypes], [int(i=='other') for i in lda_pred2]))

#KNN
knn_pred2 = knn_clf2.predict(gexpValidationc_scld)
knn_pred_prob2 = knn_clf2.predict_proba(gexpValidationc_scld)
print(pd.DataFrame(classification_report(CortexSubtypes, knn_pred2, output_dict=True, zero_division=0)))
print(confusion_matrix([int(i=='other') for i in CortexSubtypes], [int(i=='other') for i in knn_pred2]))

# RF
rf_pred2 = rf_clf2.predict(gexpValidationc_scld)
rf_pred_prob2 = rf_clf2.predict_proba(gexpValidationc_scld)
print(pd.DataFrame(classification_report(CortexSubtypes, rf_pred2, output_dict=True, zero_division=0)))
print(confusion_matrix([int(i=='other') for i in CortexSubtypes], [int(i=='other') for i in rf_pred2]))

# # CNB
# cnb_pred2 = cnb_clf2.predict(gexpValidations_scld_mm)
# cnb_pred_prob2 = cnb_clf2.predict_proba(gexpValidations_scld_mm)
# print(pd.DataFrame(classification_report(CortexSubtypes, cnb_pred2, output_dict=True, zero_division=0)))
# print(confusion_matrix([int(i=='other') for i in CortexSubtypes], [int(i=='other') for i in cnb_pred2]))

#MLP
mlp_pred2 = mlp_clf2.predict(gexpValidationc_scld)
mlp_pred_prob2 = mlp_clf2.predict_proba(gexpValidationc_scld)
print(pd.DataFrame(classification_report(CortexSubtypes, mlp_pred2, output_dict=True, zero_division=0)))
print(confusion_matrix([int(i=='other') for i in CortexSubtypes], [int(i=='other') for i in mlp_pred2]))


##################################################################
#########          Plot ROC curve for validation         #########
##################################################################
    
fig, ax = plt.subplots(figsize=(6, 6))
RocCurveDisplay.from_predictions([int(i=='other') for i in CortexSubtypes], knn_pred_prob2[:,1], name='KNN', ax=ax)
RocCurveDisplay.from_predictions([int(i=='other') for i in CortexSubtypes], svm_pred_prob2[:,1], name='SVM', ax=ax)
RocCurveDisplay.from_predictions([int(i=='other') for i in CortexSubtypes], lda_pred_prob2[:,1], name='LDA', ax=ax)
RocCurveDisplay.from_predictions([int(i=='other') for i in CortexSubtypes], rf_pred_prob2[:,1], name='RF', ax=ax)
RocCurveDisplay.from_predictions([int(i=='other') for i in CortexSubtypes], mlp_pred_prob2[:,1], name='MLP', ax=ax)
plt.savefig("ROC_pt2_Spinal_Cortex-600dpi.svg", format='svg', dpi = 600)
plt.savefig("ROC_pt2_Spinal_Cortex-300dpi.svg", format='svg', dpi = 300)
plt.show()

###############################################################################
###############################################################################


#################################
#######      Part 3       #######
#################################

## Train on novaseq data, validate on hiseq data

#Identify samples sequenced with NovaSeq platform
## Cortex:
nova_cortex = cortex_info.loc[cortex_info.Platform == 'NovaSeq'] #355 samples
## Spinal:
nova_spinal = spinal_info.loc[spinal_info.Platform == 'NovaSeq'] #334 samples

# Combine novaseq samples + hiseq samples
Nova = nova_cortex.append(nova_spinal) #689 samples
NovaPhenos = Nova.set_index('Subject')

#Identify samples sequenced with HiSeq platform
hi_cortex = cortex_info.loc[cortex_info.Platform == 'HiSeq'] #230 samples
## Spinal:
hi_spinal = spinal_info.loc[spinal_info.Platform == 'HiSeq'] #185 samples

# Combine novaseq samples + hiseq samples
Hi = hi_cortex.append(hi_spinal) #415 samples
HiPhenos = Hi.set_index('Subject')

# review phenotype data: number of ox vs notox subtypes
print(NovaPhenos['PLSSubtype'].value_counts()) #notox = 456, ox = 233
print(HiPhenos['PLSSubtype'].value_counts()) #notox = 271, ox = 144

###############################################################################
##### Gene expression data remains same:

gexp3_cortex = pd.read_csv("C:/Users/trpen/OneDrive/Desktop/ALS Classification/ML/HiSeqHoldouot_RPKM_Expression_CortexSamples_3Gene.csv",header=0, index_col=0) 
gexp3_spinal = pd.read_csv("C:/Users/trpen/OneDrive/Desktop/ALS Classification/ML/HiSeqHoldouot_RPKM_Expression_SpinalSamples_3Gene.csv",header=0, index_col=0) 


##combine cortex and spinal cord gene expression data
###gexp3_merge = gexp3_cortex.append(gexp3_spinal) #1104 samples, need to remove extra 3 rows
gexp3 = pd.concat([gexp3_cortex, gexp3_spinal], axis =1 ) #should be 1104 total samples (nova + hi)

###subset to NovaSeq samples
gexpNova = gexp3[NovaPhenos.loc[NovaPhenos['Platform']=='NovaSeq'].index] #689

###subset gexp 3 to HiSeq Samples
gexpHi = gexp3[HiPhenos.loc[HiPhenos['Platform']=='HiSeq'].index] #415

#label NovaSeq as training data, label HiSeq as validation data
gexpTrain = gexpNova
gexpValidation = gexpHi

###############################################################################
#####                   TRAINING THE CLASSIFIERS RND3                     #####
###############################################################################

#set data to train classifiers on
gexpTrain3 = gexpTrain.T
gexpTrain3_scld = StandardScaler().fit_transform(gexpTrain3)
gexpTrain3_scld_mm = MinMaxScaler().fit_transform(gexpTrain3_scld)

## Generate labels that are in the set ['Ox', 'NotOx']
# This will bulid a classifier that can discriminate Ox from the other subtypes
NovaSubtypes = pd.Series(['OX' if i=='OX' else 'other' for i in NovaPhenos['PLSSubtype']])
print(NovaSubtypes.value_counts()) #456 other, 233 OX (should be same as before)

## Set up classification comparison for 100 fold cross-validation, store results in dict.
kFolds = 100
performanceResults3 = {'KNN':[],'SVM':[],'LDA':[],'RF':[],'MLP':[]}
for i in range(kFolds):
    ## Split up data into train and test sets using random integer seed
    gexpData_train, gexpData_test, disState_train, disState_test = train_test_split(gexpTrain3_scld, NovaSubtypes, test_size=0.2, random_state=i) #20% of data is going into test set (test size)
    gexpData_train_mm, gexpData_test_mm, disState_train_mm, disState_test_mm = train_test_split(gexpTrain3_scld_mm, NovaSubtypes, test_size=0.2, random_state=i)
    ## Train and test the classifiers
    # k-nearest neighbors (KNN)
    knn_clf = KNeighborsClassifier(8)
    knn_clf.fit(gexpData_train, disState_train)
    knn_pred = knn_clf.predict(gexpData_test)
    performanceResults3['KNN'].append(classification_report(disState_test, knn_pred, output_dict=True, zero_division=0))
    # Support vector machines (SVM)
    svm_clf = SVC(kernel='linear', C = 0.025)
    svm_clf.fit(gexpData_train, disState_train)
    svm_pred = svm_clf.predict(gexpData_test)
    performanceResults3['SVM'].append(classification_report(disState_test, svm_pred, output_dict=True, zero_division=0))
    # Linear discriminant analysis (LDA)
    lda_clf = LinearDiscriminantAnalysis()
    lda_clf.fit(gexpData_train, disState_train)
    lda_pred = lda_clf.predict(gexpData_test)
    performanceResults3['LDA'].append(classification_report(disState_test, lda_pred, output_dict=True, zero_division=0))
    # Random Forest ensemble (RF)
    rf_clf = RandomForestClassifier()
    rf_clf.fit(gexpData_train, disState_train)
    rf_pred = rf_clf.predict(gexpData_test)
    performanceResults3['RF'].append(classification_report(disState_test, rf_pred, output_dict=True, zero_division=0))
    # # Complement Naive Bayes (CNB)
    # cnb_clf = ComplementNB()
    # cnb_clf.fit(gexpData_train_mm, disState_train)
    # cnb_pred = cnb_clf.predict(gexpData_test_mm)
    # performanceResults3['CNB'].append(classification_report(disState_test_mm, cnb_pred, output_dict=True, zero_division=0))
    #Multilayer perceptron (MLP)
    mlp_clf = MLPClassifier(hidden_layer_sizes=(100,100,100), alpha=1e-5, learning_rate_init=(0.0001))
    mlp_clf.fit(gexpData_train, disState_train)
    mlp_pred = mlp_clf.predict(gexpData_test)
    performanceResults3['MLP'].append(classification_report(disState_test, mlp_pred, output_dict=True, zero_division=0))


## Build figure to describe third set of classifiers: using 'f1-score'
f1_scores3 = {'KNN':{},'SVM':{},'LDA':{},'RF':{},'MLP':{}}
plotMe3 = pd.DataFrame(columns=['f1_score','method','class'])
for disState in list(set(NovaSubtypes))+['weighted avg']: # ignoring macro average
    for clf in performanceResults3.keys():
        f1_scores3[clf][disState] = sum([performanceResults3[clf][i][disState]['f1-score'] for i in range(kFolds)])/kFolds
        if not disState=='weighted avg':
            plotMe3 = pd.concat([plotMe3,pd.DataFrame({'f1_score':[performanceResults3[clf][i][disState]['f1-score'] for i in range(kFolds)], 'method':[clf for i in range(kFolds)], 'class':[disState for i in range(kFolds)]})])

# Print out the matrix of f1-scores for classifiers
print(pd.DataFrame(f1_scores3))

# Make a plot of the f1-scores to compare across all 100 cross-validated f1-scores
sns.boxplot(x='class', y='f1_score', hue='method', data=plotMe3)
plt.savefig("f1_pt3_Nova_Hi-600dpi.svg", format='svg', dpi = 600)
plt.savefig("f1_pt3_Nova_Hi-300dpi.svg", format='svg', dpi = 300)
plt.show()


## Calculate the positive predictive value (PPV) and negative predictive value (NPV)
ppvNpv4 = {'KNN':{},'SVM':{},'LDA':{},'RF':{}, 'MLP':{}}
for disState in list(set(NovaSubtypes))+['weighted avg']: # ignoring macro average
    for clf in performanceResults3.keys():
        ppvNpv4[clf][disState] = sum([performanceResults3[clf][i][disState]['precision'] for i in range(kFolds)])/kFolds

# Print out the matrix of postive and negative predictive power for classifiers
print(pd.DataFrame(ppvNpv4))


###SVM
## Build final classifiers by training on all combined spinal cord data
svm_clf3 = SVC(kernel='linear', C = 0.025, probability=True)
svm_clf3.fit(gexpTrain3_scld, NovaSubtypes)

##LDA
lda_clf3 = LinearDiscriminantAnalysis()
lda_clf3.fit(gexpTrain3_scld, NovaSubtypes)

#KNN
knn_clf3 = KNeighborsClassifier(8)
knn_clf3.fit(gexpTrain3_scld, NovaSubtypes)

# #CNB
# cnb_clf3 = ComplementNB()
# cnb_clf3.fit(gexpTrain3_scld_mm, NovaSubtypes)

#RF
rf_clf3 = RandomForestClassifier()
rf_clf3.fit(gexpTrain3_scld, NovaSubtypes)

#MLP
mlp_clf3 = MLPClassifier(hidden_layer_sizes=(100,100,100), alpha=1e-5, learning_rate_init=(0.0001))
mlp_clf3.fit(gexpTrain3_scld, NovaSubtypes)

##############
## Validate ##
##############

# Set up data for validation dataset
gexpValidation3 = gexpValidation.T
gexpValidation3_scld = StandardScaler().fit_transform(gexpValidation3)
gexpValidation3_scld_mm = MinMaxScaler().fit_transform(gexpValidation3_scld)

# Set up disease states for validation dataset
## Generate labels that are in the set ['Ox', 'NotOx']
# This will bulid a classifier that can discriminate Ox from the other subtypes
HiSubtypes = pd.Series(['OX' if i=='OX' else 'other' for i in HiPhenos['PLSSubtype']])
print(HiSubtypes.value_counts()) #271 other, 144 OX

# Get results
#SVM
svm_pred3 = svm_clf3.predict(gexpValidation3_scld)
svm_pred_prob3 = svm_clf3.predict_proba(gexpValidation3_scld)
print(pd.DataFrame(classification_report(HiSubtypes, svm_pred3, output_dict=True, zero_division=0)))
print(confusion_matrix([int(i=='other') for i in HiSubtypes], [int(i=='other') for i in svm_pred3]))

#LDA
lda_pred3 = lda_clf3.predict(gexpValidation3_scld)
lda_pred_prob3 = lda_clf3.predict_proba(gexpValidation3_scld)
print(pd.DataFrame(classification_report(HiSubtypes, lda_pred3, output_dict=True, zero_division=0)))
print(confusion_matrix([int(i=='other') for i in HiSubtypes], [int(i=='other') for i in lda_pred3]))

#KNN
knn_pred3 = knn_clf3.predict(gexpValidation3_scld)
knn_pred_prob3 = knn_clf3.predict_proba(gexpValidation3_scld)
print(pd.DataFrame(classification_report(HiSubtypes, knn_pred3, output_dict=True, zero_division=0)))
print(confusion_matrix([int(i=='other') for i in HiSubtypes], [int(i=='other') for i in knn_pred3]))

# RF
rf_pred3 = rf_clf3.predict(gexpValidation3_scld)
rf_pred_prob3 = rf_clf3.predict_proba(gexpValidation3_scld)
print(pd.DataFrame(classification_report(HiSubtypes, rf_pred3, output_dict=True, zero_division=0)))
print(confusion_matrix([int(i=='other') for i in HiSubtypes], [int(i=='other') for i in rf_pred3]))

# # CNB
# cnb_pred3= cnb_clf3.predict(gexpValidation3_scld_mm)
# cnb_pred_prob3 = cnb_clf3.predict_proba(gexpValidation3_scld_mm)
# print(pd.DataFrame(classification_report(HiSubtypes, cnb_pred3, output_dict=True, zero_division=0)))
# print(confusion_matrix([int(i=='other') for i in HiSubtypes], [int(i=='other') for i in cnb_pred3]))

#MLP
mlp_pred3 = mlp_clf3.predict(gexpValidation3_scld)
mlp_pred_prob3 = mlp_clf3.predict_proba(gexpValidation3_scld)
print(pd.DataFrame(classification_report(HiSubtypes, mlp_pred3, output_dict=True, zero_division=0)))
print(confusion_matrix([int(i=='other') for i in HiSubtypes], [int(i=='other') for i in mlp_pred3]))


##################################################################
#########          Plot ROC curve for validation         #########
##################################################################
    
fig, ax = plt.subplots(figsize=(6, 6))
RocCurveDisplay.from_predictions([int(i=='other') for i in HiSubtypes], knn_pred_prob3[:,1], name='KNN', ax=ax)
RocCurveDisplay.from_predictions([int(i=='other') for i in HiSubtypes], svm_pred_prob3[:,1], name='SVM', ax=ax)
RocCurveDisplay.from_predictions([int(i=='other') for i in HiSubtypes], lda_pred_prob3[:,1], name='LDA', ax=ax)
RocCurveDisplay.from_predictions([int(i=='other') for i in HiSubtypes], rf_pred_prob3[:,1], name='RF', ax=ax)
RocCurveDisplay.from_predictions([int(i=='other') for i in HiSubtypes], mlp_pred_prob3[:,1], name='MLP', ax=ax)
plt.savefig("ROC_pt3_Nova_Hi-600dpi.svg", format='svg', dpi = 600)
plt.savefig("ROC_pt3_Nova_Hi-300dpi.svg", format='svg', dpi = 300)
plt.show()




