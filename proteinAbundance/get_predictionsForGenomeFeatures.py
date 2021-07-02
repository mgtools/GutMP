"""
script that will use the feature space tables of two phenotypes in comparison and will train a model
"""

import pandas as pd
import numpy as np
from sklearn import svm
from collections import Counter
import os
import sys

in_dir = sys.argv[1]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/phenotype_featureSpaces/'

phenotype_1 = sys.argv[2]#'healthy'
phenotype_2 = sys.argv[3]#'CD'

min_logOdds = float(sys.argv[4])#2.0
top_n = int(sys.argv[5])#100
min_samples = int(sys.argv[6])#3
min_features = int(sys.argv[7])#1


data_points_f = f'{in_dir + phenotype_1}_vs_{phenotype_2}/{phenotype_1}_vs_{phenotype_2}_featureSpaces_minLogOdds_{min_logOdds}_top_n_{top_n}_minSamples_{min_samples}.tsv'
data_points_df = pd.read_csv(data_points_f, sep = '\t')
data_points_df = data_points_df.sample(frac = 1).reset_index(drop = True)

features = list(data_points_df.columns)[1:-1]
data_points_filtered_df = pd.DataFrame(columns = list(data_points_df.columns))
data_points_filtered_out_df = pd.DataFrame(columns = list(data_points_df.columns))

out_dir = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/prediction_plots/'

if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
out_dir = out_dir + f'{phenotype_1}_vs_{phenotype_2}/'
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

row_count = 0
filtered_out_row_count = 0
for i, row in data_points_df.iterrows():
    if len(set(list(row)[1:-1])) >= min_features:
        data_points_filtered_df.loc[row_count] = row
        row_count += 1
    else:
        data_points_filtered_out_df.loc[filtered_out_row_count] = row
        filtered_out_row_count += 1
print(f'the shape of the input matrix {data_points_df.shape}')
print(f'the shape of the input matrix after keeping data points with a minimum of {min_features} non zero values for the features, {data_points_filtered_df.shape}')

data_points_filtered_values = data_points_filtered_df.values

features = data_points_filtered_values[: , 1:-1]
labels = data_points_filtered_values[: , -1]

label_counts = Counter(labels)

from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(features, labels, test_size=0.1,random_state=109) # 70% training and 30% test

#this step is needed if your y_train and y_test are saved as type object, you need to do an explicit type casting.
y_train = y_train.astype('int')
y_test = y_test.astype('int')


X = features.astype('float')
y = labels.astype('int')

#Import svm model
from sklearn import svm
import matplotlib.pyplot as plt
from sklearn.metrics import auc
from sklearn.metrics import plot_roc_curve
from sklearn.model_selection import StratifiedKFold

cv = StratifiedKFold(n_splits = 10)
classifier = svm.SVC(kernel='linear', probability=True)

tprs = []
aucs = []
mean_fpr = np.linspace(0, 1, 100)

fig, ax = plt.subplots(figsize = (8,8))
for i, (train, test) in enumerate(cv.split(X, y)):
    classifier.fit(X[train], y[train])
    viz = plot_roc_curve(classifier, X[test], y[test],
                         name='ROC fold {}'.format(i),
                         alpha=0.3, lw=1, ax=ax)
    interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
    interp_tpr[0] = 0.0
    tprs.append(interp_tpr)
    aucs.append(viz.roc_auc)

ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
        label='Chance', alpha=.8)

mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
SVM_auc = mean_auc
std_auc = np.std(aucs)
ax.plot(mean_fpr, mean_tpr, color='b',
        label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
        lw=2, alpha=.8)

std_tpr = np.std(tprs, axis=0)
tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                label=r'$\pm$ 1 std. dev.')

ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],       
       title=f"Total Points: {len(labels)} ({phenotype_1}:{label_counts[1]}/{phenotype_2}:{label_counts[0]}) top-N DE genomes : {top_n*2} \nmin Samples {min_samples} min-features : {min_features } N-Features: {len(features[0])}")
ax.legend(loc="lower right")
plt.savefig(out_dir + f'SVM_N_fold_{phenotype_1}_vs_{phenotype_2}_logOdds_{min_logOdds}_topN_DE_genomes_{top_n*2}_min_samples_{min_samples}_minFeatures_{min_features}.png')

print('average AUC for SVM over 10 folds\n', round(mean_auc, 2), '(', data_points_filtered_values.shape[0], ')')

##############################################################################

from sklearn.ensemble import RandomForestClassifier

cv = StratifiedKFold(n_splits = 10)
classifier= RandomForestClassifier(max_depth=10, random_state=0, n_estimators = 10)

tprs = []
aucs = []
mean_fpr = np.linspace(0, 1, 100)

fig, ax = plt.subplots(figsize = (8,8))
for i, (train, test) in enumerate(cv.split(X, y)):
    classifier.fit(X[train], y[train])
    viz = plot_roc_curve(classifier, X[test], y[test],
                         name='ROC fold {}'.format(i),
                         alpha=0.3, lw=1, ax=ax)
    interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
    interp_tpr[0] = 0.0
    tprs.append(interp_tpr)
    aucs.append(viz.roc_auc)

ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
        label='Chance', alpha=.8)

mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)
ax.plot(mean_fpr, mean_tpr, color='b',
        label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
        lw=2, alpha=.8)

std_tpr = np.std(tprs, axis=0)
tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                label=r'$\pm$ 1 std. dev.')

ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
       title=f"Total Points: {len(labels)} ({phenotype_1}:{label_counts[1]}/{phenotype_2}:{label_counts[0]}) top-N DE genomes : {top_n*2} \nmin Samples {min_samples} min-features : {min_features } N-Features: {len(features[0])}")

ax.legend(loc="lower right")
plt.savefig(out_dir + f'RF_N_fold_{phenotype_1}_vs_{phenotype_2}_logOdds_{min_logOdds}_topN_DE_genomes_{top_n*2}_min_samples_{min_samples}_minFeatures_{min_features}.png')

print('average AUC for RF over 10 folds\n', round(mean_auc, 2), '(', data_points_filtered_values.shape[0], ')')
