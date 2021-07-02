"""
script that will run the different models over the different datasets (phenotype pairs), using the different feature selection methods,
each dataset pair will be repeated 100 times, with 10-fold cross validation and will report the feature sizes, the chosen features
and the total number of AUC's over each fold for each of the cases and will output this summary in a json dictionary file.
"""


from sklearn import svm
import numpy as np
from sklearn.model_selection import ShuffleSplit
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.svm import LinearSVC
from sklearn.feature_selection import SelectFromModel
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.feature_selection import RFE
import matplotlib.pyplot as plt
from sklearn.metrics import auc
from sklearn.metrics import plot_roc_curve
from sklearn.model_selection import StratifiedKFold
import json
import os

def get_predictions(X, y):   
    
    cv = ShuffleSplit(n_splits=10, test_size=0.2)    
    
    SVM_classifier = svm.SVC(kernel='rbf', probability = True)    
    SVM_scores = cross_val_score(SVM_classifier, X, y, cv = cv)
    
    RF_classifier= RandomForestClassifier(max_depth=100, random_state=0, n_estimators = 100)
    RF_scores = cross_val_score(RF_classifier, X, y, cv = cv   )
    
    return np.mean(SVM_scores), np.mean(RF_scores)

def get_predictions_2(X, y):
    
    cv = StratifiedKFold(n_splits = 10, shuffle = True)
    classifier = svm.SVC(kernel='rbf', probability=True)

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
    SVM_auc = round(mean_auc, 2)
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
           title=f"Total Points: {len(labels)} ({phen1}:{label_counts[1]}/{phen2}:{label_counts[0]})  N-Features: {X.shape[1]}")
    ax.legend(loc="lower right")
    #plt.show()
    plt.close()


    from sklearn.ensemble import RandomForestClassifier

    cv = StratifiedKFold(n_splits = 10, shuffle = True)
    classifier= RandomForestClassifier(max_depth=100, random_state=0, n_estimators = 100)

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
    RF_auc = round(mean_auc, 2)
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
           title=f"Total Points: {len(labels)} ({phen1}:{label_counts[1]}/{phen2}:{label_counts[0]})  N-Features: {X.shape[1]}")
    ax.legend(loc="lower right")
    #plt.show()
    plt.close()
    
    return SVM_auc, RF_auc




import pandas as pd
from collections import Counter

feature_space_dir = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/DESeq2_logOdds_combinedFeatureSpaces_taxa/'
selected_feature_space_dir = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/phenotype_featureSpaces/'
comparisons_f = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/comparisonLists.txt'
n_reps = 100

out_dir = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/modelComparisons_taxa/'

tax_level = 's'

min_log_odds = 2.0
min_p_val = 0.05
top_n = 50
min_samples = 5

if not os.path.isdir(out_dir):
    os.mkdir(out_dir)


n_features_dic = {'CD_subUC_vs_UC' : 80,
                  'healthy_vs_acute_leukemia': 90,
                  'healthy_vs_UC': 90,
                  'healthy_subUC_vs_UC' : 60,
                  'healthy_vs_acute_leukemia_subHealthy': 90,
                  'healthy_subT1D_vs_T1D' : 30,                  
                  'healthy_vs_CD': 100
                 }


prediction_summaries_dic = {}

with open(comparisons_f, 'r') as in_f:
    for I, line in enumerate(in_f):
        line = line.strip().split('\t')
        phen1, phen2 = line[0], line[1]

        print(f'processing {phen1} vs {phen2} ...')


        pair = f'{phen1}_vs_{phen2}'

        prediction_summaries_dic[pair] = {'all': {},
                                         'L1': {},
                                         'DT': {},
                                         'RFE': {},
                                         'SF': {}}


        #data_points_f = f'{feature_space_dir + phen1}_vs_{phen2}_combinedFeatureSpaces_minLogOdds_2.0_minPval_0.05.tsv'
        #data_points_f = f'{feature_space_dir + phen1}_vs_{phen2}_allGenomesFeatureSpaces.tsv'
        data_points_f = f'{feature_space_dir + phen1}_vs_{phen2}_{tax_level}_combinedFeatureSpaces_minLogOdds_{min_log_odds}_minPval_{min_p_val}.tsv'
        #data_points_f = f'{feature_space_dir + phen1}_vs_{phen2}_{tax_level}_allGenomesFeatureSpaces.tsv'
        data_points_df = pd.read_csv(data_points_f, sep = '\t')
        data_points_df = data_points_df.sample(frac = 1).reset_index(drop = True)

        feature_names = list(data_points_df.columns)[1:-1]

        data_points_values = data_points_df.values

        features = data_points_values[: , 1:-1]
        labels = data_points_values[: , -1]
        label_counts = Counter(labels)

        label_counts = Counter(labels)

        ## Using all features combined
        X = features.astype('float')
        y = labels.astype('int')

        print(f'begining feature space {X.shape}')
        
        
        selected_data_points_f = f'{selected_feature_space_dir + phen1}_vs_{phen2}/{phen1}_vs_{phen2}_featureSpaces_minLogOdds_{min_log_odds}_top_n_{top_n}_minSamples_{min_samples}.tsv'
        selected_data_points_df = pd.read_csv(selected_data_points_f, sep = '\t')
        selected_feature_names = list(selected_data_points_df.columns)[1:-1]
        selected_data_points_values = selected_data_points_df.values
        
        selected_feature_values = selected_data_points_values[: , 1: -1]
        selected_labels = selected_data_points_values[: , -1]
        
        X_selected = selected_feature_values.astype('float')
        y_selected = selected_labels.astype('int')
        
        ## L1 feature selection method

        SVM_auc, RF_auc = get_predictions(X, y)

        lsvc = LinearSVC(C=0.055, penalty="l1", dual=False).fit(X, y)
        L1_model = SelectFromModel(lsvc, prefit=True)
        X_L1 = L1_model.transform(X)
        L1_features = list(np.array(feature_names)[L1_model.get_support()])

        print(f'L1 feature space {X_L1.shape}')

        ## Decision tree feature selection method

        clf = ExtraTreesClassifier(n_estimators=1000)
        clf = clf.fit(X, y)
        DT_model = SelectFromModel(clf, prefit=True)
        X_DT = DT_model.transform(X)
        DT_features = list(np.array(feature_names)[DT_model.get_support()])

        print(f'DT feature space {X_DT.shape}')

        
        rfe_n_features = n_features_dic[pair]
        
        ## recursive Feature elimination

        rfe = RFE(estimator=RandomForestClassifier(), n_features_to_select = rfe_n_features, step = 1)
        rfe = rfe.fit(X, y)
        selected_features = list(rfe.support_)
        X_RFE = X[:, selected_features]
        RFE_features = list(np.array(feature_names)[selected_features])

        print(f'RFE feature space {X_RFE.shape}')
        
        
        prediction_summaries_dic[pair]['all']['features'] = feature_names
        prediction_summaries_dic[pair]['all']['n_features'] = len(feature_names)
        prediction_summaries_dic[pair]['all']['SVM_auc'] = []
        prediction_summaries_dic[pair]['all']['RF_auc'] = []

        prediction_summaries_dic[pair]['L1']['features'] = L1_features
        prediction_summaries_dic[pair]['L1']['n_features'] = len(L1_features)
        prediction_summaries_dic[pair]['L1']['SVM_auc'] = []
        prediction_summaries_dic[pair]['L1']['RF_auc'] = []

        prediction_summaries_dic[pair]['DT']['features'] = DT_features
        prediction_summaries_dic[pair]['DT']['n_features'] = len(DT_features)
        prediction_summaries_dic[pair]['DT']['SVM_auc'] = []
        prediction_summaries_dic[pair]['DT']['RF_auc'] = []

        
        prediction_summaries_dic[pair]['RFE']['features'] = RFE_features
        prediction_summaries_dic[pair]['RFE']['n_features'] = len(RFE_features)
        prediction_summaries_dic[pair]['RFE']['SVM_auc'] = []
        prediction_summaries_dic[pair]['RFE']['RF_auc'] = []
        
        prediction_summaries_dic[pair]['SF']['features'] = selected_feature_names
        prediction_summaries_dic[pair]['SF']['n_features'] = len(selected_feature_names)
        prediction_summaries_dic[pair]['SF']['SVM_auc'] = []
        prediction_summaries_dic[pair]['SF']['RF_auc'] = []
        

        for n in range(n_reps):
            
            print(f'processing [{I}] {phen1} vs {phen2} run {n} / {n_reps}')
            
            df = pd.DataFrame(data = X, columns = feature_names)
            df['class'] = y
            df_shuffled = df.sample(frac = 1).reset_index(drop = True)
            X_shuffled = df_shuffled.loc[ : , feature_names].values
            y_shuffled = df_shuffled['class']            
            all_SVM, all_RF = get_predictions_2(X_shuffled, y_shuffled)
            prediction_summaries_dic[pair]['all']['SVM_auc'].append(all_SVM)
            prediction_summaries_dic[pair]['all']['RF_auc'].append(all_RF)            
            
            L1_df = pd.DataFrame(data = X_L1, columns = L1_features)
            L1_df['class'] = y
            L1_df_shuffled = L1_df.sample(frac = 1).reset_index(drop = True)
            X_L1_shuffeled = L1_df_shuffled.loc[ : , L1_features].values
            y_L1_shuffled = L1_df_shuffled['class']
            L1_SVM, L1_RF = get_predictions_2(X_L1_shuffeled, y_L1_shuffled)
            prediction_summaries_dic[pair]['L1']['SVM_auc'].append(L1_SVM)
            prediction_summaries_dic[pair]['L1']['RF_auc'].append(L1_RF)
            
            DT_df = pd.DataFrame(data = X_DT, columns = DT_features)
            DT_df['class'] = y
            DT_df_shuffled = DT_df.sample(frac = 1).reset_index(drop = True)
            X_DT_shuffeled = DT_df_shuffled.loc[ : , DT_features].values
            y_DT_shuffled = DT_df_shuffled['class']
            DT_SVM, DT_RF = get_predictions_2(X_DT_shuffeled, y_DT_shuffled)
            prediction_summaries_dic[pair]['DT']['SVM_auc'].append(DT_SVM)
            prediction_summaries_dic[pair]['DT']['RF_auc'].append(DT_RF)
                        
            RFE_df = pd.DataFrame(data = X_RFE, columns = RFE_features)
            RFE_df['class'] = y
            RFE_df_shuffled = RFE_df.sample(frac = 1).reset_index(drop = True)
            X_RFE_shuffeled = RFE_df_shuffled.loc[ : , RFE_features].values
            y_RFE_shuffled = RFE_df_shuffled['class']
            RFE_SVM, RFE_RF = get_predictions_2(X_RFE_shuffeled, y_RFE_shuffled)
            prediction_summaries_dic[pair]['RFE']['SVM_auc'].append(RFE_SVM)
            prediction_summaries_dic[pair]['RFE']['RF_auc'].append(RFE_RF)
            
            selected_df = pd.DataFrame(data = X_selected, columns = selected_feature_names)
            selected_df['class'] = y_selected
            selected_df_shuffled = selected_df.sample(frac = 1).reset_index(drop = True)
            selected_X_shuffled = selected_df_shuffled.loc[ : , selected_feature_names].values
            selected_y_shuffled = selected_df_shuffled['class']            
            selected_SVM, selected_RF = get_predictions_2(selected_X_shuffled, selected_y_shuffled)
            prediction_summaries_dic[pair]['SF']['SVM_auc'].append(selected_SVM)
            prediction_summaries_dic[pair]['SF']['RF_auc'].append(selected_RF) 
            


with open(out_dir + 'prediction_summaries_dic.json', 'w') as out_f:
    json.dump(prediction_summaries_dic, out_f)
