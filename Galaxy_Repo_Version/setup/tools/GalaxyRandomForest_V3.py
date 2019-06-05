#!/usr/bin/env python

__author__ =     'Lukas Schmidt'
__maintainer__ = 'Lukas Schmidt'
__email__ =      'sclukas@students.uni-mainz.de'
__version__ =    '3.0'

import sys
sys.path.append('/usr/local/lib/python2.7/dist-packages')
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc
import random
import time
from math import *
import pickle
import warnings
import datetime

a = datetime.datetime.now()

warnings.filterwarnings('ignore')
time1 = (time.asctime(time.localtime(time.time())))


########################################################################################################################

def fill_list(path):
    """
    Fills a list with modified/non--modified entries of the given path
    :param path:
    :return:
    """
    my_list = []
    for i in range(len(path)):
        my_file = open(path[i], 'r')
        line = my_file.readline()
        line = my_file.readline()
        while len(line) > 3:
            my_list.append(line)
            line = my_file.readline()
        my_file.close()
    return my_list

########################################################################################################################


def random_forest(path_m1a, path_non_m1a, repetitions, splits, trees, outfile):
    """
    x repetitions of y-fold cross validation and random forest training + testing. First, the lists are randomly
    shuffled. Then, equal numbers of entries from both lists are written into the file used for training. Csv file
    is converted into pandas dataframe and type of column 'mod_type' is changed from str to int. 'm1A' is converted
    to 1, 'non-m1A' to 0 (this is required for a later step). The function then sets the predictive values (arrest-
    rate, mismatch-rate, etc.) and the target values (modification type aka m1A or nonm1A). The data is cross-
    validated using StratifiedKFold cross-validation. Data is split into training and testing sets and fed to the
    RandomForestClassifier
    :param path_m1a:
    :param path_non_m1a:
    :param repetitions:
    :param splits:
    :param trees:
    :param outfile:
    :return:
    """

    # Path to the output file comprised of a 1:1 ratio of m1A and non-m1A
    m1a_list = fill_list(path_m1a)
    non_m1a_list = fill_list(path_non_m1a)

    predictor_number = []
    for predic in predictors_in_use:
        predictor_number.append(predic)

    predictor_string = []
    for j in range(len(predictors_in_use)):
        if predictors_in_use[j] != 'pre_base':
            predictor_string.append(predictors_in_use[j])
    if pre_base:
        predictor_string.extend(['A', 'C', 'G', 'T'])
        predictor_number.extend(['A', 'C', 'G', 'T'])
        mean_feature_importance = [0] * (len(predictor_number) - 1)
    else:
        mean_feature_importance = [0] * len(predictor_number)

    # List for mean scores
    mean_sensitivity, mean_specificity, mean_ppv, mean_npv, mean_roc_auc, mean_mcc = [], [], [], [], [], []

    outfile.write('AUC' + '\t' + 'Sensitivity' + '\t' + 'Specificity' + '\t' + 'PPV' + '\t' + 'NPV' + '\t' +
                  'MCC' + '\t')

    predictors_in_use.append('mod_type')

    for j in range(repetitions):
        random.shuffle(m1a_list)
        random.shuffle(non_m1a_list)

        # Write equal numbers of m1As and non-m1As into a file
        temp_list = []
        for i in range(len(m1a_list)):
            temp_list.append(m1a_list[i].strip().split())
            temp_list.append(non_m1a_list[i].strip().split())

        # Build data pandas frame using all columns from the input file
        df = pd.DataFrame.from_records(temp_list, columns=predictor_features)
        # Remove columns that are not used
        for column in df.columns:
            if column not in predictors_in_use:
                df.drop(column, 1, inplace=True)

        # Change the modification type to numerical value
        df['mod_type'] = df['mod_type'].map({temp_list[0][-1]: 1, temp_list[1][-1]: 0})

        # Get categorical values (pre_base). This function creates 4 more columns in the pandas data frame (A, C, G, T).
        # Column 'pre_base' will be removed
        if pre_base:
            one_hot = pd.get_dummies(df['pre_base'])
            df.drop('pre_base', 1, inplace=True)
            df = df.join(one_hot)

        df_clean = df.dropna()
        df_clean.describe()

        # Use all values except for 'mod_type' as predictors
        predictors = df_clean[predictor_string]
        predictors = predictors.as_matrix()

        targets = df_clean.mod_type

        skf = StratifiedKFold(n_splits=splits, shuffle=True, random_state=None)
        forest = RandomForestClassifier(n_estimators=trees, criterion='gini', max_depth=None, max_features='sqrt',
                                        n_jobs=-1, warm_start=True, oob_score=True, random_state=None)

        splits_mean_roc, splits_sensitivity, splits_specificity, splits_ppv, splits_npv, splits_mcc = 0, 0, 0, 0, 0, 0

	if pre_base:
            temp_feature_importance = [0] * (len(predictor_number) - 1)
        else:
            temp_feature_importance = [0] * len(predictor_number)
	
	# Random forest training + testing
        for train, test in skf.split(predictors, targets):
            x_train, x_test = predictors[train], predictors[test]
            y_train, y_test = targets[train], targets[test]

            forest.fit(x_train, y_train)
            test_prediction = forest.predict(x_test)

            false_pos, true_pos, _ = roc_curve(y_test, test_prediction)
            roc_auc = auc(false_pos, true_pos)
            splits_mean_roc = splits_mean_roc + roc_auc * 100
            for k in range(len(forest.feature_importances_)):
                temp_feature_importance[k] = temp_feature_importance[k] + forest.feature_importances_[k]

            false_pos, true_pos, _ = roc_curve(y_test, test_prediction)

            # Build confusion matrix and calculate relevant values for statistical analysis
            cm = pd.crosstab(y_test, test_prediction, rownames=['Actual Class'], colnames=['Predicted Class'])
            TN = cm[0][0]
            FP = cm[0][1]
            FN = cm[1][0]
            TP = cm[1][1]
            sensitivity = (TP / (TP + FN)) * 100
            specificity = (TN / (FP + TN)) * 100
            ppv = (TP / (TP + FP)) * 100
            npv = (TN / (TN + FN)) * 100
            mcc = ((TP * TN - FP * FN) / (sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)))) * 100

            splits_sensitivity = splits_sensitivity + sensitivity
            splits_specificity = splits_specificity + specificity
            splits_ppv = splits_ppv + ppv
            splits_npv = splits_npv + npv
            splits_mcc = splits_mcc + mcc

        # Calculate the averages of n splits
        mean_sensitivity.append(splits_sensitivity / skf.n_splits)
        mean_specificity.append(splits_specificity / skf.n_splits)
        mean_ppv.append(splits_ppv / skf.n_splits)
        mean_npv.append(splits_npv / skf.n_splits)
        mean_mcc.append(splits_mcc / skf.n_splits)
        mean_roc_auc.append(splits_mean_roc / skf.n_splits)
        for l in range(len(temp_feature_importance)):
            mean_feature_importance[l] = mean_feature_importance[l] + temp_feature_importance[l] / skf.n_splits

    # Calculate the overall averages of x repetitions
    print('Sensitivity: ', sum(mean_sensitivity) / repetitions)
    print('specificity: ', sum(mean_specificity) / repetitions)
    print('Positive predicted value (PPV): ', sum(mean_ppv) / repetitions)
    print('Negative predicted value (NPV): ', sum(mean_npv) / repetitions)
    print('MCC: ', sum(mean_mcc) / repetitions)
    print('AUC: ', sum(mean_roc_auc) / repetitions)

    outfile.write(str((sum(mean_sensitivity) / repetitions)) + '\t' + str((sum(mean_specificity) / repetitions)) +
                     '\t' + str((sum(mean_ppv) / repetitions)) + '\t' + str((sum(mean_npv) / repetitions)) + '\t' +
                     str((sum(mean_mcc) / repetitions)) + '\t' + str((sum(mean_roc_auc) / repetitions)) + '\t')
    for j in range(len(mean_feature_importance)):
            outfile.write(str(mean_feature_importance[j] / repetitions) + '\t')
    outfile.write('\n')
    

    with open(sys.argv[4], 'wb') as f:
	pickle.dump(forest, f)

########################################################################################################################


m1as = [sys.argv[1]] #[sys.argv[1], sys.argv[2], sys.argv[3]]
non_m1as = [sys.argv[2]] 	#[sys.argv[4], sys.argv[5], sys.argv[6]]
outfile = open(sys.argv[3], 'w')	#open(sys.argv[7], 'w')

number_of_repetitions = int(sys.argv[5])	#int(sys.argv[8])
splits_cross_val = int(sys.argv[6])		#int(sys.argv[9])
number_of_trees = int(sys.argv[7])		#int(sys.argv[10])

outfile.write('Sensitivity' + '\t' + 'Specificity' + '\t' + 'PPV' + '\t' + 'NPV' + '\t' + 'MCC' + '\t' + 'AUC' +
                 '\t' + 'arrest_rate' + '\t' + 'mismatch_rate' + '\t' + 'cmism' + '\t' + 'gmism' + '\t' + 'tmism' + '\t'
                 + 'jump' + '\t' + 'pre_base' + '\n')


# List containing all possible features
predictor_features = ['arrest_rate', 'mism_rate', 'pre_base', 'amism', 'cmism', 'gmism', 'tmism', 'jump_rate_total']
# List containing all features to use in the RF-model
predictors_in_use = []
for i in range(8, len(sys.argv), 1):
    if sys.argv[i] == 'yes':
        predictors_in_use.append(predictor_features[i - 8])
predictor_features.append('mod_type')

print(predictors_in_use)
pre_base = True if sys.argv[10] == 'yes' else False

random_forest(m1as, non_m1as, number_of_repetitions, splits_cross_val, number_of_trees, outfile)


# forest_save = pickle._dumps(forest)

outfile.close()
b = datetime.datetime.now()
print(b-a)

# The End

