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
import warnings
import datetime

a = datetime.datetime.now()

warnings.filterwarnings("ignore")
time1 = (time.asctime(time.localtime(time.time())))


########################################################################################################################

'''Fills a list with m1A/non-m1A-entries of the given path'''
def fill_list(path):
    my_list = []
    for i in range(len(path)):
        my_file = open(path[i], "r")
        line = my_file.readline()
        header = line
        line = my_file.readline()
        while len(line) > 3:
            my_list.append(line)
            line = my_file.readline()
        my_file.close()
    return my_list, header

########################################################################################################################


''' Function implements a random forest classifier for determination of specific RNA-modifications within an
RNA-sequence
'''
def random_forest(path_m1a, path_non_m1a, repetitions, splits, trees, outfile):
    # repetitions = 100

    # Path to the output file comprised of a 1:1 ratio of m1A and non-m1A
    m1a_list = fill_list(path_m1a)[0]
    non_m1a_list, header = fill_list(path_non_m1a)[0], fill_list(path_non_m1a)[1]

    predictor_number = ((header.split())[0:-1])
    predictor_string = ''
    predictor_string = [predictor_string + predictor_number[q] for q in range(len(predictor_number))]

    # List for mean scores
    mean_roc_auc = []
    mean_feature_importance = [0] * len(predictor_number)
    
    outfile.write('AUC' + '\t')
    for k in range(len(predictor_string)):
        outfile.write(predictor_string[k] + '\t')
    outfile.write('\n')

    ''' x repetitions of y-fold cross validation and random forest training + testing. First, the lists are randomly
    shuffled. Then, equal numbers of entries from both lists are written into the file used for training. Csv file
    is converted into pandas dataframe and type of column 'mod_type' is changed from str to int. 'm1A' is converted
    to 1, 'non-m1A' to 0 (this is required for a later step). The function then sets the predictive values (arrest-
    rate, mismatch-rate, etc.) and the target values (modification type aka m1A or nonm1A). The data is cross-
    validated using StratifiedKFold cross-validation. Data is split into training and testing sets and fed to the
    RandomForestClassifier
    '''
    for j in range(repetitions):
        #if j % 50 == 0:
        #    print('Repetition: ', j)
        random.shuffle(m1a_list)
        random.shuffle(non_m1a_list)

        # Write equal numbers of m1As and non-m1As into a file
        temp_list = []
        for i in range(len(m1a_list)):
            temp_list.append(m1a_list[i].strip().split())
            temp_list.append(non_m1a_list[i].strip().split())

        df = pd.DataFrame.from_records(temp_list, columns=header.split())
        # Change the modification type to numerical value
        df['mod_type'] = df['mod_type'].map({temp_list[0][-1]: 1, temp_list[1][-1]: 0})
        df_clean = df.dropna()
        df_clean.describe()

        # Use all values except for 'mod_type' as predictors
        predictors = df_clean[predictor_string]

        predictors = predictors.as_matrix()
        targets = df_clean.mod_type

        skf = StratifiedKFold(n_splits=splits, shuffle=True, random_state=None)
        forest = RandomForestClassifier(n_estimators=trees, criterion='gini', max_depth=None, max_features='sqrt',
                                        n_jobs=-1, warm_start=True, oob_score=True, random_state=None)
        mean = 0
        temp_feature_importance = [0] * len(predictor_number)
        for train, test in skf.split(predictors, targets):
            x_train, x_test = predictors[train], predictors[test]
            y_train, y_test = targets[train], targets[test]
            forest.fit(x_train, y_train)

            # Make prediction and get the number of false positives and true positives. Calculate the AUC-value
            test_prediction = forest.predict(x_test)

            false_pos, true_pos, _ = roc_curve(y_test, test_prediction)
            roc_auc = auc(false_pos, true_pos)
            mean = mean + roc_auc

            for k in range(len(forest.feature_importances_)):
                temp_feature_importance[k] = temp_feature_importance[k] + forest.feature_importances_[k]

        mean = mean / skf.n_splits
        mean_roc_auc.append(mean)
        for l in range(len(temp_feature_importance)):
            mean_feature_importance[l] = mean_feature_importance[l] + temp_feature_importance[l] / skf.n_splits

    auc_mean = 0
    for i in range(len(mean_roc_auc)):
        auc_mean = auc_mean + mean_roc_auc[i]

    auc_final = auc_mean/repetitions
    print(auc_final)
    outfile.write(str(auc_final) + '\t')
    print(mean_feature_importance)
    for j in range(len(mean_feature_importance)):
        mean_feature_importance[j] = mean_feature_importance[j] / repetitions
        outfile.write(str(mean_feature_importance[j]) + '\t')

    outfile.write('\n')


m1as = [sys.argv[1], sys.argv[3], sys.argv[5]] #[sys.argv[1], sys.argv[2], sys.argv[3]]
non_m1as = [sys.argv[2], sys.argv[4], sys.argv[6]] 	#[sys.argv[4], sys.argv[5], sys.argv[6]]
outfile = open(sys.argv[7], 'w')	#open(sys.argv[7], 'w')

number_of_repetitions = int(sys.argv[8])	#int(sys.argv[8])
splits_cross_val = int(sys.argv[9])		#int(sys.argv[9])
number_of_trees = int(sys.argv[10])		#int(sys.argv[10])

random_forest(m1as, non_m1as, number_of_repetitions, splits_cross_val, number_of_trees, outfile)

# forest_save = pickle._dumps(forest)

outfile.close()
b = datetime.datetime.now()
print(b-a)
# The End
