#!/usr/bin/env python

__author__ =     'Lukas Schmidt'
__maintainer__ = 'Lukas Schmidt'
__email__ =      'sclukas@students.uni-mainz.de'
__version__ =    '2.0'


import sys
sys.path.append('/usr/local/lib/python2.7/dist-packages')
import pandas as pd
sys.path.append('/usr/local/lib/python3.5/dist-packages')
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc
import random
from sklearn.metrics import roc_auc_score
import time
from math import *
import pickle
import warnings
import datetime

a = datetime.datetime.now()

warnings.filterwarnings("ignore")
time1 = (time.asctime(time.localtime(time.time())))

########################################################################################################################

'''Counts the number of m1As within a given dataset (in which m1A is defined by the number "1") '''
def count_m1a(dataset):
    counter = 0
    for i in range(len(dataset)):
        if dataset[i] == 1:
            counter += 1
    return counter

########################################################################################################################

'''Fills a list with m1A/non-m1A-entries of the given path'''
def fill_list(path):
    letters = ['A', 'C', 'G', 'T']
    my_list = []
    for i in range(len(path)):
        my_file = open(path[i], 'r')
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
def random_forest(path_m1a, path_non_m1a):
    repetitions = int(sys.argv[8])

    # List for instances of m1A and non-m1A
    # if given arguments are lists
    m1a_list = fill_list(path_m1a)[0]
    non_m1a_list, header = fill_list(path_non_m1a)[0], fill_list(path_non_m1a)[1]
    predictor_number = ((header.split())[0:-1])

    predictor_string = []
    for j in range(len(header.split())):
        if header.split()[j] != 'mod_type' and header.split()[j] != 'preBase':
            predictor_string.append(header.split()[j])

    # List for mean scores
    mean_roc_auc = []
    mean_feature_importance = [0] * len(predictor_number)

    ''' x repetitions of y-fold cross validation and random forest training + testing. First, the lists are randomly
    shuffled. Then, equal numbers of entries from both lists are written into the file used for training. Csv file
    is converted into pandas dataframe and type of column 'mod_type' is changed from str to int. 'm1A' is converted
    to 1, 'non-m1A' to 0 (this is required for a later step). The function then sets the predictive values (arrest-
    rate, mismatch-rate, etc.) and the target values (modification type aka m1A or nonm1A). The data is cross-
    validated using StratifiedKFold cross-validation. Data is split into training and testing sets and fed to the
    RandomForestClassifier
    '''
    for j in range(repetitions):
        random.shuffle(m1a_list)
        random.shuffle(non_m1a_list)

        ratio_list = []
        # Write equal numbers of m1As and non-m1As into a file
        for i in range(len(m1a_list)):
            ratio_list.append(m1a_list[i].strip().split('\t'))
            ratio_list.append(non_m1a_list[i].strip().split('\t'))

        df = pd.DataFrame.from_records(ratio_list, columns=header.split())
        df['mod_type'] = df['mod_type'].map({ratio_list[0][-1]: 1, ratio_list[1][-1]: 0})
        df_clean = df.dropna()
        df_clean.describe()

        predictors = df_clean[predictor_string]

        predictors = predictors.as_matrix()
        targets = df_clean.mod_type

        skf = StratifiedKFold(n_splits=int(sys.argv[9]), shuffle=True, random_state=None)
        forest = RandomForestClassifier(n_estimators=int(sys.argv[10]), criterion='gini', max_depth=None, max_features='sqrt',
                                        n_jobs=-1, warm_start=True, oob_score=True, random_state=None)
        # max_features='sqrt'
        mean = 0
        temp_feature_importance = [0] * len(predictor_number)
        for train, test in skf.split(predictors, targets):
            x_train, x_test = predictors[train], predictors[test]
            y_train, y_test = targets[train], targets[test]

            forest.fit(x_train, y_train)
            test_prediction = forest.predict(x_test)

            roc_auc = roc_auc_score(y_test, test_prediction)
            false_pos, true_pos, _ = roc_curve(y_test, test_prediction)
            roc_auc = auc(false_pos, true_pos)
            mean = mean + roc_auc
            for k in range(len(forest.feature_importances_)):
                temp_feature_importance[k] = temp_feature_importance[k] + forest.feature_importances_[k]

        mean = mean / skf.n_splits
        mean_roc_auc.append(mean)
        for l in range(len(temp_feature_importance)):
            mean_feature_importance[l] = mean_feature_importance[l] + temp_feature_importance[l] / skf.n_splits
        # print(mean_feature_importance)

    auc_mean = 0
    for i in range(len(mean_roc_auc)):
        auc_mean = auc_mean + mean_roc_auc[i]
    print("AUC: ", auc_mean / repetitions)
    for j in range(len(mean_feature_importance)):
        mean_feature_importance[j] = mean_feature_importance[j] / repetitions
    print("Mean feature importance: ", mean_feature_importance)

    pickle.dump(forest, open(sys.argv[7], 'wb'), protocol=2)


path_m1as = [sys.argv[1], sys.argv[2], sys.argv[3]]
path_non_m1as = [sys.argv[4], sys.argv[5], sys.argv[6]]

random_forest(path_m1as, path_non_m1as)

b = datetime.datetime.now()
print(b-a)
# The End

