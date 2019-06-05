#!/usr/bin/python

__author__ = 'Lukas Schmidt'
__maintainer__ = 'Lukas Schmidt'
__email__ = 'sclukas@students.uni-mainz.de'
__version__ = '2.0'

# Steps in this order: pileup2profile -> FilterByBase -> AnnotateProfile -> PredictModifications
# This programme uses the information of a trained random forest classifier to predict modifications on a given profile.

import sys
sys.path.append('/usr/local/lib/python2.7/dist-packages')
import pandas as pd
#sys.path.append('/usr/local/lib/python3.5/dist-packages')
import pickle

########################################################################################################################

def predict_modifications(data, forest_model):
    """
    Function writes those candidates into a new file, which were predicted to be m1As (eg. entries with value 1)
    :param data: Input data frame with relevant information needed for prediction.
    :param forest_model: Trained random forest classifier.
    :return: Numpy array containing the predictions for each line of the input.
    """
    prediction = forest_model.predict(data)
    return prediction


########################################################################################################################

def write_to_outfile(splitted_line_1, splitted_line_2, mod_type):
    """
    Writes the output.
    :param splitted_line_1: Array from the annotated profile file containing all information of a given line.
    :param splitted_line_2: Array from the original profile file containing all information of a given line.
    :param mod_type: String containing the name of the modification-type.
    :return:
    """
    fout.write(
        splitted_line_2[0] + '\t' + splitted_line_2[1] + '\t' + splitted_line_2[2] + '\t' + splitted_line_2[3] + '\t' +
        splitted_line_2[4] + '\t' + splitted_line_1[0] + '\t' + splitted_line_1[1] + '\t' + splitted_line_1[3] + '\t' +
        splitted_line_1[4] + '\t' + splitted_line_1[5] + '\t' + splitted_line_1[6] + '\t' + splitted_line_1[7] + '\t' +
        mod_type + '\n'
    )


########################################################################################################################

def prepare_results(prediction):
    """
    Prepares the string containing the information on each position and the predicted modification for the output. Calls
    the 'write_to_outfile'-function.
    :param prediction: Numpy array containing the predictions on the modification for each position. The predictions are
    encoded as 0 or 1 (not modified and modified).
    """
    with open(sys.argv[1], 'r') as fin1, open(sys.argv[2], 'r') as fin2:
        fin1 = fin1.readlines()[1:]  # Read in the input file (annotated profile), leave out the header.
        fin2 = fin2.readlines()[1:]  # Original profile-file
        new_prediction = prediction.tolist()  # Convert format of the Numpy array to a list.
        for line1, mod_type, line2 in zip(fin1, new_prediction, fin2):  # Iterate through the input files and the predictions.
            split_line_1 = line1.strip().split('\t')[0:-1]
            split_line_2 = line2.strip().split('\t')
            if mod_type == 0:  # Not modified
                write_to_outfile(split_line_1, split_line_2, "'" + 'non' + sys.argv[5] + "'")  # Name of unmodified base
            elif mod_type == 1:  # Modified
                write_to_outfile(split_line_1, split_line_2, "'" + sys.argv[5] + "'")  # Name of modified base
            else:
                write_to_outfile(split_line_1, split_line_2, 'mod_type')


########################################################################################################################

with open(sys.argv[4], 'w') as fout, open(sys.argv[3], 'rb') as f:  # Outfile and RF-save-file in pkl-format
    # fout.write('arrest_rate\tmism_rate\tgmism\ttmism\tcmism\tjump_rate_total\tmod_type\n')
    fout.write('ref_seg\tpos\tref_base\tcov\tpre_base\tarrest_rate\tmism_rate\tamism\tcmism\tgmism\ttmism\t'
               'jump_rate_total\tpredicted_mod_type\n')

    #u = pickle._Unpickler(f)
    #u.encoding = 'latin1'
    random_forest = pickle.load(f)

    # List containing all possible features
    predictor_features = ['arrest_rate', 'mism_rate', 'pre_base', 'amism', 'cmism', 'gmism', 'tmism', 'jump_rate_total']
    # List containing all features to use in the RF-model
    predictors_in_use = []
    for i in range(6, len(sys.argv), 1):
        if sys.argv[i] == 'yes':
            predictors_in_use.append(predictor_features[i - 6])
    pre_base = True if sys.argv[8] == 'yes' else False

    # Read input and translate into pandas data frame
    df = pd.read_csv(sys.argv[1], delimiter='\t')

    # Remove unused columns
    for column in df.columns:
        if column not in predictors_in_use:
            df.drop(column, 1, inplace=True)

    # Get categorical values (pre_base). This function creates 4 more columns in the pandas data frame (A, C, G, T).
    # Column 'pre_base' will be removed
    if pre_base:
        one_hot = pd.get_dummies(df['pre_base'])
        df.drop('pre_base', 1, inplace=True)
        df = df.join(one_hot)

    # random_forest = pickle.load(f)  # Load the random forest save file.
    # Prepare the input data
    # d1 = pd.read_csv(sys.argv[1], delimiter='\t', nrows=1)
    # columns = d1.columns.tolist()
    # print(columns)
    # cols_to_use = columns[:len(columns)-1]
    # print(cols_to_use)
    # df = pd.read_csv(sys.argv[1], delimiter='\t', usecols=cols_to_use)
    result = predict_modifications(df, random_forest)
    prepare_results(result)

# The End

'''
/path/MH1503.profile_filteredbybase_annotated
/path/MH1503.profile_filteredbybase
/path/RT3_RF_test_no_amism.pkl
/path/MH1503.profile_results
m1A
yes     # arrest
yes     # mismatch
yes     # pre_base
no      # amism
yes     # cmism
yes     # gmism
yes     # tmism
yes     # jump
'''
