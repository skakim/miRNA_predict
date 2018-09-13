#!/usr/bin/python

import os
import sys
import time
import random
import numpy as np
import pandas as pd
import warnings

from sklearn.utils import shuffle

from sklearn import model_selection
from sklearn.metrics import confusion_matrix, f1_score, matthews_corrcoef, roc_auc_score
from sklearn.model_selection import StratifiedKFold, train_test_split
#from sklearn.ensemble import BaggingClassifier, VotingClassifier

from sklearn.neighbors import KNeighborsClassifier, NearestCentroid
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis, LinearDiscriminantAnalysis
from sklearn.neural_network import MLPClassifier
from VotingClassifier import VotingClassifier

warnings.simplefilter("ignore")

def get_uci_path():
    """
    Returns the path to the UCI datasets, depending on the 
    host name.

    :return: string
    """
    path = 'uci-datasets'

    return path


def read_uci_dataset(base_dir, dataset_idx=1):
    """
    This function returns the path to the UCI dataset

    :param base_dir: string
    The path to where the UCI datasets folder are.

    :param dataset_idx: int
    The number of the dataset

    :return: 2-tuple of ndarray
    The dataset in a Numpy array,
    the first is the data samples and the second, the labels
    """

    # Load the data
    path = os.path.join(base_dir, str(dataset_idx) + '.data')

    data = np.genfromtxt(path, delimiter=",")
    rows, cols = data.shape

    # Delete the first column of labels
    Y = data[:, 0]
    X = data[:, 1:rows]

    return X, Y

def read_mirna_dataset():

    """
    This function reads the miRNA dataset

    :return: 2-tuple of ndarray
    The dataset in a Numpy array,
    the first is the data samples and the second, the labels
    """

    # Load the data
    path = '../Analysis/attributes/final_dataset.txt'

    with open(path) as dataset:
        data = pd.read_table(dataset, sep='\t')

    print(len(data), len(data.columns))
    data = data.drop_duplicates()
    data = shuffle(data)
    print(len(data), len(data.columns))

    # convert neg and pos to 0 and 1
    data.Class = pd.Categorical(data.Class)
    data.Class = data.Class.cat.codes

    Y = data["Class"]
    X = data.loc[:, data.columns != 'Class']
    #Y = data[:, 0]
    #X = data[:, 1:rows]

    #print(X.dtypes)
    #print(Y.dtypes)
    #sys.exit()

    return X.values, Y.values

def run_ensemble(X, Y, classifiers_list = [], dry_run=False):

    K = 2
    vAcc = []
    vAUC = []
    vF1 = []
    vMCC = []
    vTime = []
    cv = StratifiedKFold(K)
    cv.get_n_splits(X, Y)
    k = 1
    for train, test in cv.split(X, Y):
        if dry_run:
            print("Fold", k)
        x_train, x_test, y_train, y_test = X[train, :], X[test, :], Y[train], Y[test]
        ensemble = VotingClassifier(classifiers=classifiers_list, z_score=False, downsample=False)

        if dry_run:
            print("Ensemble...", flush=True)
        init = time.time()
        ensemble = ensemble.fit(x_train, y_train)
        end = time.time()
        if dry_run:
            print(format(end - init, '.2f'), "s", sep='')

        # For testing x_test & y_test
        y_pred = ensemble.predict(x_test)

        cm = confusion_matrix(y_test, y_pred)
        acc = float(cm.trace()) / cm.sum()
        AUC = roc_auc_score(y_test, y_pred)
        F1 = f1_score(y_test, y_pred)
        MCC = matthews_corrcoef(y_test, y_pred)
        vAcc.append(acc)
        vAUC.append(AUC)
        vF1.append(F1)
        vMCC.append(MCC)
        vTime.append(end - init)

        if dry_run:
            print("cl_name, acc, AUC, F1, MCC")
            print(cm)
            print("Ensemble", format(acc * 100, '.2f'), format(AUC * 100, '.2f'), format(F1 * 100, '.2f'),
                  format(MCC, '.2f'), sep=',')

            for classifier in ensemble._classifiers:
                y_pred = classifier.predict(x_test)
                cm = confusion_matrix(y_test, y_pred)
                print(cm)
                acc = float(cm.trace()) / cm.sum()
                AUC = roc_auc_score(y_test, y_pred)
                F1 = f1_score(y_test, y_pred)
                MCC = matthews_corrcoef(y_test, y_pred)
                print(type(classifier).__name__, format(acc * 100, '.2f'), format(AUC * 100, '.2f'), format(MCC, '.2f'),
                      sep=',')

            print("algr, acc_avg, acc_std, auc_avg, auc_std, f1_avg, f1_std, mcc_avg, mcc_std, t_avg, t_std")

    bd_std = np.std(vAcc)
    bd_acc = np.mean(vAcc)
    auc_std = np.std(vAUC)
    auc_acc = np.mean(vAUC)
    f1_std = np.std(vF1)
    f1_acc = np.mean(vF1)
    mcc_std = np.std(vMCC)
    mcc_acc = np.mean(vMCC)
    t_std = np.std(vTime)
    t_acc = np.mean(vTime)
    if dry_run:
        print("ensemble", format(bd_acc * 100, '.2f'), format(bd_std * 100, '.2f'),
              format(auc_acc * 100, '.2f'), format(auc_std * 100, '.2f'),
              format(f1_acc * 100, '.2f'), format(f1_std * 100, '.2f'),
              format(mcc_acc, '.2f'), format(mcc_std, '.2f'),
              format(t_acc, '.2f'), format(t_std, '.2f'), sep=',')

    return ensemble, bd_acc, auc_acc, f1_acc, mcc_acc

if __name__ == '__main__':
    X, Y = read_mirna_dataset()
    run_ensemble(X, Y,
                 classifiers_list = [GaussianNB(),
                                DecisionTreeClassifier(max_depth=3, criterion='gini'),
                                DecisionTreeClassifier(max_depth=3, criterion='entropy'),
                                QuadraticDiscriminantAnalysis(),
                                MLPClassifier(hidden_layer_sizes=(1, 10), max_iter=25),
                                MLPClassifier(hidden_layer_sizes=(1, 5, 5), max_iter=25),
                                MLPClassifier(hidden_layer_sizes=(1, 3, 3, 3), max_iter=25),
                                SVC(kernel='rbf'),
                                SVC(kernel='sigmoid'),
                                KNeighborsClassifier(n_neighbors=3),
                                KNeighborsClassifier(n_neighbors=5),
                                KNeighborsClassifier(n_neighbors=7)],
                 dry_run=True)