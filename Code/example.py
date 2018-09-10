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
from sklearn.metrics import confusion_matrix, f1_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import BaggingClassifier, VotingClassifier

from sklearn.neighbors import KNeighborsClassifier, NearestCentroid
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis, LinearDiscriminantAnalysis
from sklearn.neural_network import MLPClassifier
from RotationForest import RotationForest

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

    return X.values, Y.values

if __name__ == '__main__':


    # 1. Balance
    # 2. Breast-can
    # 3. Diabetes
    # 4. Ecoli
    # 5. Iris
    # 6. Liver
    # 7. Sonar
    # 8. Soybean
    # 9. Spambase
    # 10. Waveform
    # 11. Wine
    # 12. Digit
    # 13. Hayes
    # 14. Monk1
    # 15. Monk2
    # 16. Monk3

    for i in range(1):
        if i == 4:
            continue
        print("i={}".format(i))
        #uci_path = get_uci_path()
        #X, Y = read_uci_dataset(uci_path,i)

        X, Y = read_mirna_dataset()
        
        K = 5
        vAcc = []
        vAUC = []
        vF1 = []
        vTime = []
        rAcc = []
        rAUC = []
        rF1 = []
        rTime = []
        cv = StratifiedKFold(K)
        cv.get_n_splits(X,Y)

        train_indexes = []
        test_indexes = []
        
        for train, test in cv.split(X, Y):
            train_indexes.append(train)
            test_indexes.append(test)
            x_train, x_test, y_train, y_test = X[train,:], X[test,:], Y[train], Y[test]
            #print(len(x_train), len(x_test), len(x_train)/len(X))
          
            rotfor = RotationForest(classifiers=[GaussianNB(), DecisionTreeClassifier(max_depth=5, criterion='gini'), DecisionTreeClassifier(max_depth=5, criterion='entropy'),
                                                 QuadraticDiscriminantAnalysis(), MLPClassifier(hidden_layer_sizes=(1,10), max_iter=100), SVC(),
                                                 KNeighborsClassifier(n_neighbors=1), KNeighborsClassifier(n_neighbors=3), KNeighborsClassifier(n_neighbors=5)])
            init = time.time()
            rotfor = rotfor.fit(x_train, y_train)
            end = time.time()

            # For testing x_test & y_test
            y_pred = rotfor.predict(x_test)

            cm = confusion_matrix(y_test, y_pred)
            #print(cm)
            acc = float(cm.trace())/cm.sum()
            vAcc.append(acc)
            vAUC.append(roc_auc_score(y_test, y_pred))
            vF1.append(f1_score(y_test, y_pred))
            vTime.append(end-init)

            #ranfor = BaggingClassifier(n_estimators=9)
            ranfor = VotingClassifier(estimators=[('nb', GaussianNB()),
                                                  ('dtg', DecisionTreeClassifier(max_depth=5, criterion='gini')), ('dte', DecisionTreeClassifier(max_depth=5, criterion='entropy')),
                                                  ('qda', QuadraticDiscriminantAnalysis()), ('mlp', MLPClassifier(hidden_layer_sizes=(1,10), max_iter=100)),
                                                  ('svc', SVC()), ('nc', NearestCentroid()),
                                                  ('kn1', KNeighborsClassifier(n_neighbors=1)), ('kn3', KNeighborsClassifier(n_neighbors=3)), ('kn5', KNeighborsClassifier(n_neighbors=5))])
            init = time.time()
            ranfor = ranfor.fit(x_train, y_train)
            end = time.time()

            # For testing x_test & y_test
            y_pred = ranfor.predict(x_test)
            
            cm = confusion_matrix(y_test, y_pred)
            acc = float(cm.trace())/cm.sum()
            rAcc.append(acc)
            rAUC.append(roc_auc_score(y_test, y_pred))
            rF1.append(f1_score(y_test, y_pred))
            rTime.append(end-init)

        print("algr, acc_avg, acc_std, auc_avg, auc_std, f1_avg, f1_std, t_avg, t_std")
        bd_std = np.std(vAcc)
        bd_acc = np.mean(vAcc)
        auc_std = np.std(vAUC)
        auc_acc = np.mean(vAUC)
        f1_std = np.std(vF1)
        f1_acc = np.mean(vF1)
        t_std = np.std(vTime)
        t_acc = np.mean(vTime)
        print("rotfor", format(bd_acc*100,'.2f'), format(bd_std*100,'.2f'),
              format(auc_acc * 100, '.2f'), format(auc_std * 100, '.2f'),
              format(f1_acc * 100, '.2f'), format(f1_std * 100, '.2f'),
              format(t_acc,'.2f'), format(t_std,'.2f'), sep=',')

        bd_std = np.std(rAcc)
        bd_acc = np.mean(rAcc)
        auc_std = np.std(rAUC)
        auc_acc = np.mean(rAUC)
        f1_std = np.std(rF1)
        f1_acc = np.mean(rF1)
        t_std = np.std(rTime)
        t_acc = np.mean(rTime)
        print("voting", format(bd_acc * 100, '.2f'), format(bd_std * 100, '.2f'),
              format(auc_acc * 100, '.2f'), format(auc_std * 100, '.2f'),
              format(f1_acc * 100, '.2f'), format(f1_std * 100, '.2f'),
              format(t_acc, '.2f'), format(t_std, '.2f'), sep=',')
