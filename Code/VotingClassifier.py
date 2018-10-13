"""
This module implements a custom Voting Classifier

References
----------
.. [1] Juan J. Rodriguez, et al, "Rotation Forest: A NewClassifier
          Ensemble Method", IEEE Transactions on Pattern Analysis and
          Machine Intelligence, 2006.

"""

import os
import sys
import random
import numpy as np
import pandas as pd

from sklearn.utils import shuffle
from scipy.stats import mode

# from sklearn import model_selection
from sklearn.decomposition import PCA

from sklearn.tree import DecisionTreeClassifier

np.seterr(divide='ignore', invalid='ignore')

__all__ = ["RotationForest"]


class VotingClassifier(object):
    """
    Heterogeneous Ensemble
    """

    def __init__(self, classifiers=[], z_score=False, downsample=False):
        self._rawclassifiers = classifiers
        self._classifiers = []
        self._inforotar = []
        self._std = []
        self._med = []
        self._noise = []
        self.z_score = z_score
        self.downsample = downsample

    def fit(self, X, Y):
        """
        Fit the model using X, y as training data.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape [n_samples, n_features]
        Training vectors, where n_samples is the number of samples
        and n_features is the number of features.

        y : array-like of shape [n_samples, n_outputs]
        Target values (class labels in classification, real numbers in
        regression)

        Returns
        -------
        self : object

        Returns an instance of self.
        """

        # Compute mean, std and noise for z-score
        if self.z_score:
            self._std = np.std(X, axis=0)
            self._med = np.mean(X, axis=0)
            self._noise = [random.uniform(-0.000005, 0.000005) for p in range(0, X.shape[1])]

            # Apply Z-score
            Xz = (X - self._med) / (self._std + self._noise)
        else:
            Xz = X

        for classifier in self._rawclassifiers:
            # For each classifier in the ensemble
            # Given:
            # X: the objects in the training data set (an N x n matrix)
            # Y: the labels of the training set (an N x 1 matrix)
            #print("Training", type(classifier).__name__)

            # Downsample the dataset
            if self.downsample:
                Xdown, Ydown = downsample(Xz, Y)
                # Xdown, Ydown = oversample(Xz, Y)
            else:
                Xdown, Ydown = Xz, Y

            cl = classifier
            cl.fit(Xdown, Ydown)
            self._classifiers.append(cl)

        return self

    def predict(self, X, proba=True, weights=None):
        #print(list(Y))
        #X, Y = downsample(X, Y)
        #print(list(Y))
        """
        Predict values using the model

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape [n_samples, n_features]

        Returns
        -------
        C : numpy array of shape [n_samples, n_outputs]
            Predicted values.
        """
        dim = len(self._classifiers)
        ensemble_output = np.zeros((len(X), dim))
        if self.z_score:
            # Z-score
            X = (X - self._med) / (self._std + self._noise)

        if proba:
            probas = np.asarray([clf.predict_proba(X) for clf in self._classifiers])
            y_avg = np.average(probas, axis=0, weights=weights)
            y_pred = np.argmax(y_avg, axis=1)

            """
            # generate heatmap
            classifiers_names = ["GNB",
                                 "DT (gini)",
                                 "DT (entropy)",
                                 "RF (gini)",
                                 "RF (entropy)",
                                 "QDA",
                                 "MLP 1-10",
                                 "MLP 1-5-5",
                                 "MLP 1-3-3-3",
                                 "SVC (rbf)",
                                 "SVC (sigmoid)",
                                 "3-KNN",
                                 "5-KNN",
                                 "7-KNN",
                                 "SGD",
                                 "LogisticRegression"]

            with open("results/heatmap_full.csv", 'w') as f:
                for i in range(len(self._classifiers)):
                    f.write(classifiers_names[i] + ';' +
                            ';'.join(map(str,self._classifiers[i].predict_proba(X)[:,1])) + '\n')
                f.write('EXPECTED' + ';' + ';'.join(map(str,Y)))
        """
        else:
            for i in range(0, dim):
                if proba:
                    #print(self._classifiers[i].predict_proba(X)[:,1])
                    ensemble_output[:, i] = self._classifiers[i].predict_proba(X)[:,1]
                else:
                    ensemble_output[:, i] = self._classifiers[i].predict(X)
                # print(self._classifiers[i].predict(X))
                # print(i, self._classifiers[i].predict(xrot_z))

            if proba:
                if not weights:
                    weights = [1.0]*dim
                y_avg = np.average(ensemble_output, axis=1, weights=weights)
                y_pred = np.around(y_avg)
            else:
                y_pred = mode(ensemble_output, axis=1)[0]


        return y_pred


def downsample(X, Y):
    # Convert to Dataframe
    data = np.concatenate((X, np.array([Y]).T), axis=1)
    df = pd.DataFrame(data)

    # Divide by class
    df_class_0 = df[df[34] == 0.0]
    df_class_1 = df[df[34] == 1.0]

    count_class_0 = len(df_class_0)
    count_class_1 = len(df_class_1)

    df_class_1_under = df_class_1.sample(count_class_0)
    df_test_under = shuffle(pd.concat([df_class_1_under, df_class_0], axis=0))

    #df_test_under = pd.concat([df_class_0, df_class_1], axis=0)

    Ydown = df_test_under[34]
    Xdown = df_test_under.loc[:, df_test_under.columns != 34]

    return Xdown.values, Ydown.values


def oversample(X, Y):
    # Convert to Dataframe
    data = np.concatenate((X, np.array([Y]).T), axis=1)
    df = pd.DataFrame(data)

    # Divide by class
    df_class_0 = df[df[34] == 0.0]
    df_class_1 = df[df[34] == 1.0]

    count_class_0 = len(df_class_0)
    count_class_1 = len(df_class_1)

    df_class_0_over = df_class_0.sample(count_class_1, replace=True)
    df_test_under = shuffle(pd.concat([df_class_0_over, df_class_1], axis=0))

    Ydown = df_test_under[34]
    Xdown = df_test_under.loc[:, df_test_under.columns != 34]

    return Xdown.values, Ydown.values
