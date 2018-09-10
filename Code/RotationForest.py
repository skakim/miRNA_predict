# Original Author: Borja Ayerdi [ayerdi.borja -at- gmail -dot- com]
# Copyright(c) 2016
# License: Simple BSD

"""
This module implements Rotation Forest

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

#from sklearn import model_selection
from sklearn.decomposition import PCA

from sklearn.tree import DecisionTreeClassifier

np.seterr(divide='ignore', invalid='ignore')


__all__ = ["RotationForest"]


class RotationForest(object):
    """
    Rotation Forest
    """

    def __init__(self, classifiers=[]):
        self._rawclassifiers = classifiers
        self._classifiers = []
        self._inforotar = []
        self._std = []
        self._med = []
        self._noise = []

    @staticmethod
    def _apply_pca(data, labels, n_comps=1):
        """
        Applies PCA to the data

        :param data: ndarray
        A MxN array with M samples of dimension N

        :param labels: ndarray or list
        A 1xN array with the class labels of each sample
        of data

        :return: sklearn.decomposition.PCA

        """
        # PCA
        pca = PCA(n_components=n_comps, whiten=False, copy=True)
        pca.fit(data)

        return pca


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
        self._std = np.std(X,axis=0)
        self._med = np.mean(X,axis=0)
        self._noise = [random.uniform(-0.000005, 0.000005) for p in range(0,X.shape[1])]
        
        # Apply Z-score
        Xz = (X-self._med)/(self._std+self._noise)
        
        for classifier in self._rawclassifiers:
            # For each classifier in the ensemble
            # Given:
            # X: the objects in the training data set (an N x n matrix)
            # Y: the labels of the training set (an N x 1 matrix)
            # K: the number of subsets
            # NF: the number of total features
            # {w1,w2,.., wc}: the set of class labels
            #L

            # Downsample the dataset
            #Xdown, Ydown = downsample(Xz, Y)
            #Xdown, Ydown = oversample(Xz, Y)
            Xdown, Ydown = Xz, Y
            n_samps, NF = Xdown.shape

            # Prepare the rotaton matrix R:
            # Split F (the feature set) into K subsets Fij (for j=1,..,K/4)
            # K is a random value between 1 and NF/4.
            # We want at least 1 feature for each subset.
            K = int(round(1 + NF/4*random.random()))
            
            FK = np.zeros((K,NF))
            for j in range(K):
                numSelecFeatures = int(1 + round((NF-1)*random.random()))
                rp = np.random.permutation(NF)
                v = [rp[k] for k in range(0, numSelecFeatures)]
                FK[j,v] = 1
        
            
            R = np.zeros((NF,NF))
            for l in range(K):
                # Let Xzij be the data set X for the features in Fij
                pos = np.nonzero(FK[l,:])[0]

                vpos = [pos[m] for m in range(0, len(pos))]
       
                Xzij = Xdown[:, vpos]
                pca = self._apply_pca(Xzij, Y, len(pos))

                for indI in range(0,len(pca.components_)):
                    for indJ in range(0,len(pca.components_)):
                        R[pos[indI], pos[indJ]] = pca.components_[indI,indJ]

                        
            self._inforotar.append(R)
            Xrot = Xdown.dot(R)
            
            cl = classifier
            cl.fit(Xrot, Ydown)
            self._classifiers.append(cl)

        return self

    def predict(self, X):
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
        ensemble_output = np.zeros((len(X),dim))
        
        # Z-score
        X = (X-self._med)/(self._std+self._noise)

        for i in range(0,dim):
            xrot_z = X.dot(self._inforotar[i])
            ensemble_output[:,i] = self._classifiers[i].predict(xrot_z)
            #print(i, self._classifiers[i].predict(xrot_z))

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
