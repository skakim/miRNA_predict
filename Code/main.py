from GeneticAlgorithm import GeneticAlgorithm

from sklearn.neighbors import KNeighborsClassifier, NearestCentroid
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis, LinearDiscriminantAnalysis
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import SGDClassifier, LogisticRegression

import os
import sys
import pickle
import time

base_classifiers = [GaussianNB(),
                    DecisionTreeClassifier(max_depth=5, criterion='gini'),
                    DecisionTreeClassifier(max_depth=5, criterion='entropy'),
                    RandomForestClassifier(max_depth=5, criterion='gini'),
                    RandomForestClassifier(max_depth=5, criterion='entropy'),
                    QuadraticDiscriminantAnalysis(),
                    MLPClassifier(hidden_layer_sizes=(1, 10), max_iter=50),
                    MLPClassifier(hidden_layer_sizes=(1, 5, 5), max_iter=50),
                    MLPClassifier(hidden_layer_sizes=(1, 3, 3, 3), max_iter=50),
                    SVC(kernel='rbf'),
                    SVC(kernel='sigmoid'),
                    KNeighborsClassifier(n_neighbors=3),
                    KNeighborsClassifier(n_neighbors=5),
                    KNeighborsClassifier(n_neighbors=7),
                    SGDClassifier(),
                    LogisticRegression()]


def population_size():
    k = len(base_classifiers)
    return min((5*k), ((2**k)//2))

timestr = time.strftime("%Y%m%d-%H%M%S")
os.makedirs("results/{}".format(timestr))

with open("results/{}/log.txt".format(timestr), "w") as log_file:
    pop_size = population_size()
    GA = GeneticAlgorithm(base_classifiers=base_classifiers,
                          crossover_chance=0.6,
                          mutation_chance=0.01,
                          tournament_size=pop_size//2,
                          population_size=pop_size,
                          fit_attr='AUC')

    print("0", None, ','.join(map(str,GA.fitnesses())), sep=',')
    log_file.write("0,None," + ','.join(map(str,GA.fitnesses())) + '\n')
    log_file.flush()

    i = 1
    while i <= 1000:
        GA.generate_next_population()
        print(i, GA.best_individual.phenotype, ','.join(map(str,GA.fitnesses())), sep=',')
        log_file.write(str(i) + "," + str(GA.best_individual.phenotype) + "," + ','.join(map(str, GA.fitnesses())) + '\n')
        log_file.flush()
        best_individual_filename = str(GA.fitnesses()[0]) + "-" + ''.join(map(str, GA.best_individual.phenotype))
        if not(os.path.exists("results/{}/{}.ens".format(timestr, best_individual_filename))):
            with open("results/{}/{}.ens".format(timestr, best_individual_filename), 'wb') as f:
                pickle.dump(GA.best_individual.ensemble, f, pickle.HIGHEST_PROTOCOL)
        i += 1