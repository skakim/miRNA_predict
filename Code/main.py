from GeneticAlgorithm import GeneticAlgorithm

from sklearn.neighbors import KNeighborsClassifier, NearestCentroid
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis, LinearDiscriminantAnalysis
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import SGDClassifier, LogisticRegression

import matplotlib.pyplot as plt
import numpy as np

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


def convert_binary(bitlist):
    out = 0
    for bit in bitlist:
        out = (out << 1) | bit
    return out

def create_heatmap(ga_cache, fit_attr, path):
    phenotypes = list(ga_cache.keys())
    phenotypes = sorted(phenotypes, key=lambda phenotype: getattr(ga_cache[phenotype], fit_attr))
    phenotypes_names = [''.join(map(str, phenotype)) +
                        " (" + format(getattr(ga_cache[phenotype], fit_attr), '.4f') + ")"
                        for phenotype in phenotypes]

    fitnesses = [getattr(ga_cache[phenotype], fit_attr) for phenotype in phenotypes]
    fitnesses = list(np.interp(fitnesses, (min(fitnesses), max(fitnesses)), (0.35, 0.75)))

    m = []
    for i in range(len(phenotypes)):
        phenotype = phenotypes[i]
        fitness = [fitnesses[i] if b else 1.0 for b in phenotype]
        m.append(fitness)

    #plt.matshow(m, cmap=plt.cm.hot)
    plt.imshow(m, interpolation="none", cmap=plt.cm.gist_heat)

    x_pos = np.arange(len(classifiers_names))
    plt.xticks(x_pos, classifiers_names, rotation=90)

    y_pos = np.arange(len(phenotypes))
    plt.yticks(y_pos, phenotypes_names)

    plt.tight_layout()
    plt.savefig(path)


def population_size():
    k = len(base_classifiers)
    return min((5*k), ((2**k)//2))

timestr = time.strftime("%Y%m%d-%H%M%S")
os.makedirs("results/{}".format(timestr))

with open("results/{}/log.txt".format(timestr), "w") as log_file:
    pop_size = population_size()
    fit_attr = 'AUC'
    GA = GeneticAlgorithm(base_classifiers=base_classifiers,
                          crossover_chance=0.6,
                          mutation_chance=0.01,
                          tournament_size=pop_size//10,
                          population_size=pop_size,
                          fit_attr=fit_attr)

    print("0", None, ','.join(map(str,GA.fitnesses())), sep=',')
    log_file.write("0,None," + ','.join(map(str,GA.fitnesses())) + '\n')
    log_file.flush()

    i = 1
    try:
        while i <= 30:
            GA.generate_next_population()
            print(i, GA.best_individual.phenotype, ','.join(map(str,GA.fitnesses())), sep=',')
            log_file.write(str(i) + "," + str(GA.best_individual.phenotype) + "," + ','.join(map(str, GA.fitnesses())) + '\n')
            log_file.flush()
            best_individual_filename = str(GA.fitnesses()[0]) + "-" + ''.join(map(str, GA.best_individual.phenotype))
            if not(os.path.exists("results/{}/{}.ens".format(timestr, best_individual_filename))):
                with open("results/{}/{}.ens".format(timestr, best_individual_filename), 'wb') as f:
                    pickle.dump(GA.best_individual.ensemble, f, pickle.HIGHEST_PROTOCOL)

            i += 1
    finally:
        # generate heatmap
        create_heatmap(GA.cache, fit_attr, "results/{}/heatmap.png".format(timestr))