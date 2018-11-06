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
                    #MLPClassifier(hidden_layer_sizes=(10), max_iter=50),
                    #MLPClassifier(hidden_layer_sizes=(15), max_iter=50),
                    #MLPClassifier(hidden_layer_sizes=(20), max_iter=50),
                    SVC(kernel='rbf', probability=True),
                    #SVC(kernel='sigmoid', probability=True),
                    KNeighborsClassifier(n_neighbors=3),
                    KNeighborsClassifier(n_neighbors=5),
                    KNeighborsClassifier(n_neighbors=7),
                    #SGDClassifier(loss='modified_huber'),
                    LogisticRegression()]

classifiers_names = ["GNB",
                     "DT (gini)",
                     "DT (entropy)",
                     "RF (gini)",
                     "RF (entropy)",
                     "QDA",
                     #"MLP 10",
                     #"MLP 15",
                     #"MLP 20",
                     "SVC (rbf)",
                     #"SVC (sigmoid)",
                     "3-KNN",
                     "5-KNN",
                     "7-KNN",
                     #"SGD",
                     "LogisticRegression"]


def convert_binary(bitlist):
    out = 0
    for bit in bitlist:
        out = (out << 1) | bit
    return out

def create_heatmap(ga_cache, fit_attr, path):
    #plt.figure(figsize=(654, 654))

    phenotypes = list(ga_cache.keys())
    phenotypes = sorted(phenotypes, key=lambda phenotype: getattr(ga_cache[phenotype], fit_attr))
    phenotypes_names = [''.join(map(str, phenotype)) +
                        " (" + format(getattr(ga_cache[phenotype], fit_attr), '.4f') + ")"
                        for phenotype in phenotypes]

    fitnesses = [getattr(ga_cache[phenotype], fit_attr) for phenotype in phenotypes]
    #fitnesses = list(np.interp(fitnesses, (min(fitnesses), max(fitnesses)), (0.35, 0.75)))

    m = []
    for i in range(len(phenotypes)):
        phenotype = phenotypes[i]
        fitness = [str(fitnesses[i]) if b else str(0.0) for b in phenotype]
        m.append(fitness)

    with open(path, "w") as f:
        f.write(fit_attr + "," + ",".join(classifiers_names) + "\n")
        for j in range(len(m)):
            f.write(phenotypes_names[j] + "," + ",".join(m[j]) + "\n")

    #plt.matshow(m, cmap=plt.cm.hot)
    #plt.imshow(m, interpolation="none", cmap=plt.cm.gist_heat)

    #x_pos = np.arange(len(classifiers_names))
    #plt.xticks(x_pos, classifiers_names, rotation=90)

    #y_pos = np.arange(len(phenotypes))
    #plt.yticks(y_pos, phenotypes_names)

    #plt.tight_layout()
    #plt.savefig(path)


def population_size():
    k = len(base_classifiers)
    return min((5*k), ((2**k)//2))

for j in range(0,10):
    for measure in ["AUC", "accuracy", "F1", "MCC"]:
        percentage_diversity = [0.0, 0.25, 0.5, 0.75]
        for div in percentage_diversity:
            timestr = time.strftime("%Y%m%d-%H%M%S")
            try:
                os.makedirs("results/final{}/{}+{}diversity".format(j, measure, div))
            except:
                pass

            print("============ MEASURE:", measure, div, "Diversity")

            with open("results/final{}/{}+{}diversity/log.txt".format(j, measure, div), "w") as log_file:
                pop_size = population_size()
                fit_attr = measure
                GA = GeneticAlgorithm(base_classifiers=base_classifiers,
                                      crossover_chance=0.6,
                                      mutation_chance=0.01,
                                      tournament_size=pop_size//10,
                                      population_size=pop_size,
                                      fit_attr=fit_attr,
                                      diversity_proportion=div)

                #print("0", None, ','.join(map(str,GA.fitnesses())), sep=',')
                #log_file.write("0,None," + ','.join(map(str,GA.fitnesses())) + '\n')
                #log_file.flush()

                i = 1
                try:
                    while i <= 10:
                        GA.generate_next_population()
                        print(i)
                        #log_file.write(str(i) + "," + str(GA.best_individual.phenotype) + "," + ','.join(map(str, GA.fitnesses())) + '\n')
                        log_file.write(str(i) + "," + str(GA.best_individual.accuracy) + "," +
                                       str(GA.best_individual.AUC) + "," +
                                       str(GA.best_individual.F1) + "," +
                                       str(GA.best_individual.MCC) + '\n')
                        log_file.flush()
                        #best_individual_filename = str(GA.fitnesses()[0]) + "-" + ''.join(map(str, GA.best_individual.phenotype))
                        #if not(os.path.exists("results/{}/{}.ens".format(timestr, best_individual_filename))):
                            #with open("results/{}/{}.ens".format(timestr, best_individual_filename), 'wb') as f:
                                #pickle.dump(GA.best_individual.ensemble, f, pickle.HIGHEST_PROTOCOL)

                        i += 1
                finally:
                    # generate heatmap
                    # create_heatmap(GA.cache, fit_attr, "results/{}/heatmap.csv".format(timestr))
                    pass