import sys
import numpy as np
import pandas as pd
from copy import copy, deepcopy
from operator import attrgetter
from sklearn.utils import shuffle
from sklearn.model_selection import StratifiedKFold

from Ensemble import run_ensemble

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

    #print(len(data), len(data.columns))
    data = data.drop_duplicates()
    data = shuffle(data)
    #print(len(data), len(data.columns))

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

class GeneticAlgorithm(object):
    def __init__(self, base_classifiers, crossover_chance, mutation_chance, tournament_size, population_size, fit_attr, diversity_proportion):
        self.base_classifiers = base_classifiers
        self.crossover_chance = crossover_chance
        self.mutation_chance = mutation_chance
        self.tournament_size = tournament_size
        self.population_size = population_size
        self.fit_attr = fit_attr
        self.diversity_proportion = diversity_proportion

        self.X, self.Y = read_mirna_dataset()

        self.cache = {}

        self.population = self.generate_initial_population()
        self.best_individual = None
        self.generation = 0

    def fitnesses(self):
        l = []
        for ind in self.population:
            l.append(getattr(ind, self.fit_attr))
        return sorted(l, reverse=True)

    def generate_initial_population(self):
        print("Creating initial generation")
        population = []
        for i in range(self.population_size):
            random_initial = list(np.random.randint(2, size=len(self.base_classifiers)))
            if tuple(random_initial) not in self.cache:
                ind = Individual(random_initial, self.base_classifiers, self.X, self.Y, self.fit_attr, self.diversity_proportion)
                print("Finished creating individual", i, ind.phenotype, ind.fitness)
            else:
                ind = self.cache[tuple(random_initial)]
                print("Cache used for individual", i, ind.phenotype, ind.fitness)
            population.append(ind)
            if tuple(ind.phenotype) not in self.cache:
                self.cache[tuple(ind.phenotype)] = ind
        return population

    def generate_next_population(self):
        old_population = self.population
        new_population = []

        # test if best_individual needs to be updated
        best_from_population = max(old_population, key=attrgetter('fitness'))
        if not(self.best_individual) or \
                (getattr(best_from_population, 'fitness') >= getattr(self.best_individual, 'fitness') and
                 best_from_population.phenotype != self.best_individual.phenotype):
            self.best_individual = deepcopy(best_from_population)

        # delete ensembles from cache
        # release memory (no need to maintain the ensembles in memory (the best one was already copied above)
        for phen in self.cache:
            try:
                del self.cache[phen].ensemble
            except:
                pass

        print("Generation {} completed".format(self.generation))
        self.generation += 1
        print("Creating generation {}".format(self.generation))

        # generate offspring
        # elitism
        new_population.append(self.best_individual)
        print("Copied best individual to new generation", self.best_individual.phenotype)

        # let the games begin!
        for i in range(self.population_size//2):
            # crossover
            ind1, ind2 = self.selTournament()
            phen1, phen2 = self.crossover(ind1, ind2)

            # mutation
            phen1 = self.mutation(phen1)
            phen2 = self.mutation(phen2)

            # retrieve from cache or create
            if tuple(phen1) in self.cache:
                new1 = self.cache[tuple(phen1)]
                print("Cache used for individual", (i*2) + 1, new1.phenotype, new1.fitness)
            else:
                new1 = Individual(phen1, self.base_classifiers, self.X, self.Y, self.fit_attr, self.diversity_proportion)
                print("Finished creating individual", (i*2) + 1, new1.phenotype, new1.fitness)
                if tuple(new1.phenotype) not in self.cache:
                    self.cache[tuple(new1.phenotype)] = new1

            if tuple(phen2) in self.cache:
                new2 = self.cache[tuple(phen2)]
                print("Cache used for individual", (i*2) + 2, new2.phenotype, new2.fitness)
            else:
                new2 = Individual(phen2, self.base_classifiers, self.X, self.Y, self.fit_attr, self.diversity_proportion)
                print("Finished creating individual", (i*2) + 2, new2.phenotype, new2.fitness)
                if tuple(new2.phenotype) not in self.cache:
                    self.cache[tuple(new2.phenotype)] = new2

            new_population.append(new1)
            new_population.append(new2)

        self.population = new_population

    def crossover(self, ind1, ind2):
        r = np.random.rand()
        if r > self.crossover_chance:
            return ind1.phenotype, ind2.phenotype

        phen1 = []
        phen2 = []

        # Uniform Crossover
        for i in range(len(ind1.phenotype)):
            r = np.random.rand()
            if r < 0.5:
                phen1.append(ind1.phenotype[i])
                phen2.append(ind2.phenotype[i])
            else:
                phen1.append(ind2.phenotype[i])
                phen2.append(ind1.phenotype[i])

        return phen1, phen2

    def mutation(self, phen):
        for i in range(len(phen)):
            r = np.random.rand()
            if r < self.mutation_chance:
                # binary NOT
                phen[i] = 1 - phen[i]
        return phen

    def selRandom(self, individuals, k):
        """Select *k* individuals at random from the input *individuals* with
        replacement. The list returned contains references to the input
        *individuals*.

        :param individuals: A list of individuals to select from.
        :param k: The number of individuals to select.
        :returns: A list of selected individuals.

        This function uses the :func:`~random.choice` function from the
        python base :mod:`random` module.
        """
        return np.random.choice(individuals, size=k)

    def selTournament(self, k=2):
        """Select the best individual among *tournsize* randomly chosen
        individuals, *k* times. The list returned contains
        references to the input *individuals*.

        :param self.population: A list of individuals to select from.
        :param k: The number of individuals to select.
        :param self.tournament_size: The number of individuals participating in each tournament.
        :param self.fit_attr: The attribute of individuals to use as selection criterion (accuracy, AUC, F1 or MCC)
        :returns: A list of selected individuals.
        """
        population_copy = copy(self.population)
        chosen = []
        # choose first
        aspirants = self.selRandom(population_copy, self.tournament_size)
        chosen1 = max(aspirants, key=attrgetter('fitness'))
        chosen.append(chosen1)
        # choose second
        del population_copy[population_copy.index(chosen1)]
        aspirants = self.selRandom(population_copy, self.tournament_size)
        chosen2 = max(aspirants, key=attrgetter('fitness'))
        chosen.append(chosen2)
        return chosen


class Individual(object):
    def __init__(self, phenotype, base_classifiers, X, Y, fit_attr, diversity_proportion):
        while all(v == 0 for v in phenotype): # avoid all-zero phenotypes
            phenotype = list(np.random.randint(2, size=len(base_classifiers)))
        self.phenotype = phenotype
        self.classifiers = []
        self.ensemble, self.accuracy, self.AUC, self.F1, self.MCC, self.entropy = self.calculate_measures(X, Y, base_classifiers, fit_attr)
        self.fitness = ((1.0 - diversity_proportion)*getattr(self, fit_attr)) + (diversity_proportion*self.entropy)

    def calculate_measures(self, X, Y, base_classifiers, fit_attr):
        best_ensemble = None
        vAcc = []
        vAUC = []
        vF1 = []
        vMCC = []
        vEntropy = []
        skf = StratifiedKFold(n_splits=5)
        for train_index, test_index in skf.split(X, Y):
            self.classifiers = self.generate_classifiers(base_classifiers)
            x_train, x_test, = X[train_index], X[test_index]
            y_train, y_test = Y[train_index], Y[test_index]
            ensemble, accuracy, AUC, F1, MCC, entropy = run_ensemble(x_train, x_test, y_train, y_test, classifiers_list=self.classifiers)
            vAcc.append(accuracy)
            vAUC.append(AUC)
            vF1.append(F1)
            vMCC.append(MCC)
            vEntropy.append(entropy)

            if not(best_ensemble) or eval(fit_attr) >= best_ensemble[1]:
                best_ensemble = (ensemble, eval(fit_attr))

        return best_ensemble[0], np.mean(vAcc), np.mean(vAUC), np.mean(vF1), np.mean(vMCC), np.mean(vEntropy)



    def generate_classifiers(self, base_classifiers):
        classifiers = []
        for i in range(len(self.phenotype)):
            if self.phenotype[i]:
                classifiers.append(deepcopy(base_classifiers[i]))
        return classifiers






