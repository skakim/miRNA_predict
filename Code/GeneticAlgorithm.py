import numpy as np
import pandas as pd
from copy import deepcopy
from operator import attrgetter
from sklearn.utils import shuffle

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
    def __init__(self, base_classifiers, crossover_chance, mutation_chance, tournament_size, population_size, fit_attr):
        self.base_classifiers = base_classifiers
        self.crossover_chance = crossover_chance
        self.mutation_chance = mutation_chance
        self.tournament_size = tournament_size
        self.population_size = population_size
        self.fit_attr = fit_attr

        self.X, self.Y = read_mirna_dataset()

        self.population = self.generate_initial_population()
        self.best_individual = None
        self.generation = 0

    def fitnesses(self):
        l = []
        for ind in self.population:
            l.append(getattr(ind, self.fit_attr))
        return sorted(l, reverse=True)

    def generate_initial_population(self):
        population = []
        for i in range(self.population_size):
            random_initial = list(np.random.randint(2, size=len(self.base_classifiers)))
            population.append(Individual(random_initial, self.base_classifiers, self.X, self.Y))
        return population

    def generate_next_population(self):
        old_population = self.population
        new_population = []

        # test if best_individual needs to be updated
        best_from_population = max(old_population, key=attrgetter(self.fit_attr))
        if not(self.best_individual) or getattr(best_from_population, self.fit_attr) >= getattr(self.best_individual, self.fit_attr):
            self.best_individual = best_from_population

        # generate offspring
        # elitism
        new_population.append(self.best_individual)

        # let the games begin!
        for i in range(self.population_size//2):
            # crossover
            ind1, ind2 = self.selTournament()
            phen1, phen2 = self.crossover(ind1, ind2)

            #mutation
            phen1 = self.mutation(phen1)
            phen2 = self.mutation(phen2)

            new1 = Individual(phen1, self.base_classifiers, self.X, self.Y)
            new2 = Individual(phen2, self.base_classifiers, self.X, self.Y)

            new_population.append(new1)
            new_population.append(new2)

        self.population = new_population
        self.generation += 1

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

    def selRandom(self, k):
        """Select *k* individuals at random from the input *individuals* with
        replacement. The list returned contains references to the input
        *individuals*.

        :param individuals: A list of individuals to select from.
        :param k: The number of individuals to select.
        :returns: A list of selected individuals.

        This function uses the :func:`~random.choice` function from the
        python base :mod:`random` module.
        """
        return np.random.choice(self.population, size=k)

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
        chosen = []
        for i in range(k):
            aspirants = self.selRandom(self.tournament_size)
            chosen.append(max(aspirants, key=attrgetter(self.fit_attr)))
        return chosen


class Individual(object):
    def __init__(self, phenotype, base_classifiers, X, Y):
        while all(v == 0 for v in phenotype): # avoid all-zero phenotypes
            phenotype = list(np.random.randint(2, size=len(base_classifiers)))
        self.phenotype = phenotype
        self.classifiers = self.generate_classifiers(base_classifiers)
        self.ensemble, self.accuracy, self.AUC, self.F1, self.MCC = run_ensemble(X, Y, classifiers_list=self.classifiers)

    def generate_classifiers(self, base_classifiers):
        classifiers = []
        for i in range(len(self.phenotype)):
            if self.phenotype[i]:
                classifiers.append(deepcopy(base_classifiers[i]))
        return classifiers






