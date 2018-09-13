import pickle
from sklearn.utils import shuffle
import pandas as pd
from sklearn.metrics import confusion_matrix, f1_score, matthews_corrcoef, roc_auc_score

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

with open("results/20180912-193005/0.646368614943-1010001111111001.ens", "rb") as f:
    ensemble = pickle.load(f)

X, Y = read_mirna_dataset()

# For testing x_test & y_test
y_pred = ensemble.predict(X)

cm = confusion_matrix(Y, y_pred)
acc = float(cm.trace()) / cm.sum()
AUC = roc_auc_score(Y, y_pred)
F1 = f1_score(Y, y_pred)
MCC = matthews_corrcoef(Y, y_pred)

print("cl_name, acc, AUC, F1, MCC")
print(cm)
print("Ensemble", format(acc * 100, '.2f'), format(AUC * 100, '.2f'), format(F1 * 100, '.2f'),
      format(MCC, '.2f'), sep=',')

for classifier in ensemble._classifiers:
    y_pred = classifier.predict(X)
    cm = confusion_matrix(Y, y_pred)
    print(cm)
    acc = float(cm.trace()) / cm.sum()
    AUC = roc_auc_score(Y, y_pred)
    F1 = f1_score(Y, y_pred)
    MCC = matthews_corrcoef(Y, y_pred)
    print(type(classifier).__name__, format(acc * 100, '.2f'), format(AUC * 100, '.2f'), format(MCC, '.2f'),
          sep=',')