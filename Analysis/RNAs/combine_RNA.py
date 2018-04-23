import pandas as pd

names = pd.read_csv("mart_names.txt",sep=',')
sequences = pd.read_csv("mart_sequences.txt",sep=',')

print(len(names),len(sequences))