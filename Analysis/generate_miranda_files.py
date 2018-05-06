import pandas

RNAtable = pandas.read_table("D:/TCC/miRNA_predict/Analysis/RNAs/names-sequences.txt", sep=',', doublequote=True)
miRNAtable = pandas.read_table("D:/TCC/miRNA_predict/Analysis/miRNAs/filtered_mature.fa", sep=' ', header=None)

tarbase_positive = pandas.read_table("D:/TCC/miRNA_predict/Analysis/tarbase-positive.csv", sep=' ', header=None, doublequote=True)
tarbase_negative = pandas.read_table("D:/TCC/miRNA_predict/Analysis/tarbase-negative.csv", sep=' ', header=None, doublequote=True)

#mirna_xxx.fas -> target_xxx.fas

#positive creating
i = 0
k = 0
for index, row in tarbase_positive.iterrows():
    RNA = row[tarbase_positive.columns[3]]
    miRNA = row[tarbase_positive.columns[1]]

    RNA_rows = RNAtable.loc[RNAtable["Gene.stable.ID"] == RNA]
    miRNA_row = miRNAtable.loc[miRNAtable[miRNAtable.columns[0]] == miRNA]

    if len(miRNA_row) == 1 and len(RNA_rows) > 0:
        i += 1
        k += len(RNA_rows)
        #print(">" + miRNA + "\n" + miRNA_row.iloc[0][miRNAtable.columns[5]] + "\n")
        #print("-------------------------------------------------------------")
        #j = 0
        #for index, row in RNA_rows.iterrows():
        #    j += 1
        #    print(">" + RNA + "-v" + str(j) + "\n" + row[miRNAtable.columns[4]] + "\n")

        with open("miranda_files/positive/mirna_"+str(i)+".fas","w") as mirna_file:
            mirna_file.write(">"+miRNA+"\n"+miRNA_row.iloc[0][miRNAtable.columns[5]]+"\n")

        with open("miranda_files/positive/target_" + str(i) + ".fas", "w") as rna_file:
            j = 0
            for index, row in RNA_rows.iterrows():
                j += 1
                rna_file.write(">"+RNA+"-v"+str(j)+"\n"+row[miRNAtable.columns[4]]+"\n")

print("Positive:",i,k)

#negative creating
i = 0
k = 0
for index, row in tarbase_negative.iterrows():
    RNA = row[tarbase_positive.columns[3]]
    miRNA = row[tarbase_positive.columns[1]]

    RNA_rows = RNAtable.loc[RNAtable["Gene.stable.ID"] == RNA]
    miRNA_row = miRNAtable.loc[miRNAtable[miRNAtable.columns[0]] == miRNA]

    if len(miRNA_row) == 1 and len(RNA_rows) > 0:
        i += 1
        k += len(RNA_rows)
        #print(">" + miRNA + "\n" + miRNA_row.iloc[0][miRNAtable.columns[5]] + "\n")
        #print("-------------------------------------------------------------")
        #j = 0
        #for index, row in RNA_rows.iterrows():
        #    j += 1
        #    print(">" + RNA + "-v" + str(j) + "\n" + row[miRNAtable.columns[4]] + "\n")

        with open("miranda_files/negative/mirna_"+str(i)+".fas","w") as mirna_file:
            mirna_file.write(">"+miRNA+"\n"+miRNA_row.iloc[0][miRNAtable.columns[5]]+"\n")

        with open("miranda_files/negative/target_" + str(i) + ".fas", "w") as rna_file:
            j = 0
            for index, row in RNA_rows.iterrows():
                j += 1
                rna_file.write(">"+RNA+"-v"+str(j)+"\n"+row[miRNAtable.columns[4]]+"\n")

print("Negative:",i,k)