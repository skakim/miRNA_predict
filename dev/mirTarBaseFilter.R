setwd("/home/I864741/Codigos_GabrielMoita/dev")
t = read.table("hsa_MTI.csv",sep=',',header=TRUE)

head(t,n=10)
qtd = summary(t$Support.Type)
positive = qtd["Functional MTI"]
negative = qtd["Non-Functional MTI"]

positive
negative
negative/(positive+negative)