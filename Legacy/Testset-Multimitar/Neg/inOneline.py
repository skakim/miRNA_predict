#!/usr/bin/python
import re
import sys

"""coloca todas as linhas de um fasta na mesma linha
"""
def inOneLine(inputName, outputName):
  Input = open(inputName)
  output = open(outputName, "w")
  seq1 = ""
  for line in Input.readlines():
    yesSeq = r'^A|a|G|g|C|c|T|t|R|r|H|h|K|k|D|d|E|e|S|s|N|n|Q|q|U|u|P|p|I|i|L|l|M|m|F|f|W|w|Y|y|V|v|' # agora aceita sequencia proteica
    noSeq = r'^\n'
    if line.startswith(">"):
      if re.match(yesSeq, seq1):
        output.write(seq1)
        seq1 = ""
        seqName = "\n" + line
        fileName = line.split(">")[1]
        fileName = fileName.split()[0]
        output.write(seqName)
      else:
        seq1 = ""
        seqName = "\n" + line
        fileName = line.split(">")[1]
        fileName = fileName.split()[0]
        output.write(seqName)
    elif re.match(noSeq, line):
      next
    elif re.match(yesSeq, line):
      seq = line.split()[0]
      seq1 = seq1 + seq
  output.write(seq1)
  print "It's done"

if __name__=="__main__":
  intro = "\nThis program put all fasta line only one line\n\n"
  try:
    inputA = sys.argv [1]; inputB = sys.argv[2]
  except:
    print intro, "\tUsage: \n\t\tpython",sys.argv[0], "[InputFile] [OutputFileName]\n\n\tWere:\n\t\tInputFile= Fasta File\n\t\tOutputFileName= Output file name\n"; sys.exit(1)

  inOneLine(inputA, inputB)
