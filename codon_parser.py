from __future__ import print_function
import re


"""This will parse through the Codon Table Program and will create dictionary """
CODON_TABlE_RE = re.compile("(\w+)\s+([\w?|\*?]\s+\d.\d+\s+\d+.\d)\s+\d+")   # (UGA) (* 0.21  0.7)  1345
count = 0
codonTableDict = {}
file_name = "PCC1720_codons.txt"
with open(file_name) as infile:
    for line in infile:
        linelist = line.strip().replace("(","").split(")")
        for element in linelist:
            codonline = element.strip()
            codon_match = CODON_TABlE_RE.search(codonline)
            if codon_match:
                codonValueList = CODON_TABlE_RE.search(codonline).group(2).split()
                codonKey = CODON_TABlE_RE.search(codonline).group(1).replace("U", "T")  # Convert RNA to DNA
                codonTableDict[codonKey] = codonValueList
                count += 1   # used for testing

print(count)  # used for testing
print(codonTableDict)
infile.close()

