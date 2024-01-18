import pandas as pd

def is_string_integer(s):
    try:
        num = int(s)
    except ValueError:
        return False
    return str(num) == s

data = open("chemical_chemical.links.v5.0.tsv", 'r')
fo = open("chemical_chemical.links.v5.0_processed.tab", 'w')

for line in data:
    line = line.strip().split("\t")
    if is_string_integer(line[2]):
        if int(line[2]) >= 150:
            fo.write("\t".join(line) + "\n")




