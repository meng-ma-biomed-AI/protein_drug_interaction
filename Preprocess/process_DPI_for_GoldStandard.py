import pandas as pd


def is_string_integer(s):
    try:
        num = int(s)
    except ValueError:
        return False
    return str(num) == s


data = open("anuna_protein_chemical_links_out.csv", 'r')
fo = open("anuna_protein_chemical_links_out_processed_GoldStandard.csv", 'w')
count = 0
for line in data:
    if count == 0:
        count += 1
        fo.write(line)
        continue
    else:
        line = line.strip().split(",")
        if is_string_integer(line[2]):
            if int(line[2]) >= 700:
                fo.write(",".join(line) + "\n")
    
data.close()
fo.close()







    
    





