##################################################################################################################################### 
#this script is to process the PPI data
# input: PPI data:   13 resources available 
#  BioGrid: 
#  BioPlex: 
#  Human Reference Interactome: 
#  InnateDB: 
#  Insider: 
#  Intact: 
#  Interactome3D: 
#  InWeb_IM: not available now
#  KinomeNetworkX: not available now 
#  mint: will not use it.
#  PhosphoSitePlus: not available now
#  SingnaLink:  
#  STRING:   
#  integrated_Interactions_database: 
#  UniProtKB:
# Meng, 01/09/2024 
#######################################################################################################################################

import pandas as pd
import numpy as np  
import re 

datadict = r"H:\work\DrugDiscovery\ANUNA\work\Interacome_drug_discovery\Data"

#####################################################################################################################################
# Step 1: process the HGNC to get the mapping among entrez, uniprot, ensembl. sysmbol
#####################################################################################################################################
hgnc = open(datadict + r"\HGNC\hgnc_complete_set.txt", "r", encoding='utf-8')
hgnc_mapping = open(datadict + "\HGNC\hgnc_mapping.txt", "w", encoding='utf-8')
hgnc_mapping.write("hgncid\tsymbol\tentrezid\tensemblgid\trefseqaccesssion\tccdsid\tuniprotid\n")

cur = 0
for line in hgnc:
    if cur == 0:
        cur = cur + 1 
        continue
    data = line.split("\t")
    #print(data.__len__())
    hgncid = data[0]
    symbol = data[1]
    entrezid = data[18]
    ensemblgid = data[19]
    refseqaccesssion = data[23]
    ccdsid = data[24]
    uniprotid = data[25]
    hgnc_mapping.write(hgncid + "\t" + symbol + "\t" + entrezid + "\t" + ensemblgid + "\t" + refseqaccesssion + "\t" + ccdsid + "\t" + uniprotid + "\n")

hgnc.close()
hgnc_mapping.close()

#####################################################################################################################################
# Step 2: process each PPI data sources
#####################################################################################################################################

# process 1st: Human Reference Interactome 
hri_dict = {}

hi05 = open(datadict + r"\Human_Reference_Interactome\H-I-05.tsv", 'r', encoding='utf-8')
for line in hi05:
    data = line.strip()
    hri_dict[data] = 1
hi05.close()

hiii14 = open(datadict + r"\Human_Reference_Interactome\HI-II-14.tsv", 'r', encoding='utf-8')
for line in hiii14:
    data = line.strip()
    hri_dict[data] = 1
hiii14.close()

hiunion = open(datadict + r"\Human_Reference_Interactome\HI-union.tsv", 'r', encoding='utf-8')
for line in hiunion:
    data = line.strip()
    hri_dict[data] = 1
hiunion.close()

huri = open(datadict + r"\Human_Reference_Interactome\HuRI.tsv", 'r', encoding='utf-8')
for line in huri:
    data = line.strip()
    hri_dict[data] = 1
huri.close()

litbm = open(datadict + r"\Human_Reference_Interactome\Lit-BM.tsv", 'r', encoding='utf-8')
for line in litbm:
    data = line.strip()
    hri_dict[data] = 1
litbm.close()

test = open(datadict + r"\Human_Reference_Interactome\Test_space_screens-19.tsv", 'r', encoding='utf-8')
for line in test:
    data = line.strip()
    hri_dict[data] = 1
test.close()

Venkatesan = open(datadict + r"\Human_Reference_Interactome\Venkatesan-09.tsv", 'r', encoding='utf-8')
for line in Venkatesan:
    data = line.strip()
    hri_dict[data] = 1
Venkatesan.close()

yang = open(datadict + r"\Human_Reference_Interactome\Yang-16.tsv", 'r', encoding='utf-8')
for line in yang:
    data = line.strip()
    hri_dict[data] = 1
yang.close()

yu = open(datadict + r"\Human_Reference_Interactome\Yu-11.tsv", 'r', encoding='utf-8')
for line in yu:
    data = line.strip()
    hri_dict[data] = 1
yu.close()

hri = open(datadict + r"\Human_Reference_Interactome\HUMAN_REFERENCE_INTERACTIONS.tsv", 'w', encoding='utf-8')

for key in hri_dict.keys():
    hri.write(key + "\n")
#print(hri_dict.keys().__len__()) 
hri.close()


# process 2nd: BioGrid
biogrid = open(datadict + r"\BioGrid\BIOGRID-ALL-4.4.229.mitab\BIOGRID-ALL-4.4.229.mitab.txt", 'r', encoding='utf-8')
bioppis = open(datadict + r"\BioGrid\BIOGRID-ALL-4.4.229.mitab\BioGrid_PPIs.txt", 'w', encoding='utf-8')
# load the HGNC mapping data for convert entrez id to ensembl gene id 
hgnc_mapping = open(datadict + "\HGNC\hgnc_mapping.txt", "r", encoding='utf-8')
hgnc_mapping_dict = {}
cur = 0
for line in hgnc_mapping:
    if cur == 0:
        cur = cur + 1 
        continue
    data = line.split("\t")
    hgncid = data[0]
    symbol = data[1]
    entrezid = data[2]
    ensemblgid = data[3]
    refseqaccesssion = data[4]
    ccdsid = data[5]
    # remove line break 
    uniprotid = data[6].replace("\n", "")
    if entrezid != "":
        hgnc_mapping_dict[entrezid] = ensemblgid
    if uniprotid != "":
        hgnc_mapping_dict[uniprotid] = ensemblgid
    #print(entrezid + "\t" + ensemblgid)

count = 0
for line in biogrid:
    if count == 0:
        count = count + 1
        continue
    data = line.split("\t")
    #print(data.__len__())
    entrezid1 = data[0].split(":")[1]
    entrezid2 = data[1].split(":")[1]
    #print(entrezid1 + "\t" + entrezid2)
    if entrezid1 in hgnc_mapping_dict.keys():
        ensemblgid1 = hgnc_mapping_dict[entrezid1]
    else:
        ensemblgid1 = ""    
    if entrezid2 in hgnc_mapping_dict.keys():
        ensemblgid2 = hgnc_mapping_dict[entrezid2]
    else:
        ensemblgid2 = ""
    if ensemblgid1 != "" and ensemblgid2 != "":
        bioppis.write(ensemblgid1 + "\t" + ensemblgid2 + "\n")
biogrid.close()
hgnc_mapping.close()


# process 3rd: BioPlex
bioplex = {}

data1 = open(datadict + r"\BioPlex\BioPlex_3.0_HCT116_DirectedEdges.tsv", 'r', encoding='utf-8')
count = 0
for line in data1:
    if count == 0:
        count = count + 1
        continue
    data = line.split("\t")
    #print(data.__len__())
    entrezid1 = data[0].replace("\"", "")
    entrezid2 = data[2].replace("\"", "")
    if entrezid1 in hgnc_mapping_dict.keys():
        ensemblgid1 = hgnc_mapping_dict[entrezid1]
    else:
        ensemblgid1 = ""
    if entrezid2 in hgnc_mapping_dict.keys():
        ensemblgid2 = hgnc_mapping_dict[entrezid2]
    else:
        ensemblgid2 = ""
    if ensemblgid1 != "" and ensemblgid2 != "":
        bioplex[ensemblgid1 + "\t" + ensemblgid2] = 1
data1.close()

data2 = open(datadict + r"\BioPlex\BioPlex_3.0_293T_DirectedEdges.tsv", 'r', encoding='utf-8')
count = 0
for line in data2:
    if count == 0:
        count = count + 1
        continue
    data = line.split("\t")
    #print(data.__len__())
    entrezid1 = data[0].replace("\"", "")
    entrezid2 = data[2].replace("\"", "")
    if entrezid1 in hgnc_mapping_dict.keys():
        ensemblgid1 = hgnc_mapping_dict[entrezid1]
    else:
        ensemblgid1 = ""
    if entrezid2 in hgnc_mapping_dict.keys():
        ensemblgid2 = hgnc_mapping_dict[entrezid2]
    else:
        ensemblgid2 = ""
    if ensemblgid1 != "" and ensemblgid2 != "":
        bioplex[ensemblgid1 + "\t" + ensemblgid2] = 1
data2.close()

bioplexppis = open(datadict + r"\BioPlex\BioPlex_PPIs.txt", 'w', encoding='utf-8')
for key in bioplex.keys():
    bioplexppis.write(key + "\n")
bioplexppis.close()


# process 4th: innateDB
data = open(datadict + r"\InnateDB\innatedb_ppi.mitab\innatedb_ppi.mitab", 'r', encoding='utf-8')
innateppis = open(datadict + r"\InnateDB\InnateDB_PPIs.txt", 'w', encoding='utf-8')    
count = 0
for line in data:
    if count == 0:
        count = count + 1
        continue
    ensemblgid1 = line.split("\t")[2].split(":")[1]
    ensemblgid2 = line.split("\t")[3].split(":")[1]
    if not ("ENSMUSG" in ensemblgid1 or "ENSMUSG" in ensemblgid2): # exclude the mouse data
        innateppis.write(ensemblgid1 + "\t" + ensemblgid2 + "\n")
data.close()
innateppis.close()

# process 5th: insider
data = open(datadict + r"\Insider\H_sapiens_interfacesHQ.txt", 'r', encoding='utf-8')
insiderppis = open(datadict + r"\Insider\Insider_PPIs.txt", 'w', encoding='utf-8')
count = 0
for line in data:
    if count == 0:
        count = count + 1
        continue
    uniprotid1 = line.split("\t")[0]
    uniprotid2 = line.split("\t")[1]
    if uniprotid1 in hgnc_mapping_dict.keys():
        ensemblgid1 = hgnc_mapping_dict[uniprotid1]
    else:
        ensemblgid1 = ""
    if uniprotid2 in hgnc_mapping_dict.keys():
        ensemblgid2 = hgnc_mapping_dict[uniprotid2]
    else:
        ensemblgid2 = ""
    if ensemblgid1 != "" and ensemblgid2 != "":
        insiderppis.write(ensemblgid1 + "\t" + ensemblgid2 + "\n")
data.close()
insiderppis.close()

# process 6th: intact
data = open(datadict + r"\Intact\all\psimitab\intact.txt", 'r', encoding='utf-8')
intactppis = open(datadict + "\Intact\Intact_PPIs.txt", 'w', encoding='utf-8')
count = 0
for line in data:
    if count == 0:
        count = count + 1
        continue
    #print(line.split("\t")[0])
    uniprotid1 = ""
    uniprotid2 = ""
    if "uniprotkb" in line.split("\t")[0]:
        uniprotid1 = line.split("\t")[0].split(":")[1]
    if "uniprotkb" in line.split("\t")[1]:
        uniprotid2 = line.split("\t")[1].split(":")[1]
    if uniprotid1 in hgnc_mapping_dict.keys():
        ensemblgid1 = hgnc_mapping_dict[uniprotid1]
    else:
        ensemblgid1 = ""
    if uniprotid2 in hgnc_mapping_dict.keys():
        ensemblgid2 = hgnc_mapping_dict[uniprotid2]
    else:
        ensemblgid2 = ""
    if ensemblgid1 != "" and ensemblgid2 != "":    
        intactppis.write(ensemblgid1 + "\t" + ensemblgid2 + "\n")
data.close()
intactppis.close()

# process 7th: integrated_interactions_databases

data = open(datadict + r"\integrated_interactions_dtabases\human_annotated_PPIs\human_annotated_PPIs.txt", 'r', encoding='utf-8')
iidppis = open(datadict + r"\integrated_interactions_dtabases\IID_PPIs.txt", 'w', encoding='utf-8')
count = 0
for line in data:
    if count == 0:
        count = count + 1
        continue
    #print(line.split("\t")[0])
    uniprotid1 = ""
    uniprotid2 = ""
    uniprotid1 = line.split("\t")[0]
    uniprotid2 = line.split("\t")[1]
    if uniprotid1 in hgnc_mapping_dict.keys():
        ensemblgid1 = hgnc_mapping_dict[uniprotid1]
    else:
        ensemblgid1 = ""
    if uniprotid2 in hgnc_mapping_dict.keys():
        ensemblgid2 = hgnc_mapping_dict[uniprotid2]
    else:
        ensemblgid2 = ""
    if ensemblgid1 != "" and ensemblgid2 != "":
        iidppis.write(ensemblgid1 + "\t" + ensemblgid2 + "\n")
data.close()
iidppis.close()

# process 8th: Interactome3D
data = open(datadict + r"\Interactome3D\interactions.dat", 'r', encoding='utf-8')
interactome3dppis = open(datadict + r"\Interactome3D\Interactome3D_PPIs.txt", 'w', encoding='utf-8')
count = 0
for line in data:
    if count == 0:
        count = count + 1
        continue
    uniprotid1 = ""
    uniprotid2 = ""
    uniprotid1 = line.split("\t")[0]
    uniprotid2 = line.split("\t")[1]
    if uniprotid1 in hgnc_mapping_dict.keys():
        ensemblgid1 = hgnc_mapping_dict[uniprotid1]
    else:
        ensemblgid1 = ""
    if uniprotid2 in hgnc_mapping_dict.keys():
        ensemblgid2 = hgnc_mapping_dict[uniprotid2]
    else:
        ensemblgid2 = ""
    if ensemblgid1 != "" and ensemblgid2 != "":
        interactome3dppis.write(ensemblgid1 + "\t" + ensemblgid2 + "\n")
data.close()
interactome3dppis.close()

# process 9th: SingaLink
data = open(datadict + r"\SignaLink\slk3_29accdffbdb2ef9a9224_csv\slk3_29accdffbdb2ef9a9224_csv.csv", 'r', encoding='utf-8')
slppis = open(datadict + r"\SignaLink\SignaLink_PPIs.txt", 'w', encoding='utf-8')
count = 0
for line in data:
    if count == 0:
        count = count + 1
        continue
    uniprotid1 = ""
    uniprotid2 = ""
    uniprotid1 = line.split(",")[1]
    uniprotid2 = line.split(",")[7]
    if uniprotid1 in hgnc_mapping_dict.keys():
        ensemblgid1 = hgnc_mapping_dict[uniprotid1]
    else:
        ensemblgid1 = ""
    if uniprotid2 in hgnc_mapping_dict.keys():
        ensemblgid2 = hgnc_mapping_dict[uniprotid2]
    else:
        ensemblgid2 = ""
    if ensemblgid1 != "" and ensemblgid2 != "":
        slppis.write(ensemblgid1 + "\t" + ensemblgid2 + "\n")
data.close()
slppis.close()


# process 10th: STRING
data = open(datadict + r"\STRING\9606.protein.links.v12.0.txt\9606.protein.links.v12.0.txt", 'r', encoding='utf-8')
stringppis = open(datadict + r"\STRING\STRING_PPIs_score_cut_150.txt", 'w', encoding='utf-8')

ensp_ensg = {}
string_mapping = open(datadict + r"\STRING\9606.protein.aliases.v12.0.txt\9606.protein.aliases.v12.0.txt", 'r', encoding='utf-8')
count = 0
for line in string_mapping:
    if count == 0:
        count = count + 1
        continue

    enspid = line.split("\t")[0].split(".")[1]
    ensgid = ""
    if "ENSG" in line.split("\t")[1]:
        ensgid = line.split("\t")[1]
    if ensgid != "":
        ensp_ensg[enspid] = ensgid

count = 0
for line in data:
    if count == 0:
        count = count + 1
        continue
    ensemblpid1 = line.split(' ')[0].split('.')[1]
    ensemblpid2 = line.split(' ')[1].split('.')[1]
    score = line.split(' ')[2].strip()

    ensemblgid1 = ensp_ensg[ensemblpid1]
    ensemblgid2 = ensp_ensg[ensemblpid2]
    if ensemblgid1 != "" and ensemblgid2 != "" and int(score) >= 150:
        stringppis.write(ensemblgid1 + "\t" + ensemblgid2 + "\n")
data.close()
stringppis.close()
string_mapping.close()


# process 11th: UniProtKB
data = open(datadict + r"\UniProtKB\uniprotkb_AND_model_organism_9606_2024_01_03.tsv\uniprotkb_AND_model_organism_9606_2024_01_03.tsv", 'r', encoding='utf-8')
uniprotppis = open(datadict + r"\UniProtKB\UniProtKB_PPIs.txt", 'w', encoding='utf-8')
count = 0
for line in data:
    if count == 0:
        count = count + 1
        continue
    uniprotid1 = line.split("\t")[0]
    uniprotid2s = line.split("\t")[2]
    if uniprotid1 in hgnc_mapping_dict.keys():
        ensemblgid1 = hgnc_mapping_dict[uniprotid1]
    else:
        ensemblgid1 = ""

    ensemblgid2 = ""
    if uniprotid2s != "":
        uniprotid2s = uniprotid2s.split(";")
        for uniprotid2 in uniprotid2s:
            uniprotid2 = uniprotid2.strip()
            if uniprotid2 in hgnc_mapping_dict.keys():
                ensemblgid2 = hgnc_mapping_dict[uniprotid2]
            else:
                ensemblgid2 = ""
            if ensemblgid1 != "" and ensemblgid2 != "":
                uniprotppis.write(ensemblgid1 + "\t" + ensemblgid2 + "\n")
data.close()
uniprotppis.close()


#####################################################################################################################################
# Step 3: integrate all PPI data results to generate the final PPI interactome data
#####################################################################################################################################

PPIs = {}

statistics = open(datadict + r"\PPIs_statistics.txt", 'w', encoding='utf-8')

data = open(datadict + r"\Human_Reference_Interactome\HUMAN_REFERENCE_INTERACTIONS.tsv", 'r', encoding='utf-8')
count = 0
for line in data:
    PPIs[line.strip()] = 1
    count = count + 1
data.close()
statistics.write("Human_Reference_Interactome\t" + str(count) + "\n")

data = open(datadict + r"\BioGrid\BIOGRID-ALL-4.4.229.mitab\BioGrid_PPIs.txt", 'r', encoding='utf-8')
count = 0
for line in data:
    PPIs[line.strip()] = 1
    count = count + 1
data.close()
statistics.write("BioGrid\t" + str(count) + "\n")

data = open(datadict + r"\BioPlex\BioPlex_PPIs.txt", 'r', encoding='utf-8')
count = 0
for line in data:
    PPIs[line.strip()] = 1
    count = count + 1
data.close()
statistics.write("BioPlex\t" + str(count) + "\n")

data = open(datadict + r"\InnateDB\InnateDB_PPIs.txt", 'r', encoding='utf-8')
count = 0
for line in data:
    PPIs[line.strip()] = 1
    count = count + 1    
data.close()
statistics.write("InnateDB\t" + str(count) + "\n")

data = open(datadict + r"\Insider\Insider_PPIs.txt", 'r', encoding='utf-8')
count = 0
for line in data:
    PPIs[line.strip()] = 1
    count = count + 1
data.close()
statistics.write("Insider\t" + str(count) + "\n")

data = open(datadict + r"\Intact\Intact_PPIs.txt", 'r', encoding='utf-8')
count = 0
for line in data:
    PPIs[line.strip()] = 1
    count = count + 1   
data.close()
statistics.write("Intact\t" + str(count) + "\n")

data = open(datadict + r"\integrated_interactions_dtabases\IID_PPIs.txt", 'r', encoding='utf-8')
count = 0
for line in data:
    PPIs[line.strip()] = 1
    count = count + 1
data.close()
statistics.write("IID\t" + str(count) + "\n")

data = open(datadict + r"\Interactome3D\Interactome3D_PPIs.txt", 'r', encoding='utf-8')
count = 0
for line in data:
    PPIs[line.strip()] = 1
    count = count + 1
data.close()
statistics.write("Interactome3D\t" + str(count) + "\n")

data = open(datadict + r"\SignaLink\SignaLink_PPIs.txt", 'r', encoding='utf-8')
count = 0
for line in data:
    PPIs[line.strip()] = 1
    count = count + 1
data.close()
statistics.write("SignaLink\t" + str(count) + "\n")

data = open(datadict + r"\STRING\STRING_PPIs_score_cut_150.txt", 'r', encoding='utf-8')
count = 0
for line in data:
    PPIs[line.strip()] = 1
    count = count + 1
data.close()
statistics.write("STRING\t" + str(count) + "\n")

data = open(datadict + r"\UniProtKB\UniProtKB_PPIs.txt", 'r', encoding='utf-8')
count = 0
for line in data:
    PPIs[line.strip()] = 1
    count = count + 1   
data.close()
statistics.write("UniProtKB\t" + str(count) + "\n")

ppis = open(datadict + r"\PPIs.txt", 'w', encoding='utf-8')
for key in PPIs.keys():
    ppis.write(key + "\n")
ppis.close()
statistics.close()








