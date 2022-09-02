import sys
import csv

"""
Take the scmpra object and convert to a promoter by cell matrix that can be used with MPRAnalyze
Input 1 - the csv from Siqi's script sc_crs_exp.py
Input 2 - prefix for output file
"""

scmpra = sys.argv[1]
prefix = sys.argv[2]

cellbcs = []
prombcs = []
cellbcs_promotercounts = {}
numplasmid_cellbcs_promotercounts = {}
directexp_cellbcs_promotercounts = {}
normexp_cellbcs_promotercounts = {}

print("Opening file ", scmpra, file = sys.stderr)
with open(scmpra) as scmpra_fh:
    reader = csv.DictReader(scmpra_fh)
    for line in reader:
        print(line["pBC"], line["cellBC"], line)
        cell = line['cellBC']
        promoter = line['pBC']
        if cell not in cellbcs:
            cellbcs.append(cell)
        if promoter not in prombcs:
            prombcs.append(promoter)
        if cell not in numplasmid_cellbcs_promotercounts:
            numplasmid_cellbcs_promotercounts[cell] = {}
            directexp_cellbcs_promotercounts[cell] = {}
            normexp_cellbcs_promotercounts[cell] = {}
        if promoter not in numplasmid_cellbcs_promotercounts[cell]:
            numplasmid_cellbcs_promotercounts[cell][promoter] = {}
            directexp_cellbcs_promotercounts[cell][promoter] = {}
            normexp_cellbcs_promotercounts[cell][promoter] = {}
        numplasmid_cellbcs_promotercounts[cell][promoter] = line["num_plasmid"]
        directexp_cellbcs_promotercounts[cell][promoter] = line["direct_exp"]
        normexp_cellbcs_promotercounts[cell][promoter] = line["norm_exp"]

print("Number of cells is ", len(cellbcs))
print("Number of promoters is ", len(prombcs))

def write_to_file(data, fh):
    header = ["pBC"]
    for cellbc in cellbcs:
        header.append(cellbc)
    print("\t".join(header), file = fh) #Header
    for promoter in prombcs:
        output = []
        output.append(promoter)
        for cell in cellbcs:
            #print(cell, promoter)
            if promoter in data[cell]:
                output.append(str(data[cell][promoter]))
            else:
                output.append("NA")
        print("\t".join(output), file = fh)

with open (prefix + "_numplasmid.tsv", "w") as numplasmid_fh:
    write_to_file(numplasmid_cellbcs_promotercounts, numplasmid_fh)

with open (prefix + "_directexp.tsv", "w") as directexp_fh:
    write_to_file(directexp_cellbcs_promotercounts, directexp_fh)
