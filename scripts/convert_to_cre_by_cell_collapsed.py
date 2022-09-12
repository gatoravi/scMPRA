import sys
import csv

"""
Take the scmpra object and convert to a promoter by cell matrix that can be used with MPRAnalyze
Input 1 - the csv from Siqi's script sc_crs_exp.py
Input 2 - prefix for output file
Input 3 - "replicate number" - this will be appended to all the cell names
Input 4 - The annotation file for the promoters "prom_lib_master_info.tsv"
"""

scmpra = sys.argv[1]
prefix = sys.argv[2]
rep = sys.argv[3]
annotations = sys.argv[4]

cellbcs = []
prombcs = []
promoters = []
cellbcs_promotercounts = {}
numplasmid_cellbcs_promotercounts = {}
directexp_cellbcs_promotercounts = {}
normexp_cellbcs_promotercounts = {}

pBC_promotername = {} # Dict to go from pBC to promoter name

with open(annotations) as annotations_fh:
    for line in annotations_fh:
        fields = line.split()
        promoter_name = fields[1]
        pBC = fields[-1]
        pBC_promotername[pBC] = promoter_name


print("Opening file ", scmpra, file = sys.stderr)
with open(scmpra) as scmpra_fh:
    reader = csv.DictReader(scmpra_fh)
    for line in reader:
        #print(line["pBC"], line["cellBC"], line)
        cell = line['cellBC'] + "_" + rep
        pBC = line['pBC']
        promoter = pBC
        promotername = pBC_promotername[pBC]
        cell = cell + "_" + promotername
        if cell not in cellbcs:
            cellbcs.append(cell)
        if pBC not in prombcs:
            prombcs.append(pBC)
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
print("Number of promoter barcodes is ", len(prombcs))

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

with open (prefix + "_colannot.tsv", "w") as colannot_fh:
    print("cellBC", "cell", "rep", file = colannot_fh)
    for cell in cellbcs:
        print("\t".join([cell, cell, rep]), file = colannot_fh)