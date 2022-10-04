import sys
import csv

"""
Take the scmpra object and convert to a promoter by cell matrix that can be used with MPRAnalyze
Input 1 - the csv from Siqi's script sc_crs_exp.py
Input 2 - prefix for output file
Input 3 - "replicate number" - this will be appended to all the cell names
Input 4 - The annotation file for the promoters "prom_lib_master_info.tsv"
Input 5 - the clusters file - cellbc, clusters
Input 6 - cells passing filters, usually those with the transcriptome
"""

scmpra = sys.argv[1]
prefix = sys.argv[2]
rep = sys.argv[3]
annotations = sys.argv[4]
clusters = sys.argv[5]
prefiltered_cells_file= sys.argv[6]



pref_cells = []
cellbcs = []
prombcs = []
promoters = []
cellbcs_promotercounts = {}
numplasmid_cellbcs_promotercounts = {}
directexp_cellbcs_promotercounts = {}
normexp_cellbcs_promotercounts = {}
cells_clusters = {}
pBC_promotername = {} # Dict to go from pBC to promoter name

with open(prefiltered_cells_file) as pfh:
    for line in pfh:
        pref_cells.append(line.rstrip("\n"))

with open(annotations) as annotations_fh:
    for line in annotations_fh:
        fields = line.split()
        promoter_name = fields[1]
        pBC = fields[-1]
        pBC_promotername[pBC] = promoter_name

with open(clusters) as clusters_fh:
    for line in clusters_fh:
        fields = line.split()
        cellbc = fields[0]
        cluster = fields[1]
        cells_clusters[cellbc] = cluster

print("Opening file ", scmpra, file = sys.stderr)
with open(scmpra) as scmpra_fh:
    reader = csv.DictReader(scmpra_fh)
    for line in reader:
        #print(line["pBC"], line["cellBC"], line)
        cell = line['cellBC'] + "_" + rep #Use the replicate to avoid collisions
        pBC = line['pBC']
        promotername = pBC_promotername[pBC]
        promoter = promotername #use the name here
        if cell not in cellbcs and (line['cellBC'] in pref_cells or not pref_cells):
            cellbcs.append(cell)
        if pBC not in prombcs:
            prombcs.append(pBC)
        if promotername not in promoters:
            promoters.append(promotername)
        if cell not in numplasmid_cellbcs_promotercounts:
            numplasmid_cellbcs_promotercounts[cell] = {}
            directexp_cellbcs_promotercounts[cell] = {}
            normexp_cellbcs_promotercounts[cell] = {}
        if promoter not in numplasmid_cellbcs_promotercounts[cell]:
            numplasmid_cellbcs_promotercounts[cell][promoter] = 0
            directexp_cellbcs_promotercounts[cell][promoter] = 0
            normexp_cellbcs_promotercounts[cell][promoter] = 0
        numplasmid_cellbcs_promotercounts[cell][promoter] += float(line["num_plasmid"]) #sum over all pBCs
        directexp_cellbcs_promotercounts[cell][promoter] += float(line["direct_exp"])#sum over all pBCs
        normexp_cellbcs_promotercounts[cell][promoter] += float(line["norm_exp"])#sum over all pBCs

print("Number of cells is ", len(cellbcs))
print("Number of promoter barcodes is ", len(prombcs))
print("Number of promoters is ", len(promoters))

def write_to_file(data, fh):
    header = ["promoter"]
    for cellbc in cellbcs:
        header.append(cellbc)
    print("\t".join(header), file = fh) #Header
    for promoter in promoters:
        output = []
        output.append(promoter)
        for cell in cellbcs:
            cellbc = cell.split("_")[0]
            if cellbc in pref_cells or not pref_cells: #if list of filtered cells is provided, only print those cells
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
    print("cellBC", "cell", "rep", "cluster", file = colannot_fh)
    for cell in cellbcs:
        cellbc = cell.split("_")[0]
        cluster = "unknown"
        if cellbc in cells_clusters:
            cluster = cells_clusters[cellbc]
        if cellbc in pref_cells or not pref_cells:
            print("\t".join([cell, cellbc, rep, cluster]), file = colannot_fh)