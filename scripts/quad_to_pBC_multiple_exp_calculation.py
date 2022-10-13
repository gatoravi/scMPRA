'''
This script only deals with input quint containing preassigned clusters based on 
the transcriptome data.
'''
import os,sys
from numpy.core.defchararray import count
import pandas as pd
import numpy as np
import gzip
from scipy import stats
import pickle
import argparse
import time

def return_abundant_cells(df, min_umi = 100):
    umi_dict = {}
    for _, row in df.iterrows():
        cell_bc = row['cellBC']
        if cell_bc not in umi_dict:
            umi_dict[cell_bc] = 1
        else:
            umi_dict[cell_bc] +=1
    abundant_list = []
    for key in umi_dict:
        if umi_dict[key] > min_umi:
            abundant_list.append(key)
    # Slice the dataset
    pop_df = df[df['cellBC'].isin(abundant_list)]
    return pop_df

def hamming_distance(chaine1, chaine2):
    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))

def attach_promoter_to_quad(quad_df, info_df):
    pop_list = []
    pBC_map = {}
    # Only take the variant name and the barcode
    info_df = info_df[['name', 'barcode']]
    # Extract the unique pBCs 
    pBC_list = quad_df.pBC.unique()
    pBC_to_check = info_df.barcode.unique()
    # Now attempt to error correct
    for pBC in pBC_list:
        possible_fit = []
        if pBC in pBC_to_check:
            pBC_map[pBC] = pBC
        else:
            for bg_pBC in pBC_to_check:
                temp_hamming = hamming_distance(pBC, bg_pBC)
                if temp_hamming <=2:
                    possible_fit.append(bg_pBC)
            if len(possible_fit) == 1:
                pBC_map[pBC] = possible_fit[0]
    print('finished mapping for pBC error correction')
    # Iterate through rows 
    for _, line in quad_df.iterrows():
        pBC = line['pBC']
        # fuzzy match to the info_df
        if pBC in pBC_map:
            prom_id = info_df[info_df['barcode'] == pBC_map[pBC]]['name'].values[0]
            pop_list.append([line['cellBC'], line['umi'], prom_id, pBC_map[pBC], line['rBC'], line['count'], line['cluster']])
    sexa = pd.DataFrame(pop_list,columns= ['cellBC', 'umi', 'name', 'pBC', 'rBC', 'count', 'cluster'])
    return sexa

def erorr_correct_bcs_v3(quint):
    # Here we try to use the minimum number of rBC per cell 
    start = time.time()
    cells = quint.cellBC.unique()
    pop_list = []
    trimed_pop_df = pd.DataFrame( columns = ['cellBC', 'umi', 'name', 'pBC', 'cluster', 'rBC', 'count'])
    for cell in cells:
        # Slice on cell
        data_slice = quint.loc[quint.cellBC == cell]
        # initiate a dict for holding info
        data_slice.sort_values(by = ['count'], ascending = False, inplace = True)
        data_slice = data_slice.reset_index(drop = True)

        correct_dict = {}
        pBCs = data_slice.pBC.unique()
        for pBC in pBCs:
            
            umi_to_keep = []
            rBC_to_keep = []
            pBC_slice = data_slice.loc[data_slice.pBC == pBC]
            pBC_slice = pBC_slice.values.tolist()
            for idx, row in enumerate(pBC_slice):
                umi = row[1]
                rBC = row[4]
                if idx == 0:
                    umi_to_keep.append(row[1])
                    rBC_to_keep.append(row[4])
                for umi_k in umi_to_keep:
                    if hamming_distance(umi_k, umi) <=2:
                        umi = umi_k
                        break
                for rBC_k in rBC_to_keep:
                    if hamming_distance(rBC_k, rBC) <=8:
                        rBC = rBC_k
                        break
                if rBC == row[4]:
                    rBC_to_keep.append(rBC)
                if umi == row[1]:
                    umi_to_keep.append(umi)

                info = tuple([row[0], umi, row[2], row[3], row[6]])
                if info not in correct_dict:
                     correct_dict[info] = (rBC, row[5])
                else:
                    correct_dict[info] = (correct_dict[info][0], correct_dict[info][1] + row[5])
            pop_list = []
            for key, val in correct_dict.items():
                info = list(key)
                info.append(val[0])
                info.append(val[1])
                pop_list.append(info)
        pop_df = pd.DataFrame(pop_list, columns = ['cellBC', 'umi', 'name', 'pBC', 'cluster', 'rBC', 'count'])
        # Now we remove the 1 plasmids from the pop_df
        pBCs = pop_df.pBC.unique()
        for pBC in pBCs:
            data_slice = pop_df.loc[pop_df.pBC == pBC]
            if len(data_slice) == len(data_slice.rBC.unique()) and len(data_slice.rBC.unique()) <= 1:
                zero_list = [[cell, 0, data_slice['name'].unique()[0], pBC ,data_slice['cluster'].unique()[0], 0, 0]]
                zero_df = pd.DataFrame(zero_list, columns = ['cellBC', 'umi', 'name', 'pBC', 'cluster', 'rBC', 'count'])
                trimed_pop_df = pd.concat([trimed_pop_df, zero_df])
            elif len(data_slice) == len(data_slice.rBC.unique()):
                trimed_pop_df = pd.concat([trimed_pop_df, data_slice])
            else:
                data_slice=data_slice[data_slice.rBC.duplicated(keep=False)]
                trimed_pop_df = pd.concat([trimed_pop_df, data_slice])
    end = time.time()
    print(f'The time is {end -start}')
    return  trimed_pop_df

def extract_sc_pBC_normed_exp(sexa_df):
    '''
    Input: a haxa-column def contains scMRPA information
    Output: df with [cellBC, name, pBC, norm_exp, cluster]
    '''
    pop_dict = {} 
    # Iterate through different clusters -- aka different cell type/states
    slice_df = sexa_df
    # A helper dict the records the total number of UMI for each cell
    for _,line in slice_df.iterrows():
        temp_key = line[['cellBC', 'name','pBC', 'cluster']]
        key = tuple(temp_key.values)
        rBC = str(line[['rBC']].values[0])
        print("rBC", rBC)
        if rBC == '0':
            if key not in pop_dict:
                pop_dict[key] = {}
                pop_dict[key][rBC] = 0
        else:
            if key not in pop_dict:
                pop_dict[key] = {}
            if rBC not in pop_dict[key]:
                pop_dict[key][rBC] = 1
            else:
                pop_dict[key][rBC] += 1
    pop_list = []
    print(pop_dict)
    print('finished counting')
    for key in pop_dict:
        info = list(key)
        current_cell = info[0]
        num_plasmid = len(pop_dict[key])
        # If we got an zero here 
        if any(list(pop_dict[key].keys())) == '0':
            norm_exp = 0
            direct_exp = 0
        else:
            norm_exp = sum(pop_dict[key].values())/num_plasmid
            direct_exp = sum(pop_dict[key].values())
        # The proper normalization for the expression in a cell is the total
        info.append(direct_exp)
        info.append(norm_exp)
        info.append(num_plasmid)
        pop_list.append(info)
    sc_measurement_df = pd.DataFrame(pop_list,  columns= ['cellBC', 'name', 'pBC', 'cluster', 'direct_exp', 'norm_exp', 'num_plasmid'])
    return sc_measurement_df


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    # Add the cell bc erorr corrected quint data
    parser.add_argument('--quint', help='path to quint file with clusters', required=True)
    # Add the informaiton aobut the library information that allows us to attach the promoter name to the list 
    parser.add_argument('--promlib', help='path to file named: prom_lib_master_info', required=True)
    # Add the name for the output file
    parser.add_argument('--out', help='output file name', required=True)
    parser.add_argument('--exp', help = 'name of the experiment', required = True)
    # Grab input arguments
    args= parser.parse_args()
    # Read in the quint file with columns as cellBC, umi, pBC,rBC, counts, cluster
    prom_quint = pd.read_csv(args.quint)
    # Read in the promoter lib infor
    prom_lib_info = pd.read_csv(args.promlib)
    # Add promoter id to the quint file
    prom_sext = attach_promoter_to_quad(prom_quint, prom_lib_info)
    # Save the attached promoter
    prom_sext.to_csv(args.exp + '_pBC_ec_' + args.out + '.csv',index = False)
    # Error correct the UMI and rBC, then adding 
    prom_sext_zero_padded = erorr_correct_bcs_v3(prom_sext)
    prom_sext_zero_padded.to_csv('./' + args.exp + '_pBC_ec_rBC_ec_' + args.out + '.csv',index = False)
    # Next we extract per cell expression for a barcode
    prom_exp_per_cell = extract_sc_pBC_normed_exp(prom_sext_zero_padded)
    prom_exp_per_cell.to_csv('./' + args.exp + '_pBC_sc_exp_' + args.out + '.csv',index = False)
    # 3. Return a pBC expression calculation with the min count with 2 but no min umi
    #
    # sc_exp_no_min_umi = extract_sc_pBC_normed_exp(prom_filtered)
    # sc_exp_no_min_umi.to_csv('../' + args.exp + '_rc_2_min_umi_0_pBC_exp_' + args.out + '.csv',index = False)
if __name__ == "__main__":
    main()

