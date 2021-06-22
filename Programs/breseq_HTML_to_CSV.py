#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 08:46:51 2021

@author: oliviakosterlitz
"""
import pandas as pd
import os
import argparse
import unicodedata

parser = argparse.ArgumentParser()

parser.add_argument("-p", "--plotting_csv", default = None)
parser.add_argument("-c", "--treatment_csv", default = None)
parser.add_argument("-f", "--filepairs", nargs='*')
args = parser.parse_args()
    
if len(args.filepairs) % 2:
    parser.error('filepairs arg should be pairs of values')

plot_option = args.plotting_csv
treatment_csv = args.treatment_csv
file_list = args.filepairs
file_list = zip(file_list[::2], file_list[1::2])

# for troubleshooting
#treatment_csv = "treatment.csv"
#file_list = ["EC_29_anc", "EC_29_anc/index.html", "EC_29_4", "EC_29_4/index.html"] 
#file_list = zip(file_list[::2], file_list[1::2])

def clean_dataframe (df):
    '''Cleans up the headings in the dataframe'''
    headings = []
    current_columns = df.columns.tolist()
    for i in current_columns:
        headings.append(i[1])
    df.columns = headings
    return (df)

# Merge all of the dataframes into one
def Combine_predicted_mutations (dict_dfs):
    '''Takes in a dictionary of dataframes. Combines them and adds the key of 
    the dictionary as an index'''
    dfs_keys = list(dict_dfs.keys())
    dfs_frames = []
    for i in dfs_keys:
        dfs_frames.append(dict_dfs[i])
    combined_frame = pd.concat(dfs_frames, keys = dfs_keys, sort=False)
    return (combined_frame)

dfs = {}
for i in file_list:
    sample_name = i[0]
    df = pd.read_html(i[1])
    predicted_mutations = clean_dataframe(df[1])
    predicted_mutations.to_csv((sample_name + '_predicted_mutations.csv'), encoding = "UTF-16")
    if len(df) >= 3:
        unassigned_missing_coverage = clean_dataframe(df[2])
        unassigned_missing_coverage.to_csv((sample_name + '_unassigned_missing_coverage.csv'), encoding = "UTF-16")
    if len(df) >= 4:
        unassigned_new_junction = clean_dataframe(df[3])
        unassigned_new_junction.to_csv((sample_name + '_unassigned_new_junction.csv'), encoding = "UTF-16")
    dfs[sample_name] = predicted_mutations



if plot_option != None:

    samples_combined = Combine_predicted_mutations(dfs)
    samples_combined = samples_combined.reset_index()
    samples_combined.to_csv('combined_predicted_mutations.csv', encoding = "UTF-16")
    
    df = samples_combined
    df = df.reset_index()
    df = df[["index", "level_0",  "position", "mutation", "annotation", "gene"]]
    df_dict = df.set_index('index').T.to_dict('list')
    strain_name = list(df.level_0.unique())
    strain_name_dict = {}
    
    for k in range(0, len(strain_name)):
        strain_name_dict[strain_name[k]] = k
    
    total = 0
    strain_dict = {}
    pandas_dict = {}
    for key in df_dict:
        strain = (df_dict[key][0])
        strain_number = strain_name_dict[strain]
        if treatment_csv != None:
            treatment_strain_map = pd.read_csv(treatment_csv)
            treatment = treatment_strain_map.treatment[treatment_strain_map.strains == strain].values[0]
        else: 
            treatment = strain
        position = (df_dict[key][1])
        variant_type = (df_dict[key][2])
        annotation = (df_dict[key][3])
        gene = (df_dict[key][4])
        if strain not in strain_dict:
            strain_dict[strain] = []
        if variant_type[1] == "→": #finds the SNPS
            variant = "snp"
            total += 1
            clean_gene = unicodedata.normalize("NFKD",gene).replace("←", "").replace("→","").replace(" ","").replace("/", ",")
            add_to_dict = [treatment, strain, strain_number, position, variant, 0, clean_gene]
            strain_dict[strain].append(add_to_dict)
            pandas_dict[key]= add_to_dict
        elif variant_type[0] == "Δ": #find the deletions
            variant = "del"
            deletion_size = int(variant_type.replace(",","").replace("Δ","").replace("bp",""))
            end_position = position + deletion_size
            clean_gene = unicodedata.normalize("NFKD",gene).replace("←", "").replace("→","").replace(" ","").replace("[", "").replace("]", "").replace("–", "/")
            clean_gene = clean_gene.split("/")
            add_to_dict_deletion_start = [treatment, strain, strain_number, position, variant, deletion_size, clean_gene[0]]
            add_to_dict_deletion_end = [treatment, strain, strain_number, end_position, variant, deletion_size*-1, clean_gene[-1]]
            strain_dict[strain].append(add_to_dict_deletion_start)
            strain_dict[strain].append(add_to_dict_deletion_end)
            pandas_dict[key]= add_to_dict_deletion_start
            pandas_dict[int(key)+.1]= add_to_dict_deletion_end
            total += 1
        else: #this is insertions and structural variants
            variant = "ins"
            clean_gene = unicodedata.normalize("NFKD",gene).replace("←", "").replace("→","").replace(" ","").replace("[", "").replace("]", "").replace("–", "/")
            clean_gene = clean_gene.split("/")
            add_to_dict = [treatment, strain, strain_number, position, variant, 0, clean_gene[0]]
            strain_dict[strain].append(add_to_dict)
            pandas_dict[key]= add_to_dict
            total += 1
    
    df_plot = pd.DataFrame.from_dict(pandas_dict, orient = 'index')
    df_plot = df_plot.rename(columns={0: "treatment", 1: "ID_strain", 2:"strain", 3:"base", 4:"mut.type", 5:"delta.bases",6:"gene"})
    df_plot = df_plot.reset_index()
    df_plot = df_plot.sort_values(by=['treatment', 'strain', 'base'])
    df_plot = df_plot.drop('index', axis = 1)
    df_plot.to_csv('plot_variants.csv')







