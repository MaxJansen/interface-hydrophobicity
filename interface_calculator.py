# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 13:10:21 2021

@author: Max Jansen
"""

# Import packages
import os
import argparse
import time
import numpy as np
import pandas as pd
import Bio
from Bio.PDB import *

# Time how long the whole script takes
START_TIME = time.time()


def convert_pdb(pdb_file):
    """Selects lines from pdb-file that contain information about the atoms
    in the structure and converts this to a pandas DataFrame."""

    selected = [[]]
    for line in pdb_file:
        line_list = line.split()
        if line_list[0] == 'ATOM':
            selected.append(line_list)
    del selected[0]
    pdb_atom_df = pd.DataFrame(selected)
    return pdb_atom_df


def comp_dist(pdb_atom_line, pdb_atom_line2):
    """Computes euclidian distance between two xyz-coordinates from a pdb."""

    arr1 = np.array(pdb_atom_line[6:9])
    arr1 = arr1.astype(np.float)
    arr2 = np.array(pdb_atom_line2[6:9])
    arr2 = arr2.astype(np.float)
    dist = np.linalg.norm(arr2-arr1)
    return dist


def iterate_chain(input_df, chain_1, chain_2):
    """Iterates throught two chains selected from a pdb-files and determines
    which CA atoms (from a respective residue) from either protein chains are
    within a selected distance threshold from each other."""

    atom_match_list = []
    atom_match_df = False
    protein_ca_df = input_df[input_df[2] == 'CA']
    df_1 = protein_ca_df[protein_ca_df[4] == chain_1]
    df_2 = protein_ca_df[protein_ca_df[4] == chain_2]

    for index, row in df_1.iterrows():
        for index_t, row_t in df_2.iterrows():
            dist = comp_dist(row, row_t)
            # Distance threshold
            if dist < 7:
                # Select relevant data from atom-atom distances, remove index
                    # if you want to select the whole row from the df
                    # Currently selects residue name and residue seq number
                atom_match = [row[3:6:2], row_t[3:6:2]]
                atom_match_list.append(atom_match)

    # What happens if no atoms from two chains are within distance threshold?
    if len(atom_match_list) != 0:
        atom_match_df = pd.DataFrame(atom_match_list)
        atom_match_df.columns = [chain_1, chain_2]
    else:
        pass
    return atom_match_df


def chain_hydrophobicity(df_of_aa):
    """Takes a df with two chains out of the protein as input and returns
    the chain names and their respective average hydrophobicity scores"""

    chains_list = list(df_of_aa.columns)

    interface_1 = df_of_aa[[chains_list[0]]]
    hydro_df_1 = add_hydrophobicity(interface_1, chains_list[0])
    print(hydro_df_1)
    chain_1_score = [chains_list[0], hydro_df_1['Hydrophobicity'].mean()]

    interface_2 = df_of_aa[[chains_list[1]]]
    hydro_df_2 = add_hydrophobicity(interface_2, chains_list[1])
    print(hydro_df_2)
    chain_2_score = [chains_list[1], hydro_df_2['Hydrophobicity'].mean()]

    return[chain_1_score, chain_2_score]


def add_hydrophobicity(chain_df, chain_name):
    """Takes chain from a protein as input in pd df format. Maps hydrophobicity
    scores from a dict onto each residue. Ensures that each residue is only
    counted once."""

    hydro_dict = {"PHE": 100, "ILE":	99, "TRP": 97, "LEU":	97, "VAL":	76,
                  "MET"	: 74, "TYR":	63, "CYS": 49, "ALA":	41, "THR":	13,
                  "HIS": 8, "GLY": 0, "SER": -5, "GLN": -10, "ARG": -14,
                  "LYS": -23, "ASN": -28, "GLU": -31, "PRO": -46, "ASP": -55}

    chain_df = chain_df[chain_name].apply(lambda x: x)
    chain_df.columns = ['Res_name', 'Res_seq_num']

    chain_df.drop_duplicates(subset="Res_seq_num", keep=False, inplace=True)
    chain_df['Hydrophobicity'] = chain_df['Res_name'].map(hydro_dict)

    return chain_df


def parser():
    """Retrieves the arguments from the command line."""

    ding = argparse.ArgumentParser(description='A program that computes \
                                   hydrophobicity of antibody interfaces.')
    ding.add_argument('-ab', required=True, metavar='pdb_id_and_chains',
                      dest='pdbAndChains', help='[-ab] S')
    arguments = ding.parse_args()  # takes the arguments
    return [arguments]


# ==============================================================================
#  Main commands
# ==============================================================================

# Get command line arguments, allows running the script for
# a whole list of pdb-ids.
COMMAND_LINE = parser()
FULL_PDB_ID = COMMAND_LINE[0].pdbAndChains

# Split the long pdb id into parts. Allows matching former chains to latter
# chains
PDB_ID_LIST = FULL_PDB_ID.split("_")
PDB_ID = PDB_ID_LIST[0].lower()
FORMER_CHAINS = PDB_ID_LIST[1]
LATTER_CHAINS = PDB_ID_LIST[2]

# DOWNLOAD PDB FILE
pdbl = PDBList()
prot_wd = os.getcwd()
pdb_filename = pdbl.retrieve_pdb_file(PDB_ID, pdir=prot_wd+'/data',
                                      file_format='pdb')
# Get a dataframe from pdb file
F = open('data/pdb' + PDB_ID + ".ent", "r")
X = F.readlines()
PROTEIN_DF = convert_pdb(X)

# Determine unique chains based on column in pdb, compare to filename.
UNIQUE_CHAINS = PROTEIN_DF[4].unique()

for i in FORMER_CHAINS + LATTER_CHAINS:
    if any(i in s for s in UNIQUE_CHAINS):
        pass
    else:
        print("Warning! Chain listed if in filename is not present in ATOM\
              part of pdb")
# Make a dict that you can use to write final output
d = {"Protein_pdb": [], "Total_mean_hydrophobicity": [],
     "Former_chains_hydrophobicity": []}

# Iterate through two chains at a time to find boundary between the two
chain_score_list = []
for chain_id1 in FORMER_CHAINS:
    for chain_id2 in LATTER_CHAINS:
        matched_chains_df = iterate_chain(PROTEIN_DF, chain_id1, chain_id2)
        # What to do if two chains do not match?
        if type(matched_chains_df) == pd.core.frame.DataFrame:
            chainchain_scores = chain_hydrophobicity(matched_chains_df)
            chain_score_list.append(chainchain_scores)

        else:
            pass

# Make a neat df with rows for each chain match and four columns:
    # One for former chain, former score, latter chain and latter score.
chain_score_df = pd.DataFrame(chain_score_list, columns=[
        "former_chain", "latter_chain"])
score_df1 = pd.DataFrame(chain_score_df['former_chain'].tolist(), columns=[
        'Former_chains', 'Former_scores'])
score_df2 = pd.DataFrame(chain_score_df['latter_chain'].tolist(), columns=[
        'Latter_chains', 'Latter_scores'])
full_score_df = pd.concat([score_df1, score_df2], axis=1)

# Final scores, take mean of all chain match scores to get 1 final pdb name and
# 2 summary statistics. Write to file
d["Protein_pdb"].append(FULL_PDB_ID)
d["Total_mean_hydrophobicity"].append(pd.Series([
        full_score_df.Former_scores.mean(),
        full_score_df.Latter_scores.mean()]).mean())
d["Former_chains_hydrophobicity"].append(full_score_df.Former_scores.mean())

f = open("interface_hydrophobicity_output.txt", "a")
f.write(str(d) + '\n')
f.close()

print(d)
END_TIME = time.time()
print(str(END_TIME - START_TIME) + "seconds")
