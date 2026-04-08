#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed 8 April 11:22:02 2026

@author: ahmed
"""
##############

import pandas as pd
from Bio import Entrez
from ete3 import NCBITaxa, TreeStyle, PieChartFace, faces
import argparse
import warnings
import math
##############

Entrez.email = 'drahmedelsherbini@gmail.com'


warnings.filterwarnings("ignore")

# ----------------------- ARGUMENT PARSING ----------------------- #
my_parser = argparse.ArgumentParser(description='NCBI Species Tree Builder with just a filter')
my_parser.add_argument('-i', '--input', metavar='input', type=str, required=False, help='Path to input CSV file')
args = my_parser.parse_args()


# ----------------------- MAIN VARIABLES ----------------------- #
f_name = args.input
##############

#f_name = "species.csv"
df = pd.read_csv(f_name, on_bad_lines='skip')
ncbi = NCBITaxa()
species_taxids = {}
##############
species_names = pd.Series(df["Species"])
species_names = species_names[~species_names.str.contains("sp.", regex=False)]
species_names = species_names[~species_names.str.contains("Species", regex=False)]
##############
for species_name in species_names:
    taxid = ncbi.get_name_translator([species_name])
    species_taxids[species_name] = taxid[species_name][0] if taxid else None

taxa_ids = [taxid for taxid in species_taxids.values() if taxid]
tree = ncbi.get_topology(taxa_ids)
##############
def annotate_tree_with_scientific_names(tree):
    for node in tree.traverse():
        if node.is_leaf():
            taxid = int(node.name)
            scientific_name = ncbi.get_taxid_translator([taxid]).get(taxid)
            node.name = scientific_name if scientific_name else "Unknown"
##############
annotate_tree_with_scientific_names(tree)
print(tree)
output_file = f"{f_name[:-4]}_tree.nwk"
tree.write(outfile=output_file)
print("Tree with pie charts rendered successfully!")
##############
