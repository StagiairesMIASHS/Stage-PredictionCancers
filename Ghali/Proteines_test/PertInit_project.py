from operator import index
import pandas as pd
from pandas import read_csv
import numpy as np
import argparse
import sys
import json
import os
import statistics
import matplotlib.pyplot as plt

###importation des données

pfam2go = read_csv('pfam.csv')
print(pfam2go)

###tri des données (équivalent du na.omit sur R)

"""
# On vérifie qu'il n'y a aucune donnée sur la colone, et dans ce cas, on la supprime
test = pfam2go[ 'metadata/member_databases']
print (test.shape) #on a 19179 lignes

x=float('NaN')
for i in range (19180) :
  if not (test[i] == x) :
      print ("good")

pfam2go = print(pd.read_csv('pfam.csv',keep_default_nan=False, na_values=[""]))
pfam2go[ pfam2go["metadata/go_terms"] != np.NaN ]
print (pfam2go['metadata/accession'])

def replaceUnknown(row):
    if row["notes"] == "unknown" :
        if row["game_location"] == "H" : #home
            return row["fran_id"]
        elif row["game_location"] == "A": #away
            return row["opp_fran"]
    return row["notes"]
"""

###Affichage des données sous forme de graphiques

"""
#plt.plot(pfam2go['metadata/accession'],pfam2go['metadata/type'])
#plt.show()

#plt.hist(pfam2go['metadata/type'],pfam2go['metadata/source_database'])
#plt.show()
"""

### On renomme les colonnes

"""
#pfam2go.index = ['accession' , 'name' , 'source_database' , 'type' , 'integrated' , 'member_databases', 'go_terms']
#print(pfam2go)


#Pour séparer les données, on a les catégories suivantes:
#metadata/accession , metadata/name , metadata/source_database , metadata/type , metadata/integrated , metadata/member_databases, metadata/go_terms

"""

###Programme source code

"""
if __name__ == "__main__":

    # ------------------------------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description='Highlight somatically perturbed mechanisms in putative cancer'
                                                 'driver genes.')
    parser.add_argument('--pertinint_results', type=str, help='Path to output file from PertInInt', default=None)
    parser.add_argument('--out_file', type=str, help='Results output file', default=None)
    parser.add_argument('--pfam2go_file', type=str, help='Path to Pfam2Go file (will be downloaded if not present)',
                        default='pfam2go.txt')
    parser.add_argument('--ligand_groups_file', type=str, help='Path to ligand ID -> ligand group mapping file ' +
                        '(will be downloaded if not present)', default='interacdome_ligandgrps.txt')
    args = parser.parse_args()

    # ------------------------------------------------------------------------------------------------
    # (i) confirm that we can open PertInInt results file
    if not args.pertinint_results or not os.path.isfile(args.pertinint_results):
        sys.stderr.write('Could not open PertInInt results file: '+str(args.pertinint_results)+'\\n' +
                         'Usage: python '+sys.argv[0]+' --pertinint_results <results_file> --out_file <output_file>\\n')
        sys.exit(1)

    # (ii) confirm that we can write to mechanisms output file
    if not args.out_file:
        sys.stderr.write('Could not write to specified output file: ' + str(args.out_file) + '\\n' +
                         'Usage: python '+sys.argv[0]+' --pertinint_results <results_file> --out_file <output_file>\\n')
        sys.exit(1)

    for subdir in ['/'.join(args.out_file.split('/')[:ind]) for ind in xrange(2, args.out_file.count('/')+1)]:
        if not os.path.isdir(subdir):
            if call(['mkdir', subdir]):  # any code returned other than "0"
                sys.stderr.write('Could not write to '+args.out_file+'. Exiting\\n')
                sys.exit(1)

    # ------------------------------------------------------------------------------------------------
    parse_ordered_genes(args.pertinint_results,
                        args.out_file,
                        args.pfam2go_file,
                        args.ligand_groups_file)
    sys.stderr.write('Wrote tab-delimited mechanisms to: '+args.out_file+'\\n')

"""