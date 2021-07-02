"""                                                                                                                                                                                                        
script to calculate the differentially expressed protein groups between two phenotypes at protein level
In this case I will take three things into consideration later when selecting candidate protein groups as features

I- differential gene expression coefficients                                                                                                                                                              
II- number of samples the protein group is epxressed in (representativeness)
III- Average expression level, in case for choosing protein groups that are not shared between the two phenotypes

this script will output differentially expressed list of protein groups the protein groups that are shared between the two
phenotyeps, it will also output information about the proitein groups that are exlusive for just one of the phenotype
alongside with the number of samples they are expressed in                                                                                      
"""

import numpy as np
import pandas as pd
import math
from collections import Counter
import json
import os
import sys

in_dir = sys.argv[1]#/data/mstambou/proteome_landscapes/cluster_abundances_HAPiID_6k/proteinGroup_NSAFs_phenotypes/

phenotype_1 = sys.argv[2]#'healthy'
phenotype_2 = sys.argv[3]#'CD'

phen1_df = pd.read_csv(f'{in_dir}/{phenotype_1}_proteinGroup2NSAF.tsv', sep = '\t')
phen2_df = pd.read_csv(f'{in_dir}/{phenotype_2}_proteinGroup2NSAF.tsv', sep = '\t')

out_dir = sys.argv[4]#'/data/mstambou/proteome_landscapes/cluster_abundances_HAPiID_6k/proteinGroups_logOdds_phenotypes/'

if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
    
out_dir = f'{out_dir + phenotype_1}_vs_{phenotype_2}/'
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
    
print(f'{phenotype_1} expressed protein groups', phen1_df.shape)
print(f'{phenotype_2} expressed protein groups', phen2_df.shape)

common_expressed_proteinGroups = list(set(phen1_df['clusterID']).intersection(set(phen2_df['clusterID'])))
print('common expressed protein groups', len(common_expressed_proteinGroups))
phen1_only_proteinGroups = list(set(phen1_df['clusterID']).difference(set(common_expressed_proteinGroups)))
phen2_only_proteinGroups = list(set(phen2_df['clusterID']).difference(set(common_expressed_proteinGroups)))

phen1_phen2_df = pd.DataFrame(columns = ['clusterID', f'{phenotype_1}_rel_NSAF', f'{phenotype_2}_rel_NSAF', f'{phenotype_1}_n_samples', f'{phenotype_2}_n_samples'])

for i, clusterID in enumerate(common_expressed_proteinGroups):
    phen1_rel_NSAF = float(phen1_df[phen1_df['clusterID'] == clusterID]['NSAF_%'].values[0])
    phen2_rel_NSAF = float(phen2_df[phen2_df['clusterID'] == clusterID]['NSAF_%'].values[0])
    phen1_n_samples = int(phen1_df[phen1_df['clusterID'] == clusterID]['n_samples'].values[0])
    phen2_n_samples = int(phen2_df[phen2_df['clusterID'] == clusterID]['n_samples'].values[0])
    phen1_phen2_df.loc[i] = [clusterID, phen1_rel_NSAF, phen2_rel_NSAF, phen1_n_samples, phen2_n_samples]
    
log_odds = []

print("calculating log odds ...")
for i, row in phen1_phen2_df.iterrows():
    ratio = row[f'{phenotype_1}_rel_NSAF'] / row[f'{phenotype_2}_rel_NSAF']
    log_odds.append(math.log2(ratio))
    
phen1_phen2_df['log_odds_ratio'] = log_odds
phen1_phen2_df.sort_values(by = 'log_odds_ratio', inplace = True)
phen1_phen2_df.reset_index(inplace = True, drop = True)

print('writing to file ...')
phen1_phen2_df.to_csv(out_dir + f'{phenotype_1}_{phenotype_2}_logOdds.tsv', sep = '\t', index = False)

phen1_only_proteinGroups_df = phen1_df[(phen1_df['clusterID'].isin(phen1_only_proteinGroups))].sort_values(by = 'n_samples', ascending = False)
phen2_only_proteinGroups_df = phen2_df[(phen2_df['clusterID'].isin(phen2_only_proteinGroups))].sort_values(by = 'n_samples', ascending = False)

phen1_only_proteinGroups_df.to_csv(out_dir + f'{phenotype_1}_onlyProteinGroups.tsv', sep = '\t', index = False)
phen2_only_proteinGroups_df.to_csv(out_dir + f'{phenotype_2}_onlyProteinGroups.tsv', sep = '\t', index = False)
