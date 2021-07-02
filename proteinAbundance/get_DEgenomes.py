"""
script to calculate the differentially expressed genomes between two phenotypes at the genome (i.e. the strain level)
In this case I will take three things into consideration later when selecting candidate genomes as features

I- differential gene expression coefficients
II- number of samples the genome is epxressed in (representativeness)
III- Average expression level, in case for choosing genomes that are not shared between the two phenotypes

this script will output differentially expressed list of genomes for the genomes that are shared between the two
phenotyeps, it will also output information about genomes that are exlusive for just one of the phenotype
alongside with the number of samples they are expressed in
"""

import numpy as np
import pandas as pd
import math
from collections import Counter
import json
import sys
import os

in_dir = sys.argv[1]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/genome_spectralCount_phenotypes/'

phenotype_1 = sys.argv[2]#'healthy'
phenotype_2 = sys.argv[3]#'CD'

out_dir = sys.argv[4]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/genome_logOdds_phenotypes/'

in_f1 = f'{in_dir + phenotype_1}_genome2NormalizedNspectra.tsv'
in_f2 = f'{in_dir + phenotype_2}_genome2NormalizedNspectra.tsv'


phen1_df =  pd.read_csv(in_f1, sep = '\t')
phen2_df = pd.read_csv(in_f2, sep = '\t')


if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
    
out_dir = out_dir + f'{phenotype_1}_vs_{phenotype_2}/'
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

print(f'{phenotype_1} Expressed genomes: {len(phen1_df)}')
print(f'{phenotype_2} Expressed genomes: {len(phen2_df)}')

common_expressed_genomes = list(set(phen1_df['genome']).intersection(set(phen2_df['genome'])))
print(f'Common Expressed genomes: {len(common_expressed_genomes)}')
phen1_only_genomes = list(set(phen1_df['genome']).difference(set(common_expressed_genomes)))
phen2_only_genomes = list(set(phen2_df['genome']).difference(set(common_expressed_genomes)))

phen1_phen2_df = pd.DataFrame(columns = ['genome', f'{phenotype_1}_rel_specs_normalized', f'{phenotype_2}_rel_specs_normalized', f'{phenotype_1}_n_samples', f'{phenotype_2}_n_samples' ])

for i, genome in enumerate(common_expressed_genomes):
    phen1_rel_spec = phen1_df[phen1_df['genome'] == genome]['rel_n_spectra_%_normalized'].values[0]
    phen2_rel_spec = phen2_df[phen2_df['genome'] == genome]['rel_n_spectra_%_normalized'].values[0]
    phen1_n_samples = phen1_df[phen1_df['genome'] == genome]['n_samples'].values[0]
    phen2_n_samples = phen2_df[phen2_df['genome'] == genome]['n_samples'].values[0]    
    phen1_phen2_df.loc[i] = [genome, phen1_rel_spec, phen2_rel_spec, phen1_n_samples, phen2_n_samples]
    
log_odds = []

for i, row in phen1_phen2_df.iterrows():
    ratio = row[f'{phenotype_1}_rel_specs_normalized']/row[f'{phenotype_2}_rel_specs_normalized']
    log_odds.append(math.log2(ratio))
    
phen1_phen2_df['log_odds_ratio'] = log_odds
phen1_phen2_df.sort_values(by = 'log_odds_ratio', inplace = True)
phen1_phen2_df.reset_index(inplace = True, drop = True)

phen1_phen2_df.to_csv(out_dir + f'{phenotype_1}_{phenotype_2}_logOdds.tsv', sep = '\t', index = False)


phen1_only_genomes_df = phen1_df[(phen1_df['genome'].isin(phen1_only_genomes))].sort_values(by = 'n_samples', ascending = False)
phen2_only_genomes_df = phen2_df[(phen2_df['genome'].isin(phen2_only_genomes))].sort_values(by = 'n_samples', ascending = False)

phen1_only_genomes_df.to_csv(out_dir + f'{phenotype_1}_onlyGenomes.tsv', sep = '\t', index = False)
phen2_only_genomes_df.to_csv(out_dir + f'{phenotype_2}_onlyGenomes.tsv', sep = '\t', index = False)
