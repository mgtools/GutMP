"""
script to assign most specific taxonomy for each of the genomes
and then extract genomes that have the same species assignments but differentially expressed
at the strain levels
"""

import pandas as pd
from collections import Counter
import os
import copy
import sys
import json

in_dir = sys.argv[1]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/genome_logOdds_phenotypes/'
min_logOdds = float(sys.argv[2])#2.0
tax_level = sys.argv[3]#'s'

phenotype_1 = sys.argv[4]#'healthy'
phenotype_2 = sys.argv[5]#'CD'

allBin2taxon_dic_f = '/data/mstambou/proteome_landscapes/auxilary/allBin2taxon_dic.json'
with open(allBin2taxon_dic_f, 'r') as in_f:
    allBin2taxon_dic = json.load(in_f)
    
def getMST(genome, allBin2taxon_dic = allBin2taxon_dic, l = tax_level):
    if allBin2taxon_dic[genome][l] != '':
        return l + '__' + allBin2taxon_dic[genome][l]
    else:
        for c in 'sgfocpd':
            taxon = allBin2taxon_dic[genome][c]
            if taxon != '' and taxon != 'GCF_' and taxon != 'GCA_':
                return c+'__'+taxon
            
in_f = f'{in_dir}{phenotype_1}_vs_{phenotype_2}/{phenotype_1}_{phenotype_2}_logOdds.tsv'
phen1_f = f'{in_dir}{phenotype_1}_vs_{phenotype_2}/{phenotype_1}_onlyGenomes.tsv'
phen2_f = f'{in_dir}{phenotype_1}_vs_{phenotype_2}/{phenotype_2}_onlyGenomes.tsv'

out_dir = f'{in_dir}{phenotype_1}_vs_{phenotype_2}/'

in_df = pd.read_csv(in_f, sep = '\t')
phen1_df = pd.read_csv(phen1_f, sep = '\t')
phen2_df = pd.read_csv(phen2_f, sep = '\t')

in_df['taxa'] = [getMST(item) for item in in_df['genome']]
phen1_df['taxa'] = [getMST(item) for item in phen1_df['genome']]
phen2_df['taxa'] = [getMST(item) for item in phen2_df['genome']]

filtered_logODds_df = copy.deepcopy(in_df[(in_df['log_odds_ratio'] >= min_logOdds) | (in_df['log_odds_ratio'] <= -min_logOdds)])

DE_species = list({k:v for k,v in Counter(filtered_logODds_df['taxa']).items() if v > 1}.keys())

DE_species_df = copy.deepcopy(filtered_logODds_df[filtered_logODds_df['taxa'].isin(DE_species)])
DE_species_df.sort_values(by = 'taxa', inplace = True)
DE_species_df.reset_index(drop = True, inplace = True)
DE_taxa = []
for taxa in list(set(DE_species_df['taxa'])):
    taxa_rows = DE_species_df[DE_species_df['taxa'] == taxa]
    if (min(taxa_rows['log_odds_ratio']) < 0) & (max(taxa_rows['log_odds_ratio']) > 0):
        DE_taxa.append(taxa)

DE_species_df = DE_species_df[DE_species_df['taxa'].isin(DE_taxa)]
DE_species_df.to_csv(out_dir + f'{phenotype_1}_vs_{phenotype_2}_DE_common_genome_species.tsv', sep = '\t', index = False)
sameSpeciesDiffStrains = list(set(phen1_df['taxa']).intersection(set(phen2_df['taxa'])))

phen1_commomSpecies_df = copy.deepcopy(phen1_df[phen1_df['taxa'].isin(sameSpeciesDiffStrains)])
phen1_commomSpecies_df['phenotype'] = [phenotype_1]*len(phen1_commomSpecies_df)
phen2_commomSpecies_df = copy.deepcopy(phen2_df[phen2_df['taxa'].isin(sameSpeciesDiffStrains)])
phen2_commomSpecies_df['phenotype'] = [phenotype_2]*len(phen2_commomSpecies_df)
phen1_phen2_commonSpecies_combined_df = pd.concat([phen1_commomSpecies_df, phen2_commomSpecies_df])
phen1_phen2_commonSpecies_combined_df.sort_values(by = 'taxa', inplace = True)
phen1_phen2_commonSpecies_combined_df.to_csv(out_dir + f'{phenotype_1}_{phenotype_2}_commonSpecies_combined.tsv', sep = '\t', index = False)
