"""
script to assign most specific taxonomy for each of the genomes
and then extract genomes that are differentially expressed
at the strain levels, starting from the set of genomes suggested by DESeq2 program
"""

import pandas as pd
from collections import Counter
import os
import copy
import sys
import json

in_dir = sys.argv[1]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/DE_genomes/'
min_p_val = float(sys.argv[2])#0.05
tax_level = sys.argv[3]#'s'

phenotype_1 = sys.argv[4]#'healthy'
phenotype_2 = sys.argv[5]#'CD'

#in_dir = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/DE_genomes/'
#min_p_val = 0.05
#tax_level = 's'
#phenotype_1 = 'healthy_subT1D'
#phenotype_2 = 'T1D'

reverse_list = ['healthy_subT1D_vs_T1D', 'healthy_subUC_vs_UC', 'healthy_vs_UC', 'CD_subUC_vs_UC']

phen2samples_dic_f = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/samplesSubSamples2phenotypes_dic.json'
with open(phen2samples_dic_f, 'r') as in_f:
    phen2samples_dic = json.load(in_f)
    
    
genome2sample2points_dic_f = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/genome2sample2points_dic.json'
with open(genome2sample2points_dic_f, 'r') as in_f:
    genome2samples_dic = json.load(in_f)
    
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

phen1_samples = set(phen2samples_dic[phenotype_1])
phen2_samples = set(phen2samples_dic[phenotype_2])
            
in_f = f'{in_dir}{phenotype_1}_vs_{phenotype_2}_genomeProfiles_DE_genomes.csv'

out_dir = f'{in_dir}{phenotype_1}_vs_{phenotype_2}/'

if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

in_df = pd.read_csv(in_f, sep = '\t')

in_df['genome'] = list(in_df.index)
in_df['taxa'] = [getMST(item) for item in list(in_df.index)]

strain2phen_df = pd.DataFrame(columns = ['genome', 'taxa', 'phenotype'])
for i, row in in_df.iterrows():
    genome = row['genome']
    taxa = row['taxa']
    genome_samples = set(genome2samples_dic[genome].keys())
    genome_phen1_samples = genome_samples.intersection(phen1_samples)
    genome_phen2_samples = genome_samples.intersection(phen2_samples)
    
    if len(genome_phen1_samples) > 0 and len(genome_phen2_samples) > 0:
        phen = 'both'
    elif len(genome_phen1_samples) > 0 and len(genome_phen2_samples) == 0:
        phen = phenotype_1
    elif len(genome_phen1_samples) == 0 and len(genome_phen2_samples) > 0:
        phen = phenotype_2
        
    strain2phen_df.loc[i] = [genome, taxa, phen]    

#now after this point everything is DE
filtered_logOdds_df = copy.deepcopy(in_df[in_df['padj'] <= min_p_val]) 


phen1_DE_strains_df = pd.DataFrame(columns = list(filtered_logOdds_df.columns))
phen2_DE_strains_df = pd.DataFrame(columns = list(filtered_logOdds_df.columns))
for i, row in filtered_logOdds_df.iterrows():
    if f'{phenotype_1}_vs_{phenotype_2}' in reverse_list:
        if row['log2FoldChange'] <= 0:
            phen1_DE_strains_df.loc[len(phen1_DE_strains_df)] = list(row)
        else:
            phen2_DE_strains_df.loc[len(phen2_DE_strains_df)] = list(row)
    else:
        if row['log2FoldChange'] >= 0:
            phen1_DE_strains_df.loc[len(phen1_DE_strains_df)] = list(row)
        else:
            phen2_DE_strains_df.loc[len(phen2_DE_strains_df)] = list(row)

phen1_DE_strains_df.sort_values(by = ['taxa'], inplace = True)
phen2_DE_strains_df.sort_values(by = ['taxa'], inplace = True)
strain2phen_df.sort_values(by = ['taxa'], inplace = True)
phen1_DE_strains_df.to_csv(out_dir + f'{phenotype_1}_DE_strains.tsv', sep = '\t', index = False)
phen2_DE_strains_df.to_csv(out_dir + f'{phenotype_2}_DE_strains.tsv', sep = '\t', index = False)
strain2phen_df.to_csv(out_dir + f'{phenotype_1}_vs_{phenotype_2}_strain2phen.tsv', sep = '\t', index = False)
