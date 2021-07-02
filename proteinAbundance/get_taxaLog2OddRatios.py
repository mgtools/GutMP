"""
script to compute the log 2 odds rations for genome expression profiles between two pairs of samples belonging to different phenotypes
"""

import numpy as np
import pandas as pd
import math
from collections import Counter
import json
import sys
import os

phenotype_pairs_f = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/comparisonLists.txt'
samplesSubsamples2phenotypes_dic_f = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/samplesSubSamples2phenotypes_dic.json'

genome_profile_dir = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/phenotype_taxa_relNorm_profiles_scaled/'
out_dir = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/log2_ratios_taxa_scaled/'

tax_level = 's'

if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

with open(samplesSubsamples2phenotypes_dic_f, 'r') as in_f:
    samplesSubsamples2phenotypes_dic = json.load(in_f)
    with open(phenotype_pairs_f, 'r') as in_f:
        for line in in_f:
            line = line.strip().split('\t')
            phen1, phen2 = line[0], line[1]
            print(f'processing {phen1} vs {phen2} ...')
            phen1_samples = samplesSubsamples2phenotypes_dic[phen1]        
            phen2_samples = samplesSubsamples2phenotypes_dic[phen2]
            profile_df = pd.read_csv(genome_profile_dir + f'{phen1}_vs_{phen2}_{tax_level}_taxaProfiles.tsv', sep = '\t')
            with open(out_dir + f'{phen1}_vs_{phen2}_{tax_level}_taxaLog2oddsRatios.tsv', 'w') as out_f:
                out_f.write('taxa\tlog2_odds\n')                
                for i, row in profile_df.iterrows():
                    taxa = row['taxa']
                    taxa_phen1_profile = list(row.loc[phen1_samples])
                    taxa_phen2_profile = list(row.loc[phen2_samples])
                    taxa_phen1_exp = sum(taxa_phen1_profile)
                    taxa_phen2_exp = sum(taxa_phen2_profile)
                    ratio = taxa_phen1_exp / taxa_phen2_exp
                    log_odds = math.log2(ratio)                                        
                    out_f.write(f'{taxa}\t{log_odds}\n')
                    

                    
