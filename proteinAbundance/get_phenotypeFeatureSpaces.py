"""
script that will be used for feature extraction from the log odds ratio tables, I will basically get the top most differentially 
expressed genomes accoss two classes of phenotypes and will use them as features to train an SVM later
"""

import pandas as pd
import json
import copy
import os
import sys

in_dir = sys.argv[1]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/genome_logOdds_phenotypes/'

phenotype_1 = sys.argv[2]#'healthy'
phenotype_2 = sys.argv[3]#'CD'

samples_dir = sys.argv[4]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/genome_spectralCount/'
out_dir = sys.argv[5]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/phenotype_featureSpaces/'
samples2phenotypes_dic_f = sys.argv[6]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/samples2phenotypes_dic.json'

min_logOdds = float(sys.argv[7])#2.0
top_n = int(sys.argv[8])#50
min_samples = int(sys.argv[9])#5

print(f'constructing a feature space for {phenotype_1} and {phenotype_2} samples, with minimum log odds of {min_logOdds}, \nand\
selecting top {top_n} most differentially expressed for each phenotype, and selecting genomes that are only expressed in \n\
one phenotype with a minimum of {min_samples} samples within the respective phenotype samples')

phenotypes = {f'{phenotype_1}' : 1, f'{phenotype_2}': 0}

logOdds_f = f'{in_dir}{phenotype_1}_vs_{phenotype_2}/{phenotype_1}_{phenotype_2}_logOdds.tsv'
logOdds_df = pd.read_csv(logOdds_f, sep = '\t')

phenotype1_only_f = f'{in_dir}{phenotype_1}_vs_{phenotype_2}/{phenotype_1}_onlyGenomes.tsv'
phenotype2_only_f = f'{in_dir}{phenotype_1}_vs_{phenotype_2}/{phenotype_2}_onlyGenomes.tsv'

phenotype1_only_df = pd.read_csv(phenotype1_only_f, sep = '\t')
phenotype2_only_df = pd.read_csv(phenotype2_only_f, sep = '\t')

in_suffix = '_genome2NspectraNormalized.tsv'


print(f'extracting features and creating a feature space for {phenotype_1} samples and {phenotype_2} samples')

if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
    
out_dir = f'{out_dir+phenotype_1}_vs_{phenotype_2}/'
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

with open(samples2phenotypes_dic_f, 'r') as in_f:
    samples2phenotypes_dic = json.load(in_f)
    
print('number of differentially expressed genomes', logOdds_df.shape[0])

filtered_logODds_df = logOdds_df[(logOdds_df['log_odds_ratio'] >= min_logOdds) | (logOdds_df['log_odds_ratio'] <= -min_logOdds)]
print('number of differentially expressed genomes with absolute log2 odds of greater than', min_logOdds, ':', filtered_logODds_df.shape[0])
print(f'number of genomes only expressed in {phenotype_1}, {len(phenotype1_only_df)}')
print(f'number of genomes only expressed in {phenotype_2}, {len(phenotype2_only_df)}')

phen1_filtered_logOdds_df = copy.deepcopy(filtered_logODds_df[filtered_logODds_df['log_odds_ratio'] > 0])
phen2_filtered_logOdds_df = copy.deepcopy(filtered_logODds_df[filtered_logODds_df['log_odds_ratio'] < 0])

phen1_filtered_logOdds_df.sort_values(by = f'{phenotype_1}_n_samples', inplace = True, ascending = False)
phen2_filtered_logOdds_df.sort_values(by = f'{phenotype_2}_n_samples', inplace = True, ascending = False)
phen1_filtered_logOdds_df.reset_index(drop = True, inplace = True)
phen2_filtered_logOdds_df.reset_index(drop = True, inplace = True)

featured_genomes = list(phen1_filtered_logOdds_df['genome'].loc[:top_n])
featured_genomes.extend(list(phen2_filtered_logOdds_df['genome'].loc[:top_n]))
featured_genomes.extend(list(phenotype1_only_df[phenotype1_only_df['n_samples'] >= min_samples]['genome']))
featured_genomes.extend(list(phenotype2_only_df[phenotype2_only_df['n_samples'] >= min_samples]['genome']))

feature_space = copy.deepcopy(featured_genomes)
feature_space.insert(0, 'sample_name')
feature_space.append('class')

feature_space_idx = {item:i for i, item in enumerate(feature_space)}

featured_df = pd.DataFrame(columns = feature_space)
row_counter = 0
for phenotype in phenotypes:
    phenotype_samples = samples2phenotypes_dic[phenotype]
    for sample in phenotype_samples:
        sample_df = pd.read_csv(samples_dir + sample + in_suffix, sep = '\t')
        sample_row = [0.0]*len(feature_space_idx)
        for i, row in sample_df.iterrows():
            genome = row['genome']
            if genome in featured_genomes:
                rel_n_spectra = row['rel_n_spectra_%_normalized']
                sample_row[feature_space_idx[genome]] = rel_n_spectra
        sample_row[feature_space_idx['class']] = phenotypes[phenotype]
        sample_row[feature_space_idx['sample_name']] = sample
        featured_df.loc[row_counter] = sample_row
        row_counter += 1
print(f'feature space dimenions: {featured_df.shape[0]} data-points and {featured_df.shape[1] -2} features')
featured_df.to_csv(out_dir + f'{phenotype_1}_vs_{phenotype_2}_featureSpaces_minLogOdds_{min_logOdds}_top_n_{top_n}_minSamples_{min_samples}.tsv', sep = '\t', index = False)
