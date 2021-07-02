"""
script that will be used for feature extraction from the log odds ratio tables, I will basically get the top most differentially 
expressed genomes accoss two classes of phenotypes and will use them as features to train an SVM later
"""

import pandas as pd
import json
import copy
import sys
import os 

in_dir = sys.argv[1]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/genome_logOdds_phenotypes/'

phenotype_1 = sys.argv[2]#'healthy'
phenotype_2 = sys.argv[3]#'CD'

samples_dir = sys.argv[4]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/genome_spectralCount/'
out_dir = sys.argv[5]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/phenotype_featureSpaces/'
samples2phenotypes_dic_f = sys.argv[6]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/samples2phenotypes_dic.json'

min_logOdds = float(sys.argv[7])#2.0
top_n = int(sys.argv[8])#100
min_samples = int(sys.argv[9])#3

allBin2taxon_dic_f = sys.argv[10]#'/data/mstambou/proteome_landscapes/auxilary/allBin2taxon_dic.json'
tax_level = sys.argv[11]#'s'

print(f'constructing a feature space for {phenotype_1} and {phenotype_2} samples, with minimum log odds of {min_logOdds}, \nand\
selecting top {top_n} most differentially expressed for each phenotype, and selecting genomes at the taxonomic level {tax_level} that are only expressed in \n\
one phenotype with a minimum of {min_samples} samples within the respective phenotype samples')

in_suffix = '_genome2NspectraNormalized.tsv'

phenotypes = {f'{phenotype_1}' : 1, f'{phenotype_2}': 0}

in_dir = in_dir + f'{phenotype_1}_vs_{phenotype_2}_{tax_level}_taxa/'
out_dir = out_dir + f'{phenotype_1}_vs_{phenotype_2}_{tax_level}_taxa/'

print(f'extracting features and creating a feature space for {phenotype_1} samples and {phenotype_2} samples')

if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

logOdds_f = in_dir + f'{phenotype_1}_{phenotype_2}_{tax_level}_logOdds.tsv'
logOdds_df = pd.read_csv(logOdds_f, sep = '\t')

phenotype1_only_f = f'{in_dir + phenotype_1}_only_{tax_level}.tsv'
phenotype2_only_f = f'{in_dir + phenotype_2}_only_{tax_level}.tsv'

phenotype1_only_df = pd.read_csv(phenotype1_only_f, sep = '\t')
phenotype2_only_df = pd.read_csv(phenotype2_only_f, sep = '\t')

print('there are', len(phenotype1_only_df), tax_level, 'taxonimic level genomes only in ', phenotype_1)
print('there are', len(phenotype2_only_df), tax_level, 'taxonimic level genomes only in ', phenotype_2)



with open(samples2phenotypes_dic_f, 'r') as in_f:
    samples2phenotypes_dic = json.load(in_f)
    
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

print('number of differentially expressed genomes', logOdds_df.shape[0])
            
filtered_logODds_df = logOdds_df[(logOdds_df['log_odds_ratio'] >= min_logOdds) | (logOdds_df['log_odds_ratio'] <= -min_logOdds)]
print('number of differentially expressed genomes with absolute log2 odds of greater than', min_logOdds, ':', filtered_logODds_df.shape[0])
print(f'number of genomes only expressed in {phenotype_1}, {len(phenotype1_only_df)}')
print(f'number of genomes only expressed in {phenotype_2}, {len(phenotype2_only_df)}')

phen1_filtered_logOdds_df = copy.deepcopy(filtered_logODds_df[filtered_logODds_df['log_odds_ratio'] > 0])
phen2_filtered_logOdds_df = copy.deepcopy(filtered_logODds_df[filtered_logODds_df['log_odds_ratio'] < 0])

print('there are ', len(phen1_filtered_logOdds_df), 'differentially enriched ', tax_level, 'genomes enriched in', phenotype_1,'with minimum log odds of ', min_logOdds)
print('there are ', len(phen2_filtered_logOdds_df), 'differentially enriched ', tax_level, 'genomes enriched in', phenotype_2,'with minimum log odds of ', min_logOdds)

phen1_filtered_logOdds_df.sort_values(by = f'{phenotype_1}_n_samples', inplace = True, ascending = False)
phen2_filtered_logOdds_df.sort_values(by = f'{phenotype_2}_n_samples', inplace = True, ascending = False)
phen1_filtered_logOdds_df.reset_index(drop = True, inplace = True)
phen2_filtered_logOdds_df.reset_index(drop = True, inplace = True)

featured_taxa = list(phen1_filtered_logOdds_df['taxa'].loc[:top_n])
featured_taxa.extend(list(phen1_filtered_logOdds_df['taxa'].loc[:top_n]))
featured_taxa.extend(list(phenotype1_only_df[phenotype1_only_df['n_samples'] >= min_samples]['taxa']))
featured_taxa.extend(list(phenotype2_only_df[phenotype2_only_df['n_samples'] >= min_samples]['taxa']))

featured_taxa = list(set(featured_taxa))

feature_space = copy.deepcopy(featured_taxa)
feature_space.insert(0, 'sample_name')
feature_space.append('class')


feature_space_idx = {item:i for i, item in enumerate(feature_space)}

print('feature space size', len(feature_space) -2)

featured_df = pd.DataFrame(columns = feature_space)
row_counter = 0
for phenotype in phenotypes:
    phenotype_samples = samples2phenotypes_dic[phenotype]
    for sample in phenotype_samples:
        sample_df = pd.read_csv(samples_dir + sample + in_suffix, sep = '\t')
        sample_row = [0.0]*len(feature_space_idx)
        for i, row in sample_df.iterrows():
            genome = row['genome']
            taxa = getMST(genome)
            if taxa in featured_taxa:
                rel_n_spectra = row['rel_n_spectra_%_normalized']
                sample_row[feature_space_idx[taxa]] = rel_n_spectra
        sample_row[feature_space_idx['class']] = phenotypes[phenotype]
        sample_row[feature_space_idx['sample_name']] = sample
        featured_df.loc[row_counter] = sample_row
        row_counter += 1
print(f'feature space dimenions: {featured_df.shape[0]} data-points and {featured_df.shape[1] -2} features')
featured_df.to_csv(out_dir + f'{phenotype_1}_vs_{phenotype_2}_{tax_level}_featureSpaces_minLogOdds_{min_logOdds}_top_n_{top_n}_minSamples_{min_samples}.tsv', sep = '\t', index = False)
