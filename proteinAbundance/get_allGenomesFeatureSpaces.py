"""
script that will make combined feature spaces from the DESeq2 suggested genomes and the log odds ratio suggested genomes
"""
import copy
import pandas as pd
from multiprocessing import Pool
import os
import sys
import json


n_threads = 7

phenotype_pairs_f = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/comparisonLists.txt'

log2_ratios_dir = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/log2_ratios_genomes_scaled/'
samplesSubsamples2phenotypes_dic_f = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/samplesSubSamples2phenotypes_dic.json'
samples_dir = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/genome_spectralCount/'

out_dir = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/allGenomes_featureSpaces/'

in_suffix = '_genome2NspectraNormalized.tsv'

if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

with open(samplesSubsamples2phenotypes_dic_f, 'r') as in_f:
    samplesSubsamples2phenotypes_dic = json.load(in_f)

phens = []
with open(phenotype_pairs_f, 'r') as in_f:
     for line in in_f:
            line = line.strip().split('\t')
            phen1, phen2 = line[0], line[1]
            phens.append((phen1, phen2))

phen1, phen2 = 'healthy', 'CD'

def get_phenotypeFeatureSpace(phen):
    
    phen1 = phen[0]
    phen2 = phen[1]
    print(f'processing {phen1} vs {phen2} ...')
    phenotypes = {f'{phen1}' : 1, f'{phen2}': 0}

    log2_ratios_df = pd.read_csv(log2_ratios_dir + f'{phen1}_vs_{phen2}_genomeLog2oddsRatios.tsv', sep = '\t')

    featured_genomes = list(log2_ratios_df['genome_id'])

    feature_space = copy.deepcopy(featured_genomes)
    feature_space.insert(0, 'sample_name')
    feature_space.append('class')

    feature_space_idx = {item:i for i, item in enumerate(feature_space)}

    featured_df = pd.DataFrame(columns = feature_space)
    row_counter = 0
    for phenotype in phenotypes:
        phenotype_samples = samplesSubsamples2phenotypes_dic[phenotype]
        for sample in phenotype_samples:
            sample_df = pd.read_csv(samples_dir + sample + in_suffix, sep = '\t')
            sample_row = [0.0]*len(feature_space_idx)
            for i, row in sample_df.iterrows():
                genome = row['genome']
                if genome in featured_genomes:
                    rel_n_spectra = row['rel_n_spectra_%_normalized']
                    sample_row[feature_space_idx[genome]] += rel_n_spectra
            sample_row[feature_space_idx['class']] = phenotypes[phenotype]
            sample_row[feature_space_idx['sample_name']] = sample
            featured_df.loc[row_counter] = sample_row
            row_counter += 1

    print(f'{phen1} vs {phen2} feature space dimenions: {featured_df.shape[0]} data-points and {featured_df.shape[1] -2} features')
    featured_df.to_csv(out_dir + f'{phen1}_vs_{phen2}_allGenomesFeatureSpaces.tsv', sep = '\t', index = False)

    
with Pool(n_threads) as p:
    p.map(get_phenotypeFeatureSpace, phens)
