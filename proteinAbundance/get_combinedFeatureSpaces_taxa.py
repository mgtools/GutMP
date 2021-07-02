"""
script that will make combined feature spaces from the DESeq2 suggested genomes and the log odds ratio suggested genomes
"""
import copy
import pandas as pd
from multiprocessing import Pool
import os
import sys
import json

min_log_odds = 2.0
min_p_val = 0.05

n_threads = 7

phenotype_pairs_f = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/comparisonLists.txt'

DESeq2_dir = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/DE_taxa/'
log2_ratios_dir = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/log2_ratios_taxa_scaled/'
samplesSubsamples2phenotypes_dic_f = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/samplesSubSamples2phenotypes_dic.json'
samples_dir = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/genome_spectralCount/'
allBin2taxon_dic_f = '/data/mstambou/proteome_landscapes/auxilary/allBin2taxon_dic.json'

out_dir = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/DESeq2_logOdds_combinedFeatureSpaces_taxa/'

in_suffix = '_genome2NspectraNormalized.tsv'

tax_level = 's'

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

phen1, phen2 = 'healthy', 'CD'

def get_phenotypeFeatureSpace(phen):
    
    phen1 = phen[0]
    phen2 = phen[1]
    print(f'processing {phen1} vs {phen2} ...')
    phenotypes = {f'{phen1}' : 1, f'{phen2}': 0}

    DESeq2_df = pd.read_csv(DESeq2_dir + f'{phen1}_vs_{phen2}_{tax_level}_taxaProfiles_DE_taxa.csv', sep = '\t')
    log2_ratios_df = pd.read_csv(log2_ratios_dir + f'{phen1}_vs_{phen2}_{tax_level}_taxaLog2oddsRatios.tsv', sep = '\t')

    DESeq_taxa = list(DESeq2_df[DESeq2_df['padj'] <= min_p_val].index)
    log_2_ratios_taxa = list(log2_ratios_df[(log2_ratios_df['log2_odds'] <= -min_log_odds) | (log2_ratios_df['log2_odds'] >= min_log_odds)]['taxa'])

    featured_taxa = list(set(DESeq_taxa).union(set(log_2_ratios_taxa)))

    feature_space = copy.deepcopy(featured_taxa)
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
                taxa = getMST(genome)
                if taxa in featured_taxa:
                    rel_n_spectra = row['rel_n_spectra_%_normalized']
                    sample_row[feature_space_idx[taxa]] += rel_n_spectra
            sample_row[feature_space_idx['class']] = phenotypes[phenotype]
            sample_row[feature_space_idx['sample_name']] = sample
            featured_df.loc[row_counter] = sample_row
            row_counter += 1

    print(f'{phen1} vs {phen2} feature space dimenions: {featured_df.shape[0]} data-points and {featured_df.shape[1] -2} features')
    featured_df.to_csv(out_dir + f'{phen1}_vs_{phen2}_{tax_level}_combinedFeatureSpaces_minLogOdds_{min_log_odds}_minPval_{min_p_val}.tsv', sep = '\t', index = False)

    
with Pool(n_threads) as p:
    p.map(get_phenotypeFeatureSpace, phens)
