"""
script to create protein group profile matrices for log2 folds ratios calculation
"""

from multiprocessing import Pool
import pandas as pd
import json
import os

phenotype_pairs_f = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/comparisonLists.txt'
phen2samples_dic_f = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/samplesSubSamples2phenotypes_dic.json'
clusterAbundance_dir = '/data/mstambou/proteome_landscapes/cluster_abundances_HAPiID_6k/cluster_NSAFs/'

out_dir = '/data/mstambou/proteome_landscapes/cluster_abundances_HAPiID_6k/phenotype_proteinGroup_abundance_profiles_log2folds/'

suffix = '_covering_80_percent_ribP_elonF_spectra_NSAF.tsv'

#scaling_factor = 1000
n_pool = 7
phen2proteinGroups_dic_f = '/data/mstambou/proteome_landscapes/cluster_abundances_HAPiID_6k/phen2proteinGroups_dic.json'
with open(phen2proteinGroups_dic_f, 'r') as in_f:
    phen2proteinGroups_dic = json.load(in_f)

if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

with open(phen2samples_dic_f, 'r') as in_f:
    phen2samples_dic = json.load(in_f)

with open(phenotype_pairs_f, 'r') as in_f:
    phenotype_pairs = []
    for line in in_f:
        line = line.strip().split('\t')
        phen1, phen2 = line[0], line[1]
        phenotype_pairs.append((phen1, phen2))
        
def get_phenPair2protGroupProfiles(phen_pair, out_dir = out_dir, phen2proteinGroups_dic = phen2proteinGroups_dic):
    phen1, phen2 = phen_pair[0], phen_pair[1]

    print(f'processing {phen1} vs {phen2} ...')

    phen1_samples = phen2samples_dic[phen1]
    phen2_samples = phen2samples_dic[phen2]

    pair_df = pd.DataFrame()
    pair_clusters = phen2proteinGroups_dic[phen1]
    pair_clusters.extend(phen2proteinGroups_dic[phen2])
    pair_clusters = list(set(pair_clusters))

    cluster2idx_dic = {clstr:i for i, clstr in enumerate(pair_clusters)}

    pair_df['cluster_ID'] = pair_clusters
    metadata_df = pd.DataFrame()
    all_samples, all_phenotypes = [], []
    for sample in phen1_samples:
        all_samples.append(sample)
        all_phenotypes.append(phen1)
        sample_df = pd.read_csv(clusterAbundance_dir + sample + suffix, sep = '\t')
        sample_col = [0.0001]*len(cluster2idx_dic)
        for i, row in sample_df.iterrows():
            clusterID = row['clusterID']
            n_specs = row['NSAF_%']
            if n_specs == 0:
                n_specs = 1
            sample_col[cluster2idx_dic[clusterID]] = n_specs
        pair_df[sample] = sample_col

    for sample in phen2_samples:
        all_samples.append(sample)
        all_phenotypes.append(phen2)
        sample_df = pd.read_csv(clusterAbundance_dir + sample + suffix, sep = '\t')
        sample_col = [0.0001]*len(cluster2idx_dic)
        for i, row in sample_df.iterrows():
            clusterID = row['clusterID']
            n_specs = row['NSAF_%']
            if n_specs == 0:
                n_specs = 1
            sample_col[cluster2idx_dic[clusterID]] = n_specs
        pair_df[sample] = sample_col
    metadata_df['sample_ID'] = all_samples
    metadata_df['phenotype'] = all_phenotypes

    pair_df.to_csv(out_dir + f'{phen1}_vs_{phen2}_proteinGroupProfiles.tsv', sep = '\t', index = False)
        
with Pool(n_pool) as p:
    p.map(get_phenPair2protGroupProfiles, phenotype_pairs)
