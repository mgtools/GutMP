"""
script to create phenotype comparison tables for taxonomies, at specified level, protein expression profiles, and metadata files accompaniying them.
These tables will be used as inputs for the DESeq package in R for differential genome analysis and extracting meaningfull
differentially expressed taxonomic groups from samples.
"""

import pandas as pd
import json
import os

phenotype_pairs_f = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/comparisonLists.txt'
samplesSubsamples2phenotypes_dic_f = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/samplesSubSamples2phenotypes_dic.json'
allBin2taxon_dic_f = '/data/mstambou/proteome_landscapes/auxilary/allBin2taxon_dic.json'

samples_summary_f = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/samples_1_genomes_1000_peps_withPhenotypes.tsv'
genome_abundance_dir = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/genome_spectralCount/'

out_dir = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/phenotype_taxa_abundance_profiles_DESeq2/'

tax_level = 's'

if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

with open(samplesSubsamples2phenotypes_dic_f, 'r') as in_f:
    samplesSubsamples2phenotypes_dic = json.load(in_f)
    
samples_summary_df = pd.read_csv(samples_summary_f, sep = '\t')

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
    
with open(phenotype_pairs_f, 'r') as in_f:
    for line in in_f:
        line = line.strip().split('\t')
        phen1, phen2 = line[0], line[1]
        print(f'processing {phen1} vs {phen2} ...')
        phen1_samples = samplesSubsamples2phenotypes_dic[phen1]
        phen1_df = samples_summary_df[samples_summary_df['sample_name'].isin(phen1_samples)]
        pair_genomes = [item for sublist in list(phen1_df['expanded_genomes']) for item in sublist.split('|')]
        phen2_samples = samplesSubsamples2phenotypes_dic[phen2]
        phen2_df = samples_summary_df[samples_summary_df['sample_name'].isin(phen2_samples)]
        pair_genomes.extend([item for sublist in list(phen2_df['expanded_genomes']) for item in sublist.split('|')])
        
        pair_genomes = list(set(pair_genomes))
        
        pair_taxa = list(set([getMST(genome) for genome in pair_genomes]))
        
        taxa2idx_dic = {taxa:i for i, taxa in enumerate(pair_taxa)}
        pair_df = pd.DataFrame()
        pair_df['taxa'] = pair_taxa
        pair_metadata_df = pd.DataFrame(columns = ['sample_ID', 'phenotype'])
        
        for sample in phen1_samples:
            sample_df = pd.read_csv(f'{genome_abundance_dir + sample}_genome2Nspectra.tsv', sep = '\t')
            sample_vec = [1]*len(pair_taxa)
            pair_metadata_df.loc[len(pair_metadata_df)] = [sample, phen1]
            for i, row in sample_df.iterrows():
                genome_id, n_spectra = row['genome'], int(round(row['n_spectra']))
                taxa_id = getMST(genome_id)
                sample_vec[taxa2idx_dic[taxa_id]] += n_spectra
            pair_df[sample] = sample_vec
            
        for sample in phen2_samples:
            sample_df = pd.read_csv(f'{genome_abundance_dir + sample}_genome2Nspectra.tsv', sep = '\t')
            sample_vec = [1]*len(pair_taxa)
            pair_metadata_df.loc[len(pair_metadata_df)] = [sample, phen2]
            for i, row in sample_df.iterrows():
                genome_id, n_spectra = row['genome'], int(round(row['n_spectra']))
                taxa_id = getMST(genome_id)
                sample_vec[taxa2idx_dic[taxa_id]] += n_spectra
            pair_df[sample] = sample_vec
            
        pair_df.to_csv(out_dir + phen1 + '_vs_' + phen2 + f'_{tax_level}_taxaProfiles.tsv', sep = '\t', index = False)
        pair_metadata_df.to_csv(out_dir + phen1 + '_vs_' + phen2 + '_metadata.tsv', sep = '\t', index = False)
