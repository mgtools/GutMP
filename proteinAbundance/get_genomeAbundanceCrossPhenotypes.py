"""
script that will go over the abundant genome quantifications over each sample
and will quantify genome abundances from all the samples that belong to a certain phenotype
in parallel, the script will normalize genome abundances by number of samples, sample depth and proteome size
"""

import json
import pandas as pd
import os
from multiprocessing import Pool
import sys

samples2phenotypes_dic_f = sys.argv[1]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/samples2phenotypes_dic.json'

with open(samples2phenotypes_dic_f, 'r') as in_f:
    samples2phenotypes_dic = json.load(in_f)
    
summary_f = sys.argv[2]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/samples_1_genomes_1000_peps_withPhenotypes.tsv'
summary_df = pd.read_csv(summary_f, sep = '\t')

phenotypes = list(set(summary_df['diagnosis']))

sample_quantifications_dir = sys.argv[3]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/genome_spectralCount/'
out_dir = sys.argv[4]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/genome_spectralCount_phenotypes/'

n_threads = int(sys.argv[5])#20

in_suffix = '_genome2NspectraNormalized.tsv'
out_suffix = '_genome2NormalizedNspectra.tsv'


def getPhenotye2NormalizedGenomeAbundances(phenotype, samples2phenotypes_dic = samples2phenotypes_dic, sample_quantifications_dir = sample_quantifications_dir):
    """
    function that given a phenotype it will calculate normalized genome abundances for each phenotype
    It will normalize the values by sample depth, proteome sizes and number of samples
    """
    print('processing phenotype', phenotype, '...')
    phen_samples = samples2phenotypes_dic[phenotype]
    phen_df = pd.DataFrame(columns = ['genome', 'n_spectra', 'rel_n_spectra_%', 'n_spectra_normalized', 'rel_n_spectra_%_normalized', 'n_samples'])
    phen_dic = {}
    n_samples = len(phen_samples)
    genome2Nsamples_dic = {}
    for sample in phen_samples:
        sample_df = pd.read_csv(sample_quantifications_dir + sample + in_suffix, sep = '\t')
        for i, row in sample_df.iterrows():
            genome = row['genome']
            n_spectra = row['n_spectra']
            rel_n_spectra = row['rel_n_spectra_%']
            n_spectra_normalized = row['n_spectra_normalized']
            rel_n_spectra_normalized = row['rel_n_spectra_%_normalized']
            if genome not in phen_dic:
                phen_dic[genome] = {'n_spectra' : 0, 'rel_n_spectra_%' : 0, 'n_spectra_normalized' : 0, 'rel_n_spectra_%_normalized' : 0}
            phen_dic[genome]['n_spectra'] += n_spectra
            phen_dic[genome]['rel_n_spectra_%'] += rel_n_spectra
            phen_dic[genome]['n_spectra_normalized'] += n_spectra_normalized
            phen_dic[genome]['rel_n_spectra_%_normalized'] += rel_n_spectra_normalized
            if genome not in genome2Nsamples_dic:
                genome2Nsamples_dic[genome] = 0
            genome2Nsamples_dic[genome] += 1
            
    for i, genome in enumerate(phen_dic):
        phen_df.loc[i] = [genome, phen_dic[genome]['n_spectra'], phen_dic[genome]['rel_n_spectra_%']/n_samples, phen_dic[genome]['n_spectra_normalized'], phen_dic[genome]['rel_n_spectra_%_normalized']/n_samples, genome2Nsamples_dic[genome]]
    phen_df.to_csv(out_dir + phenotype + out_suffix, sep = '\t', index = False)
        
with Pool(n_threads) as p:
    p.map(getPhenotye2NormalizedGenomeAbundances, phenotypes)
