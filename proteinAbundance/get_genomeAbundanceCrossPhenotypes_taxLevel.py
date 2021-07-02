"""
script that will go over the abundant genome quantifications over each sample
and will quantify taxonomic abundances (at a specified level) from all the samples that belong to a certain phenotype
this is a generic script that will process all phenotypes
"""

import json
import pandas as pd
import os
import sys

samples2phenotypes_dic_f = sys.argv[1]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/samples2phenotypes_dic.json'

summary_f = sys.argv[2]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/samples_1_genomes_1000_peps_withPhenotypes.tsv'

sample_quantifications_dir = sys.argv[3]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/genome_spectralCount/'
out_dir = sys.argv[4]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/genome_spectralCount_phenotypes/'

tax_level = sys.argv[5]#'s'

with open(samples2phenotypes_dic_f, 'r') as in_f:
    samples2phenotypes_dic = json.load(in_f)
    
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

summary_df = pd.read_csv(summary_f, sep = '\t')


phenotypes = list(set(summary_df['diagnosis']))


if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
out_dir = out_dir + tax_level + '_taxa/'
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

in_suffix = '_genome2NspectraNormalized.tsv'
out_suffix = 'genome2NormalizedNspectra.tsv'

for phenotype in phenotypes:
    print('processing phenotye', phenotype, '...')
    phen_samples = samples2phenotypes_dic[phenotype]
    phen_df = pd.DataFrame(columns = ['taxa', 'n_spectra', 'rel_n_spectra_%', 'n_spectra_normalized', 'rel_n_spectra_%_normalized', 'n_samples'])
    phen_dic = {}
    n_samples = len(phen_samples)
    taxa2Nsamples_dic = {}
    for sample in phen_samples:
        sample_df = pd.read_csv(sample_quantifications_dir + sample + in_suffix, sep = '\t')
        sample_taxon = set()
        for i, row in sample_df.iterrows():
            genome = row['genome']
            taxon = getMST(genome)
            n_spectra = row['n_spectra']
            rel_n_spectra = row['rel_n_spectra_%']
            n_spectra_normalized = row['n_spectra_normalized']
            rel_n_spectra_normalized = row['rel_n_spectra_%_normalized']
            if taxon not in taxa2Nsamples_dic:
                taxa2Nsamples_dic[taxon] = 0
            if taxon not in sample_taxon:                
                taxa2Nsamples_dic[taxon] += 1
                sample_taxon.add(taxon)
            if taxon not in phen_dic:
                phen_dic[taxon] = {'n_spectra' : 0, 'rel_n_spectra_%' : 0, 'n_spectra_normalized' : 0, 'rel_n_spectra_%_normalized' : 0}
            phen_dic[taxon]['n_spectra'] += n_spectra
            phen_dic[taxon]['rel_n_spectra_%'] += rel_n_spectra
            phen_dic[taxon]['n_spectra_normalized'] += n_spectra_normalized
            phen_dic[taxon]['rel_n_spectra_%_normalized'] += rel_n_spectra_normalized
    for i, taxon in enumerate(phen_dic):
        phen_df.loc[i] = [taxon, phen_dic[taxon]['n_spectra'], phen_dic[taxon]['rel_n_spectra_%']/n_samples, phen_dic[taxon]['n_spectra_normalized'], phen_dic[taxon]['rel_n_spectra_%_normalized']/n_samples, taxa2Nsamples_dic[taxon]]
    phen_df.to_csv(out_dir + phenotype + '_' +tax_level + '_' + out_suffix, sep = '\t', index = False)
        
