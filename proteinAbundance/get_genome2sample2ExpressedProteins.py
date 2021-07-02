"""
script to study proportion of proteomes covered as as function of expressed samples
"""

import os
import pandas as pd
import json

samples_dir = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/samples_1_genomes_1000_peps/'
summary_f = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/samples_1_genomes_1000_peps.tsv'
protein2genome_dic_f = '/data/mstambou/proteome_landscapes/auxilary/proteins2genome_dic.json'

with open(protein2genome_dic_f, 'r') as in_f:
    proteins2genomes_dic = json.load(in_f)
    
summary_df = pd.read_csv(summary_f, sep = '\t')
genome2sample2ExpressedProteins_dic = dict()
genome2sample2Proteins2Nspecs_dic = dict()

for i, row in summary_df.iterrows():

    sample = row['sample_name']
    expanded_genomes = row['expanded_genomes'].split('|')
    
    for genome in expanded_genomes:
        if genome not in genome2sample2ExpressedProteins_dic:
            genome2sample2ExpressedProteins_dic[genome] = {}
            genome2sample2Proteins2Nspecs_dic[genome] = {}
        if sample not in genome2sample2ExpressedProteins_dic[genome]:
            genome2sample2ExpressedProteins_dic[genome][sample] = []
            genome2sample2Proteins2Nspecs_dic[genome][sample] = {}
    print('processing sample', i, sample)

    sample_df = pd.read_csv(samples_dir + sample + '_covering_80_percent_ribP_elonF_spectra.tsv.0.01.tsv', sep = '\t')

    for j, sample_row in sample_df.iterrows():
        proteins = [item.split('(pre')[0] for item in sample_row['Protein'].split(';') if item.startswith('XXX_') == False]            
        if proteins:        
            for protein in proteins:
                genome = proteins2genomes_dic[protein]
                if genome in expanded_genomes:
                    if genome not in genome2sample2ExpressedProteins_dic:
                        genome2sample2ExpressedProteins_dic[genome] = {}
                        genome2sample2Proteins2Nspecs_dic[genome] = {}
                    if sample not in genome2sample2ExpressedProteins_dic[genome]:
                        genome2sample2ExpressedProteins_dic[genome][sample] = []
                        genome2sample2Proteins2Nspecs_dic[genome][sample] = {}
                    genome2sample2ExpressedProteins_dic[genome][sample].append(protein)
                    if protein not in genome2sample2Proteins2Nspecs_dic[genome][sample]:
                        genome2sample2Proteins2Nspecs_dic[genome][sample][protein] = 1
                    elif protein in genome2sample2Proteins2Nspecs_dic[genome][sample]:
                        genome2sample2Proteins2Nspecs_dic[genome][sample][protein] += 1
                        
    for genome in expanded_genomes:
        genome2sample2ExpressedProteins_dic[genome][sample] = list(set(genome2sample2ExpressedProteins_dic[genome][sample]))

with open('genome2sample2ExpressedProteins_dic.json', 'w') as out_f:
    json.dump(genome2sample2ExpressedProteins_dic, out_f)
    
with open('genome2sample2Proteins2Nspecs_dic.json', 'w') as out_f:
    json.dump(genome2sample2Proteins2Nspecs_dic, out_f)
