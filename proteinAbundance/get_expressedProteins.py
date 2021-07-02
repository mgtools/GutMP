"""
script to get the exprssed proteins of the highly expressed genomes
accross all choses samples
"""

import pandas as pd
import json

samples_dir = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/samples_1_genomes_1000_peps/'
summary_f = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/samples_1_genomes_1000_peps.tsv'
protein2genome_dic_f = '/data/mstambou/proteome_landscapes/auxilary/proteins2genome_dic.json'

genome2expressedProteins_dic = {}

with open(protein2genome_dic_f, 'r') as in_f:
    protein2genome_dic = json.load(in_f)

summary_df = pd.read_csv(summary_f, sep = '\t')

row = summary_df.loc[0]

for i, row in summary_df.iterrows():
    
    sample, genomes = row['sample_name'], row['expanded_genomes'].split('|')
    print('processing sample', i, sample, '...')
    sample_df = pd.read_csv(samples_dir + sample + '_covering_80_percent_ribP_elonF_spectra.tsv.0.01.tsv', sep = '\t')

    for j, sample_row in sample_df.iterrows():
        protein_seqs = [item.split('(pre')[0] for item in sample_row['Protein'].split(';') if item.startswith('XXX_') == False]
        if protein_seqs:
            for protein in protein_seqs:
                genome = protein2genome_dic[protein]
                if genome not in genome2expressedProteins_dic:
                    genome2expressedProteins_dic[genome] = []
                genome2expressedProteins_dic[genome].append(protein)
            
            
genome2expressedProteins_dic = {genome:list(set(proteins)) for genome, proteins in genome2expressedProteins_dic.items()}
with open('genome2expressedProteins_dic.json', 'w') as out_f:
    json.dump(genome2expressedProteins_dic, out_f)
