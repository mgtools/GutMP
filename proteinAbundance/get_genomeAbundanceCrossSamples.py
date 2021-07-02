"""
script to quantify genome abundance accorss the different samples so far
"""

import json
import pandas as pd
import os

allBin2taxon_dic_f = '/data/mstambou/proteome_landscapes/auxilary/allBin2taxon_dic.json'
spectral_counts_dir = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/genome_spectralCount/'
suffix = '_genome2NspectraNormalized.tsv'

genome2proteome_dic_f = '/data/mstambou/proteome_landscapes/auxilary/genome2protein_dic.json'
genome2expressedProteins_dic_f = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/genome2expressedProteins_dic.json'

out_dir = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/'

genome2totalPoints_dic, genome2sample2points_dic = {}, {}
genome2absoluteTotalPoints_dic = {}

genome2totalPointsNormalized_dic, genome2sample2pointsNormalized_dic = {}, {}
genome2absoluteTotalPointsNormalized_dic = {}

with open(allBin2taxon_dic_f, 'r') as in_f:
    allBin2taxon_dic = json.load(in_f)

for fname in [item for item in os.listdir(spectral_counts_dir) if item.endswith(suffix)]:
    df = pd.read_csv(spectral_counts_dir + fname, sep = '\t')
    sample = fname.replace(suffix, '')
    for i, row in df.iterrows():
        genome, rel_spec = row['genome'], row['rel_n_spectra_%']
        abs_spec = row['n_spectra']
        rel_spec_normalized = row['rel_n_spectra_%_normalized']
        abs_spec_normalized = row['n_spectra_normalized']        
        if genome not in genome2totalPoints_dic:
            genome2totalPoints_dic[genome] = 0
            genome2sample2points_dic[genome] = {}
            genome2absoluteTotalPoints_dic[genome] = 0
            genome2totalPointsNormalized_dic[genome] = 0
            genome2sample2pointsNormalized_dic[genome] = {}
            genome2absoluteTotalPointsNormalized_dic[genome] = 0
        genome2totalPoints_dic[genome] += rel_spec
        genome2sample2points_dic[genome][sample] = rel_spec
        genome2absoluteTotalPoints_dic[genome] += abs_spec
        genome2totalPointsNormalized_dic[genome] += rel_spec_normalized
        genome2sample2pointsNormalized_dic[genome][sample] = rel_spec_normalized
        genome2absoluteTotalPointsNormalized_dic[genome] += abs_spec_normalized
      
with open(out_dir + 'genome2totalPoints_dic.json', 'w') as out_f:
    json.dump(genome2totalPoints_dic, out_f)
    
with open(out_dir + 'genome2sample2points_dic.json', 'w') as out_f:
    json.dump(genome2sample2points_dic, out_f)
    
with open(out_dir + 'genome2absoluteTotalPoints_dic.json', 'w') as out_f:
    json.dump(genome2absoluteTotalPoints_dic, out_f)
        
with open(out_dir + 'genome2totalPointsNormalized_dic.json', 'w') as out_f:
    json.dump(genome2totalPointsNormalized_dic, out_f)
    
with open(out_dir + 'genome2sample2pointsNormalized_dic.json', 'w') as out_f:
    json.dump(genome2sample2pointsNormalized_dic, out_f)
    
with open(out_dir + 'genome2absoluteTotalPointsNormalized_dic.json', 'w') as out_f:
    json.dump(genome2absoluteTotalPointsNormalized_dic, out_f)

def getMST(genome, allBin2taxon_dic = allBin2taxon_dic):
    for c in 'sgfocpd':
        taxon = allBin2taxon_dic[genome][c]
        if taxon != '' and taxon != 'GCF_' and taxon != 'GCA_':
            return c+'__'+taxon
genome2totalPoints_df = pd.DataFrame(columns = ['genome_ID', 'taxonomy', 'total_points', 'absolute_n_spectra', 'total_points_normalized', 'absolute_n_spectra_normalized'])

print('calculating total spectra per genome ...')
for i, genome in enumerate(genome2totalPoints_dic):
    totalPoints = genome2totalPoints_dic[genome]
    taxon = getMST(genome)
    n_spectra = genome2absoluteTotalPoints_dic[genome]
    
    totalPointsNormalized = genome2totalPointsNormalized_dic[genome]    
    n_spectraNormalized = genome2absoluteTotalPointsNormalized_dic[genome]
    genome2totalPoints_df.loc[i] = [genome, taxon, totalPoints, n_spectra, totalPointsNormalized, n_spectraNormalized]
    
print('retrieving total number of proteomes and expressed proteomes for each genome ...')    

with open(genome2proteome_dic_f, 'r') as in_f:
    genome2proteome_dic = json.load(in_f)
with open(genome2expressedProteins_dic_f, 'r') as in_f:
    genome2expressedProteins_dic = json.load(in_f)
    

genome2Nproteome_dic = {k:len(set(v)) for k,v in genome2proteome_dic.items()}
genome2NexpressedProteins_dic = {k:len(set(v)) for k, v in genome2expressedProteins_dic.items()}

print('calculating the proportions of expressed genomes as a percentage ...')
genome2totalPoints_df['%_expressed_proteomes'] = [(genome2NexpressedProteins_dic[item] / genome2Nproteome_dic[item])*100 for item in list(genome2totalPoints_df['genome_ID'])]

    
genome2totalPoints_df.sort_values(by = 'total_points', ascending = False, inplace = True)

genome2totalPoints_df.to_csv('/data/mstambou/proteome_landscapes/highly_abundant_genomes/genome2totalPoints.tsv', sep = '\t', index = False)
