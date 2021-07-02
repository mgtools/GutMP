"""                                                                                                                                
script that will quantify genomes per sample                                                                                       
"""

import pandas as pd
import json
import copy
import time
import os
import sys
from multiprocessing import Pool

time_start = time.time()

summary_f = sys.argv[1]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/samples_1_genomes_1000_peps.tsv'               
identifications_dir = sys.argv[2]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/samples_1_genomes_1000_peps/'        
summary_df = pd.read_csv(summary_f, sep = '\t')
proteins2genomes_dic_f = sys.argv[3]#'/data/mstambou/proteome_landscapes/HAPiID_6k/data/proteins2genome_dic.json'                  

genome2ProteomeSize_dic_f = sys.argv[4] #/data/mstambou/proteome_landscapes/auxilary/genome2ProteomeSize_dic.json

out_dir = sys.argv[5]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/genome_spectralCount/'                           

n_threads = int(sys.argv[6])#4


if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

with open(proteins2genomes_dic_f, 'r') as in_f:
    proteins2genomes_dic = json.load(in_f)
    
with open(genome2ProteomeSize_dic_f, 'r') as in_f:
    genome2ProteomeSize_dic = json.load(in_f)
    

genome2proteomeSizePM_dic = {k:(v/10**6) for k,v in genome2ProteomeSize_dic.items()}

def getSpectra2genome(sample_df, expanded_genomes):
    """
    Function that goes over the FDR filtered PSM identifications 
    and will return spectra to genome maps that will be used for 
    identifying uniquely and multi-mapped spectra to genomes
    """
    spectra2protein_dic, spectra2genome_dic = {}, {}
    for j, sample_row in sample_df.iterrows():
        proteins = [item.split('(pre')[0] for item in sample_row['Protein'].split(';') if item.startswith('XXX_') == False]
        
        if proteins:
            spectrum = sample_row['SpecID']
            spectra2protein_dic[spectrum] = proteins
            spectra2genome_dic[spectrum] = []
            for protein in proteins:
                genome = proteins2genomes_dic[protein]
                if genome in expanded_genomes:
                    if genome not in spectra2genome_dic[spectrum]: #make sure you only mention each genome once per spectrum at most
                        spectra2genome_dic[spectrum].append(genome)
    return spectra2genome_dic

def getGenome2uniqueSpecs(spectra2genome_dic, expanded_genomes):
    """
    Function that will return genome to unique spectra maps 
    """
    genome2uniqueSpec_dic = {genome:0 for genome in expanded_genomes}
    multiMappedSpecs2genomes_dic = {}
    for spec in spectra2genome_dic:
        genomes = spectra2genome_dic[spec]
        if len(genomes) == 1:
            genome2uniqueSpec_dic[genomes[0]] += 1
        else:
            multiMappedSpecs2genomes_dic[spec] = genomes    
    return genome2uniqueSpec_dic, multiMappedSpecs2genomes_dic

def getGenomeCoefficients(genomes, genome2uniqueSpec_dic):
    total_specs = sum([genome2uniqueSpec_dic[genome] for genome in genomes])
    return {genome:genome2uniqueSpec_dic[genome]/total_specs for genome in genomes}

def getGenome2Nspectra(row, genome2proteomeSizePM_dic = genome2proteomeSizePM_dic):

    #sample = row['sample_name']                                                                                                   
    sample = row[0]
    print('processing sample', sample, '...')

    #expanded_genomes = set(row['expanded_genomes'].split('|'))                                                                    
    expanded_genomes = set(row[-1].split('|'))    
    
    sample_df = pd.read_csv(identifications_dir + sample + '_covering_80_percent_ribP_elonF_spectra.tsv.0.01.tsv', sep = '\t')

    spectra2genome_dic = getSpectra2genome(sample_df, expanded_genomes)
    genome2uniqueSpec_dic, multiMappedSpecs2genomes_dic = getGenome2uniqueSpecs(spectra2genome_dic, expanded_genomes)

    #just in case I have genomes with no unique spectral support (should not be considerede further)                               
    genome2uniqueSpec_dic = {k:v for k,v in genome2uniqueSpec_dic.items() if v != 0}
    multiMappedSpecs2genomes_dic = {k:list(set(v).intersection(set(genome2uniqueSpec_dic.keys()))) for k,v in multiMappedSpecs2genomes_dic.items()}

    genome2Nspecs_dic = copy.deepcopy(genome2uniqueSpec_dic)    

    #fraction assignments for the multimapped reads shared by 2 or more genomes within a sample
    #I use the coeficients based on unique spectra mapped between the competing genomes
    for spec in multiMappedSpecs2genomes_dic:
        genomes = multiMappedSpecs2genomes_dic[spec]
        genome2coefficients_dic = getGenomeCoefficients(genomes, genome2uniqueSpec_dic)        
        for genome in genomes:
            genome2Nspecs_dic[genome] += genome2coefficients_dic[genome]
    
    genome2NormalizedNspecs_dic = {k:(v/genome2proteomeSizePM_dic[k] ) for k,v in genome2Nspecs_dic.items()}
    
    #print(genome2Nspecs_dic)
    #print(genome2NormalizedNspecs_dic)
    df = pd.DataFrame(columns=['genome', 'n_spectra', 'rel_n_spectra_%', 'n_spectra_normalized', 'rel_n_spectra_%_normalized'])

    #total_spectra = len(spectra2genome_dic)
    total_spectra = sum(genome2Nspecs_dic.values())
    total_normalized_spectra = sum(genome2NormalizedNspecs_dic.values())    
    
    print('total:',total_spectra, 'total Normalized:',total_normalized_spectra)

    for i, genome in enumerate(genome2Nspecs_dic):
        n_specs = round(genome2Nspecs_dic[genome], 3)
        rel_n_specs = round((n_specs/total_spectra)*100, 3)
        n_specs_normalized = round(genome2NormalizedNspecs_dic[genome], 3)
        rel_n_specs_normalized = round((n_specs_normalized/total_spectra)*100, 3)
        df.loc[i] = [genome, n_specs, rel_n_specs, n_specs_normalized, rel_n_specs_normalized]
    df.sort_values(by = 'n_spectra', ascending = False, inplace = True)
    df.to_csv(out_dir + sample + '_genome2NspectraNormalized.tsv', sep = '\t', index = False)

rows = [list(row) for row in summary_df.values]
#rows = [rows[0]]
#getGenome2Nspectra(row)                                                                                                           

with Pool(n_threads) as p:
    p.map(getGenome2Nspectra, rows)

time_end = time.time()
print(f'processed {len(rows)} in {time_end - time_start} seconds ...')
