"""
script that will go over the abundant protein group quantifications over each sample
and will quantify protein group abundances from all the samples that belong to a certain phenotype
"""

import json
import pandas as pd
import os
import time
from multiprocessing import Pool
import sys

n_pool = 20

#this is basically the same dictionary mapping so you can keep this here the same
samples2phenotypes_dic_f = sys.argv[1]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/samples2phenotypes_dic.json'
proteinGroup2Annotation_dic_f = sys.argv[2]#'/data/mstambou/proteome_landscapes/cluster_abundances_HAPiID_6k/proteinGroup2Annotation_dic.json'
proteinGroup2COG_dic_f = sys.argv[3]#'/data/mstambou/proteome_landscapes/cluster_abundances_HAPiID_6k/proteinGroup2COGs_dic.json'

summary_f = sys.argv[4]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/samples_1_genomes_1000_peps_withPhenotypes.tsv'
sample_quantifications_dir = sys.argv[5]#'/data/mstambou/proteome_landscapes/cluster_abundances_HAPiID_6k/cluster_NSAFs/'
out_dir = sys.argv[6]#'/data/mstambou/proteome_landscapes/cluster_abundances_HAPiID_6k/proteinGroup_NSAFs_phenotypes/'

with open(proteinGroup2COG_dic_f, 'r') as in_f:
    proteinGroup2COG_dic = json.load(in_f)

with open(samples2phenotypes_dic_f, 'r') as in_f:
    samples2phenotypes_dic = json.load(in_f)    

with open(proteinGroup2Annotation_dic_f, 'r') as in_f:
    proteinGroup2Annotation_dic = json.load(in_f)

summary_df = pd.read_csv(summary_f, sep = '\t')

phenotypes = list(set(summary_df['diagnosis']))


if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

in_suffix = '_covering_80_percent_ribP_elonF_spectra_NSAF.tsv'
out_suffix = '_proteinGroup2NSAF.tsv'
out_suffix_2 = '_PSSM2NSAF.tsv'
out_suffix_3 = '_COG2NSAF.tsv'

def get_phenotype2NSAFs(phenotype, proteinGroup2Annotation_dic = proteinGroup2Annotation_dic, proteinGroup2COG_dic = proteinGroup2COG_dic):
    time_start = time.time()
    print('processing phenotype', phenotype)
    phen_samples = samples2phenotypes_dic[phenotype]
    phen_df = pd.DataFrame(columns = ['clusterID', 'NSAF', 'NSAF_%', 'n_samples'])
    PSSM_phen_df = pd.DataFrame(columns = ['PSSM_ID', 'NSAF', 'NSAF_%', 'n_samples'])
    COG_phen_df = pd.DataFrame(columns = ['COG', 'NSAF', 'NSAF_%', 'n_samples'])
    phen_dic = {}
    PSSM_phen_dic = {}
    COG_phen_dic = {}
    cluster2Nsamples_dic = {}
    PSSM2Nsamples_dic = {}
    COG2Nsamples_dic = {}
    sample2PSSM_dic = {}
    sample2COG_dic = {}
    n_samples = len(phen_samples)
    for sample in phen_samples:
        sample2PSSM_dic[sample] = set()
        sample2COG_dic[sample] = set()
        sample_df = pd.read_csv(sample_quantifications_dir + sample + in_suffix, sep = '\t')
        for i, row in sample_df.iterrows():
            clusterID = row['clusterID']
            if clusterID not in cluster2Nsamples_dic:
                cluster2Nsamples_dic[clusterID] = 0            
            cluster2Nsamples_dic[clusterID] += 1
            
            annotation = proteinGroup2Annotation_dic[clusterID]
            COG_annotation = proteinGroup2COG_dic[clusterID]
            if annotation != 'NA':
                PSSMs = [item for sublist in proteinGroup2Annotation_dic[clusterID]["PSSM_IDs"] for item in sublist.split('|')]
            else:
                PSSMs = [clusterID]
            
            COGs = proteinGroup2COG_dic[clusterID]
            if len(COGs) == 1 and COGs[0] == 'NA':
                COGs = [clusterID]
            
            NSAF = row['NSAF']
            rel_NSAF = row['NSAF_%']
            if clusterID not in phen_dic:
                phen_dic[clusterID] = {'NSAF' : 0, 'NSAF_%' : 0}                
                
            phen_dic[clusterID]['NSAF'] += NSAF
            phen_dic[clusterID]['NSAF_%'] += rel_NSAF
            
            for PSSM_ID in PSSMs:
                if PSSM_ID not in PSSM2Nsamples_dic:
                    PSSM2Nsamples_dic[PSSM_ID] = 0
                if PSSM_ID not in sample2PSSM_dic[sample]:
                    PSSM2Nsamples_dic[PSSM_ID] += 1
                    sample2PSSM_dic[sample].add(PSSM_ID)
                
                if PSSM_ID not in PSSM_phen_dic:
                    PSSM_phen_dic[PSSM_ID] = {'NSAF' : 0, 'NSAF_%' : 0}
                PSSM_phen_dic[PSSM_ID]['NSAF'] += NSAF
                PSSM_phen_dic[PSSM_ID]['NSAF_%'] += rel_NSAF            
            
            for COG in COGs:
                if COG not in COG2Nsamples_dic:
                    COG2Nsamples_dic[COG] = 0
                if COG not in sample2COG_dic[sample]:
                    COG2Nsamples_dic[COG] += 1
                    sample2COG_dic[sample].add(COG)
                
                if COG not in COG_phen_dic:
                    COG_phen_dic[COG] = {'NSAF' : 0, 'NSAF_%' : 0}
                COG_phen_dic[COG]['NSAF'] += NSAF
                COG_phen_dic[COG]['NSAF_%'] += rel_NSAF            
            
    for i, clusterID in enumerate(phen_dic):
        phen_df.loc[i] = [clusterID, phen_dic[clusterID]['NSAF'], phen_dic[clusterID]['NSAF_%']/n_samples, cluster2Nsamples_dic[clusterID]]
    phen_df.to_csv(out_dir + phenotype + out_suffix, sep = '\t', index = False)
    
    for i, PSSM_ID in enumerate(PSSM_phen_dic):
        PSSM_phen_df.loc[i] = [PSSM_ID, PSSM_phen_dic[PSSM_ID]['NSAF'], PSSM_phen_dic[PSSM_ID]['NSAF_%']/n_samples, PSSM2Nsamples_dic[PSSM_ID]]
    PSSM_phen_df.to_csv(out_dir + phenotype + out_suffix_2, sep = '\t', index = False)
    
    for i, COG in enumerate(COG_phen_dic):
        COG_phen_df.loc[i] = [COG, COG_phen_dic[COG]['NSAF'], COG_phen_dic[COG]['NSAF_%']/n_samples, COG2Nsamples_dic[COG]]
    COG_phen_df.to_csv(out_dir + phenotype + out_suffix_3, sep = '\t', index = False)
    print('finished processing phenotype', phenotype, 'in', time.time() - time_start, 'seconds ...')

        
with Pool(n_pool) as p:
    p.map(get_phenotype2NSAFs, phenotypes)
