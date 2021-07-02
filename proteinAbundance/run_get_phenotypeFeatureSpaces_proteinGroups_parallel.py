"""
script to parellel run the scirpt get_phenotypeFeatureSpaces_proteinGroups.py across all pairs of samples comparisons specified by a list
"""
from multiprocessing import Pool
import sys



import os
in_f = sys.argv[1]#'/data/mstambou/proteome_landscapes/cluster_abundances_HAPiID_6k/comparisonLists.txt'

min_log = float(sys.argv[2])#2.5
top_n = int(sys.argv[3])#50
min_n_samples = int(sys.argv[4])#15
n_pool = int(sys.argv[5])#3

comparisons = []
with open(in_f, 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        phen1, phen2 = line[0], line[1]
        comparisons.append((phen1, phen2))
        
def run_script(comparison):
    phen1 = comparison[0]
    phen2 = comparison[1]
    print(f'processing {phen1} vs {phen2}')
    cmd = f'python3 get_phenotypeFeatureSpaces_proteinGroups.py proteinGroups_logOdds_phenotypes/ {phen1} {phen2} cluster_NSAFs/ phenotype_featureSpaces/ /data/mstambou/proteome_landscapes/highly_abundant_genomes/samplesSubSamples2phenotypes_dic.json {min_log} {top_n} {min_n_samples}'
    print(cmd)
    os.system(cmd)

with Pool(n_pool) as p:
    p.map(run_script, comparisons)
