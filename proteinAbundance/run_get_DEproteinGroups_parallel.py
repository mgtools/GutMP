"""
script to parellel run the scirpt get_DEproteinGroups.py across all pairs of samples comparisons specified by a list
"""
from multiprocessing import Pool
n_pool = 20

import os
in_f = '/data/mstambou/proteome_landscapes/cluster_abundances_HAPiID_6k/comparisonLists.txt'

comparisons = []
with open(in_f, 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        phen1, phen2 = line[0], line[1]
        comparisons.append((phen1, phen2))
        
def run_DEproteinGroups(comparison):
    phen1 = comparison[0]
    phen2 = comparison[1]
    print(f'processing {phen1} vs {phen2}')
    cmd = f'python3 get_DEproteinGroups.py /data/mstambou/proteome_landscapes/cluster_abundances_HAPiID_6k/proteinGroup_NSAFs_phenotypes/ {phen1} {phen2} /data/mstambou/proteome_landscapes/cluster_abundances_HAPiID_6k/proteinGroups_logOdds_phenotypes/'
    print(cmd)
    os.system(cmd)

with Pool(n_pool) as p:
    p.map(run_DEproteinGroups, comparisons)
