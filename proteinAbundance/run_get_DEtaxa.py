"""                                                                                                                                                                                                        
script to parellel run the scirpt get_predictionsForGenomeFeatures.py across all pairs of samples comparisons specified by a list                                                                  
"""
from multiprocessing import Pool
import sys
import os

in_f = sys.argv[1]#'/data/mstambou/proteome_landscapes/cluster_abundances_HAPiID_6k/comparisonLists.txt'                                                                                                   
in_dir = sys.argv[2]

out_dir = sys.argv[3]

if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

n_pool = int(sys.argv[4])#3                                                                                                                                                                               

tax_level = sys.argv[5]

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
    cmd = f'Rscript get_DEtaxa.R {in_dir + phen1}_vs_{phen2}_{tax_level}_taxaProfiles.tsv {in_dir + phen1}_vs_{phen2}_metadata.tsv {out_dir}'
    print(cmd)
    os.system(cmd)

with Pool(n_pool) as p:
    p.map(run_script, comparisons)

