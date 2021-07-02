"""                                                                                                                                                                                                       
script to parellel run the scirpt get_predictionsForGenomeFeatures.py across all pairs of samples comparisons specified by a list                                                                  
"""
from multiprocessing import Pool
import sys
import os

in_f = sys.argv[1]#'/data/mstambou/proteome_landscapes/cluster_abundances_HAPiID_6k/comparisonLists.txt'                                                                                                   

min_spec = int(sys.argv[2])#10                                                                                                                            
max_pid = float(sys.argv[3])#70                                                                                                                                                                                

comparisons = []
with open(in_f, 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        phen1, phen2 = line[0], line[1]
        print(f'currently processing {phen1} vs {phen2} ... ')
        cmd = f'python3 get_strainSpecificHighlyExpressedProteins.py /data/mstambou/proteome_landscapes/highly_abundant_genomes/genome_logOdds_phenotypes/ {phen1} {phen2} {min_spec} {max_pid}'
        os.system(cmd)

