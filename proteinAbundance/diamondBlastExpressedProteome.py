
"""
script where I will extract for a genome, all of it's proteins that have spectra support in a list of samples
where that genome is expressed for a particular phenotype, diamond blast these proteins against the COG
database, extract the COGs for the protein hits and then maintain a list of COG to spectra and a list of 
proteins with no COG hits and their respective spectra as well
"""

import json
from Bio import SeqIO
import os
import copy
import sys
import pandas as pd

genome_id = sys.argv[1]#'12718_7_61'
phenotype = sys.argv[2]#'healthy'
out_dir = sys.argv[3]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/genome_logOdds_phenotypes/healthy_vs_CD/'

genome2sample2Proteins2Nspecs_dic_f = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/genome2sample2Proteins2Nspecs_dic.json'
samples2phenotypes_dic_f = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/samples2phenotypes_dic.json'
samples2phenotypesSubsampled_dic_f = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/samples2phenotypesSubsampled_dic.json'
proteomes_dir = '/data/mstambou/proteome_landscapes/HAPiID_6k/data/proteomes/'
seqID2COGs_dic_f = '/data/mstambou/proteome_landscapes/COG_db/seqID2COGs_dic.json'

with open(seqID2COGs_dic_f, 'r') as in_f:
    seqID2COGs_dic = json.load(in_f)

with open(genome2sample2Proteins2Nspecs_dic_f, 'r') as in_f:
    genome2sample2Proteins2Nspecs_dic = json.load(in_f)

with open(samples2phenotypes_dic_f, 'r') as in_f:
    samples2phenotypes_dic = json.load(in_f)
    
with open(samples2phenotypesSubsampled_dic_f, 'r') as in_f:
    samples2phenotypesSubsampled_dic = json.load(in_f)
    
allSamples2phenotypes_dic = copy.deepcopy(samples2phenotypes_dic)
allSamples2phenotypes_dic.update(samples2phenotypesSubsampled_dic)
    
out_dir = f'{out_dir + phenotype}/'
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

expressed_samples = set(genome2sample2Proteins2Nspecs_dic[genome_id].keys()).intersection(set(allSamples2phenotypes_dic[phenotype]))

protein2spectra_dic = {}

for sample in expressed_samples:
    proteins = genome2sample2Proteins2Nspecs_dic[genome_id][sample]
    for protein in proteins:
        if protein not in protein2spectra_dic:
            protein2spectra_dic[protein] = genome2sample2Proteins2Nspecs_dic[genome_id][sample][protein]
        else:
            protein2spectra_dic[protein] += genome2sample2Proteins2Nspecs_dic[genome_id][sample][protein]
            
proteome = SeqIO.to_dict(SeqIO.parse(proteomes_dir + genome_id + '.fna.FGS.faa', 'fasta'))

with open(out_dir + genome_id+'_protein2spectra_dic.json', 'w') as out_f:
    json.dump(protein2spectra_dic, out_f)
    
with open(out_dir + genome_id+'_expressedProteinSeqs.faa', 'w') as out_f:
    for pid in protein2spectra_dic:
        seq = str(proteome[pid].seq)
        out_f.write('>'+pid+'\n')
        out_f.write(seq+'\n')

cmd = f'diamond blastp --more-sensitive --query {out_dir}{genome_id}_expressedProteinSeqs.faa --db /data/mstambou/proteome_landscapes/COG_db/cog_diamond_db/cog-20_COGsMerged.dmnd --out {out_dir}{genome_id}_diamondBlast.txt --threads 60 --outfmt 6'
os.system(cmd)

expressed_proteins2spectra_dic_f = out_dir + genome_id+'_protein2spectra_dic.json'
diamond_out_f = f'{out_dir}{genome_id}_diamondBlast.txt'

with open(expressed_proteins2spectra_dic_f, 'r') as in_f:
    expressed_proteins2spectra_dic = json.load(in_f)
    
diamond_out_df = pd.read_csv(diamond_out_f, sep = '\t')
diamond_out_df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
diamond_out_best_df = diamond_out_df.drop_duplicates(subset = 'qseqid', keep = 'first')

def check_overlap(sstart, send, cog_start, cog_end, rate = 0.5):
    cog_len = cog_end - cog_start + 1
    if (sstart <= cog_start) and (send > cog_start):
        if (send - cog_start) >= (cog_len*rate):
            return True
    elif (sstart >= cog_start) and (send <= cog_end):
        if (send - sstart + 1) >= (cog_len*rate):
            return True
    elif (sstart <= cog_end) and (send >= cog_end):
        if (cog_end - sstart) >= (cog_len*rate):
            return True
    elif (sstart <= cog_start) and (send >= cog_end):
        return True
    return False

row = diamond_out_best_df.iloc[1]
proteome2COG_dic = {}

for i, row in diamond_out_best_df.iterrows():    
    q_seq = row['qseqid']
    ref_seq = row['sseqid']
    sstart, send = row['sstart'], row['send']
    cogs = seqID2COGs_dic[ref_seq]    
    if q_seq not in proteome2COG_dic:
        proteome2COG_dic[q_seq] = []   
    
    for i, cog in enumerate(cogs['cog']):
        coords = cogs['coord'][i]
        if '=' in coords:
            sub_coords = coords.split('=')
            for sub_coord in sub_coords:
                cog_start, cog_end = int(sub_coord.split('-')[0]), int(sub_coord.split('-')[1])
                overlap = check_overlap(sstart, send, cog_start, cog_end)
                if overlap:                    
                    if cog not in proteome2COG_dic[q_seq]:
                        proteome2COG_dic[q_seq].append(cog)                
        else:
            cog_start, cog_end = int(coords.split('-')[0]), int(coords.split('-')[1])
            overlap = check_overlap(sstart, send, cog_start, cog_end)
            if overlap:    
                if cog not in proteome2COG_dic[q_seq]:
                    proteome2COG_dic[q_seq].append(cog)

unknonwFunctions_dic = {}
cogs2spectra_dic = {}
for protein in expressed_proteins2spectra_dic:
    if protein not in proteome2COG_dic:
        unknonwFunctions_dic[protein] = expressed_proteins2spectra_dic[protein]
    elif protein in proteome2COG_dic:
        if proteome2COG_dic[protein] == []:
            unknonwFunctions_dic[protein] = expressed_proteins2spectra_dic[protein]
        else:
            cogs = proteome2COG_dic[protein]
            n_spectra = expressed_proteins2spectra_dic[protein]
            for cog in cogs:
                if cog not in cogs2spectra_dic:
                    cogs2spectra_dic[cog] = n_spectra
                else:
                    cogs2spectra_dic[cog] += n_spectra
                    
                    
with open(out_dir + genome_id+'_cogs2pectra_dic.json', 'w') as out_f:
    json.dump(cogs2spectra_dic, out_f)
    
with open(out_dir + genome_id+'_proteinsUnknonwCogs_dic.json', 'w') as out_f:
    json.dump(unknonwFunctions_dic, out_f)

with open(out_dir + genome_id+'_proteome2COG_dic.json', 'w') as out_f:
    json.dump(proteome2COG_dic, out_f)
