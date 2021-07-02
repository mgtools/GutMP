#!/usr/bin/python3
"""
script that will get the cogs with spectral support, expressed protein sequences and expressed proteins with spectral support
and will blast these proteins against the proteomes of the other strain and vice versa, and will map the COGs to the respective
proteins if any found and will also map the spectral support as another column to this.
"""
import json
import pandas as pd
import os
from Bio import SeqIO
import copy 
import sys


in_dir = sys.argv[1]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/genome_logOdds_phenotypes/'
phen_1 = sys.argv[2]#'healthy'
phen_2 = sys.argv[3]#'CD'

genome1 = sys.argv[4]#'12718_7_61'
genome2 = sys.argv[5]#'14207_7_24'

min_spec = 10
max_pid = 70

in_dir = f'{in_dir}{phen_1}_vs_{phen_2}/'

out_dir = in_dir + f'{genome1}_vs_{genome2}/'
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

allBin2taxa_dic_f = '/data/mstambou/proteome_landscapes/auxilary/allBin2taxon_dic.json'
cog2cat_dic_f = '/data/mstambou/proteome_landscapes/COG_db/cog2cat_dic.json'
proteomes_dir = '/data/mstambou/proteome_landscapes/HAPiID_6k/data/proteomes/'
with open(cog2cat_dic_f, 'r') as in_f:
    cog2cat_dic = json.load(in_f)
    
with open(allBin2taxa_dic_f, 'r') as in_f:
    allBin2taxa_dic = json.load(in_f)
    
def get_MRCA(genome1, genome2, allBin2taxa_dic = allBin2taxa_dic):
    for t in 'sgfocpk':
        if allBin2taxa_dic[genome1][t] == allBin2taxa_dic[genome2][t]:
            return f'{t}__{allBin2taxa_dic[genome1][t]}'
        
   
phen1_cogs_f = f'{in_dir}{phen_1}/{genome1}_cogs2pectra_dic.json'
phen1_unknown_f = f'{in_dir}{phen_1}/{genome1}_proteinsUnknonwCogs_dic.json'
phen1_expressed_protSeqs_f = f'{in_dir}{phen_1}/{genome1}_expressedProteinSeqs.faa'
phen1_prot2spectra_f = f'{in_dir}{phen_1}/{genome1}_protein2spectra_dic.json'
phen1_prot2cog_f = f'{in_dir}{phen_1}/{genome1}_proteome2COG_dic.json'

phen2_cogs_f = f'{in_dir}{phen_2}/{genome2}_cogs2pectra_dic.json'
phen2_unknown_f = f'{in_dir}{phen_2}/{genome2}_proteinsUnknonwCogs_dic.json'
phen2_expressed_protSeqs_f = f'{in_dir}{phen_2}/{genome2}_expressedProteinSeqs.faa'
phen2_prot2spectra_f = f'{in_dir}{phen_2}/{genome2}_protein2spectra_dic.json'
phen2_prot2cog_f = f'{in_dir}{phen_2}/{genome2}_proteome2COG_dic.json'

with open(phen1_cogs_f, 'r') as in_f:
    phen1_cogs_dic = json.load(in_f)    
with open(phen1_unknown_f, 'r') as in_f:
    phen1_unknown_dic = json.load(in_f)
with open(phen1_prot2spectra_f, 'r') as in_f:
    phen1_prot2spectra_dic = json.load(in_f)
with open(phen1_prot2cog_f, 'r') as in_f:
    phen1_prot2cog_dic = json.load(in_f)
    
with open(phen2_cogs_f, 'r') as in_f:
    phen2_cogs_dic = json.load(in_f)
with open(phen2_unknown_f, 'r') as in_f:
    phen2_unknown_dic = json.load(in_f)
with open(phen2_prot2spectra_f, 'r') as in_f:
    phen2_prot2spectra_dic = json.load(in_f)
with open(phen2_prot2cog_f, 'r') as in_f:
    phen2_prot2cog_dic = json.load(in_f)
    
def get_cat2spec(phen_cogs_dic, phen_unknown_dic, cog2cat_dic = cog2cat_dic):
    phen_cat2spec_dic = {}

    for cog in phen_cogs_dic:
        cat = cog2cat_dic[cog]
        if len(cat) == 1:
            if cat not in phen_cat2spec_dic:
                phen_cat2spec_dic[cat] = phen_cogs_dic[cog]
            else:
                phen_cat2spec_dic[cat] += phen_cogs_dic[cog]
        else:
            for c in cat:
                if c not in phen_cat2spec_dic:
                    phen_cat2spec_dic[c] = phen_cogs_dic[cog]
                else:
                    phen_cat2spec_dic[c] += phen_cogs_dic[cog]

    for prot in phen_unknown_dic:
        phen_cat2spec_dic['S'] += phen_unknown_dic[prot]

    return phen_cat2spec_dic

phen1_cat2spec_dic = get_cat2spec(phen1_cogs_dic, phen1_unknown_dic)
phen1_total_spec = sum(phen1_cat2spec_dic.values())
phen1_normalizedCat2spec_dic = {k:(v*100)/phen1_total_spec for k,v in phen1_cat2spec_dic.items()}
phen2_cat2spec_dic = get_cat2spec(phen2_cogs_dic, phen2_unknown_dic)
phen2_total_spec = sum(phen2_cat2spec_dic.values())
phen2_normalizedCat2spec_dic = {k:(v*100)/phen2_total_spec for k,v in phen2_cat2spec_dic.items()}

phen2cats_df = pd.DataFrame(columns = ['COG_category', 'rel_spectra', 'phenotype'])

phen2cogs_df = pd.DataFrame(columns = ['COG', f'{phen_1}_spectra', f'{phen_2}_spectra'])

phen1_totalSpec, phen2_totalSpec = 0, 0
for i, cog in enumerate(set(phen1_cogs_dic.keys()).union(set(phen2_cogs_dic.keys()))):
    phen1_spec, phen2_spec = 0, 0
    if cog in phen1_cogs_dic:
        phen1_spec = phen1_cogs_dic[cog]
        phen1_totalSpec += phen1_cogs_dic[cog]
    if cog in phen2_cogs_dic:
        phen2_spec = phen2_cogs_dic[cog]
        phen2_totalSpec += phen2_cogs_dic[cog]
    phen2cogs_df.loc[i] = [cog, phen1_spec, phen2_spec]    
phen2cogs_df[f'{phen_1}_rel_spectra'] = [(item/phen1_totalSpec)*100 for item in phen2cogs_df[f'{phen_1}_spectra']]
phen2cogs_df[f'{phen_2}_rel_spectra'] = [(item/phen2_totalSpec)*100 for item in phen2cogs_df[f'{phen_2}_spectra']]

phen2cogs_df.sort_values(by = f'{phen_1}_spectra', ascending = False, inplace = True)
phen2cogs_df.to_csv(out_dir + f'{genome1}_vs_{genome2}_cogExpressions.tsv', sep = '\t', index = False)

r = 0
for cat in phen1_normalizedCat2spec_dic:
    n_spec = phen1_normalizedCat2spec_dic[cat]
    phen2cats_df.loc[r] = [cat, n_spec, phen_1]
    r += 1
    
for cat in phen2_normalizedCat2spec_dic:
    n_spec = phen2_normalizedCat2spec_dic[cat]
    phen2cats_df.loc[r] = [cat, n_spec, phen_2]
    r += 1

phen2cats_df.to_csv(out_dir + f'{genome1}_vs_{genome2}_cogCatExpressions.tsv', sep = '\t', index = False)

mrca = get_MRCA(genome1, genome2)
    
import plotly.express as px
fig = px.bar(phen2cats_df, x = 'COG_category', y = 'rel_spectra', color = 'phenotype', barmode = 'group', title = f'genus:{genome1} vs genus:{genome2} <br>Taxonomy {mrca}')

fig.write_image(out_dir + f'{genome1}_vs_{genome2}_cogCategories.png')



cmd1 = f'diamond makedb --in {proteomes_dir + genome1}.fna.FGS.faa -d {out_dir + genome1}.dmnd'
cmd2 = f'diamond makedb --in {proteomes_dir + genome2}.fna.FGS.faa -d {out_dir + genome2}.dmnd'
os.system(cmd1)
os.system(cmd2)

cmd1 = f'diamond blastp --more-sensitive --query {phen1_expressed_protSeqs_f} --db {out_dir + genome2}.dmnd --out {out_dir}{genome1}_Query_{genome2}_db_diamondBlast.txt --threads 60 --outfmt 6'
cmd2 = f'diamond blastp --more-sensitive --query {phen2_expressed_protSeqs_f} --db {out_dir + genome1}.dmnd --out {out_dir}{genome2}_Query_{genome1}_db_diamondBlast.txt --threads 60 --outfmt 6'
os.system(cmd1)
os.system(cmd2)

out_fmt_6_cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

genome1_diamond_out_df = pd.read_csv(f'{out_dir}{genome1}_Query_{genome2}_db_diamondBlast.txt', sep = '\t')
genome2_diamond_out_df = pd.read_csv(f'{out_dir}{genome2}_Query_{genome1}_db_diamondBlast.txt', sep = '\t')

genome1_diamond_out_df.columns = out_fmt_6_cols
genome1_diamond_out_best_df = genome1_diamond_out_df.drop_duplicates(subset = 'qseqid', keep = 'first')
genome1_diamond_out_best_df.reset_index(drop = True, inplace = True)

genome2_diamond_out_df.columns = out_fmt_6_cols
genome2_diamond_out_best_df = genome2_diamond_out_df.drop_duplicates(subset = 'qseqid', keep = 'first')
genome2_diamond_out_best_df.reset_index(drop = True, inplace = True)

def add_spectralSupport(df, prot2spectra_dic, prot2cog_dic):
    specs = []
    for proteinID in list(df['qseqid']):
        n_specs = prot2spectra_dic[proteinID]
        specs.append(n_specs)
    df = copy.deepcopy(df)
    df['n_specs'] = specs

    for proteinID in set(prot2spectra_dic.keys()).difference(set(df['qseqid'])):
        n_specs = prot2spectra_dic[proteinID]
        row = [proteinID, 'NA', -1, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', n_specs]
        df.loc[len(df)] = row
    df.sort_values(by = 'n_specs', ascending = False, inplace = True)
    
    cogs = ['|'.join(prot2cog_dic[item]) if item in prot2cog_dic else 'NA' for item in df['qseqid']]
    cogs = ['NA' if item == '' else item for item in cogs]
    df['COG'] = cogs
    
    return df

genome1_diamond_out_best_df = add_spectralSupport(genome1_diamond_out_best_df, phen1_prot2spectra_dic, phen1_prot2cog_dic) 
genome1_diamond_out_best_df.to_csv(f'{out_dir}{genome1}_Query_{genome2}_db_diamondBlastWithSpectralSupport.tsv', sep = '\t', index = False)
genome1_specific_df = genome1_diamond_out_best_df[(genome1_diamond_out_best_df['pident'] < max_pid) & (genome1_diamond_out_best_df['n_specs'] >= min_spec )]
genome1_specific_df.to_csv(f'{out_dir}{genome1}_Specific_db_diamondBlastWithSpectralSupport.tsv', sep = '\t', index = False)

genome2_diamond_out_best_df = add_spectralSupport(genome2_diamond_out_best_df, phen2_prot2spectra_dic, phen2_prot2cog_dic) 
genome2_diamond_out_best_df.to_csv(f'{out_dir}{genome2}_Query_{genome1}_db_diamondBlastWithSpectralSupport.tsv', sep = '\t', index = False)
genome2_specific_df = genome2_diamond_out_best_df[(genome2_diamond_out_best_df['pident'] < max_pid) & (genome2_diamond_out_best_df['n_specs'] >= min_spec )]
genome2_specific_df.to_csv(f'{out_dir}{genome2}_Specific_db_diamondBlastWithSpectralSupport.tsv', sep = '\t', index = False)
