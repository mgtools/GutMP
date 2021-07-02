"""
script that for each taxonomic group that is common between two given phenotypes, it will go over the strains that are preferred 
to each phenotype, for each taxonomic group it will maintain two lists of strains, one for each phenotype, then it will extract the 
proteomes of each strain group for each phenotype, and also extract information about protein to number of expressed spectra
spectra
starting from DESeq2 suggested genomes
"""
import json
import pandas as pd
import os
from Bio import SeqIO
import copy 
import sys

in_dir = sys.argv[1]#'/data/mstambou/proteome_landscapes/highly_abundant_genomes/DE_genomes/'

phen_1 = sys.argv[2]#'healthy'
phen_2 = sys.argv[3]#'CD'

min_spec = int(sys.argv[4])#10
max_pid = float(sys.argv[5])#70


#in_dir = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/DE_genomes/'

#phen_1 = 'healthy_subT1D'
#phen_2 = 'T1D'

#min_spec = 10
#max_pid = 70

reverse_list = ['healthy_subT1D_vs_T1D', 'healthy_subUC_vs_UC', 'healthy_vs_UC', 'CD_subUC_vs_UC']


phen1_DE_strains_f = f'{in_dir}{phen_1}_vs_{phen_2}/{phen_1}_DE_strains.tsv'
phen2_DE_strains_f = f'{in_dir}{phen_1}_vs_{phen_2}/{phen_2}_DE_strains.tsv'
strain2phen_f = f'{in_dir}{phen_1}_vs_{phen_2}/{phen_1}_vs_{phen_2}_strain2phen.tsv'

phen1_DE_strains_df = pd.read_csv(phen1_DE_strains_f, sep = '\t')
phen2_DE_strains_df = pd.read_csv(phen2_DE_strains_f, sep = '\t')
strain2phen_df = pd.read_csv(strain2phen_f, sep = '\t')

out_dir = in_dir + f'{phen_1}_vs_{phen_2}/DE_strains/'
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

phen2samples_dic_f = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/samplesSubSamples2phenotypes_dic.json'
genome2sample2Proteins2Nspecs_dic_f = '/data/mstambou/proteome_landscapes/highly_abundant_genomes/genome2sample2Proteins2Nspecs_dic.json'
allBin2taxa_dic_f = '/data/mstambou/proteome_landscapes/auxilary/allBin2taxon_dic.json'
proteome2COG_dic_f = '/data/mstambou/proteome_landscapes/auxilary/allProteome2COG_moreSensitive_dic.json'
proteomes_dir = '/data/mstambou/proteome_landscapes/HAPiID_6k/data/proteomes/'

print('loading annotations into memory, please be patient ...')
with open(proteome2COG_dic_f, 'r') as in_f:
    proteome2COG_dic = json.load(in_f)    
with open(allBin2taxa_dic_f, 'r') as in_f:
    allBin2taxa_dic = json.load(in_f)    
with open(phen2samples_dic_f, 'r') as in_f:
    phen2samples_dic = json.load(in_f)
with open(genome2sample2Proteins2Nspecs_dic_f, 'r') as in_f:
    genome2sample2Proteins2Nspecs_dic = json.load(in_f)
    

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


phen1_dir = out_dir + f'{phen_1}_strainSpecific/'

if not os.path.isdir(phen1_dir):
    os.mkdir(phen1_dir)

for taxa in set(phen1_DE_strains_df['taxa']):
    phen1_genomes = list(phen1_DE_strains_df[phen1_DE_strains_df['taxa'] == taxa]['genome'])
    taxa_phens = strain2phen_df[(strain2phen_df['taxa'] == taxa) & (strain2phen_df['phenotype'].isin([f'{phen_2}', 'both']) ) ]
    sister_strains = set(taxa_phens['genome']).difference(phen1_genomes)
    if sister_strains:        
        taxa_name = '_'.join(taxa.split(' '))
        taxa_dir = phen1_dir + taxa_name +'/'

        if not os.path.isdir(taxa_dir):
            os.mkdir(taxa_dir)

        phen1_samples = phen2samples_dic[f'{phen_1}']

        phen1_expressed_proteins2Nspecs_dic = {}

        phen1_expressed_proteinSeqs_dic = {}

        sister_strains_proteomes_dic = {}

        for strain in phen1_genomes:    
            samples = set(phen1_samples).intersection(set(genome2sample2Proteins2Nspecs_dic[strain].keys()))
            strain_proteins = set()
            strain_proteome = SeqIO.to_dict(SeqIO.parse(proteomes_dir + strain +'.fna.FGS.faa', 'fasta'))
            for sample in samples:
                proteins2Nspecs = genome2sample2Proteins2Nspecs_dic[strain][sample]
                for protein in proteins2Nspecs:
                    if protein not in phen1_expressed_proteins2Nspecs_dic:
                        phen1_expressed_proteins2Nspecs_dic[protein] = proteins2Nspecs[protein]
                    else:
                        phen1_expressed_proteins2Nspecs_dic[protein] += proteins2Nspecs[protein]
                    strain_proteins.add(protein)
            phen1_expressed_proteinSeqs_dic.update({k:v for k,v in strain_proteome.items() if k in strain_proteins})

        for sister_strain in sister_strains:
            sister_strain_proteome = SeqIO.to_dict(SeqIO.parse(proteomes_dir + sister_strain +'.fna.FGS.faa', 'fasta'))
            sister_strains_proteomes_dic.update({k:v for k,v in sister_strain_proteome.items() })

        phen1_expressed_proteinSeqs_f = f'{taxa_dir + taxa_name}_expressedProtSeqs.faa'
        sister_strains_proteomes_f = f'{taxa_dir + taxa_name}_sisterStrainsProteomes.faa'

        with open(phen1_expressed_proteinSeqs_f, 'w') as out_f:
            {out_f.write(f'>{k}\n{str(phen1_expressed_proteinSeqs_dic[k].seq)}\n') for k in phen1_expressed_proteinSeqs_dic}

        with open(sister_strains_proteomes_f, 'w') as out_f:
            {out_f.write(f'>{k}\n{str(sister_strains_proteomes_dic[k].seq)}\n') for k in sister_strains_proteomes_dic}

        cmd1 = f'diamond makedb --in {sister_strains_proteomes_f} -d {taxa_dir + taxa_name}_sisterStrainsDB.dmnd'
        cmd2 = f'diamond blastp --query {phen1_expressed_proteinSeqs_f} --db {taxa_dir + taxa_name}_sisterStrainsDB.dmnd --out {taxa_dir}{taxa_name}_{phen_1}_Query_{phen_2}_db_diamondBlast.txt --threads 60 --outfmt 6'
        os.system(cmd1)
        os.system(cmd2)


        out_fmt_6_cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

        phen1_diamond_out_df = pd.read_csv(f'{taxa_dir}{taxa_name}_{phen_1}_Query_{phen_2}_db_diamondBlast.txt', sep = '\t')
        phen1_diamond_out_df.columns = out_fmt_6_cols
        phen1_diamond_out_best_df = phen1_diamond_out_df.drop_duplicates(subset = 'qseqid', keep = 'first')
        phen1_diamond_out_best_df.reset_index(drop = True, inplace = True)

        phen1_diamond_out_best_df = add_spectralSupport(phen1_diamond_out_best_df, phen1_expressed_proteins2Nspecs_dic, proteome2COG_dic) 
        phen1_diamond_out_best_df.to_csv(f'{taxa_dir}{taxa_name}_{phen_1}_Query_{phen_2}_db_diamondBlastWithSpectralSupport.tsv', sep = '\t', index = False)
        phen1_specific_df = phen1_diamond_out_best_df[(phen1_diamond_out_best_df['pident'] < max_pid) & (phen1_diamond_out_best_df['n_specs'] >= min_spec )]
        phen1_specific_df.to_csv(f'{taxa_dir}{taxa_name}_{phen_1}_Query_{phen_2}_Specific_db_diamondBlastWithSpectralSupport.tsv', sep = '\t', index = False)

        
phen2_dir = out_dir + f'{phen_2}_strainSpecific/'

if not os.path.isdir(phen2_dir):
    os.mkdir(phen2_dir)

for taxa in set(phen2_DE_strains_df['taxa']):
    phen2_genomes = list(phen2_DE_strains_df[phen2_DE_strains_df['taxa'] == taxa]['genome'])
    taxa_phens = strain2phen_df[(strain2phen_df['taxa'] == taxa) & (strain2phen_df['phenotype'].isin([f'{phen_1}', 'both']) ) ]
    sister_strains = set(taxa_phens['genome']).difference(phen2_genomes)
    if sister_strains:        
        taxa_name = '_'.join(taxa.split(' '))
        taxa_dir = phen2_dir + taxa_name +'/'

        if not os.path.isdir(taxa_dir):
            os.mkdir(taxa_dir)

        phen2_samples = phen2samples_dic[f'{phen_2}']

        phen2_expressed_proteins2Nspecs_dic = {}

        phen2_expressed_proteinSeqs_dic = {}

        sister_strains_proteomes_dic = {}

        for strain in phen2_genomes:    
            samples = set(phen2_samples).intersection(set(genome2sample2Proteins2Nspecs_dic[strain].keys()))
            strain_proteins = set()
            strain_proteome = SeqIO.to_dict(SeqIO.parse(proteomes_dir + strain +'.fna.FGS.faa', 'fasta'))
            for sample in samples:
                proteins2Nspecs = genome2sample2Proteins2Nspecs_dic[strain][sample]
                for protein in proteins2Nspecs:
                    if protein not in phen2_expressed_proteins2Nspecs_dic:
                        phen2_expressed_proteins2Nspecs_dic[protein] = proteins2Nspecs[protein]
                    else:
                        phen2_expressed_proteins2Nspecs_dic[protein] += proteins2Nspecs[protein]
                    strain_proteins.add(protein)
            phen2_expressed_proteinSeqs_dic.update({k:v for k,v in strain_proteome.items() if k in strain_proteins})

        for sister_strain in sister_strains:
            sister_strain_proteome = SeqIO.to_dict(SeqIO.parse(proteomes_dir + sister_strain +'.fna.FGS.faa', 'fasta'))
            sister_strains_proteomes_dic.update({k:v for k,v in sister_strain_proteome.items() })

        phen2_expressed_proteinSeqs_f = f'{taxa_dir + taxa_name}_expressedProtSeqs.faa'
        sister_strains_proteomes_f = f'{taxa_dir + taxa_name}_sisterStrainsProteomes.faa'

        with open(phen2_expressed_proteinSeqs_f, 'w') as out_f:
            {out_f.write(f'>{k}\n{str(phen2_expressed_proteinSeqs_dic[k].seq)}\n') for k in phen2_expressed_proteinSeqs_dic}

        with open(sister_strains_proteomes_f, 'w') as out_f:
            {out_f.write(f'>{k}\n{str(sister_strains_proteomes_dic[k].seq)}\n') for k in sister_strains_proteomes_dic}

        cmd1 = f'diamond makedb --in {sister_strains_proteomes_f} -d {taxa_dir + taxa_name}_sisterStrainsDB.dmnd'
        cmd2 = f'diamond blastp --query {phen2_expressed_proteinSeqs_f} --db {taxa_dir + taxa_name}_sisterStrainsDB.dmnd --out {taxa_dir}{taxa_name}_{phen_2}_Query_{phen_1}_db_diamondBlast.txt --threads 60 --outfmt 6'
        os.system(cmd1)
        os.system(cmd2)


        out_fmt_6_cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

        phen2_diamond_out_df = pd.read_csv(f'{taxa_dir}{taxa_name}_{phen_2}_Query_{phen_1}_db_diamondBlast.txt', sep = '\t')
        phen2_diamond_out_df.columns = out_fmt_6_cols
        phen2_diamond_out_best_df = phen2_diamond_out_df.drop_duplicates(subset = 'qseqid', keep = 'first')
        phen2_diamond_out_best_df.reset_index(drop = True, inplace = True)

        phen2_diamond_out_best_df = add_spectralSupport(phen2_diamond_out_best_df, phen2_expressed_proteins2Nspecs_dic, proteome2COG_dic) 
        phen2_diamond_out_best_df.to_csv(f'{taxa_dir}{taxa_name}_{phen_2}_Query_{phen_1}_db_diamondBlastWithSpectralSupport.tsv', sep = '\t', index = False)
        phen2_specific_df = phen2_diamond_out_best_df[(phen2_diamond_out_best_df['pident'] < max_pid) & (phen2_diamond_out_best_df['n_specs'] >= min_spec )]
        phen2_specific_df.to_csv(f'{taxa_dir}{taxa_name}_{phen_2}_Query_{phen_1}_Specific_db_diamondBlastWithSpectralSupport.tsv', sep = '\t', index = False)
