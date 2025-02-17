#!/usr/bin/env python
# -*- coding: utf-8 -*-
# python count_miRNA.py -i input1.txt input2.txt ... -o output.txt

import argparse
import os
import pandas as pd

def read_parsed_bowtie(file, qscore=0):
    """Reads parsed Bowtie results and counts miRNA occurrences."""
    bowtie = {}
    with open(file, 'r') as f:
        l = f.readline()
        if l[0] == '#':
            l = f.readline()
        while l:
            data = l.rstrip('\n').split('\t')
            if MM_Qscore_pass(data[1], qscore):
                mir_name = data[0] + '#' + remove_Qscore(data[1])
                if mir_name in bowtie:
                    bowtie[mir_name][0] += 1
                else:
                    bowtie[mir_name] = [1]
            l = f.readline()
    return bowtie

def MM_Qscore_pass(mut, qscore=0):
    """Checks if the mutation passes the given Q-score threshold."""
    if mut == 'PM':
        return True
    else:
        for mm in mut.split(','):
            if not int(mm.split(':')[-1]) >= qscore:
                return False
        return True

def remove_Qscore(mut):
    """Removes Q-score values from mutation annotations for sorting."""
    if mut == 'PM':
        return '0'
    else:
        return ','.join([':'.join(x.split(':')[:-1]) for x in mut.split(',')])

# Argument Parsing
parser = argparse.ArgumentParser(usage='python count_miRNA.py [ -i input1.txt input2.txt ... / -d input_dir ] -o output.txt')
parser.add_argument('-i', metavar='Input files', nargs='*', help='Parsed bowtie result files.')
parser.add_argument('-d', metavar='Input directory', help='Directory with parsed bowtie result files. (suffix= ".parsed.txt")')
parser.add_argument('-o', metavar='Output file', required=True, help='Name for output file.')
parser.add_argument('-q', metavar='Q-score', choices=range(42), type=int, default=38,
                    help='[optional] Q-score filter threshold for mismatch base. (higher or equal, Q>=38)')
p = vars(parser.parse_args())

# Get input files
in_files = p['i']
if p['d']:
    in_files = []
    if p['d'][-1] != '/':
        p['d'] += '/'
    for file in os.listdir(p['d']):
        if '.parsed.txt' in file:
            in_files.append(p['d'] + file)

# Check for duplicate files
if len(in_files) != len(set(in_files)):
    print("Warning: duplicate files given. Please check the file names.\n")
    quit()

# Merge files into a table
final_df = pd.DataFrame({})
for file in in_files:
    bowtie_dict = read_parsed_bowtie(file, p['q'])
    miR_list = sorted(bowtie_dict.keys())
    bowtie_df = pd.DataFrame({file: [bowtie_dict[x][0] for x in miR_list]}, index=miR_list)
    final_df = pd.concat([final_df, bowtie_df], axis=1, sort=False).fillna(0)

# Process miRNA names and mutations
final_df['miRNA name#pos:mut'] = final_df.index
final_df[['miRNA name', 'pos:mut']] = final_df.pop('miRNA name#pos:mut').str.split("#", expand=True)

# Split "pos:mut" column and handle cases where only one column is returned
tmp_df = final_df['pos:mut'].str.split(",", expand=True).fillna(0)
if tmp_df.shape[1] == 1:
    tmp_df.columns = ['1']
    tmp_df['2'] = '0'  # Assign a default value to the second column
elif tmp_df.shape[1] >= 2:
    tmp_df = tmp_df.iloc[:, :2]
    tmp_df.columns = ['1', '2']
final_df = final_df.join(tmp_df)

# Further split columns into position and mutation
for col in list(tmp_df.columns):
    split_df = final_df.pop(col).astype(str).str.split(":", expand=True).fillna('0')

    # If only one column exists after splitting, add a second column with default '0'
    if split_df.shape[1] == 1:
        split_df['1'] = '0'

    split_df.columns = ['pos' + col, 'mut' + col]
    split_df['pos' + col] = split_df['pos' + col].astype(int)
    final_df = final_df.join(split_df)

# Count mismatches
final_df['mm_numb'] = final_df['pos:mut'].str.count(':')
# Sort by miRNA name, mismatch number, and mutation positions
final_df = final_df.sort_values(['miRNA name', 'mm_numb', 'pos1', 'mut1', 'pos2', 'mut2'])
# Reorder columns
final_df = final_df[list(final_df.columns)[-7:-5] + list(final_df.columns)[:-7]]
# Replace '0' with 'PM'
final_df['pos:mut'] = final_df['pos:mut'].replace('0', 'PM')
# Remove file paths and suffixes from column names
final_df.columns = [os.path.basename(x).split('.')[0] if i > 1 else x for i, x in enumerate(final_df.columns)]

# Add total (PM+1MM+2MM) column for downstream calculations
tmp_sum = final_df.groupby(['miRNA name']).sum()
tmp_sum.columns = [x + ' (PM+1MM+2MM)' for x in tmp_sum.columns]
final_df = pd.merge(final_df, tmp_sum, left_on='miRNA name', right_index=True)

# Final output
final_df.to_csv(p['o'], sep='\t', index=False)

