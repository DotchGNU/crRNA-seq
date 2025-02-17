#!/usr/bin/env python

# ./Mapping_miRNA_with_bowtie.py -i Input_folder/ -r miRNA_sequence.fa -o Output_folder/

# bowtie parameter: --norc -l 7 -n 2 -a -f
# All scripts must be present in same directory!
# Scripts used: FastQcollapse.py, parse_bowtie_result.py, count_miRNA.py, extract_bowtie_unmapped_reads.py


def get_input_file(file_list, dir, format = ['.fq', '.fastq']):
	out_list = []
	if file_list:
		dir = ''
	else:
		if dir[-1] != '/': dir += '/'
		file_list = os.listdir(dir)
	for file in file_list:
		suffix = os.path.splitext(file)[-1]
		if suffix in format:
			fq_file = dir + file
			out_list.append(fq_file)
	return out_list

def run_bowtie_and_parse((fq_file, miR_file, out_name)):
	##### STEP 1) Collapse fastq file #####
	fq_col = out_name + '.col.fa'
	run = subprocess.check_output( "python {0}/FastQcollapse.py {1} > {2}".format(script_loc, fq_file, fq_col), shell=True )

	##### STEP 2) Construct Bowtie DB #####
	discard = open(os.devnull, 'w')	# discard warning messages from bowtie-build, "Warning: Encountered reference sequence with only gaps"
	run = subprocess.check_output( 'bowtie-build -q {0} {1}'.format(fq_col, out_name), stderr=discard, shell=True )

	##### STEP 3) Run bowtie with "--norc -l 7 -n 2 -a -f" #####
	fq_mapped = out_name + '.bowtie_mapped.txt'
	run = subprocess.check_output( 'bowtie {0} --norc -l 7 -n 2 -a -q -f {1} > {2}'.format(out_name, miR_file, fq_mapped), stderr=discard, shell=True )

	##### STEP 4) Select uniquely mapped reads #####
	if os.path.getsize(fq_mapped) != 0:
		out = fq_mapped.replace('.txt', '.parsed.txt')
		run = subprocess.check_output( 'python {0}/parse_bowtie_result.py -i {1} -fq {2} -o {3}'.format(script_loc, fq_mapped, fq_file, out), shell=True )

		### Get unmapped reads (optional) ###
		if p['un']:
			fq_unmapped = out_name + '.bowtie_unmapped.fastq'
			run = subprocess.check_output( 'python {0}/extract_bowtie_unmapped_reads.py -mapped {1} -fq {2} -out {3}'.format(script_loc, fq_mapped, fq_file, fq_unmapped), shell=True)


def get_miRNA_count(mapped_list, Qscore, dir):	# count only mismatch with its Q-score above threshold
	count_file = dir + 'miRNA_count.Q' + str(Qscore) + '.txt'
	run = subprocess.check_output( 'python {0}/count_miRNA.py -i {1} -q {2} -o {3}'.format(script_loc, ' '.join(mapped_list), str(Qscore), count_file), shell=True )
	return count_mapped_read(count_file)


def count_total_read(fq_file):
	return subprocess.check_output( 'wc -l ' + fq_file + ' | awk \'{printf $1/4}\'', shell=True)

def count_mapped_read(count_file):
	df = pd.read_table(count_file)			# read-in count file
	df_col = [ c for c in df.columns if c.split()[-1] != '(PM+1MM+2MM)' ]
	df = df[df_col]					# remove colums with PM+1MM+2MM count
	df['mm_numb'] = df['pos:mut'].str.count(':')	# count no. of mismatch
	df = df.groupby(['mm_numb']).sum()
	df = df.rename(index={0: 'PM', 1: '1MM', 2: '2MM'})
	return df.to_dict()	# { 'file1': {'PM':10, '1MM':2, '2MM':1}, 'file2': {'PM':12, '1MM':2, '2MM':1} }


import subprocess, os, sys, argparse
import pandas as pd
# from multiprocessing import Pool


parser=argparse.ArgumentParser(usage='./Mapping_miRNA_with_bowtie.py [options] -r miRNA_sequence.fa', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-r', metavar='miR_seq.fa', required=True, help='\t\t\t(REQUIRED)')
parser.add_argument('-i', metavar='input.fq', nargs='*', help='\tdefault: None \t(optional)')
parser.add_argument('-d', metavar='input_dir/', default='./', help='\tdefault: ./ \t(optional)')
parser.add_argument('-o', metavar='output_dir/', default='./', help='\tdefault: ./ \t(optional)')
parser.add_argument('-p', metavar='int', choices=range(10), default=4, type=int, help='\tdefault: 4 \t(no. of worker threads)')
parser.add_argument('-q', metavar='int', choices=range(42), default=38, type=int, help='\tdefault: 38 \t(mismatch\'s Q-score must be higher or equal (>=) to this value)')
parser.add_argument('-un', action='store_true', help='\tdefault: False \t(save unmapped reads to fastq file)')
p=vars(parser.parse_args())

### This script's location ###
script_loc = '/'.join( os.path.realpath(__file__).split('/')[:-1] )

fq_file_list = get_input_file(p['i'], p['d'])
out_dir = p['o']
miR_file = p['r']
q_threshold = p['q']	# Q-score threshold value. (only at mutated bases)
threads = p['p']

### Set output file names ###
if out_dir[-1] != '/': out_dir += '/'
out_file_list = [ out_dir + os.path.basename(f).split('.')[0] for f in fq_file_list ]

### Mapping miRNAs with bowtie & parse output ###
# pool = Pool( processes = threads )	# No. of workers
if __name__ == '__main__':
#	processes = pool.map( run_bowtie_and_parse, [(f, miR_file, out_file_list[i]) for i,f in enumerate(fq_file_list)] )
    for i, f in enumerate(fq_file_list):
        run_bowtie_and_parse((f, miR_file, out_file_list[i]))
### Count mapped miRNAs ###
bowtie_file_list = [ f+'.bowtie_mapped.parsed.txt' for f in out_file_list if os.path.getsize(f+'.bowtie_mapped.txt') != 0 ]
mapped_count = get_miRNA_count(bowtie_file_list, q_threshold, out_dir)

### FINAL OUTPUT ###
summary_file = out_dir + 'Mapping_summary.Q' + str(q_threshold) + '.txt'
header = [ 'File name', 'Total read', 'PM read', '%', '1MM read', '%', '2MM read', '%' ]
with open(summary_file, 'w') as f:
	f.writelines( '\t'.join(header) + '\n' )
	for file in fq_file_list:
		total = count_total_read(file)
		file_name = os.path.basename(file).split('.')[0]
		print_out = [file_name, total]
		if mapped_count.has_key(file_name):
			for mm in ['PM', '1MM', '2MM']:
				mapped = mapped_count[file_name][mm]
				print_out.append( str(mapped) )					# mapped read count
				print_out.append( str( round(mapped / float(total) * 100, 2) ) )	# mapping %
		else:	# files with no mapped read
			print_out += [ '0', '0', '0', '0', '0', '0' ]
		f.writelines( '\t'.join(print_out) + '\n' )
