# python parse_bowtie_result.py -i bowtie_result.txt -fq bowtie_input.fastq -o output.txt

""" * Alignment criteria
	1) PM > 1MM > 2MM
 	2) Counts "concatemer" as 1. (where, "concatemer" = >=2 identical miRNA sequence (including mismatch) present in a read
    - In case of concatemer, only the first-encountered miRNA is considered. (for Q-score)
"""


def read_bowtie_result(file):
	mapped = {}
	with open(file,'r') as f:
		l=f.readline()
		while l:
			data = l.rstrip('\n').split('\t')	# ['mir-1', '+', 'GTATGCAGGCTT#11', '0', 'GTATGC', 'IIIIII', '2', '']
			if data[1] == '+':			# only forward mapped read (equivalent of bowtie option: --norc)
				miR_name = data[0].split()[0]
				fq_seq = data[2]
				fq_count = 1
				if '#' in data[2]:		# if fastq file is collapsed before bowtie mapping
					fq_seq = data[2].split('#')[0]
					fq_count = data[2].split('#')[-1]
				mm_count = data[7].count('>')				# No. of mismatch (for PM, its 0)
				end = str( int(data[3]) + len(data[4]) )		# index for alignment end
				data = [ data[3], end, miR_name, data[7], fq_count ]	# ['0', '23', 'mir-1', '', '11']

				if mapped.has_key( fq_seq ):
					if mapped[ fq_seq ].has_key( mm_count ):
						mapped[ fq_seq ][ mm_count ].append(data)
					else:
						mapped[ fq_seq ][ mm_count ] = [data]
				else:
					mapped[ fq_seq ] = { mm_count : [data] }

				# mapped = { 'read_1': { 0: [ ['0','23','mir-1','','11'], ... ]
				#			,1: [ [], ... ]
				#			,2: [ [], ... ]	}
				#	    ,'read_2': {}, ... }

			l=f.readline()
        return mapped


# bowtie_dict = { 'read_1': { 0: [ ['0','23','mir-1','','11'], ... ] }, 'read_2': {} }
def select_uniquely_mapped_read(bowtie_dict):
	uniq_mapped = {}
	multi_mapped = {}
        multi_mm_count_list = []
	for fq_seq, score_dict in bowtie_dict.iteritems():	# score_dict = { 0: [ ['0','23','mir-1','','11'], ... ] }
                mm_count = min(score_dict.keys())

		# Uniquely mapped read (with the least mismatch; PM > 1MM > 2MM)
		if len( score_dict[mm_count] ) == 1:
			uniq_mapped[fq_seq] = score_dict[mm_count][0]

		# Concatamer (>1 "identical" miRNAs mapped to a read)
		elif len( set([ tuple(x[2:-1]) for x in score_dict[mm_count] ]) ) == 1:
			start_index_list = [ int(x[0]) for x in score_dict[mm_count] ]	# [0, 4, 25 ...] = list of start indices for concatamer
			start_index = start_index_list.index(min(start_index_list))
			uniq_mapped[fq_seq] = score_dict[mm_count][start_index]		# include an alignment with lowest start index

		# Unambiguously mapped reads (>1 "different" miRNAs mapped to a read)
		else:
			multi_mapped[fq_seq] = score_dict[mm_count]
                        multi_mm_count_list.append(mm_count)

                # if multiple PM reads is exist : warning
                if 0 in multi_mm_count_list:
                    print ("Warning: Some perfect match reads is discarded")
	return uniq_mapped, multi_mapped


def get_qscore(file, seq_list=[]):
	res = {}
        with open(file,'r') as r:
                l=r.readline()	# first line of fastq
		if l[0]=='@':
	                while l:
                                l=r.readline()	# sequence
                                seq=l.rstrip('\n\r\a')
				l=r.readline()	# strand
				l=r.readline()	# Q-score
				qscore=l.rstrip('\n\r\a')
				if res.has_key(seq):
					res[seq].append(qscore)
				else:
					res[seq] = [qscore]
	                        l=r.readline()
		else:
			print 'Warning! "' + file + '" invalid fastq format. (does not start with @)\n'
			quit()
	if seq_list != []:
		return { seq : res[seq] for seq in seq_list }
	else:
		return res


#q_ref = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"  # Illumina 1.8+ Q-score
q_ref = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
q_dict = { x:i for i,x in enumerate(q_ref) }

def correct_bowtie_mut_Qscore(mut, start, qscore):
        if mut == '':
                return 'PM', qscore
        else:
                mismatch = []
                for mm in mut.split(','):
                        tmp = mm.split(':')
                        mm_pos = int(tmp[0])			# position according to miRNA
			qscore_pos = mm_pos + int(start)	# position according to Q-score
                        mm_qscore = q_dict[ qscore[qscore_pos] ]			# get mismatch Q-score
			new_qscore = qscore[:qscore_pos] + ' ' + qscore[qscore_pos+1:]	# exclude mismatch Q-score
                        mismatch.append( str(mm_pos+1) + ':' + tmp[1][-1] + tmp[1][0] + ':' + str(mm_qscore) )
                return ','.join(mismatch), new_qscore


import argparse

parser=argparse.ArgumentParser(usage='python parse_bowtie_result.py [options] -i bowtie_result.txt -fq bowtie_input.fastq -o output.txt')
parser.add_argument('-i', metavar='bowtie result', required=True, help='bowtie_result file in its default table format')
parser.add_argument('-o', metavar='file name', help='By default, ".parsed" will be added to the input file\'s name.')
parser.add_argument('-fq', metavar='fastq file', required=True, help='fastq file used for bowtie mapping')
p=vars(parser.parse_args())


if p['o']==None:
	p['o'] = p['i'] + '.parsed'

mapped_reads, discard = select_uniquely_mapped_read( read_bowtie_result(p['i']) )
q_score_dict = get_qscore( p['fq'], mapped_reads.keys() )

heading = ['# miRNA_name', 'pos:mut:Q-score', 'Q-score', '5p_sequence', '5p_Q-score', '3p_sequence', '3p_Q-score']
with open(p['o'],'w') as f:
	f.writelines('\t'.join(heading)+'\n')

	for fq_read, data in mapped_reads.iteritems():	# data = ['0', '23', 'mir-1', '', '11']
		start = int(data[0])
		end = int(data[1])
		for q_score in q_score_dict[ fq_read ]:
			mut, q_score = correct_bowtie_mut_Qscore( data[-2], start, q_score )
			f.writelines( '\t'.join([ data[2], mut, q_score[start:end], fq_read[:start], q_score[:start], fq_read[end:], q_score[end:] ]) + '\n' )
