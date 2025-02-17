# python FastQcollapse.py "input.fastq" > "output.fasta"

import sys, operator

in_in=sys.argv[1]

with open(in_in,'r') as r:
	l=r.readline()
	res={}
	while l:
		if l[0]=='@':
			l=r.readline()
			if '@' not in l:
				if l:
		                        seq=l.rstrip('\n\r\a')
				else:			# if '@' sign in Q-score is present at EOF.
					continue	# Skip the while loop
			else:
				l=r.readline()
				seq=l.rstrip('\n\r\a')
			if res.has_key(seq):
				res[seq]+=1
			else:
				res[seq]=1
		l=r.readline()

for data in sorted(res.items(), key=operator.itemgetter(1), reverse=True):
	print '>'+data[0]+'#'+str(data[1])
	print data[0]
