from Bio import SeqIO
import argparse

# Function to get mapped reads
def bowtie_mapped(file):
    mapped = set()
    with open(file, 'r') as f:
        for line in f:
            read_id = line.split('\t')[2].split('#')[0]
            mapped.add(read_id)
    return mapped

# Function to extract unmapped reads
def extract_unmapped_reads(fastq_file, mapped):
    unmapped_reads = []
    for record in SeqIO.parse(fastq_file, "fastq"):
        if record.seq not in mapped:
            unmapped_reads.append(record)
    return unmapped_reads

# Argument parser
parser = argparse.ArgumentParser(description='Extract unmapped reads from a FastQ file.')
parser.add_argument('-fq', required=True, help='Input FastQ file')
parser.add_argument('-mapped', required=True, help='File containing mapped reads')
parser.add_argument('-out', required=True, help='Output FastQ file for unmapped reads')

# Parse arguments
args = parser.parse_args()

# Process files
mapped_reads = bowtie_mapped(args.mapped)

unmapped_reads = extract_unmapped_reads(args.fq, mapped_reads)

# Write unmapped reads to output file
with open(args.out, 'w') as output_handle:
    SeqIO.write(unmapped_reads, output_handle, "fastq")

