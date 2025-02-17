import argparse

def read_fastq(file):
    with open(file, 'r') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            sequence = f.readline().strip()
            plus = f.readline().strip()
            quality = f.readline().strip()
            yield header, sequence, plus, quality

def filter_reads_by_length(input_file, output_file, min_length, max_length):
    with open(output_file, 'w') as out:
        for header, sequence, plus, quality in read_fastq(input_file):
            if min_length <= len(sequence) <= max_length: # 이상 & 이하
                out.write(f"{header}\n{sequence}\n{plus}\n{quality}\n")

def main():
    parser = argparse.ArgumentParser(description='Filter FASTQ reads by length.')
    parser.add_argument('-i', '--input', required=True, help='Input FASTQ file')
    parser.add_argument('-o', '--output', required=True, help='Output FASTQ file')
    parser.add_argument('-m', '--min', type=int, required=True, help='Minimum length of reads')
    parser.add_argument('-M', '--max', type=int, required=True, help='Maximum length of reads')

    args = parser.parse_args()

    filter_reads_by_length(args.input, args.output, args.min, args.max)

if __name__ == "__main__":
    main()

