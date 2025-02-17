import argparse

def read_fastq(file):
    """read fastq -> header, sequence, plus, quality"""
    with open(file, 'r') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            sequence = f.readline().strip()
            plus = f.readline().strip()
            quality = f.readline().strip()
            yield header, sequence, plus, quality

def filter_reads(input_file, output_file, sequence_to_filter):
    """python string search : if str in str"""
    with open(output_file, 'w') as out:
        for header, sequence, plus, quality in read_fastq(input_file):
            if sequence_to_filter in sequence:
                out.write(f"{header}\n{sequence}\n{plus}\n{quality}\n")

def main():
    parser = argparse.ArgumentParser(description='Filter FASTQ reads containing a specific sequence.')
    parser.add_argument('-i', '--input', required=True, help='Input FASTQ file')
    parser.add_argument('-o', '--output', required=True, help='Output FASTQ file')
    parser.add_argument('-s', '--sequence', required=True, help='Sequence to filter')

    args = parser.parse_args()

    filter_reads(args.input, args.output, args.sequence)

if __name__ == "__main__":
    main()

