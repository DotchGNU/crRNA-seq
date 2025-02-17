import argparse
import re
from Bio import SeqIO

def modify_pattern(pattern, replace_index):
    """
    Modify the pattern by replacing the character at replace_index with '.'
    replace_index is 1-based.
    """
    if replace_index is not None and 1 <= replace_index <= len(pattern):
        # Convert to 0-based index and replace character at index with '.'
        pattern_list = list(pattern)
        pattern_list[replace_index - 1] = '.'
        pattern = ''.join(pattern_list)
    return pattern

def read_fastq(input_file):
    """
    Read FASTQ file and return a list of sequences.
    """
    return [str(record.seq) for record in SeqIO.parse(input_file, "fastq")]

def filter_sequences(sequences, pattern):
    """
    Filter sequences that do not match the given pattern.
    The pattern includes a start, any character within 5nt, the pattern, and the end.
    Returns the count of matched sequences and a list of unmatched sequences.
    """
    pattern_regex = re.compile(f"^{pattern}$")
    filtered_sequences = [seq for seq in sequences if not pattern_regex.match(seq)]
    match_count = len(sequences) - len(filtered_sequences)
    return match_count, filtered_sequences

def main(input_file, initial_pattern, replace_index):
    pattern = modify_pattern(initial_pattern, replace_index)
    current_sequences = read_fastq(input_file)

    print("pattern\tlength\tmatch_count\tremained_reads")
    while len(pattern) > 0:
        match_count, current_sequences = filter_sequences(current_sequences, pattern)
        print(f"{pattern}\t{len(pattern)}\t{match_count}\t{len(current_sequences)}")

        if pattern[0:2] == ".?":
            pattern = pattern[2:]
        else :
            pattern = pattern[1:]

        #if args.l :
        #    pattern = pattern[:-1] # count length from 5'-end
        #else:
        #    pattern = pattern[1:] # count length from 3'-end / default

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find patterns in FASTQ sequences using Biopython.")
    parser.add_argument("-i", "--input_file", required=True, help="Path to the input FASTQ file")
    parser.add_argument("-s", "--pattern", required=True, help="The search pattern")
    parser.add_argument("-r", type=int, help="The position in the pattern to replace with any character", required=False)
    parser.add_argument('-l', action='store_true', help='count length from 3-end / Default = 3-end')


    args = parser.parse_args()
    main(args.input_file, args.pattern, args.r)

