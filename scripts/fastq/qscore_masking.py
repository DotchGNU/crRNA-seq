import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def replace_low_quality_bases(record, threshold):
    new_seq = []
    for base, score in zip(record.seq, record.letter_annotations["phred_quality"]):
        if score <= threshold:
            new_seq.append('N')
        else:
            new_seq.append(base)
    new_record = SeqRecord(Seq(''.join(new_seq)), id=record.id, description=record.description)
    new_record.letter_annotations["phred_quality"] = record.letter_annotations["phred_quality"]
    return new_record


def process_fastq(input_file, output_file, threshold):
    with open(input_file, "r") as input_handle, open(output_file, "w") as output_handle:
        records = SeqIO.parse(input_handle, "fastq")
        modified_records = (replace_low_quality_bases(record, threshold) for record in records)
        SeqIO.write(modified_records, output_handle, "fastq")

def main():
    parser = argparse.ArgumentParser(description="Replace low quality bases in FASTQ files.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ file path")
    parser.add_argument("-o", "--output", required=True, help="Output FASTQ file path")
    parser.add_argument("-q", "--quality", type=int, default=30, help="Quality score threshold")
    args = parser.parse_args()

    process_fastq(args.input, args.output, args.quality)

if __name__ == "__main__":
    main()

