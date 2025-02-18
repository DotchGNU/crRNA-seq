#!/bin/bash
# main.sh: Main pipeline script for processing crRNA-seq data.
# This script executes the entire pipeline while ensuring each step is run
# in its designated conda environment using "conda run".
# All outputs are stored in the "result" directory.

# Set directories
FASTQ_DIR="./fastq"
RESULT_DIR="./result"
mkdir -p "$RESULT_DIR"

# Function: Adapter trimming using cutadapt in the crRNA-seq environment
adapter_trimming() {
    echo "### Adapter trimming ###"
    #ln -s "$FASTQ_DIR" fastq
    mkdir -p "$RESULT_DIR/cutadapt"
    for gz in fastq/*R1*.fastq.gz; do
        sample=$(basename "$gz" | cut -d'_' -f1)
        echo "Processing sample: $sample"
        conda run -n crRNA-seq cutadapt -a TGGAATTCTCGGGTGCCAAGG \
            -o "$RESULT_DIR/cutadapt/${sample}.cutadapt.fastq" "$gz" \
            &> "$RESULT_DIR/cutadapt/${sample}.cutadapt_log.txt"
    done
}

# Function: Generate adapter trimming summary
adapter_trimming_summary() {
    echo ""
    echo "### Summary: Adapter trimming ###"
    echo -e "Sample_name\tTotal\tWith_adapters(%)"
    for log in "$RESULT_DIR/cutadapt"/*.cutadapt_log.txt; do
        echo -n "$(basename "$log" | cut -d'.' -f1) "
        awk -F: '
            /Total reads processed:/ {gsub(/,/, "", $2); total=$2}
            /Reads with adapters:/ {gsub(/,/, "", $2); with_adapters=$2}
            END {print total, with_adapters}' "$log"
    done
}

# Function: Search for crRNA repeats using a Python script in the crRNA-seq environment
search_crRNA_repeats() {
    echo ""
    echo "### Search crRNA repeats ###"
    mkdir -p "$RESULT_DIR/crRNA_repeat_16nt"
    for fq in "$RESULT_DIR/cutadapt"/*.fastq; do
        sample=$(basename "$fq" | cut -d'.' -f1-2)
        echo "Searching repeats in sample: $sample"
        conda run -n crRNA-seq python ./scripts/fastq/filter_by_sequence_fastq.py \
            -i "$fq" \
            -o "$RESULT_DIR/crRNA_repeat_16nt/${sample}.crRNA_repeat.fastq" \
            -s GTTTTAGAGCTATGCT
    done
}

# Function: Trim crRNA repeats using cutadapt in the crRNA-seq environment
trim_crRNA_repeat() {
    echo ""
    echo "### Trim crRNA repeat ###"
    for fq in "$RESULT_DIR/crRNA_repeat_16nt"/*.fastq; do
        echo "Trimming file: $fq"
        conda run -n crRNA-seq cutadapt -a GTTTTAGAGCTATGCT "$fq" \
            > "${fq/.fastq/.trim.fq}" \
            2> "${fq/.fastq/.cutadapt_log.txt}"
    done
}

# Function: Check trimmed crRNA repeat summary
check_trimmed_crRNA_repeat() {
    echo ""
    echo "### Summary: Repeat trimming ###"
    echo -e "Sample_name\treads_with_crRNA_repeat\ttrimmed(%)"
    for log in "$RESULT_DIR/crRNA_repeat_16nt"/*.cutadapt_log.txt; do
        echo -n "$(basename "$log" | cut -d'.' -f1) "
        awk -F: '
            /Total reads processed:/ {gsub(/,/, "", $2); total=$2}
            /Reads with adapters:/ {gsub(/,/, "", $2); with_adapters=$2}
            END {print total, with_adapters}' "$log"
    done
}

# Function: Spacer length processing
function spacer_length_processing {
    echo ""
    echo "########################################################"
    echo "######### Spacer Length ################################"
    echo "########################################################"
    mkdir -p "$RESULT_DIR/spacer_length"
    for fq in "$RESULT_DIR/crRNA_repeat_16nt"/*.trim.fq; do
        for seq in TGCGCTGGTTGATTTCTTCTTGCGCTTTTT TTATATGAACATAACTCAATTTGTAAAAAA AGGAATATCCGCAATAATTAATTGCGCTCT AGTGCCGAGGAAAAATTAGGTGCGCTTGGC TAAATTTGTTTAGCAGGTAAACCGTGCTTT TTCAGCACACTGAGACTTGTTGAGTTCCAT; do
            name=$(basename "$fq")
            echo "Processing sample:$name / Spacer:$seq"
            conda run -n crRNA-seq python ./scripts/grep_count_by_length_using_regex.py \
                -i "$fq" \
                -s "$seq" \
                >> "$RESULT_DIR/spacer_length/${name/.fq/.spacer.regex.txt}"
        done
    done
}

# Function: Merge spacer length results using a Python script in the crRNA-seq environment
merge_spacer_length_results() {
    conda run -n crRNA-seq python ./scripts/merge_and_lookup_table.spacer_regex.py \
        -i "$RESULT_DIR/spacer_length"/*.spacer*.regex.txt \
        -o "$RESULT_DIR/SF370.Spacer_length_regex.txt"
}

# Function: Filter reads by length (5-30 nt) using a Python script in the crRNA-seq environment
length_cutoff() {
    mkdir -p "$RESULT_DIR/crRNA_repeat_16nt/5-30nt"
    for fq in "$RESULT_DIR/crRNA_repeat_16nt"/*.trim.fq; do
        sample=$(basename "$fq")
        echo "Filtering file by length: $sample"
        conda run -n crRNA-seq python ./scripts/fastq/filter_fastq_by_read_length.py \
            -i "$fq" \
            -o "$RESULT_DIR/crRNA_repeat_16nt/5-30nt/${sample/.fq/.5-30.fq}" \
            -m 5 -M 30
    done
}

# Function: Reverse mapping (20nt crRNA) in the reverse_mapping environment
reverse_mapping() {
    echo ""
    echo "######### Reverse Mapping: 20nt crRNA #########"
    mkdir -p "$RESULT_DIR/reverse_mapping/spacer_20nt"
    conda run -n reverse_mapping python -W ignore::FutureWarning ./scripts/Reverse_mapping/Reverse_mapping_with_bowtie.py \
        -d "$RESULT_DIR/crRNA_repeat_16nt/5-30nt" \
        -r ./data/SF370_spacer_20nt.fa \
        -o "$RESULT_DIR/reverse_mapping/spacer_20nt" \
        -q 30
}

# Function: Extract exact 20nt aligned reads using awk
extract_exact_20nt_aligned_reads() {
    mkdir -p "$RESULT_DIR/reverse_mapping/spacer_20nt/20nt_parsed"
    for txt in "$RESULT_DIR/reverse_mapping/spacer_20nt"/*.parsed.txt; do
        sample=$(basename "$txt")
        echo "Extracting exact 20nt aligned reads for: $sample"
        awk -F'\t' 'NR==1 {print $0} NR>1 && $4 == "" && $6 == "" {print $0}' \
            "$txt" > "$RESULT_DIR/reverse_mapping/spacer_20nt/20nt_parsed/${sample/.parsed.txt/.20nt.parsed.txt}"
    done
}

# Function: Parse 20nt aligned reads using reverse_mapping in the reverse_mapping environment
parse_20nt_aligned_reads() {
    conda run -n reverse_mapping python -W ignore::FutureWarning ./scripts/Reverse_mapping/count_RNA.py \
        -d "$RESULT_DIR/reverse_mapping/spacer_20nt/20nt_parsed" \
        -o "$RESULT_DIR/SF370_targeted.20nt.count.Q30.txt" \
        -q 30
}

# Function: Extract unmapped reads in the reverse_mapping environment
extract_unmapped_reads() {
    for txt in "$RESULT_DIR/reverse_mapping/spacer_20nt"/*mapped.txt; do
        sample=$(basename "$txt")
        echo "Extracting unmapped reads for: $sample"
        conda run -n crRNA-seq python ./scripts/Reverse_mapping/extract_bowtie_unmapped_reads.biopython.py \
            -fq "$RESULT_DIR/crRNA_repeat_16nt/5-30nt/${sample/.bowtie_mapped.txt/.cutadapt.crRNA_repeat.trim.5-30.fq}" \
            -mapped "$txt" \
            -out "${txt/.bowtie_mapped.txt/.20nt_un.fq}"
    done
}

# Function: Reverse mapping with deletion in the reverse_mapping environment
reverse_mapping_with_deletion() {
    echo ""
    echo "######### Reverse Mapping with Deletion: 20nt crRNA #########"
    mkdir -p "$RESULT_DIR/reverse_mapping/spacer_20nt_with_deletion"
    conda run -n reverse_mapping python -W ignore::FutureWarning ./scripts/Reverse_mapping/Reverse_mapping_with_bowtie.py \
        -i "$RESULT_DIR/reverse_mapping/spacer_20nt"/*.20nt_un.fq \
        -r ./data/SF370_spacer_20nt_with_deletion.fa \
        -o "$RESULT_DIR/reverse_mapping/spacer_20nt_with_deletion" \
        -q 30
}

# Function: Extract exact 20nt aligned reads with deletion using awk
extract_exact_20nt_aligned_reads_with_deletion() {
    mkdir -p "$RESULT_DIR/reverse_mapping/spacer_20nt_with_deletion/20nt_parsed"
    for txt in "$RESULT_DIR/reverse_mapping/spacer_20nt_with_deletion"/*.parsed.txt; do
        sample=$(basename "$txt")
        echo "Extracting exact 20nt aligned reads with deletion for: $sample"
        awk -F'\t' 'NR==1 {print $0} NR>1 && $4 == "" && $6 == "" {print $0}' \
            "$txt" > "$RESULT_DIR/reverse_mapping/spacer_20nt_with_deletion/20nt_parsed/${sample/.parsed.txt/.20nt.parsed.txt}"
    done
}

# Function: Parse 20nt aligned reads with deletion in the reverse_mapping environment
parse_20nt_aligned_reads_with_deletion() {
    conda run -n reverse_mapping python ./scripts/Reverse_mapping/count_RNA.py \
        -d "$RESULT_DIR/reverse_mapping/spacer_20nt_with_deletion/20nt_parsed" \
        -o "$RESULT_DIR/SF370_targeted.20nt_1deletion.count.Q30.txt" \
        -q 30
}

# Function: Generate 5' sequence distribution by calling an external shell script
generate_5p_distribution() {
    echo ""
    echo "######### Generate 5' Sequence Distribution #########"
    pushd "$RESULT_DIR" > /dev/null
    bash ../scripts/collect_5p_sequence_and_calculate_histogram_for_q30.SF370.sh
    popd > /dev/null
}

# Execute pipeline steps sequentially
adapter_trimming
adapter_trimming_summary
search_crRNA_repeats
trim_crRNA_repeat
check_trimmed_crRNA_repeat
spacer_length_processing
merge_spacer_length_results
length_cutoff
reverse_mapping
extract_exact_20nt_aligned_reads
parse_20nt_aligned_reads
extract_unmapped_reads
reverse_mapping_with_deletion
extract_exact_20nt_aligned_reads_with_deletion
parse_20nt_aligned_reads_with_deletion
generate_5p_distribution

echo "Pipeline completed successfully."

