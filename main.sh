#!/bin/bash
# main.sh: Main pipeline script for mapping crRNA with Bowtie and regex.
# This version uses a relative directory structure and stores all output in result/

# Activate conda (adjust if necessary)
if [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
else
    echo "Error: conda.sh not found. Please check your Anaconda/Miniconda installation."
    exit 1
fi

# Set the fastq directory (now using the local "./fastq" directory)
FASTQ_DIR="./fastq"

# Define result directory and create it if it doesn't exist
RESULT_DIR="./result"
mkdir -p "$RESULT_DIR"

# Function: Adapter trimming
function adapter_trimming {
    echo "### Adapter trimming ###"
    ln -s "$FASTQ_DIR" fastq  # Create a symbolic link to the fastq directory
    mkdir -p "$RESULT_DIR/cutadapt"
    for gz in fastq/*R1*.fastq.gz; do
        out=$(basename "$gz" | cut -d'_' -f1)
        echo "$out"
        cutadapt -a TGGAATTCTCGGGTGCCAAGG \
            -o "$RESULT_DIR/cutadapt/${out}.cutadapt.fastq" "$gz" \
            &> "$RESULT_DIR/cutadapt/${out}.cutadapt_log.txt"
    done
}

# Function: Adapter trimming summary
function adapter_trimming_summary {
    echo ""
    echo "### Summary : Adapter trimming ###"
    echo -e "Sample_name\tTotal\tWith_adapters(%)"
    for log in "$RESULT_DIR/cutadapt"/*.cutadapt_log.txt; do
        echo -n "$(basename "$log" | cut -d'.' -f1) "
        awk -F: '
            /Total reads processed:/ {gsub(/,/, "", $2); total=$2}
            /Reads with adapters:/ {gsub(/,/, "", $2); with_adapters=$2}
            END {print total, with_adapters}' "$log"
    done
}

# Function: Search crRNA repeats
function search_crRNA_repeats {
    echo ""
    echo "### Search crRNA repeats ###"
    mkdir -p "$RESULT_DIR/crRNA_repeat_16nt"
    conda activate GWK
    for fq in "$RESULT_DIR/cutadapt"/*.fastq; do
        out=$(basename "$fq" | cut -d'.' -f1-2)
        echo "$out"
        python ./scripts/fastq/filter_by_sequence_fastq.py \
            -i "$fq" \
            -o "$RESULT_DIR/crRNA_repeat_16nt/${out}.crRNA_repeat.fastq" \
            -s GTTTTAGAGCTATGCT
    done
}

# Function: Check crRNA repeat count
function check_crRNA_repeat_count {
    echo ""
    echo "### Reads with crRNA repeats (16nt) ###"
    echo -e "Sample_name\treads_with_crRNA_repeat\ttrimmed(%)"
    for log in "$RESULT_DIR/crRNA_repeat_16nt"/*.cutadapt_log.txt; do
        echo -n "$(basename "$log" | cut -d'.' -f1) "
        awk -F: '
            /Total reads processed:/ {gsub(/,/, "", $2); total=$2}
            /Reads with adapters:/ {gsub(/,/, "", $2); with_adapters=$2}
            END {print total, with_adapters}' "$log"
    done
}

# Function: Trim crRNA repeats
function trim_crRNA_repeat {
    echo ""
    echo "### Trim crRNA repeat ###"
    for fq in "$RESULT_DIR/crRNA_repeat_16nt"/*.fastq; do
        echo "$fq"
        cutadapt -a GTTTTAGAGCTATGCT "$fq" \
            > "${fq/.fastq/.trim.fq}" \
            2> "${fq/.fastq/.cutadapt_log.txt}"
    done
}

# Function: Check trimmed crRNA repeats
function check_trimmed_crRNA_repeat {
    echo ""
    echo "### Summary : Repeat trimming ###"
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
    conda activate GWK
    for fq in "$RESULT_DIR/crRNA_repeat_16nt"/*.trim.fq; do
        for seq in TGCGCTGGTTGATTTCTTCTTGCGCTTTTT TTATATGAACATAACTCAATTTGTAAAAAA AGGAATATCCGCAATAATTAATTGCGCTCT AGTGCCGAGGAAAAATTAGGTGCGCTTGGC TAAATTTGTTTAGCAGGTAAACCGTGCTTT TTCAGCACACTGAGACTTGTTGAGTTCCAT; do
            name=$(basename "$fq")
            echo "$name"
            python ./scripts/grep_count_by_length_using_regex.py \
                -i "$fq" \
                -s "$seq" \
                >> "$RESULT_DIR/spacer_length/${name/.fq/.spacer.regex.txt}"
        done
    done
}

# Function: Merge spacer length results
function merge_spacer_length_results {
    python ./scripts/merge_and_lookup_table.spacer_regex.py \
        -i "$RESULT_DIR/spacer_length"/*.spacer.regex.txt \
        -o "$RESULT_DIR/SF370.Spacer_length_regex.txt"
}

# Function: Length cutoff
function length_cutoff {
    mkdir -p "$RESULT_DIR/crRNA_repeat_16nt/5-30nt"
    for fq in "$RESULT_DIR/crRNA_repeat_16nt"/*.trim.fq; do
        name=$(basename "$fq")
        echo "$name"
        python ./scripts/fastq/filter_fastq_by_read_length.py \
            -i "$fq" \
            -o "$RESULT_DIR/crRNA_repeat_16nt/5-30nt/${name/.fq/.5-30.fq}" \
            -m 5 -M 30
    done
}

# Function: Reverse mapping (20nt crRNA)
function reverse_mapping {
    echo ""
    echo "########################################################"
    echo "######### Reverse mapping : 20nt crRNA #################"
    echo "########################################################"
    conda activate reverse_mapping
    mkdir -p "$RESULT_DIR/reverse_mapping/spacer_20nt"
    python ./scripts/Mapping_miRNA/Mapping_miRNA_with_bowtie.py \
        -d "$RESULT_DIR/crRNA_repeat_16nt/5-30nt" \
        -r ./data/SF370_spacer_20nt.fa \
        -o "$RESULT_DIR/reverse_mapping/spacer_20nt" \
        -q 30
}

# Function: Extract exact 20nt aligned reads
function extract_exact_20nt_aligned_reads {
    mkdir -p "$RESULT_DIR/reverse_mapping/spacer_20nt/20nt_parsed"
    for txt in "$RESULT_DIR/reverse_mapping/spacer_20nt"/*.parsed.txt; do
        name=$(basename "$txt")
        echo "$name"
        awk -F'\t' 'NR==1 {print $0} NR>1 && $4 == "" && $6 == "" {print $0}' \
            "$txt" > "$RESULT_DIR/reverse_mapping/spacer_20nt/20nt_parsed/${name/.parsed.txt/.20nt.parsed.txt}"
    done
}

# Function: Parse 20nt aligned reads
function parse_20nt_aligned_reads {
    python ./scripts/Mapping_miRNA/count_miRNA.py \
        -d "$RESULT_DIR/reverse_mapping/spacer_20nt/20nt_parsed" \
        -o "$RESULT_DIR/SF370_targeted.20nt.count.Q30.txt" \
        -q 30
}

# Function: Extract unmapped reads
function extract_unmapped_reads {
    conda activate gwk
    for txt in "$RESULT_DIR/reverse_mapping/spacer_20nt"/*mapped.txt; do
        name=$(basename "$txt")
        echo "$name"
        python ./scripts/Mapping_miRNA/extract_bowtie_unmapped_reads.biopython.py \
            -fq "$RESULT_DIR/crRNA_repeat_16nt/5-30nt/${name/.bowtie_mapped.txt/.cutadapt.crRNA_repeat.trim.5-30.fq}" \
            -mapped "$txt" \
            -out "${txt/.bowtie_mapped.txt/.20nt_un.fq}"
    done
}

# Function: Reverse mapping with deletion
function reverse_mapping_with_deletion {
    echo ""
    echo "#########################################################################"
    echo "######### Reverse mapping : 20nt crRNA with deletion ####################"
    echo "#########################################################################"
    conda activate reverse_mapping
    mkdir -p "$RESULT_DIR/reverse_mapping/spacer_20nt_with_deletion"
    python ./scripts/Mapping_miRNA/Mapping_miRNA_with_bowtie.py \
        -i "$RESULT_DIR/reverse_mapping/spacer_20nt"/*.20nt_un.fq \
        -r ./data/SF370_spacer_20nt_with_deletion.fa \
        -o "$RESULT_DIR/reverse_mapping/spacer_20nt_with_deletion" \
        -q 30
}

# Function: Extract exact 20nt aligned reads with deletion
function extract_exact_20nt_aligned_reads_with_deletion {
    mkdir -p "$RESULT_DIR/reverse_mapping/spacer_20nt_with_deletion/20nt_parsed"
    conda activate GWK
    for txt in "$RESULT_DIR/reverse_mapping/spacer_20nt_with_deletion"/*.parsed.txt; do
        name=$(basename "$txt")
        echo "$name"
        awk -F'\t' 'NR==1 {print $0} NR>1 && $4 == "" && $6 == "" {print $0}' \
            "$txt" > "$RESULT_DIR/reverse_mapping/spacer_20nt_with_deletion/20nt_parsed/${name/.parsed.txt/.20nt.parsed.txt}"
    done
}

# Function: Parse 20nt aligned reads with deletion
function parse_20nt_aligned_reads_with_deletion {
    conda activate reverse_mapping
    python ./scripts/Mapping_miRNA/count_miRNA.py \
        -d "$RESULT_DIR/reverse_mapping/spacer_20nt_with_deletion/20nt_parsed" \
        -o "$RESULT_DIR/SF370_targeted.20nt_1deletion.count.Q30.txt" \
        -q 30
}

# Function: Generate 5p distribution
function generate_5p_distribution {
    echo ""
    echo "####################################################################"
    echo "######### Gather 5p sequence (Q30)  ################################"
    echo "####################################################################"
    # Run the script in the result directory so that its output is stored there
    pushd "$RESULT_DIR" > /dev/null
    bash ../scripts/collect_5p_sequence_and_calculate_histogram_for_q30.SF370.sh
    popd > /dev/null
}

# Run all pipeline steps sequentially
adapter_trimming
adapter_trimming_summary
search_crRNA_repeats
check_crRNA_repeat_count
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

