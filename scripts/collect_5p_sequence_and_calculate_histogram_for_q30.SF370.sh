#!/bin/bash
# collect_5p_sequence_and_calculate_histogram_for_q30.SF370.sh
# This script collects 5' sequence data from reverse mapping results,
# applies Q-score masking, calculates the 5' distribution, and merges the results.

# Set input and output directories
input_dir="reverse_mapping/spacer_20nt/"
output_dir="5p_sequence_dist_Q30"
mkdir -p "$output_dir"

# Initialize conda for the current shell and activate the environment
eval "$(conda shell.bash hook)"  # For bash
conda activate GWK

# Collect 5' sequences only for PM reads (excluding redundant 3' sequence)
for file in ${input_dir}*.bowtie_mapped.parsed.txt; do
    txt=$(echo "${file}" | cut -d"/" -f3)
    echo "$txt"
    for spacer in SF370_spacer1 SF370_spacer2 SF370_spacer3 SF370_spacer4 SF370_spacer5 SF370_spacer6; do
        output="${output_dir}/${txt%.txt}.${spacer}.PM.5p.fastq"
        echo "Collecting 5' sequence for: $txt, Spacer: $spacer"
        awk -F'\t' -v spacer="$spacer" 'NR>1 && $1 == spacer && $2 == "PM" && $6 == "" {print "@"NR "\n" $4 "\n+\n"$5}' "$file" > "$output"
    done
done

# Q-score masking: retain reads with Q-score <= 30
for fq in ${output_dir}/*.PM.5p.fastq; do
    echo "Applying Q<=30 masking on: $fq"
    conda run -n crRNA-seq python ../scripts/fastq/qscore_masking.py -q 30 -i "$fq" -o "${fq/.fastq/.Q30.fastq}"
done

# Calculate distribution of 5' sequences (for reads with Q > 30)
for fq in ${output_dir}/*.PM.5p.Q30.fastq; do
    echo "Calculating 5' distribution for: $fq"
    awk 'NR%4==2 { if ($0 == "") print "no_5p_tail"; else print}' "$fq" | sort | uniq -c | sort -nr | awk -F' ' '{print $2"\t"$1}' > "${fq/.fastq/}.5p_distribution.txt"
done

# Merge all distributions for each spacer
for spacer in SF370_spacer1 SF370_spacer2 SF370_spacer3 SF370_spacer4 SF370_spacer5 SF370_spacer6; do
    echo "Merging 5' distribution for spacer: $spacer"
    conda run -n crRNA-seq python ../scripts/merge_and_lookup_table.manual.py -i ${output_dir}/*${spacer}*.5p_distribution.txt -o ${output_dir}/${spacer}.Q30.5p_distribution.merged.txt
done

