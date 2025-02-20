# crRNA-seq

## Overview

This pipeline processes crRNA-Seq using Bowtie, regular expressions, and several Python/shell scripts. The main steps include adapter trimming, crRNA repeat search, spacer length analysis, reverse mapping, and 5' sequence distribution analysis. To ensure reproducibility, the required conda environment YAML files are provided. All scripts use relative paths according to the project structure.

## Directory Structure

```
crRNA_pipeline/
├── LICENSE
├── README.md					# This file
├── main.sh                                     # Main pipeline bash script 
├── data					# Reference files of crRNA spacer
│   ├── SF370_spacer_20nt.fa
│   ├── SF370_spacer_20nt_with_deletion.fa
│   └── SF370_spacer_30nt.fa
├── env						
│   ├── crRNA-seq.yml				# Conda environment YAML file for crRNA-seq
│   └── reverse_mapping.yml			# Conda environment YAML file for reverse_mapping
├── fastq
│   └── readme.txt				# You can get example fastq file using wget
├── output_example.tar.gz			# You can check out example results 
└── scripts
    ├── Reverse_mapping
    │   ├── FastQcollapse.py
    │   ├── Reverse_mapping_with_bowtie.py
    │   ├── count_RNA.py
    │   ├── extract_bowtie_unmapped_reads.biopython.py
    │   └── parse_bowtie_result.py
    ├── collect_5p_sequence_and_calculate_histogram_for_q30.SF370.sh
    ├── fastq
    │   ├── filter_by_sequence_fastq.py
    │   ├── filter_fastq_by_read_length.py
    │   └── qscore_masking.py
    ├── grep_count_by_length_using_regex.py
    ├── merge_and_lookup_table.manual.py
    └── merge_and_lookup_table.spacer_regex.py
```

## Requirements

- **Operating System:** Linux/Mac (Bash shell)
- **Conda & Python:** Python 2.7.14 (reverse_mapping) & 3.7.12 (crRNA-seq) with Conda installed
- **Tools:** cutadapt 4.1, bowtie 1.2.3, and other tools as required (should be included in the conda environments)
- **Dependencies:** Python libraries and other dependencies are managed within the conda environments.

## Installation & Setup

1. **Clone the Repository and Navigate to the Project Directory (~1min)**

   ```bash
   git clone https://github.com/DotchGNU/crRNA-seq.git
   cd crRNA-seq
   ```

2. **Create Conda Environments (~10min)**

   Use the provided YAML files to create the necessary conda environments:

   ```bash
   conda env create -f env/crRNA-seq.yml
   conda env create -f env/reverse_mapping.yml
   ```

   > **Note:** Adjust or update the environment files as needed.


3. **Prepare Data**

   The pipeline expects FASTQ files to be located in the `./fastq` directory. Place your FASTQ files in `./fastq` or update the FASTQ directory path in `main.sh` accordingly.

   **Example Input File:**
   If you need a sample FASTQ file for testing, you can download one using `wget`:
   ```bash
   mkdir -p fastq
   cd fastq
   wget http://ago.korea.ac.kr/crispr_abasic/crRNA_sequencing_SF370/SF370_1a.R1.fastq.gz
   cd ..
   ```

## Usage

**Run the Entire Pipeline**

   From the project root directory, execute the main pipeline script:

   ```bash
   bash main.sh
   ```

   You can check out the example results in `output_example.tar.gz` 

**Expected Run Time**

    This project was tested on a Mac Pro with the following specs:

    - OS: macOS High Sierra (Darwin 17.5.0)
    - CPU: 12-Core Intel Xeon E5 @ 2.7 GHz
    - Memory: 64 GB

When running the `example_input` file under these conditions, the total execution time is **approximately 5 minutes**.  

## Additional Information

- **Relative Paths:**\
  All scripts use relative paths based on the project root (`crRNA-seq/`). Ensure the folder structure remains unchanged or update paths as necessary.

- **Dependency Files:**\
  Files such as `qscore_masking.py` and `merge_and_lookup_table.manual.py` must exist in the appropriate locations (e.g., `./scripts/fastq/` or `./scripts/`). Verify and adjust the paths in the scripts if needed.

## Troubleshooting

- **File Path Errors:**\
  Ensure that all required files exist in the specified directories. Verify the relative paths in the scripts if you encounter errors.

- **Dependency Problems:**\
  Confirm that all required Python packages are installed within the respective conda environments. Adjust the YAML files if necessary.

## License & Contact

- **License:**\
This software is licensed under the MIT License. See the LICENSE file for details.

- **Contact:**\
  For questions or bug reports, please refer to [gwpia0409@korea.ac.kr].


