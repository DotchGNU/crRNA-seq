import argparse
import pandas as pd

def merge_and_lookup(input_files, output_file):
    # Initialize an empty DataFrame to store the merged data
    merged_data = pd.DataFrame()

    # Iterate over the input files
    for file in input_files:
        # Extract the file name without extension
        file_name = file.split('.')[0]

        # Read the data from the current file
        data = pd.read_csv(file, sep='\s+',  names=['value', file_name])
        #print(file, data.columns)

        # Merge the current data with the existing merged data
        if merged_data.empty:
            merged_data = data
        else:
            merged_data = pd.merge(merged_data, data, on='value', how='outer')

    # Save the merged data to the output file
    merged_data.fillna(0, inplace=True) # remove Nan value
    merged_data.to_csv(output_file, index=False, sep="\t")

    # Display the merged data
    #print("Merged data:")
    #print(merged_data)

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Merge and lookup data from multiple files.")
    parser.add_argument("-i", "--input_files", nargs='+', help="List of input files", required=True)
    parser.add_argument("-o", "--output_file", help="Output file name", required=True)
    args = parser.parse_args()

    # Call the merge_and_lookup function with the provided arguments
    merge_and_lookup(args.input_files, args.output_file)

