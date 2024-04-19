import os
import sys
import argparse

def process_files(file_paths):
    # Dictionary to store aggregated data by sample name
    sample_data = {}

    # Process each file
    for file_path in file_paths:
        if not os.path.exists(file_path):
            continue

        # Extract sample name from file basename
        filename = os.path.basename(file_path)
        if filename.startswith('ribodetector.') and filename.endswith('.log'):
            sample_name = filename[len('ribodetector.'):-len('.log')]

            # Initialize counts for non-rRNA and rRNA
            counts = {'non-rRNA': 0, 'rRNA': 0}

            # Read through lines in the file
            with open(file_path, 'r') as file:
                for line in file:
                    # example lines of interest
                    # 2024-04-18 13:45:51 : INFO  Detected 87139369 non-rRNA sequences
                    # 2024-04-18 13:45:51 : INFO  Detected 1060179 rRNA sequences
                    if 'Detected' in line:
                        parts = line.split('Detected')[1]
                        # print(parts)
                        # 87139369 non-rRNA sequences
                        
                        parts = parts.split()
                        # print(parts)
                        # ['87139369', 'non-rRNA', 'sequences']
                        
                        if len(parts) >= 2:
                            sequence_type = parts[1]
                            
                            # have had colour codes added to strings, need to remove^[[1m^[[96m87139369^[[0m non-rRNA sequences '\x1b[1m\x1b[96m87139369\x1b[0m'
                            parsed_count = parts[0].split("96m")[1].rstrip("\x1b[0m")
                            # print(f"Parsed count - {parsed_count}")
                            count = int(parsed_count)
                            if sequence_type in counts:
                                counts[sequence_type] += count

            # Store counts in the dictionary
            sample_data[sample_name] = counts

    # Calculate total and fractions
    for sample_name, counts in sample_data.items():
        total = counts['non-rRNA'] + counts['rRNA']
        fraction_non_rRNA = counts['non-rRNA'] / total if total > 0 else 0
        fraction_rRNA = counts['rRNA'] / total if total > 0 else 0
        counts.update({
            'total': total,
            'fraction_non-rRNA': fraction_non_rRNA,
            'fraction_rRNA': fraction_rRNA
        })

    return sample_data

def write_tsv(output_file, sample_data):
    # Write aggregated data to a TSV file
    with open(output_file, 'w') as file:
        # Write header
        file.write("sample_name\tnon-rRNA\trRNA\ttotal\tfraction_non-rRNA\tfraction_rRNA\n")

        # Write data rows
        for sample_name, counts in sample_data.items():
            file.write(f"{sample_name}\t{counts['non-rRNA']}\t{counts['rRNA']}\t{counts['total']}\t{counts['fraction_non-rRNA']}\t{counts['fraction_rRNA']}\n")


def main(file_list, output_file):

    # extract counts from log file - also calculates fractions of total for each sample
    log_read_counts = process_files(file_list)

    # write to TSV
    write_tsv(output_file, log_read_counts)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate non-rRNA and rRNA read counts from series of Ribodetector log files")

    # Positional argument for file paths
    parser.add_argument("file_list", nargs="+", help="List of file paths containing log files.")

    # Optional argument for output file
    parser.add_argument("-o", "--output_file", required=True, help="Output file name (TSV format).")

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    # Parse arguments from the command line
    args = parser.parse_args()

    # Call the main function with the parsed arguments
    main(args.file_list, args.output_file)



