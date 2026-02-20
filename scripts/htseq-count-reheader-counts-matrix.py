#!/usr/bin/env python3
"""
Reformat htseq-count output headers to replace empty fields with attribute names
and simplify BAM file paths to sample names.
"""

import argparse
import sys
from pathlib import Path


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Reformat htseq-count output headers',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  %(prog)s --idattr gene_id --additional-attr gene_name --additional-attr gene_type \\
           --bam-suffix .chr8.bam input.txt -o output.txt
  
  %(prog)s --idattr gene_id --additional-attr gene_name --additional-attr gene_type \\
           --add-chromosome-info --bam-suffix .chr8.bam input.txt -o output.txt
        """
    )
    
    parser.add_argument(
        'input',
        help='Input htseq-count file (use - for stdin)'
    )
    
    parser.add_argument(
        '-o', '--output',
        default='-',
        help='Output file (default: stdout)'
    )
    
    parser.add_argument(
        '--idattr',
        required=True,
        help='ID attribute name (e.g., gene_id)'
    )
    
    parser.add_argument(
        '--additional-attr',
        action='append',
        default=[],
        dest='additional_attrs',
        help='Additional attribute names in order (can be specified multiple times)'
    )
    
    parser.add_argument(
        '--bam-suffix',
        default='',
        help='Suffix to remove from BAM filenames (e.g., .chr8.bam)'
    )
    
    parser.add_argument(
        '--add-chromosome-info',
        action='store_true',
        help='Indicates htseq-count was run with --add-chromosome-info (adds chromosome column)'
    )
    
    return parser.parse_args()


def reformat_header(header_line, idattr, additional_attrs, bam_suffix, add_chromosome_info):
    """
    Reformat the header line from htseq-count output.
    
    Args:
        header_line: Original header line (tab-separated)
        idattr: ID attribute name
        additional_attrs: List of additional attribute names
        bam_suffix: Suffix to remove from BAM filenames
        add_chromosome_info: Whether chromosome info column is present
    
    Returns:
        Reformatted header line
    """
    fields = header_line.rstrip('\n').split('\t')
    
    # Calculate expected number of empty fields
    expected_empty = 1 + len(additional_attrs)  # idattr + additional attrs
    if add_chromosome_info:
        expected_empty += 1  # chromosome column
    
    # Count actual empty fields at the beginning
    empty_count = 0
    for field in fields:
        if field == '':
            empty_count += 1
        else:
            break
    
    # Validate that we have the right number of empty fields
    if empty_count != expected_empty:
        if add_chromosome_info:
            raise ValueError(
                f"Expected {expected_empty} empty fields (1 idattr + "
                f"{len(additional_attrs)} additional attrs + 1 chromosome) but found {empty_count}"
            )
        else:
            raise ValueError(
                f"Expected {expected_empty} empty fields (1 idattr + "
                f"{len(additional_attrs)} additional attrs) but found {empty_count}"
            )
    
    # Build new header: attribute names + chromosome (if applicable) + reformatted sample names
    new_fields = [idattr] + additional_attrs
    
    if add_chromosome_info:
        new_fields.append('chromosome')
    
    # Process BAM file paths
    for bam_path in fields[empty_count:]:
        # Extract basename
        basename = Path(bam_path).name
        
        # Remove suffix if provided
        if bam_suffix:
            if not basename.endswith(bam_suffix):
                raise ValueError(
                    f"Basename '{basename}' does not end with the provided "
                    f"bam_suffix '{bam_suffix}'"
                )
            sample_name = basename[:-len(bam_suffix)]
        else:
            sample_name = basename
        
        new_fields.append(sample_name)
    
    return '\t'.join(new_fields) + '\n'


def main():
    """Main function."""
    args = parse_args()
    
    # Open input file
    if args.input == '-':
        input_file = sys.stdin
    else:
        try:
            input_file = open(args.input, 'r')
        except IOError as e:
            sys.stderr.write(f"Error opening input file: {e}\n")
            sys.exit(1)
    
    # Open output file
    if args.output == '-':
        output_file = sys.stdout
    else:
        try:
            output_file = open(args.output, 'w')
        except IOError as e:
            sys.stderr.write(f"Error opening output file: {e}\n")
            if input_file != sys.stdin:
                input_file.close()
            sys.exit(1)
    
    try:
        # Process first line (header)
        first_line = input_file.readline()
        if not first_line:
            sys.stderr.write("Error: Input file is empty\n")
            sys.exit(1)
        
        try:
            new_header = reformat_header(
                first_line,
                args.idattr,
                args.additional_attrs,
                args.bam_suffix,
                args.add_chromosome_info
            )
            output_file.write(new_header)
        except ValueError as e:
            sys.stderr.write(f"Error processing header: {e}\n")
            sys.exit(1)
        
        # Copy remaining lines unchanged
        for line in input_file:
            output_file.write(line)
    
    finally:
        # Close files if not stdin/stdout
        if input_file != sys.stdin:
            input_file.close()
        if output_file != sys.stdout:
            output_file.close()


if __name__ == '__main__':
    main()