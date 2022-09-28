#!/usr/bin/env python3

import pandas as pd
import os
import sys

'''
Quick and dirty script to combine extracted 'SN' sections of samtools stats output into a single file
'''

def main(in_dir, suffix, out_file):

    # find file names that match suffix
    samples2path = {f.replace(suffix, ""): os.path.join(in_dir, f)
                    for f in os.listdir(in_dir) if f.endswith(suffix)
                    }

    print(samples2path)

    assert len(samples2path) != 0, f"No files found in {in_dir} with suffix {out_file}"

    # Read in TSVs into dataframes
    dfs = {sample: pd.read_csv(path, sep="\t", header=None, names=["metric", "value", "comment"])
           for sample, path in samples2path.items()}

    # Concat into a single df, appending the sample name as a column
    df = pd.concat(dfs.values(),
                   keys=dfs.keys(),
                   names=["sample_name", None]).reset_index("sample_name").reset_index(drop=True)

    df.sort_values(by="sample_name", inplace=True)

    df.to_csv(out_file, sep="\t", header=True, index=False, na_rep="NA")

if __name__ == '__main__':

    help = """python combine_stats_sn_tables.py INPUT_DIR FILE_SUFFIX OUT_FILE
              INPUT_DIR - path to input directory containing 'SN' samtools stats sections to combine
              FILE_SUFFIX - string suffix to identify files containing 'SN' samtools stats output for a given file (i.e. output extracted with grep ^SN <samtools_stats.txt> | cut -f 2-)
              OUT_FILE - name of output table of combined SN outputs for all identified files
             """

    if len(sys.argv) == 1:
        print(help)
        sys.exit(0)

    elif "-h" in sys.argv or "--help" in sys.argv:
        print(help)
        sys.exit(0)

    main(sys.argv[1], sys.argv[2], sys.argv[3])
