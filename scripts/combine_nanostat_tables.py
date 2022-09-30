#!/usr/bin/env python3

import pandas as pd
import os
import sys

'''
Quick and dirty script to combine NanoStat output TSVs into a single file, adding a sample identifier as a column
'''

def main(in_dir, suffix, out_file):

    # find file names that match suffix
    samples2path = {f.replace(suffix, ""): os.path.join(in_dir, f)
                    for f in os.listdir(in_dir) if f.endswith(suffix)
                    }

    print(samples2path)

    assert len(samples2path) != 0, f"No files found in {in_dir} with suffix {suffix}"

    # Read in TSVs into dataframes
    dfs = {sample: pd.read_csv(path, sep="\t", header=None)
           for sample, path in samples2path.items()}

    # Concat into a single df, appending the sample name as a column
    df = pd.concat(dfs.values(),
                   keys=dfs.keys(),
                   names=["sample_name", None]).reset_index("sample_name").reset_index(drop=True)

    df.sort_values(by="sample_name", inplace=True)

    df.to_csv(out_file, sep="\t", header=True, index=False, na_rep="NA")

if __name__ == '__main__':

    help = """python combine_nanostat_tables.py INPUT_DIR FILE_SUFFIX OUT_FILE
              INPUT_DIR - path to input directory containing per-sample NanoStat summary TSVs to combine
              FILE_SUFFIX - string suffix to identify files containing NanoStat summary TSVs' for a given file (e.g. <sample_name><FILE_SUFFIX>)
              OUT_FILE - name of output table of combined NanoStat outputs for all identified files
             """

    if len(sys.argv) == 1:
        print(help)
        sys.exit(0)

    elif "-h" in sys.argv or "--help" in sys.argv:
        print(help)
        sys.exit(0)

    main(sys.argv[1], sys.argv[2], sys.argv[3])
