#!/usr/bin/env python3
import pandas as pd
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", required=True,
        help="Input PDB file")
parser.add_argument("-r", "--reference", required=True,
        help="Input PDB file")
parser.add_argument("-v", "--variant_hash", required=False,
        help="Alphanumeric character STRING linked to a specific amino acid sequence")
parser.add_argument("-t", "--timestamp", required=False,
        help="System TIMESTAMP applied to distinct analytical applications")
args = parser.parse_args()
print(args.input_file)
input=(args.input_file)

def prepare_(test_file):
    with open(test_file) as parsed_file:
        parsed_file = pd.read_table(parsed_file,engine="python",delimiter="\s+",header=None)
        #print(parsed_file)
        df = pd.DataFrame(data=parsed_file)
        df2=df[[1,2,3,4,5,6,7,8]]
        #"Residue-1","Residue-2","Residue-Id","Residue-Id","Missing atoms(T/F)","RMSD(residue-wise by asl)","Difference in average b-factor (residue-wise)"])
    return df2

run_prep_df = prepare_(input)
variant_hash = args.variant_hash
timestamp  = args.timestamp
reference = args.reference

if reference:
    run_prep_df['reference'] = reference
if variant_hash:
    run_prep_df['vhash'] = variant_hash
if timestamp:
    run_prep_df['timestamp'] = pd.to_datetime(timestamp)
out_file=f"/scicomp/groups/OID/NCIRD/ID/VSDB/GAT/cdp_archive/{timestamp}/{variant_hash}"+".per_res_rmsd.tsv"
run_prep_df.to_csv(out_file,sep ='\t',header=False,index=False)

