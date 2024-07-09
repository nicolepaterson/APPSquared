#!/usr/bin/env python3
import pandas as pd
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", required=True,
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
        parsed_file = pd.read_table(parsed_file,engine="python",delimiter="\s+",header=0)
        print(parsed_file)
        rosetta_df = pd.DataFrame(data=parsed_file).drop(columns=["pose_id","fa_atr","fa_rep","fa_sol","fa_intra_rep","fa_intra_sol_xover4","lk_ball_wtd","fa_elec","pro_close","fa_dun","p_aa_pp","yhh_planarity","ref","rama_prepro","description"])
            #columns=["pose_id","resi1","pdbid1","restype1","resi2","pdbid2","restype2","fa_atr","fa_rep","fa_sol","fa_intra_rep","fa_intra_sol_xover4","lk_ball_wtd","fa_elec","pro_close","hbond_sr_bb","hbond_lr_bb","hbond_bb_sc","hbond_sc","dslf_fa13","omega","fa_dun","p_aa_pp","yhh_planarity","ref","rama_prepro","total","description"]
        print(input)
    return rosetta_df

run_prep_df = prepare_(input)
variant_hash = args.variant_hash
timestamp  = args.timestamp

if variant_hash:
    run_prep_df['variant_hash'] = variant_hash
if timestamp:
    run_prep_df['timestamp'] = pd.to_datetime(timestamp)
run_prep_df.info(verbose=True)
out_file=f"/scicomp/groups/OID/NCIRD/ID/VSDB/GAT/cdp_archive/{timestamp}/{variant_hash}"+".per_res_scorefile.tsv"
run_prep_df.to_csv(out_file,sep ='\t',header=False,index=False)
print(run_prep_df)
