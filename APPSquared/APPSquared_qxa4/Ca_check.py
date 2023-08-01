#!bin/env/python
import argparse
import os
import pandas as pd
import re

parser = argparse.ArgumentParser()

parser.add_argument("-i", dest="input", help="input file generated by running contact_diffs.py", required=True)
parser.add_argument("-o", dest="output", help="output file name", required=True)
args = parser.parse_args()

def compare_files(test_file):
    with open(test_file) as parsed_file:
        Ca_sites = []
        for line in parsed_file:
            #match = re.search(r'501',line)
            match = re.search(r'134|135|137|138|139|140|141|142|150|166|167|168|169|170|177|203|204|205|215|221|222|235|236|237|239|240|246|248|249',line)
            if match:
                print(line)
                Ca_sites.append(line)
    return Ca_sites

Ca_adj = compare_files(args.input)
#print(rbs_adj)
diff_df = pd.DataFrame(Ca_adj)
diff_df.to_csv(args.output)
