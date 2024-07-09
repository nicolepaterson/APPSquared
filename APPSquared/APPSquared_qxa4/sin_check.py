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
        rbs_sites = []
        for line in parsed_file:
            #match = re.search(r'501',line)
            match = re.search(r'98|153|155|183|190|194|226',line)
            #rbs = [114,115,116,117,118,122,123,164,165,166,167,168,169,234]
            #for i in rbs:
                #if item in parsed_file.read():
            if match:
                print(line)
                rbs_sites.append(line)
    return rbs_sites

rbs_adj = compare_files(args.input)
#print(rbs_adj)
diff_df = pd.DataFrame(rbs_adj)
diff_df.to_csv(args.output)
