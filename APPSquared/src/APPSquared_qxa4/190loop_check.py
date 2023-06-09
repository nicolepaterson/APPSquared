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
        _190loop_sites = []
        for line in parsed_file:
            match = re.search(r'193|194|195|196|197|198|199|200|201',line)
            #rbs = [98|134|135|136|137|138|153|155|183|190|194|221|222|223|224|225|226|227|228]
            #for i in rbs:
                #if item in parsed_file.read():
            if match:
                print(line)
                _190loop_sites.append(line)
    return _190loop_sites

_190loop_adj = compare_files(args.input)
#print(rbs_adj)
diff_df = pd.DataFrame(_190loop_adj)
diff_df.to_csv(args.output)

