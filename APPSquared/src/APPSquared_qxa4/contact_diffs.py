#!bin/env/python

import argparse
import os
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument("-i", dest="input", help="input frequency file generated by running getcontacts", required=True)
parser.add_argument("-r", dest="reference", help="reference frequency file generated by running getcontacts",required=True)
parser.add_argument("-o", dest="output", help="output file name", required=True)
args = parser.parse_args()

def compare_files(reference, test_file):
    with open(test_file) as parsed_file:
        diffs = []
        for line in parsed_file:
            reference_file = open(reference)
            if line in reference_file.read():
                continue
            else:
                print(line)
                diffs.append(line)
    return diffs

differences = compare_files(args.reference, args.input)
print(differences)
diff_df = pd.DataFrame(differences)
diff_df.to_csv(args.output, index=False)
