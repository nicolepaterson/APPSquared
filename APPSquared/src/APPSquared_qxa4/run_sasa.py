#!bin/env python
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-r", "--refpdb", dest="refpdb", required=True,help="ref for calculating delta SASA")
parser.add_argument("-p", "--testpdb",dest="testpdb", required=True,help="test for generating delta SASA")
parser.add_argument("-d", "--deltaSASA",dest="dSASA", required=False,help="delta SASA value file")
args = parser.parse_args()
p = PDBParser(QUIET=1)
struct = p.get_structure(args.refpdb,args.refpdb)
struct2 = p.get_structure(args.testpdb,args.testpdb)
sr = ShrakeRupley()
sr.compute(struct, level="S")
sr.compute(struct2, level="S")
calc = round(struct.sasa, 2)
calc2 = round(struct2.sasa, 2)
print(calc-calc2)
dSASA = calc-calc2
sasa1 = round(struct.sasa, 2)
sasa2 = round(struct2.sasa, 2)
list = f'isolate {args.dSASA},SASA1 {sasa1},SASA2 {sasa2},dSASA {dSASA}'
f = open(args.dSASA,"w")
print(list)
f.write(list)
