import warnings
from Bio import BiopythonWarning
from Bio.PDB import PDBParser, Selection
import pandas as pd
import argparse
from pathlib import Path

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_pdb_file", required=True,help="Input PDB file")

    parser.add_argument("-v", "--variant_hash", required=True,
        help="Alphanumeric character STRING linked to a specific amino acid sequence")

    parser.add_argument("-t", "--timestamp_version", required=True,
        help="System TIMESTAMP applied to distinct analytical applications")

    parser.add_argument("-o", "--output_csv_file", required=True,
        help="""Output csv file. If no argument given, then csv file with the same name as 
              the input is written to same directory""")

    return parser

def ignore_pdb_warnings():
    """
    Supress the biopython pdb import warnings
    """
    warnings.simplefilter("ignore", BiopythonWarning)


args = get_parser().parse_args()
pdb_file            = args.input_pdb_file
pdb_path            = Path(args.input_pdb_file)
variant_hash        = args.variant_hash
timestamp_version   = args.timestamp_version
if args.output_csv_file:
    output_csv_file = Path(args.output_csv_file)
else:
    output_csv_file = Path(pdb_file.parent.joinpath(str(pdb_file)+"raw_pdb.tsv"))

def make_PDB_df(structure):
    """
    Create a dataframe of residue properties contained with the PDB structure input
    """
    structure_list = []
    for atom in Selection.unfold_entities(structure, "A"):
        atom_name = str(atom.id)
        if atom_name == "CA":
            resnum = str(atom.get_parent().id[1])
            resname = str(atom.get_parent().resname)
            bfactor = float(atom.bfactor)
            if atom.get_parent().id[2] ==  ":":
                insertion_code = " "
            else:
                insertion_code = atom.get_parent().id[2]
            chain = str(atom.get_parent().get_parent().id)
            node_id = chain + ":" + resnum + ":" + insertion_code + ":" + resname
            structure_list.append([node_id, atom.coord[0], atom.coord[1], atom.coord[2], bfactor])

    structure_df = pd.DataFrame(structure_list, columns=["node_id", "x_coordinate", "y_coordinate", "z_coordinate", "b_factor"])
    structure_df["chain"] = structure_df.node_id.str.split(":", expand=True)[0]
    structure_df["site"] = structure_df.node_id.str.split(":", expand=True)[1].astype(int)
    structure_df["insertion_code"] = structure_df.node_id.str.split(":", expand=True)[2]
    structure_df["residue"] = structure_df.node_id.str.split(":", expand=True)[3]

    return structure_df

ignore_pdb_warnings()
parser = PDBParser()
pdb_file            = args.input_pdb_file
structure = parser.get_structure(pdb_file,pdb_path)
structure_df = make_PDB_df(structure)
print(structure_df)

if variant_hash:
    structure_df["variant_hash"] = variant_hash
if timestamp_version:
    structure_df["version"] = pd.to_datetime(timestamp_version)

structure_df.to_csv(output_csv_file, sep ="\t",header=False,index=False)
print(structure_df.dtypes)
