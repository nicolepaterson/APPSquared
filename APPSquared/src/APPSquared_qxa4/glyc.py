from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
import pandas as pd
from pypackage.struct.utils import nodeId_column
import sys
from sklearn.neighbors import radius_neighbors_graph

print(sys.argv[1])
input=(sys.argv[1])
parser = PDBParser()
structure= parser.get_structure(input, sys.argv[1])
io=PDBIO(structure)
io.set_structure(structure)

def make_PDB_df(structure):
    from Bio.PDB import Selection
    """
    Create BioPython dataframe of atom properties
    """
    structure_list = []
    for atom in Selection.unfold_entities(structure, 'A'):
        atom_name = str(atom.id)
        if atom_name == 'CA':
            resnum = str(atom.get_parent().id[1])
            resname = str(atom.get_parent().resname)
            bfactor = float(atom.bfactor)
            if atom.get_parent().id[2] == ' ':
                insertion_code = '_'
            else:
                insertion_code = atom.get_parent().id[2]
            chain = str(atom.get_parent().get_parent().id)
            nodeid = chain + ':' + resnum + ':' + insertion_code + ':' + resname
            structure_list.append([nodeid, atom.coord[0], atom.coord[1], atom.coord[2], bfactor])

    structure_df = pd.DataFrame(structure_list, columns=['NodeId', 'x', 'y', 'z', 'B_factor'])
    structure_df['Chain'] = structure_df.NodeId.str.split(':', expand=True)[0]
    structure_df['ResNum'] = structure_df.NodeId.str.split(':', expand=True)[1].astype(int)
    structure_df['InsertionCode'] = structure_df.NodeId.str.split(':', expand=True)[2]
    structure_df['Residue'] = structure_df.NodeId.str.split(':', expand=True)[3]
    return structure_df
structure_dff = make_PDB_df(structure)
print(structure_dff)
print(structure_dff['Residue'])

def add_glycosylation_site_col(df, column, suffix=None):
    if suffix:
        suffix = ' '+suffix
    else:
        suffix = ''
    df.loc[(df[column] == "ASN") &
           (df[column].shift(-1) != "PRO") &
           (df[column].shift(-2).isin(["SER", "THR"])), f'Glycosylation Site{suffix}'] = 1
    df.loc[(df[column] != "PRO") &
           (df[column].shift(1) == "ASN") &
           (df[column].shift(-1).isin(["SER", "THR"])), f'Glycosylation Site{suffix}'] = 2
    df.loc[(df[column].isin(["SER", "THR"])) &
           (df[column].shift(1) != "PRO") &
           (df[column].shift(2) == "ASN"), f'Glycosylation Site{suffix}'] = 3
    df.fillna(False, inplace=True)
    print(df)
glyc_df = add_glycosylation_site_col(structure_dff,"Residue")
print(glyc_df)
#output.to_csv(glyc_df)

def calculate_distance_to_nearest_glycosylation(df, gly_col='Glycosylation Site', output_col='Distance to Glycosylation', radius=20):
    model_coords = df[['x', 'y', 'z']].values
    gly_sites = df.loc[(df[gly_col] == 1)].copy()
    gly_sites['GlySite'] = 'Distance_to_N' + gly_sites['ResNum'].astype(str) + ':' + gly_sites['Chain']
    cols_rename_dict = gly_sites['GlySite'].to_dict()
    m = radius_neighbors_graph(model_coords, 20, mode='distance', include_self=True)
    distance_matrix = pd.DataFrame(m.toarray()).replace(0, 20)
    distance_matrix.rename(columns=cols_rename_dict, inplace=True)
    distance_matrix = distance_matrix[list(cols_rename_dict.values())].copy()
    df = df.merge(distance_matrix, left_index=True, right_index=True)
    glycosylation_cols = [col for col in df if 'Distance_to_N' in col]
    for gly_col_name in glycosylation_cols:
        df.loc[(df.ResNum == int(gly_col_name.split(':')[0][13:])) & (df.Chain == gly_col_name[-1]), gly_col_name] = 0.0
    df[output_col] = df[glycosylation_cols].min(axis=1)
    return df
    print(df)
gly_dist2 = calculate_distance_to_nearest_glycosylation(structure_dff)
#gly_dist = calculate_distance_to_nearest_glycosylation(glyc_df)
print(gly_dist2)

def add_dssp_cols(structure,file,structure_dff):
    from Bio.PDB import DSSP
    from Bio.PDB.DSSP import residue_max_acc
    dssp_data = DSSP(structure[0], file) # , dssp=dssp_loc
    dssp_df = pd.DataFrame([tuple(x[0]) + tuple(x[1]) + y for x, y in list(zip(dssp_data.keys(), list(dssp_data)))], columns=['Chain', 'HetFlag', 'ResNum', 'InsertionCode','Index', 'Residue', 'SS', 'RelASA', 'Phi', 'Psi','NH-->O_1_relidx', 'NH-->O_1_energy','O-->NH_1_relidx', 'O-->NH_1_energy','NH-->O_2_relidx', 'NH-->O_2_energy','O-->NH_2_relidx', 'O-->NH_2_energy'])
    aa_mapper = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS','I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN','G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP','A': 'ALA', 'V': 'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}
    insertioncode_mapper = {' ': '_'}
    dssp_df = dssp_df.replace({"InsertionCode": insertioncode_mapper})
    dssp_df = dssp_df.replace({"Residue": aa_mapper})
    # quick change, not sure if function works
    dssp_df['NodeId'] = nodeId_column(dssp_df)
    #dssp_df['NodeId'] = dssp_df['Chain'] + ':' + dssp_df['ResNum'].astype('str') + ':' + \
    #                            dssp_df['InsertionCode'] + ':' + dssp_df['Residue']
    dssp_df['maxASA'] = dssp_df['Residue'].map(residue_max_acc['Sander'])
    dssp_df['ASA'] = dssp_df['maxASA'] * dssp_df['RelASA']
    dssp_df['ResNum'] = dssp_df['ResNum'].astype('int')
    dssp_df = dssp_df[['NodeId', 'SS', 'RelASA', 'ASA']]
    structure_df = pd.merge(structure_dff, dssp_df, on='NodeId', how='left')
    return structure_df
ASA = add_dssp_cols(structure,input,gly_dist2)
output1= pd.DataFrame(ASA)
#output2= pd.DataFrame(ASA)
filename1 = input +"gly_ASA_dist.csv"
#filename2 = input +"ASA_dist.csv"
output1.to_csv(filename1)
#output2.to_csv(filename2)