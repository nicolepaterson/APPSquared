# APPSquared 
## Antigenic Prediction from Protein Sequences Pipeline 
![]
()

This pipeline generates biochemical analysis from HA and NA proteins relative to a reference vaccine strain


```usage for running H1 subtype analysis:
		bash run_H1_BCP.sh -r <RMSD ref> -d <directory> -n <analysis_name>
```

```usage for running H3 subtype analysis:
`
bash run_H3_BCP.sh -r <RMSD ref> -d <directory> -n <analysis_name>
```
```usage for running N2 subtype analysis:`
`
bash run_N2_BCP.sh -r <RMSD ref> -d <directory> -n <analysis_name>
```

```optional flags
-q turn off antibody dockings
-p turn off nearest neighbor calls 
-o turn off rosetta energy scores
```

`Set up getcontacts library git clone https://github.com/getcontacts/getcontacts.git 
echo "export PATH=`pwd`/getcontacts:\$PATH" >> ~/.bashrc source ~/.bashrc`

`dependencies, can be installed with conda.
numpy==1.23.5
pandas==1.5.3
biopython==1.78
dssp==3.0.0
scikit-learn==1.2.2
scipy==1.10.1
matplotlib=3.6.2
seaborn==0.12.2
cython==0.29.32
pip install ticc==0.1.4 
`
