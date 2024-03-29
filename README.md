# APPSquared 
## Antigenic Prediction from Protein Sequences Pipeline 
![APPSquared](https://github.com/nicolepaterson/APPSquared/blob/main/pipeline.png)

This pipeline generates biochemical analysis from HA and NA proteins relative to a reference vaccine strain


```usage for running H1 subtype analysis:
		bash run_H1_BCP.sh -r <RMSD ref> -d <directory> -n <analysis_name> -v <variant hash>
```

```usage for running H3 subtype analysis:
		bash run_H3_BCP.sh -r <RMSD ref> -d <directory> -n <analysis_name> -v <variant hash>
```
```usage for running N2 subtype analysis:
		bash run_N2_BCP.sh -r <RMSD ref> -d <directory> -n <analysis_name> -v <variant hash>
```
```usage for running Bvic subtype analysis:
		bash run_Bvic_BCP.sh -r <RMSD ref> -d <directory> -n <analysis_name> -v <variant hash>
```
```
set up necessary conda  environments: 
conda env create --name glyc --file=glyc.yaml
conda env create --name getcontacts --file=getcontacts.yaml
```
```
Set up getcontacts library:

git clone https://github.com/getcontacts/getcontacts.git 
echo "export PATH=`pwd`/getcontacts:\$PATH" >> ~/.bashrc source ~/.bashrc
```
