#!usr/bin/bash
# @nicolemariepaterson
# NA pipeline
usage() { echo "Usage: $0 -r H3 reference structure for rmsd in .mae format, -d output directory, -n name of new output directory "; exit 1; }
while getopts ":r:d:n:" arg; 
do
	case "${arg}" in
                r)
                        reference=$OPTARG
                        ;;
                d)
                        directory=$OPTARG
                        ;;
                n)
                        name=$OPTARG
                        ;;
                v)
                        variant_hash=$OPTARG
                        ;;
	esac
done

BCPhome=/scicomp/groups/OID/NCIRD/ID/VSDB/GAT
#schrodingerhome=/scicomp/groups/OID/NCIRD/ID/VSDB/GAT/schrodinger

#ml conda
#make the directory for storing the data
mkdir $name
#cd $directory/$name

ml schrodinger

if ! [ -f $filename.mae ] ; then
    for file in $directory/prep*.pdb
    do
        filename=$(basename $file);
        run pdbconvert -ipdb $filename -omae $name/$filename.mae
        #run utilities/prepwizard $filename prep_$filename
    done
fi

#Runs AlphaFold on sequence fasta file
#run Rosetta energy score for ddG

module load rosetta/3.13-cluster
#bash run_rosetta_scores.sh -d $reference -n $name
#for file in $directory/prep*.pdb
#do
#    bash run_rosetta_scores.sh -d $directory -n $name
#done

# Compare SNPs to DMS data

#Generate a control file for antibody docks
#if ! [ $antibody_dock_off ] ; then
#    for file in $name/*.mae
#    do
#        run -FROM psp piper.py -jobname $filename_dock -receptor prep_1G01.pdb -receptor_chain J,K -poses 10 -rotations 70000 -OMPI 1 -JOBID -antibody -ligand $filename -ligand_chain A -use_nonstandard_residue yes -HOST localhost:1 -TMPLAUNCHDIR
#    done
#fi
#Computes the glycosylation distances for each site from pdb file
ml conda
conda activate glyc
bash run_glyc.sh -i $directory -v $variant_hash -s N2 -p NA
#bash run_glyc_diff.sh -i $directory -r $reference
#python glyc_score.py -i $directory -r $reference

conda activate getcontacts
for file in $directory/prep*.pdb
do
    filename=$(basename "$file")
    echo $filename
    #sudo /home/qxa4/schrodinger/utilities/prepwizard $filename prep_$filename
    python get_static_contacts.py --structure $filename --output $filename.NAstatic_contacts.csv --itypes all
    python get_contact_frequencies.py --input_files $filename.NAstatic_contacts.csv --output_file $filename.NAfreq.tsv --itypes all
    python contact_diffs.py -i $filename.NAstatic_contacts.csv -r $reference.NAstatic_contacts.csv -o $filename.contactdiffs.NAstatic.csv
    python get_contact_flare.py --input $filename.NAstatic_contacts.csv --out $filename.NAcontact_flare.json
done

for file in $name/*.pdb
do
    filename=$(basename "$file")
    python run_sasa.py -r $reference -p $filename -d $name.sasa.csv
done
