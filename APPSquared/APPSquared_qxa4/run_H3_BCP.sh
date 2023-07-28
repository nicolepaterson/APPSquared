#!usr/bin/bash
# @nicolemariepaterson
# H3 pipeline


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
                o)
                        name=$rosetta
                        ;;
                p)
                        name=$nearest_neighbor
                        ;;
                q)
                        name=$antibody_dock_off
                        ;;
        esac
done

BCPhome=/scicomp/groups/OID/NCIRD/ID/VSDB/GAT
#schrodingerhome=/scicomp/groups/OID/NCIRD/ID/VSDB/GAT/schrodinger

#ml conda
#make the directory for storing the data
mkdir $directory/$name


#run Rosetta energy score for ddG

module load rosetta/3.13-cluster
if [ -e "-o" ]; then
    bash run_rosetta_scores.sh -d $reference -n $directory/$name/$name
#for file in $directory/prep*.pdb
#bash run_rosetta_scores.sh -d $directory -n $name
fi
#Computes the glycosylation distances for each site from pdb file
ml conda
    conda activate glyc
    bash run_glyc.sh -i $directory
#bash run_glyc_diff.sh -i $directory -r $reference
#python glyc_score.py -i $directory -r $reference

ml conda
conda activate getcontacts
for file in $directory/prep*.pdb
do
    filename=$(basename "$file")
    echo $filename
    #sudo /home/qxa4/schrodinger/utilities/prepwizard $filename prep_$filename
    python get_static_contacts.py --structure $filename --output $directory/$name/$filename.H3static_contacts.csv --itypes all
    sleep 5
    python get_contact_frequencies.py --input_files $directory/$name/$filename.H3static_contacts.csv --output_file $directory/$name/$filename.H3freq.tsv --itypes all
    sleep 5
    python contact_diffs.py -i $directory/$name/$filename.H3static_contacts.csv -r $reference.H3static_contacts.csv -o $directory/$name/$filename.H3contactdiffs.static.csv
    sleep 5
    python get_contact_flare.py --input $directory/$name/$filename.H3static_contacts.csv --out $directory/$name/$filename.H3contact_flare.json
    sleep 5
    python H3_rbs_check.py -i $directory/$name/$filename.H3contactdiffs.static.csv -o $directory/$name/$filename.H3rbs.csv
    sleep 5
    python sin_check.py -i $directory/$name/$filename.H3contactdiffs.static.csv -o $directory/$name/$filename.H3_SIN.csv
    python siteA_check.py -i $directory/$name/$filename.H3contactdiffs.static.csv -o $directory/$name/$filename.H3siteA.csv
    python siteB_check.py -i $directory/$name/$filename.H3contactdiffs.static.csv -o $directory/$name/$filename.H3siteB.csv
    python siteC_check.py -i $directory/$name/$filename.H3contactdiffs.static.csv -o $directory/$name/$filename.H3siteC.csv
    python siteD_check.py -i $directory/$name/$filename.H3contactdiffs.static.csv -o $directory/$name/$filename.H3siteD.csv
    python siteE_check.py -i $directory/$name/$filename.H3contactdiffs.static.csv -o $directory/$name/$filename.H3siteE.csv
    python H3_rbs_check.py -i $directory/$name/$filename.H3siteA.csv -o $directory/$name/$filename.H3rbsnet_antA.csv
    python H3_rbs_check.py -i $directory/$name/$filename.H3siteB.csv -o $directory/$name/$filename.H3rbsnet_antB.csv
    python H3_rbs_check.py -i $directory/$name/$filename.H3siteC.csv -o $directory/$name/$filename.H3rbsnet_antC.csv
    python H3_rbs_check.py -i $directory/$name/$filename.H3siteD.csv -o $directory/$name/$filename.H3rbsnet_antD.csv
    python H3_rbs_check.py -i $directory/$name/$filename.H3siteE.csv -o $directory/$name/$filename.H3rbsnet_antE.csv
    #sed -i "s/:/\,/g" $filename.contactdiffs.static.csv
    sed -i "s/\"//g" $directory/$name/$filename.H3contactdiffs.static.csv
    tail -n +2 $directory/$name/$filename.H3contactdiffs.static.csv
    sed '1i zero,bond_type,chain1,chain2' $directory/$name/$filename.H3contactdiffs.static.csv
    sed -i "s/\"//g" $directory/$name/$filename.H3rbs.csv
    tail -n +2 $directory/$name/$filename.H3rbs.csv
    sed '1i zero,bond_type,chain1,chain2' $directory/$name/$filename.H3rbs.csv
    sed -i "s/\"//g" $directory/$name/$filename.H3siteA.csv
    tail -n +2 $directory/$name/$filename.H3siteA.csv
    sed '1i zero,bond_type,chain1,chain2' $directory/$name/$filename.H3siteA.csv
    sed -i "s/\"//g" $directory/$name/$filename.H3siteB.csv
    tail -n +2 $directory/$name/$filename.H3siteB.csv
    sed '1i zero,bond_type,chain1,chain2' $directory/$name/$filename.H3siteB.csv
    sed -i "s/\"//g" $directory/$name/$filename.H3siteC.csv
    tail -n +2 $directory/$name/$filename.H3siteC.csv
    sed '1i zero,bond_type,chain1,chain2' $directory/$name/$filename.H3siteC.csv
    sed -i "s/\"//g" $directory/$name/$filename.H3siteD.csv
    tail -n +2 $directory/$name/$filename.H3siteD.csv
    sed '1i zero,bond_type,chain1,chain2' $directory/$name/$filename.H3siteD.csv
    sed -i "s/\"//g" $directory/$name/$filename.H3siteE.csv
    tail -n +2 $directory/$name/$filename.H3siteE.csv
    sed '1i zero,bond_type,chain1,chain2' $directory/$name/$filename.H3siteE.csv
done
cp /scicomp/home-pure/qxa4/BCP_V2/count_contacts_compile.sh $directory/$name/count_contacts_compile.sh
bash $directory/$name/count_contacts_compile.sh -d $directory/$name -n $filename"_contacts_compile.txt"
#do
#    filename=$(basename "$file")
#    python run_sasa.py -r $reference -p $filename -d $directory/$name/$filename.sasa.csv
#done
ml schrodinger
#if [ -e "-q" ]; then
for file in $directory/*.pdb
do
    filename=$(basename $file);
    if [ ! -f  $filename.mae ] ; then
        run pdbconvert -ipdb $filename -omae $directory/$name/$filename.mae
        while ! [ -f $directory/$name/$filename.mae ]
        do
            sleep 10
        done
        run utilities/prepwizard $directory/$name/$filename.mae $name/$filename.mae
    fi
done
#Generates a control file for running grid file for Glide
#while ! [ -f $name/$filename.mae ]
#    do
#        sleep 10
#    done

for file in $name/prep*.mae
do
    filename=$(basename "$file");
    gridtable="FORCEFIELD   OPLS_2005\nGRID_CENTER_ASL residue.num 153\nGRIDFILE    $directory/$name/"$filename"grid.zip\nINNERBOX   10, 10, 10\nOUTERBOX   30, 30, 30\nRECEP_FILE   "$file""
    echo -e $gridtable > $directory/$name/$filename"grid.inp"
#       bash run_grid_gen.sh -g $directory/$name
done

#Generate a control file for running Glide sialic acid docks from grid file
for file in $directory/$name/prep*.mae
do
    filename=$(basename "$file");
    #docktable="FORCEFIELD   OPLS_2005\nGRIDFILE   "$filename".mae_grid.zip\nLIGANDFILE   /scicomp/home-pure/qxa4/BCP_V2/ligprep_3-sialyl-out.maegz\nNREPORT   2\nPOSTDOCK_XP_DELE   0.5\nPRECISION   XP\nWRITE_XP_DESC   True\nWRITE_RES_INTERACTION   True\nWRITE_CSV       True"
    docktable="FORCEFIELD   OPLS_2005\nGRIDFILE   $directory/$name/"$filename"grid.zip\nLIGANDFILE   /scicomp/home-pure/qxa4/BCP_V2/ligprep_6-sialyl-out.maegz\nNREPORT   2\nPOSTDOCK_XP_DELE   0.5\nPRECISION   XP\nWRITE_XP_DESC   True\nWRITE_RES_INTERACTION   True\nWRITE_CSV       True"
    echo -e $docktable > $directory/$name/$filename.dock.inp
#    bash run_ligdock.sh -l $directory
done

#Runs the sialic acid docking from the control file
for file in $directory/$name/prep*.pdb.maegrid.inp
do
    filename=$(basename "$file");
    run glide $directory/$name/$filename -OVERWRITE -HOST localhost -TMPLAUNCHDIR -LOCAL
    while ! [ -f $directory/$name/$filename_grid.zip ]
    do
        sleep 20
    done
done
for file in $directory/$name/prep*.pdb.mae.dock.inp
do
    filename=$(basename "$file");
    run glide $directory/$name/$filename -OVERWRITE -HOST localhost -TMPLAUNCHDIR -LOCAL
done
#fi
#Generate a control file for antibody docks
#if ! [ $antibody_dock_off ] ; then
#if [ -e "-q" ]; then
#    for file in $directory/$name/prep*.mae
#    do
#        run -FROM psp piper.py -jobname $filename_dock -receptor prep_7TZ5.pdb -receptor_chain B,C -poses 10 -rotations 70000 -OMPI 1 -JOBID -antibody -ligand $filename -ligand_chain A -use_nonstandard_residue yes -HOST localhost:1 -TMPLAUNCHDIR
#    done
#fi
