#!/bin/bash
# @nicolemariepaterson 
usage() { }
while getopts ":d:" arg;
do
        case "${arg}" in
                d)
                        directory=$OPTARG
                        ;;
        esac
done
for file in $directory/prep*pdb
do
    filename=$(basename "$file")
    echo $filename
    #sudo /home/qxa4/schrodinger/utilities/prepwizard $filename prep_$filename
    python H3_rbs_check.py -i $filename.siteA.csv -o $filename.rbsnet_antA.csv
    python H3_rbs_check.py -i $filename.siteB.csv -o $filename.rbsnet_antB.csv
    python H3_rbs_check.py -i $filename.siteC.csv -o $filename.rbsnet_antC.csv
    python H3_rbs_check.py -i $filename.siteD.csv -o $filename.rbsnet_antD.csv
    python H3_rbs_check.py -i $filename.siteE.csv -o $filename.rbsnet_antE.csv
    python rbs_check.py -i $filename.Ca.csv -o $filename.rbsnet_Ca.csv
    python rbs_check.py -i $filename.Cb.csv -o $filename.rbsnet_Cb.csv
    python rbs_check.py -i $filename.KM.csv -o $filename.rbsnet_KM.csv
    python rbs_check.py -i $filename.Sa.csv -o $filename.rbsnet_Sa.csv
    python rbs_check.py -i $filename.Sb.csv -o $filename.rbsnet_Sb.csv
    #python get_contact_flare.py --input $filename.rbs.csv --out $filename.rbs.contact_flare.json
done


