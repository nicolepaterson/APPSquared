#!/bin/bash
# @nicolemariepaterson 

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
    python get_static_contacts.py --structure $filename --output $filename.static_contacts.csv --itypes all
    python get_contact_frequencies.py --input_files $filename.static_contacts.csv --output_file $filename.freq.tsv --itypes all
    python contact_diffs.py -i $filename.static_contacts.csv -r prep_Darwin_6_2021_trimer.pdb.contactdiffs.static.csv -o $filename.contactdiffs.static.csv
    python get_contact_flare.py --input $filename.static_contacts.csv --out $filename.contact_flare.json
    python rbs_check.py -i $filename.contactdiffs.static.csv -o $filename.rbs.csv
    #python get_contact_flare.py --input $filename.rbs.csv --out $filename.rbs.contact_flare.json
done
