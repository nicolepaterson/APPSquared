#!/bin/bash
# @nicolemariepaterson 

while getopts ":d:r:" arg;
do
    case "${arg}" in
        d)
                directory=$OPTARG
                ;;
	r)
		reference=$OPTARG
		;;
    esac
done
for file in $directory/*.pdb
do
    filename=$(basename "$file")
    echo $filename
    #sudo /home/qxa4/schrodinger/utilities/prepwizard $filename prep_$filename
    python get_static_contacts.py --structure $filename --output $filename.static_contacts.csv --itypes all
    python get_contact_frequencies.py --input_files $filename.static_contacts.csv --output_file $filename.freq.tsv --itypes all
    python contact_diffs.py -i $filename.static_contacts.csv -r $reference.static_contacts.csv -o $filename.contactdiffs.static.csv
    python get_contact_flare.py --input $filename.static_contacts.csv --out $filename.contact_flare.json
    python rbs_check.py -i $filename.contactdiffs.static.csv -o $filename.H1rbs.csv
    python H3_rbs_check.py -i $filename.contactdiffs.static.csv -o $filename.H3rbs.csv
    python siteA_check.py -i $filename.contactdiffs.static.csv -o $filename.siteA.csv
    python siteB_check.py -i $filename.contactdiffs.static.csv -o $filename.siteB.csv
    python siteC_check.py -i $filename.contactdiffs.static.csv -o $filename.siteC.csv
    python siteD_check.py -i $filename.contactdiffs.static.csv -o $filename.siteD.csv
    python siteE_check.py -i $filename.contactdiffs.static.csv -o $filename.siteE.csv
    python Ca_check.py -i $filename.contactdiffs.static.csv -o $filename.Ca.csv
    python Cb_check.py -i $filename.contactdiffs.static.csv -o $filename.Cb.csv
    python KM_check.py -i $filename.contactdiffs.static.csv -o $filename.KM.csv
    python Sa_check.py -i $filename.contactdiffs.static.csv -o $filename.Sa.csv
    python Sb_check.py -i $filename.contactdiffs.static.csv -o $filename.Sb.csv
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

