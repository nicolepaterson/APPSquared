#!/bin/bash
# @nicolemariepaterson 

while getopts "d:n:" arg; 
do
	case "${arg}" in
		d)
           	directory=$OPTARG
           	;;
       	n)
           	name=$OPTARG
           	;;
	esac
done


out_home=/scicomp/groups/OID/NCIRD/ID/VSDB/GAT
#schrodingerhome=/scicomp/groups/OID/NCIRD/ID/VSDB/GAT/schrodinger

module load rosetta/3.13-cluster


mkdir $out_home/$name

for file in $directory/*.pdb
do
    filename=$(basename "$file")
    score_jd2.default.linuxgccrelease -in:file:s $file -out:file:scorefile $out_home/$name/$filename.scorefile.txt
    per_residue_energies.linuxgccrelease -in:file:s $file -out:file:scorefile $out_home/$name/$filename.res_energy.txt @flag_per_residue
    residue_energy_breakdown.linuxgccrelease -in:file:s $file -out:file:scorefile $out_home/$name/$filename.res_breakdown.txt @flag_residue_energy_breakdown
done

