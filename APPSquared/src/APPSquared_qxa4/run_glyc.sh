#!bin/bash
while getopts ":i:r:" i;
do
    case "${i}" in
        i)
                input_dir=$OPTARG
        ;;
    esac
done

for file in $input_dir/*.pdb
do
    python glyc.py $file
done


