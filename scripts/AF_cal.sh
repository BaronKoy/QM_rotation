#!/usr/bin/env bash

path=/home/baron/Documents/rotation_2/QM_rotation/scripts/outputs/SLiM_AFs/

infile=$path/genstart.txt

AF_file=$path/start_AFs.txt
unfreqsAFs=$path/un_freqs.txt

fin_file=$path/final_file.txt

grep 'm1' $infile > $AF_file
sed -i 's/ /\t/g' $AF_file
awk '{print$9}' $AF_file > $unfreqsAFs
while read -r line
    do
        echo "$(("$line"/100))"
    done < $unfreqsAFs > freqAFs.txt