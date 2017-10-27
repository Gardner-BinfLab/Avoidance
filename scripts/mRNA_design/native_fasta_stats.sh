#!/bin/bash

# $1 A FASTA file of a selected mRNA. (This is a protein coding gene)
# $2 core protein coding genes. (This is for calculating codon index)


if [ $# -ge 2 ] && [ -f $1 ] && [ -f $2 ]; then

native_cai=$(./CAI_calculator.py -mrna $2 -s $1 | cut -f2)
native_avoidance_mfe=$(cat $1 | gawk '!/>/{print substr($0,1,21)}' | while read i; do grep -i $i avoidance_table.tsv | cut -f2; done)
native_folding_mfe=$(cat $1 | gawk '!/>/{print substr($0,1,37)}' | RNAfold --noPS | gawk 'match($0,/ \((.*)\)$/,m) {print m[1]}')

cat $1 | gawk -v c="$native_cai" -v a="$native_avoidance_mfe" -v f="$native_folding_mfe" '{if(/>/) print $0" avoidance_mfe: "a" folding_mfe: "f" cai: "c; else print}'

fi;
