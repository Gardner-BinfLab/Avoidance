#!/bin/bash

# $1 A FASTA file of a selected mRNA. (This is a protein coding gene)
# $2 core protein coding genes. (This is for calculating codon index)
# $3 core ncRNA genes of a selected species. (This is for calculating avoidance)

if [ $# -ge 3 ] && [ -f $1 ] && [ -f $2 ] && [ -f $3 ]; then



cat $1 | gawk '{if(/>/) print; else  printf $0;}END{printf "\n"}' > $1.tmp;


cat $1.tmp | gawk '{if(/>/) print; else print substr($0,1,21)}' > $1.1.21.fasta #slice out first 21 nucleotides

./RNAup_avoidance_calculator.py -m $1.1.21.fasta -n $3 #calculate avoidance files

filename=$(cat $1.1.21.fasta | gawk 'match($0,/>(.*)/,m){print m[1]}') #detect the filename



native_cai=$(./CAI_calculator.py -mrna $2 -s $1.tmp | cut -f2)
native_avoidance_mfe=$(cat "$filename".result | gawk 'match($0,/.*\s\(([-]*[0-9]+.[0-9]+)\s=\s.*/,m){avoid+=m[1]}END{print avoid}')
native_folding_mfe=$(cat $1.tmp | gawk '!/>/{print substr($0,1,37)}' | RNAfold --noPS | gawk 'match($0,/ \((.*)\)$/,m) {print m[1]}')

cat $1.tmp | gawk -v c="$native_cai" -v a="$native_avoidance_mfe" -v f="$native_folding_mfe" '{if(/>/) print $0" avoidance_mfe: "a" folding_mfe: "f" cai: "c; else print}'


#clean-up
rm $1.tmp
rm $1.1.21.fasta
rm "$filename".result
rm "$filename".seq


fi;
