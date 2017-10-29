#!/bin/bash

# $1 A FASTA file of a selected mRNA. (This is a protein coding gene)
# $2 core protein coding genes. (This is for calculating codon index)
# $3 core ncRNA genes of a selected species. (This is for calculating avoidance)

# $4 core protein coding genes or protein coding genes of a selected species. (This is for RNA codon distribution. Not necessarily core genes)
# $5 the number of sequences to be created. (If not selected, the default is 100 per the line in avoidance table.)



if [ $# -ge 4 ] && [ -f $1 ] && [ -f $2 ] && [ -f $3 ] && [ -f $4 ]; then

cat $1 | gawk '{if(/>/) print; else print substr($0,1,21)}' > $1.1.21.fasta # slice out the avoidance region

./mRNA_sample.py -mrna $1.1.21.fasta > possible_leading_sequences.fasta #sample from the avoidance region

./CAI_calculator.py -mrna $4 -c > codon_frequency_table.csv #calculate the frequency table of codons



#calculate the avoidance table, this may take time. increase the number of CPUs, the default is 2
mkdir -p ./avoidance_scores
cd ./avoidance_scores

../RNAup_avoidance_calculator.py -mrna ../possible_leading_sequences.fasta -ncrna ../$3 -parallel 2 #calculate the avoidance scores based on ncRNAs


cat ../possible_leading_sequences.fasta | gawk 'match($0,/>(.*)/,m){print m[1]}' | while read i;
do 
avoidance_mfe=$(cat $i.result | gawk 'match($0,/.*\s\(([-]*[0-9]+.[0-9]+)\s=\s.*/,m){avoid+=m[1]}END{print avoid}') #combine the MFEs
sequence=$(grep -A1 $i ../possible_leading_sequences.fasta | grep -v ">") #find the sequence

echo -e $sequence"\t"$avoidance_mfe

done > ../avoidance_table.tsv # create the avoidance table

cd ..

#NOW sample using the available data. This part can be run isolated if all the tables are available.
./mRNA_sample.py -m $1 -f codon_frequency_table.csv -a avoidance_table.tsv -t $2 -n ${5:-100} > designed_mRNAs.fasta


fi;







