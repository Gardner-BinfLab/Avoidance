#!/usr/bin/env python2.7


'''
Created on 16/12/2016

@author: suu13
'''


from __future__ import print_function

from Bio.Data.CodonTable import standard_dna_table
from Bio.Seq import Seq
from string import ascii_lowercase
from numpy import random
from itertools import product
import subprocess
from Bio.SeqUtils import CodonUsage
from Bio import SeqIO
import argparse

from hashlib import md5





def codon_weights(amino_acid):


    #initialize for the first use
    if 'calculated_weights' not in vars():
        global calculated_weights
        calculated_weights={}
        with open(args.frequencytable) as freq_table:
            global codon_frequency_table
            codon_frequency_table=freq_table.readlines()

    if amino_acid in calculated_weights.keys():
        total_counts=sum([count_of_codon for codon_type,count_of_codon in calculated_weights[amino_acid]])
        probabilities=[(codon_type,count_of_codon/total_counts) for codon_type,count_of_codon in calculated_weights[amino_acid]]
        return probabilities

    else:
        for line in codon_frequency_table:
            try:
                calculated_weights[line.split()[1]].append((line.split()[0],float(line.split()[2])))
            except:
                try:
                    calculated_weights[line.split()[1]]=[(line.split()[0],float(line.split()[2]))]
                except:
                    pass
        total_counts=sum([count_of_codon for codon_type,count_of_codon in calculated_weights[amino_acid]])
        probabilities=[(codon_type,count_of_codon/total_counts) for codon_type,count_of_codon in calculated_weights[amino_acid]]
        return probabilities



def create_leading_region(mRNA_piece):
    slice_to_codons=map(''.join,zip(*[iter(mRNA_piece)]*3)) #slice into codons eg. ['ATG', 'AAG', 'CAG', 'TCA', 'TCT', 'CGA', 'ATC', 'AAT', 'ATT']
    amino_acids=map(lambda x: standard_dna_table.forward_table[x], slice_to_codons) #find amino acids eg.['M', 'K', 'Q', 'S', 'S', 'R', 'I', 'N', 'I']
    possible_codons=map(lambda x: [codon for codon, aa in standard_dna_table.forward_table.iteritems() if aa==x],amino_acids)


    map(lambda x: print(">%s\n%s" % (md5(''.join(x).lower()).hexdigest(),''.join(x))),list(product(*possible_codons)))
    return


def shuffle_global_with_initial_seq_with_frequencies(original_mRNA,sequence_piece,codon_frequency_table):
    mRNA_sequence=list(str(original_mRNA.seq))

    mRNA_sequence[0:len(sequence_piece)]=list(sequence_piece) #place the initial sequence
    for i in range(len(sequence_piece),len(mRNA_sequence)-3,3): #shuffle the rest except the stop codon
        original_codon=''.join(mRNA_sequence[i:i+3]).upper()
        codon_amino_acid=standard_dna_table.forward_table[original_codon]
        alternative_codons_with_frequencies=codon_weights(codon_amino_acid)
        replacement_codon=random.choice([codon for codon,p in alternative_codons_with_frequencies],p=[p for codon,p in alternative_codons_with_frequencies])
        if replacement_codon != original_codon: #to mark changes
            mRNA_sequence[i:i+3]=list(replacement_codon)
    mRNA_sequence=(''.join(mRNA_sequence))

    if str(Seq(str(original_mRNA.seq)).translate()) == str(Seq(mRNA_sequence).translate()): #just check the final translation
        return mRNA_sequence


def calculate_5prime_folding(mRNA_mutated):

    shell_command="echo %s | RNAfold --noPS | gawk 'match($0,/ \((.*)\)$/,m) {print m[1]}'" % (mRNA_mutated[0:37]) #first 37 nucleotides
    #shell_command="echo %s | RNAfold --noPS" % (mRNA_mutated[0:37]) #first 37 nucleotides
    mfe_folding=subprocess.check_output(shell_command,shell=True).strip()

    return float(mfe_folding)



def main():
    original_mRNA=SeqIO.parse(args.mrna,"fasta").next() #read only the first sequence, the rest is not read.


    if args.frequencytable is not None and args.avoidancetable is not None and args.trainingmrnas is not None:
        training_cai_index=CodonUsage.CodonAdaptationIndex() #init CAI object
        training_cai_index.generate_index(args.trainingmrnas) #read mRNA file and create CAI index

        with open(args.avoidancetable) as avoidance_file, open(args.frequencytable) as frequency_file:
            avoidance_table=avoidance_file.readlines()
            avoidance_table=filter(None, map(lambda x: x.strip(), avoidance_table))

        avoidance_dictionary={}
        for i in avoidance_table:
            avoidance_dictionary[i.split()[0]]=float(i.split()[1])

        for sequence_piece in sorted(avoidance_dictionary.keys()):
            for _ in xrange(0,args.n):
                mRNA_mutated=shuffle_global_with_initial_seq_with_frequencies(original_mRNA,sequence_piece,args.frequencytable)

                assigned_id=md5(mRNA_mutated.lower()).hexdigest()
                print(">%s avoidance_mfe: %f folding_mfe: %f cai: %f\n%s" %(assigned_id,
                                                                      avoidance_dictionary[sequence_piece],
                                                                      calculate_5prime_folding(mRNA_mutated),
                                                                      training_cai_index.cai_for_gene(str(mRNA_mutated).lower()),mRNA_mutated))
        return

    else:
        create_leading_region(str(original_mRNA.seq).upper())
        return



if __name__ == '__main__':
    Argument_Parser=argparse.ArgumentParser(prog="mRNA_sample.py")
    Argument_Parser.add_argument('-mrna',type=str,help="A FASTA file of single mRNA to mutate. "
                                                       "It must be in nucleotide format. If not other arguments supplied, create all possible combination of input sequence."
                                                       "Be careful about the input sequence size, it may take time.",required=True)
    Argument_Parser.add_argument('-trainingmrnas',type=str,help="Reference mRNAs file to calculate CAI, (training dataset).")
    Argument_Parser.add_argument('-n',type=int,default=100,help="Number of sequences to create per leading sequence. default 100")
    Argument_Parser.add_argument('-frequencytable',type=str,help="Frequency table of codons which is used to select more biologically relevant codons.")
    Argument_Parser.add_argument('-avoidancetable',type=str,help="Leading sequences to merge with GFPs and report avoidance. e.g. 21 nucleotides leading")

    args=Argument_Parser.parse_args()
    main()