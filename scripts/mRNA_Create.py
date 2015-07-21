#!/usr/bin/env python2.7


'''
Created on 30/06/2014

@author: suu13
'''

from Bio.SeqUtils import CodonUsage
from Bio.Data import CodonTable
from Bio.Data.CodonTable import standard_dna_table

from hashlib import sha1
from Bio.Seq import Seq
from Bio import SeqIO
import argparse
import time
from random import sample
from random import choice

from numpy import random
from re import sub
import subprocess


def codon_weights(codon_frequency_table,amino_acid):
    if 'calculated_weights' not in vars():
        global calculated_weights
        calculated_weights={}
        
    if amino_acid in calculated_weights.keys():
        total_counts=sum([Count for Codon,Count in calculated_weights[amino_acid]])
        probabilities=[(Codon,Count/total_counts) for Codon,Count in calculated_weights[amino_acid]]
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
        total_counts=sum([Count for Codon,Count in calculated_weights[amino_acid]])
        probabilities=[(Codon,Count/total_counts) for Codon,Count in calculated_weights[amino_acid]]
        return probabilities
                
            
            
        
    

class mRNAmutate:

    def __init__(self,original_mRNA):
        self.original_mRNA_translate=Seq(original_mRNA).translate()
        self.original_mRNA=original_mRNA
    def translation_print(self):
        print self.original_mRNA_translate
    def shuffle_global(self): #shuffle global
        mRNA_seq=list(self.original_mRNA)
        for i in range(0,len(self.original_mRNA)-3,3):
            OriginalCodon=''.join(mRNA_seq[i:i+3]).upper()
            CodonAminoAcid=standard_dna_table.forward_table[OriginalCodon]
            AlternativeCodons=[Codon for Codon,aa in standard_dna_table.forward_table.iteritems() if aa==CodonAminoAcid] #bu ilginc oldu
            ReplacementCodon=sample(AlternativeCodons,1)[0]
            if ReplacementCodon != OriginalCodon and choice([False,True]) is True: #yazi tura atip degistirme karari veriyor
                mRNA_seq[i:i+3]=list(ReplacementCodon)
            else:
                pass

        mRNA_seq=(''.join(mRNA_seq))
        if str(self.original_mRNA_translate) == str(Seq(mRNA_seq).translate()):
            return mRNA_seq

    def shuffle_global_with_inital_seq(self,initial_sequence_piece):
        mRNA_seq=list(self.original_mRNA)
        mRNA_seq[0:len(initial_sequence_piece)]=list(initial_sequence_piece)
        for i in range(len(initial_sequence_piece),len(self.original_mRNA)-3,3):
            OriginalCodon=''.join(mRNA_seq[i:i+3]).upper()
            CodonAminoAcid=standard_dna_table.forward_table[OriginalCodon]
            AlternativeCodons=[Codon for Codon,aa in standard_dna_table.forward_table.iteritems() if aa==CodonAminoAcid] #bu ilginc oldu
            ReplacementCodon=sample(AlternativeCodons,1)[0]
            if ReplacementCodon != OriginalCodon and choice([False,True]) is True: #yazi tura atip degistirme karari veriyor
                mRNA_seq[i:i+3]=list(ReplacementCodon)
            else:
                pass

        mRNA_seq=(''.join(mRNA_seq))
        if str(self.original_mRNA_translate) == str(Seq(mRNA_seq).translate()):
            return mRNA_seq
    def shuffle_global_with_initial_seq_with_frequencies(self,initial_sequence_piece,codon_frequency_table):
        mRNA_seq=list(self.original_mRNA)
        mRNA_seq[0:len(initial_sequence_piece)]=list(initial_sequence_piece) #place initial sequence
        for i in range(len(initial_sequence_piece),len(self.original_mRNA)-3,3): #shuffle rest
            OriginalCodon=''.join(mRNA_seq[i:i+3]).upper()
            CodonAminoAcid=standard_dna_table.forward_table[OriginalCodon]
            Alternative_Codons_with_Frequencies=codon_weights(codon_frequency_table,CodonAminoAcid)
            ReplacementCodon=random.choice([codon for codon,p in Alternative_Codons_with_Frequencies],p=[p for codon,p in Alternative_Codons_with_Frequencies])
            if ReplacementCodon != OriginalCodon: #to mark changes
                mRNA_seq[i:i+3]=list(ReplacementCodon)
        mRNA_seq=(''.join(mRNA_seq))
        if str(self.original_mRNA_translate) == str(Seq(mRNA_seq).translate()):
            return mRNA_seq

    def shuffle_window(self,start,end):
        mRNA_seq=list(self.original_mRNA)
        for i in range(start,end,3):
            OriginalCodon=''.join(mRNA_seq[i:i+3]).upper()
            CodonAminoAcid=standard_dna_table.forward_table[OriginalCodon]
            AlternativeCodons=[Codon for Codon,aa in standard_dna_table.forward_table.iteritems() if aa==CodonAminoAcid] #bu ilginc oldu
            ReplacementCodon=sample(AlternativeCodons,1)[0]
            if ReplacementCodon != OriginalCodon and  choice([False,True]) is True:
                mRNA_seq[i:i+3]=list(ReplacementCodon)
            else:
                pass

        mRNA_seq=(''.join(mRNA_seq))
        if str(self.original_mRNA_translate) == str(Seq(mRNA_seq).translate()):
            return mRNA_seq



def main():
    mRNA=SeqIO.parse(args.mRNAfasta,"fasta").next()
    Mutation_List=[]
    Initial_Sequences=[]
    if args.fastaprev is not None:
        for i in SeqIO.parse(args.fastaprev,"fasta"):
            Mutation_List.append(sha1(str(i.seq).lower()).hexdigest()) #sekanslari lower yap ki emin olalim
    if args.initialsequences is not None:
        with open(args.initialsequences) as initial_seqs_file:
            Initial_Sequences=initial_seqs_file.read().splitlines()

    
    if args.window is None and args.initialsequences is None:
        for _ in range(0,args.n):
            mRNA_Mutated=mRNAmutate(str(mRNA.seq).lower()).shuffle_global()
            if mRNA_Mutated not in Mutation_List: #make sequence and first 50 nucleotides uniq
                GFPid=sha1(mRNA_Mutated.lower()).hexdigest() #yine sekansin kucultulmus halini kullan
                print ">%s\n%s" % (GFPid,mRNA_Mutated)
            else:
                continue
            Mutation_List.append(GFPid)
    
    
    elif args.window is None and args.frequencytable is None and args.initialsequences is not None:
        for Sequence_Piece in Initial_Sequences:
            for _ in range(0,args.n):
                mRNA_Mutated=mRNAmutate(str(mRNA.seq).lower()).shuffle_global_with_inital_seq(Sequence_Piece)
                if mRNA_Mutated not in Mutation_List: #make sequence and first 50 nucleotides uniq
                    GFPid=sha1(mRNA_Mutated.lower()).hexdigest() #yine sekansin kucultulmus halini kullan
                    print ">SeqID:%s-InitSeqID:%s\n%s" % (GFPid,sha1(Sequence_Piece).hexdigest(),mRNA_Mutated)
                else:
                    continue
                Mutation_List.append(GFPid)

    elif args.window is None and args.frequencytable is not None and args.initialsequences is not None:
        with open(args.frequencytable) as freq_table:
            codon_frequency_table=freq_table.readlines()
        for Sequence_Piece in Initial_Sequences:
            for _ in range(0,args.n):
                mRNA_Mutated=mRNAmutate(str(mRNA.seq).lower()).shuffle_global_with_initial_seq_with_frequencies(Sequence_Piece,codon_frequency_table)
                if mRNA_Mutated not in Mutation_List: #make sequence and first 50 nucleotides uniq
                    GFPid=sha1(mRNA_Mutated.lower()).hexdigest() #yine sekansin kucultulmus halini kullan
                    print ">SeqID:%s-InitSeqID:%s\n%s" % (GFPid,sha1(Sequence_Piece).hexdigest(),mRNA_Mutated)
                else:
                    continue
                Mutation_List.append(GFPid)
            
    
    else:
        for _ in range(0,args.n):
            mRNA_Mutated=mRNAmutate(str(mRNA.seq).lower()).shuffle_window(args.window[0],args.window[1])
            if mRNA_Mutated not in Mutation_List: #make sequence and first 50 nucleotides uniq
                GFPid=sha1(mRNA_Mutated.lower()).hexdigest() #yine sekansin kucultulmus halini kullan
                print ">SeqID:%s\n%s" % (GFPid,mRNA_Mutated)
            else:
                continue
            Mutation_List.append(GFPid)



if __name__ == '__main__':
    Argument_Parser=argparse.ArgumentParser(prog="mRNA_Create.py")
    Argument_Parser.add_argument('-mRNAfasta',type=str,help="FASTA file of single mRNA to mutate",required=True)
    Argument_Parser.add_argument('-n',type=int,help="Number of sequences to create",required=True)
    Argument_Parser.add_argument('-window',type=int,nargs=2,help="Window of mutation")
    Argument_Parser.add_argument('-fastaprev',type=str,help="Previously created FASTA file to continue")
    Argument_Parser.add_argument('-initialsequences',type=str,help="Initial sequences to merge with GFPs")
    Argument_Parser.add_argument('-frequencytable',type=str,help="Frequency table of Codons")
    args=Argument_Parser.parse_args()
    main()