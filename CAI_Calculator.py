#!/usr/bin/env python2.7

'''
Created on 18/11/2013

@author: suu13
'''
from collections import Counter
from Bio.SeqUtils import CodonUsage
from Bio import SeqIO
import argparse
import re
from random import sample

from Bio.Data.CodonTable import standard_dna_table


from os import remove


def codon_count(sequence):
    codon_list=[sequence[i:i+3].upper() for i in range(0,len(sequence),3)]
    return Counter(codon_list)
        
    




def Nucleic_Acid_Detect(strg, search=re.compile(r'[^a|t|g|c|A|T|G|C]').search):
    return not bool(search(strg))


def Exception_Fixer(file_argument):
    mRNA_dict={}
    for mRNA in SeqIO.parse(file_argument,"fasta"):
        if(len(str(mRNA.seq))%3 == 0) and Nucleic_Acid_Detect(str(mRNA.seq))==True:
            #print len(str(mRNA.seq))
            mRNA_dict[str(mRNA.description)]=str(mRNA.seq)
    output_tmp_file_name=''.join((sample('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ',8)+list(".txt")))
    with open(output_tmp_file_name,"w") as tmp:
        for seq_record in mRNA_dict:
            tmp.write('>' + seq_record +'\n' + mRNA_dict[seq_record] +'\n')
    return output_tmp_file_name

def CAI_print(fasta_file,CDS_CAI_Index):
    FASTA_iterator_obj=SeqIO.parse(fasta_file,"fasta")
    for seq_record in FASTA_iterator_obj:
        print "%s\t%f" % (str(seq_record.description),CDS_CAI_Index.cai_for_gene(str(seq_record.seq).lower()))
    


def main():
    
    
    if args.codoncount is True:
        Codons=Counter()
        for sequence in SeqIO.parse(args.mRNAfasta,"fasta"):
            Codons=Codons + codon_count(str(sequence.seq))
        for key,value in Codons.iteritems():
            try:
                print "%s\t%s\t%d" % (key,standard_dna_table.forward_table[key.upper()],value)
            except:
                print "%s\t%s\t%d" % (key,'Stop',value)
        return

            

    try:
        CDS_CAI_Index=CodonUsage.CodonAdaptationIndex() #init CAI object
        CDS_CAI_Index.generate_index(args.mRNAfasta) #read mRNA file and create CAI index
        #CDS_CAI_Index.print_index()
        if args.othersfasta is not None:
            try:
                CAI_print(args.othersfasta,CDS_CAI_Index)
            except TypeError:
                output_tmp_file_name=Exception_Fixer(args.othersfasta)
                CAI_print(output_tmp_file_name,CDS_CAI_Index)
                #print "Exception in othersfasta file which probably have wrong codons or codon numbers..."
                #raise
            except:
                print "Unexpected error"
        else:
            #CDS_CAI_Index.print_index()
            print "\nCodon\tAA\tFrequency"
            for key,value in CDS_CAI_Index.codon_count.iteritems():
                try:
                    print "%s\t%s\t%d" % (key,standard_dna_table.forward_table[key.upper()],value)
                except:
                    print "%s\t%s\t%d" % (key,'Stop',value)
            CAI_print(args.mRNAfasta,CDS_CAI_Index)
    
    except:
        output_tmp_file_name=Exception_Fixer(args.mRNAfasta)
        CDS_CAI_Index=CodonUsage.CodonAdaptationIndex() #init CAI object
        CDS_CAI_Index.generate_index(output_tmp_file_name) #read mRNA file and create CAI index
        if args.othersfasta is not None:
            try:
                CAI_print(args.othersfasta,CDS_CAI_Index)
            except TypeError:
                output_tmp_file_name=Exception_Fixer(args.othersfasta)
                CAI_print(output_tmp_file_name,CDS_CAI_Index)
                #print "Exception in othersfasta file which probably have wrong codons or codon numbers..."
                #raise
            except:
                print "Unexpected error"
        else:
            CDS_CAI_Index.print_index()
            print "\nCodon\tAA\tFrequency"
            for key,value in CDS_CAI_Index.codon_count.iteritems():
                try:
                    print "%s\t%s\t%d" % (key,standard_dna_table.forward_table[key],value)
                except:
                    print "%s\t%s\t%d" % (key,'Stop',value)
            CAI_print(output_tmp_file_name,CDS_CAI_Index)
        remove(output_tmp_file_name)

    

            

if __name__ == '__main__':
    Argument_Parser=argparse.ArgumentParser(prog="CAI_Calculator.py")
    Argument_Parser.add_argument('-mRNAfasta',type=str,help="Reference mRNAs file to calculate CAI, (training dataset)",required=True)
    Argument_Parser.add_argument('-othersfasta',type=str,help="FASTA of other RNAs")
    Argument_Parser.add_argument('-codoncount',action='store_true',help="Only codon count")
    args=Argument_Parser.parse_args()
    main()
    pass