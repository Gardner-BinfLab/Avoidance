#!/usr/bin/env python2.7
'''
Created on 18/04/2013

@author: suu13

Bioinformatics alignment format converter.


Format name     Reads     Writes     Notes
clustal     1.46     1.46     The alignment format of Clustal X and Clustal W.
emboss     1.46     No     The EMBOSS simple/pairs alignment format.
fasta     1.46     1.48     This refers to the input file format introduced for Bill Pearson's FASTA tool, where each record starts with a ">" line. Note that storing more than one alignment in this format is ambiguous. Writing FASTA files with AlignIO failed prior to release 1.48 (Bug 2557).
fasta-m10     1.46     No     This refers to the pairwise alignment output from Bill Pearson's FASTA tools, specifically the machine readable version when the -m 10 command line option is used. The default free format text output from the FASTA tools is not supported.
ig     1.47     No     The refers to the IntelliGenetics file format often used for ordinary un-aligned sequences. The tool MASE also appears to use the same file format for alignments, hence its inclusion in this table. See MASE format.
maf     TBD     TBD    Multiple Alignment Format (MAF) produced by Multiz. Used to store whole-genome alignments, such as the 30-way alignments available from the UCSC genome browser.
nexus     1.46     1.48     Also known as PAUP format. Uses Bio.Nexus internally. Only one alignment per file is supported.
phylip     1.46     1.46     This is a strict interpretation of the interlaced PHYLIP format which truncates names at 10 characters.
phylip-sequential     1.59     1.59     This is a strict interpretation of the sequential PHYLIP format which truncates names at 10 characters.
phylip-relaxed     1.58     1.58     This is a relaxed interpretation of the PHYLIP format which allows long names.
stockholm     1.46     1.46     Also known as PFAM format, this file format supports rich annotation.

'''
from Bio import AlignIO,SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
import argparse



def Sequence_Format_Changer(InputFile,Input_Format,Output_Format):
    Sequences=SeqIO.parse(InputFile,Input_Format)
    
    with open(InputFile+"."+Output_Format,"w") as Output_File:
        Total_Seq=SeqIO.write(Sequences, Output_File, Output_Format)
    print "Converted %i records... DONE..." % Total_Seq
    return




def Alignment_Reverse_Complementer(InputFile,Input_Format,Output_Format):
    with open(InputFile) as Input:
        AlignmentPrev=AlignIO.parse(Input,Input_Format).next()
        New_Alignment_List=[]
        for Sequence in AlignmentPrev:
            Rev_Comp=Sequence.reverse_complement()
            Rev_Comp.id=Sequence.id
            Rev_Comp.name=Sequence.name
            Rev_Comp.description=Sequence.description
            New_Alignment_List.append(Rev_Comp)
        New_Alignment=MultipleSeqAlignment(New_Alignment_List)
        print New_Alignment.format(Output_Format)
    


def Alignment_Format_Changer(InputFile,Input_Format,Output_Format,Size):
    if (Size ==None):
#        try:
        with open(InputFile) as Input:
            AlignmentPrev=AlignIO.parse(Input,Input_Format).next()
            print AlignmentPrev.format(Output_Format)
#        except:
#            print "File read or format error.\n"
#            #Sequence_Format_Changer(InputFile,Input_Format,Output_Format)
#            #return
    else:#ulan bu kisim ne ise yariyor.
        try:
            with open(InputFile) as Input:
                AlignmentPrev=AlignIO.parse(Input,Input_Format).next()
                subalignment=AlignmentPrev[:Size,:] #virgulden onceki alignmenti slice ediyor, virgulden sonrasi sekanslari
                print subalignment.format(Output_Format)
                    
        except:
            print "File read or format error.\n"

def Sequence_Reverse_Complementer(InputFile):

    Sequences=SeqIO.parse(InputFile,"fasta")
    for Sequence in Sequences:
        print ">"+str(Sequence.description)+"_Reverse_Complement\n"+str(Sequence.seq.reverse_complement())


def Alignment_to_Ungapped_FASTA(InputFile,Input_Format):
    with open(InputFile) as Input:
        AlignmentPrev=AlignIO.parse(Input,Input_Format).next()
        for Sequence in AlignmentPrev:
            print ">"+str(Sequence.id)+"\n"+str(Sequence.seq.ungap('-').ungap('.'))

def RNA_to_DNA(InputFile):
    with open(InputFile) as Input:
        AlignmentPrev=AlignIO.parse(Input,"fasta").next()
        for Sequence in AlignmentPrev:
            print ">"+str(Sequence.id)+"\n"+str(Sequence.seq.back_transcribe()).lower()


def Collapse_Alignment(InputFile,Input_Format):
    with open(InputFile) as Input:
        AlignmentPrev=AlignIO.parse(Input,Input_Format).next() #read alignment file
    Collapsed_Alignment=list()
    for i in range(0,AlignmentPrev.get_alignment_length()): #iterate over alignment file
        Collapsed_Region=list()
        for a in AlignmentPrev[:,i:i+1]:
            Collapsed_Region.append(str(a.seq).upper())
        if(len(set(Collapsed_Region)) == 1): #eger uzunlugu 1 ise zaten hepsi ayni direkt al
            Collapsed_Alignment.append(set(Collapsed_Region).pop())
        else:
            Collapsed_Alignment.append([x for x in Collapsed_Region if (x != '-')][0]) #eger uzunlugu 1 den fazla ise

    print ">%s\n%s" % (InputFile,''.join(Collapsed_Alignment))


def main():
    if args.reversecomplement==True:
        Alignment_Reverse_Complementer(args.file,args.inputformat,args.outputformat)
    elif args.inputformat!=None and args.outputformat !=None:
        Alignment_Format_Changer(args.file,args.inputformat,args.outputformat,args.size)
    elif args.ungappedfasta is True:
        Alignment_to_Ungapped_FASTA(args.file,args.inputformat)
    elif args.inputformat is not None and args.file is not None and args.collapse is False:
        Sequence_Reverse_Complementer(args.file)
    elif args.collapse is True:
        Collapse_Alignment(args.file,args.inputformat)
    else:
        RNA_to_DNA(args.file)
        
    




if __name__ == '__main__':
    Argument_Parser=argparse.ArgumentParser(prog="biof_converter.py")
    Argument_Parser.add_argument('-file',type=str,help="Input file of alignment",required=True)
    Argument_Parser.add_argument('-inputformat',type=str,help="Input format")
    Argument_Parser.add_argument('-outputformat',type=str,help="Output format")
    Argument_Parser.add_argument('-size',type=int,help="The number of output sequences")
    Argument_Parser.add_argument('-reversecomplement',action='store_true',help="Take reverse complement of the alignment")
    Argument_Parser.add_argument('-ungappedfasta',action='store_true',help="Remove gaps of the alignment and write as FASTA file")
    Argument_Parser.add_argument('-collapse',action='store_true',help="Collapse Alignment")
    args=Argument_Parser.parse_args()
    main()
    pass
