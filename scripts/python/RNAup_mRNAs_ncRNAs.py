#!/usr/bin/env python2.7


'''
Created on 03/07/2014

@author: suu13
'''

from Bio.SeqUtils import CodonUsage
from Bio.Data import CodonTable
from Bio.Data.CodonTable import standard_dna_table

from Bio.Seq import Seq
from Bio import SeqIO
import argparse
import time
from random import sample
from random import choice

from re import sub
import subprocess
from joblib import Parallel,delayed


def RNAup_Executer(Single_RNA,Multiple_RNA,Interaction_Region): #Single RNA vs ncRNAs
    RNA_to_Interact=SeqIO.parse(Single_RNA,"fasta")
    RNAs=list(SeqIO.parse(Multiple_RNA,"fasta"))


    for single_seq in RNA_to_Interact:
        Sample_Name=str(single_seq.description)
        with open(Sample_Name+'.seq',"w") as tmp_file:
            tmp_file.write('>' + sub('\s','$',str(single_seq.description)) +'\n' + str(single_seq.seq)[Interaction_Region[0]:Interaction_Region[1]] +'\n')
            for sequence in RNAs:
                tmp_file.write('>' + sub('\s','$',str(sequence.description)) +'\n' + str(sequence.seq) +'\n')
        Shell_Command="RNAup -b -o --interaction_first < '%s' > '%s' " % (Sample_Name+'.seq',Sample_Name+'.result')
        print str(single_seq.description)
        subprocess.check_call(Shell_Command,shell=True)
    return

def RNAup_Executer_Parallel(Single_RNA,Multiple_RNA,Interaction_Region,n_jobs):
    RNA_to_Interact=SeqIO.parse(Single_RNA,"fasta")
    RNAs=list(SeqIO.parse(Multiple_RNA,"fasta"))

    Parallel(n_jobs=n_jobs,verbose=5)(delayed(Parallel_Function)(single_seq,RNAs,Interaction_Region) for single_seq in RNA_to_Interact)
    return

def Parallel_Function(single_seq,RNAs,Interaction_Region):
    Sample_Name=str(single_seq.description)
    with open(Sample_Name+'.seq',"w") as tmp_file:
        tmp_file.write('>' + sub('\s','$',str(single_seq.description)) +'\n' + str(single_seq.seq)[Interaction_Region[0]:Interaction_Region[1]] +'\n')
        for sequence in RNAs:
            tmp_file.write('>' + sub('\s','$',str(sequence.description)) +'\n' + str(sequence.seq) +'\n')
    Shell_Command="RNAup -b -o --interaction_first < '%s' > '%s' " % (Sample_Name+'.seq',Sample_Name+'.result') #-b yi ekledim
    print str(single_seq.description)
    subprocess.check_call(Shell_Command,shell=True)
    return


def main():
    if args.parallel is None:
        if args.window is None:
            Interaction_Region=(None,None)
            RNAup_Executer(args.singlerna,args.multiplerna,Interaction_Region)
        else:
            Interaction_Region=(args.window[0],args.window[1])
            RNAup_Executer(args.singlerna,args.multiplerna,Interaction_Region)
    else:
        if args.window is None:
            Interaction_Region=(None,None)
            RNAup_Executer_Parallel(args.singlerna,args.multiplerna,Interaction_Region,args.parallel)
        else:
            Interaction_Region=(args.window[0],args.window[1])
            RNAup_Executer_Parallel(args.singlerna,args.multiplerna,Interaction_Region,args.parallel)

if __name__ == '__main__':
    Argument_Parser=argparse.ArgumentParser(prog="Binding_Energy_Minimize_RNAup.py")
    Argument_Parser.add_argument('-singlerna',type=str,help="Single RNA to interact (mRNAs)",required=True)
    Argument_Parser.add_argument('-multiplerna',type=str,help="Other RNAs to interact (ncRNAs)",required=True)
    Argument_Parser.add_argument('-window',type=int,nargs=2,help="mRNAs interaction region")
    Argument_Parser.add_argument('-parallel',type=int,help="Use parallel computing")
    args=Argument_Parser.parse_args()
    main()

