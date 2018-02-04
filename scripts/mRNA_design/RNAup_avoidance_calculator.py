#!/usr/bin/env python2.7


'''
Created on 17/12/2016

@author: suu13
'''


from Bio import SeqIO
import argparse
from re import sub
import subprocess
from joblib import Parallel,delayed


def RNAup_executer_parallel(mRNA_file,ncRNA_file,n_jobs):
    mRNA_list=SeqIO.parse(mRNA_file,"fasta")
    ncRNA_list=list(SeqIO.parse(ncRNA_file,"fasta"))

    Parallel(n_jobs=n_jobs,verbose=5)(delayed(compute_function)(a_single_mRNA,ncRNA_list) for a_single_mRNA in mRNA_list)
    return

def compute_function(a_single_mRNA,ncRNA_list):
    sample_name=str(a_single_mRNA.description).strip()
    with open(sample_name+'.seq',"w") as tmp_file:
        tmp_file.write('>' + sub('\s','$',str(a_single_mRNA.description)) +'\n' + str(a_single_mRNA.seq) +'\n')
        for sequence in ncRNA_list:
            tmp_file.write('>' + sub('\s','$',str(sequence.description)) +'\n' + str(sequence.seq) +'\n')
    shell_command="RNAup -b -o --interaction_first < '%s' > '%s' " % (sample_name+'.seq',sample_name+'.result') #-b
    print str(a_single_mRNA.description)
    subprocess.check_call(shell_command,shell=True)
    return


def main():
    RNAup_executer_parallel(args.mrna,args.ncrna,args.parallel)
    return


if __name__ == '__main__':
    Argument_Parser=argparse.ArgumentParser(prog="RNAup_avoidance_calculator.py")
    Argument_Parser.add_argument('-mrna',type=str,help="A list of mRNAs in FASTA format.",required=True)
    Argument_Parser.add_argument('-ncrna',type=str,help="A list of ncRNAs in FASTA format.",required=True)
    Argument_Parser.add_argument('-parallel',type=int,default=1,help="The number of parallel process.")
    args=Argument_Parser.parse_args()
    main()

