How to use;  

./mRNA_design.sh trial.fasta AM946981.2.HMM_Core.fasta AM946981.2-rfam-filtered.fasta 10  


To install required packages:

pip install biopython  #this is biopython, you need this to calculate codon index and read bio-files etc.
pip install joblib  #this is a parallelization library, you need this to run RNAup jobs in parallel.

Optional:
Add 'mRNA_design' directory to your path, in your ~/.bashrc file:

export PATH=$PATH:/path/to/mRNA_design

