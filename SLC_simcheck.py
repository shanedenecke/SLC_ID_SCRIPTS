#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 13:57:09 2019


SLC simcheck.py


1st argument SLC marked Fasta file
2nd argument SLC dict

Objective:
    Check whether each candidate has at least 20% similarity to at least 1 family member

@author: shanedenecke
"""
import os
import pandas as pd
from Bio import SeqIO
import sys 


os.path.dirname(os.path.abspath(__file__))
#os.chdir('/home/shanedenecke/Documents/SLC_id/iterative_search/iterative_search_BemTab')


#sys.argv=['','./simcheck/proteome_SLC_mark_simcheck.fa','./final_output/SLC_final_output.csv','../../proteomes/BemTab_unigene.faa']


##define argument variables
if len(sys.argv)>2:
    mark_fa=list(SeqIO.parse(str(sys.argv[1]), 'fasta'))
    slc_table=pd.read_csv(sys.argv[2],sep=',',header=0)
    ori_fa=list(SeqIO.parse(str(sys.argv[3]), 'fasta'))
 
## run for loop
    
for index in range(0,slc_table.shape[0]):
    row=slc_table.iloc[index,]
    code=row['code']
    family=row['name']
    sequence=[x for x in ori_fa if code in x.description]
    SeqIO.write(sequence,'test.fa','fasta')
    bst_cmd='blastp -query test.fa -db '+str(sys.argv[1])+' -outfmt "6 qseqid sseqid pident evalue qcovs" -max_target_seqs 5 > blast_output.tsv'
    os.system(bst_cmd)
