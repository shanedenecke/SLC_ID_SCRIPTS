#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 10 22:26:07 2019

## download and clean SLC_ID genomes

@author: shanedenecke
"""

##Import Libraries
import os
import pandas as pd
from sh import gunzip 
import re

## Set WD

os.chdir(os.path.expanduser('~')+'/Documents/SLC_id')
if 'proteomes' not in os.listdir(): 
    os.mkdir('./proteomes')
os.chdir('./proteomes')
## Import Data
spin=pd.read_csv('../general_reference/Pest_species_list.csv')

## Remove all current data from directory
os.system('rm -f *')

## Loop to download and clean data
for index in range(0,spin.shape[0]):
    
    ##slice row and print information
    row=spin.iloc[index,]
    print(row['Genus']+' '+row['Species'])
    print(index)
    
    ## download genome
    os.system('wget '+row.loc['link'])
    
    ## Unzip and rename downloaded genome
    files=os.listdir()
    #raw=[x for x in files if '.gz' in files or len(x)!=18]
    for x in files:
        if '.gz' in x:
            raw=x
        elif len(x)!=18:
            raw=x
        elif '.faa' not in x:
            raw=x
    if '.gz' in raw:
        gunzip(raw)
        raw=re.sub('.gz','',raw)
        os.system('rm *.gz')
    os.rename(src=raw,dst='raw_fasta.faa')
    
    ## If needed modify in place fasta file by appending __ to each geneid
    if row.loc['sed']=='yes':
        os.system("perl -p -i -e 's/(^>.+$)/$1__/g' raw_fasta.faa")
    
    ##Perform cleaning command from table
    bash='cat raw_fasta.faa '+row.loc['Command']+' '+row.loc['Final_file']
    os.system(bash)
    
    ## remove junk fasta file
    os.remove('raw_fasta.faa')
    
        
        
        
