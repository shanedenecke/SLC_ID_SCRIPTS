#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 20:45:09 2019

Summarize counts of fasta tables

@author: shanedenecke
"""


## import modules
import os
import pandas as pd


##working directory
os.chdir('/home/shanedenecke/Documents/SLC_id')
os.mkdir('SLC_family_counts')


##import refrence data
dros_slc=pd.read_csv('/home/shanedenecke/Documents/SLC_id/Dm_Final_Database/SLC_dict.csv')
human_slc=pd.read_csv('/home/shanedenecke/Documents/SLC_id/Human_HMM_SLC/SLC_dict.csv')
slc_fams=human_slc['name'].str.split('_',expand=True)[[0,1]].apply(lambda x: '_'.join(x), axis=1).unique().tolist()

##read dataframes into python list
for i in os.listdir('./iterative_search/'):
    slc_table=pd.read_csv('./iterative_search/'+i+'/SLC_dict.csv')
    abbreviation= i.split('_')[-1]

    
    ## clean up data frame
    family=slc_table['name'].str.split('_',expand=True)[[0,1]].apply(lambda x: '_'.join(x), axis=1)
    ids=slc_table.iloc[:,0] ## pd.Series.to_frame(
    clean=pd.DataFrame(dict(ids=ids,family=family))
        
    ##summarize data frame
    summary_count=clean.groupby(['family'])[['ids']].count()
    for j in slc_fams:
        if j not in summary_count.index.tolist():
            add=pd.DataFrame.from_dict({j:0},orient='index',columns=['ids'])
            summary_count=summary_count.append(pd.DataFrame(data=add))
    summary_count=summary_count.sort_index()
    summary_count.index.name='family'   #count_tables.append(summary_count)
    
    ##Write to csv
    summary_count.to_csv('./SLC_family_counts/'+i+'_SLC_counts.csv')

