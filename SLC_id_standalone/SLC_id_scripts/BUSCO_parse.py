#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 10:12:36 2020

BUSCO output parser

@author: shanedenecke
"""
### Packages
import os
import pandas as pd
import argparse
import re
import sys


def list_to_df(x):
    k=[]
    v=[]
    namedict={'C':'Completeness','S':'Single_copy','D':'Duplicated','F':'Fragmented','M':'Missing','Species':'Species'}
    for j in x:
        sublist=j.split(':')
        k.append(namedict[sublist[0]])
        v.append([sublist[1]])
    final_data=pd.DataFrame(data=dict(zip(k,v)))
    return(final_data)



################################################## Read in ARGS
CLI=argparse.ArgumentParser()
CLI.add_argument("-file",type=str,help='Summary output file from BUSCO')
CLI.add_argument("-dir",type=str,help='Directory containing BUSCO outputs')
CLI.add_argument("-thresh",type=int,default=0,help='threshold for inclusion. e.g. 75 would remove species with below this value')
args = CLI.parse_args()

#args.dir='./clean_summary/'

if args.dir:
    os.chdir(args.dir)
    blank=pd.DataFrame()
    log=[]
    for i in [x for x in os.listdir() if 'txt' in x]:
        spname=i[0:6]
        try:
            raw=[line.rstrip('\n') for line in open(i)][7]
        except:
            log.append('Cannot read in '+i+' File')
            pass
        values=re.findall(r'[A-Z].[0-9]+\.[0-9]',raw)
        values.append('Species:'+spname)
        sp_dataframe=list_to_df(values)
        blank=blank.append(sp_dataframe,ignore_index=True)

    #blank.columns=['Completeness','Single_copy','Duplicated','Fragmented','Missing','Species']
    blank['Completeness']=pd.to_numeric(blank['Completeness'])
    blank=blank[blank.Completeness>args.thresh]
    blank.to_csv(sys.stdout,index=False,sep='\t')

elif args.file:
    i=args.file
    b=os.path.basename(i)
    blank=pd.DataFrame()
    spname=b[0:6]
    try:
        raw=[line.rstrip('\n') for line in open(i)][7]
    except:
        print('Cannot read in '+i+' File')
        pass 
    values=re.findall(r'[A-Z].[0-9]+\.[0-9]',raw)
    values.append('Species:'+spname)
    sp_dataframe=list_to_df(values)
    blank=blank.append(sp_dataframe,ignore_index=True)
    blank.to_csv(sys.stdout,header=None,index=False,sep='\t')
