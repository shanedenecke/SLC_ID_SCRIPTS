#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 14:47:12 2019

@author: shanedenecke
"""

import pandas as pd
import os


a=pd.read_csv('/home/shanedenecke/Documents/SLC_id/general_reference/odb10v0_species.tab',sep='\t')
b=pd.read_csv('/home/shanedenecke/Documents/SLC_id/general_reference/raul_species.csv',sep=',')
d=pd.read_csv('/home/shanedenecke/Documents/SLC_id/general_reference/Pest_species_list.csv',sep=',')

d["Full"]=d["Genus"]+' '+d["Species"]

orsp=a.iloc[:,1].tolist()
rasp=b.loc[:,'Species name'].tolist()
mysp=d.loc[:,'Full'].tolist()

c=[x for x in rasp if x in orsp]
with open('/home/shanedenecke/Documents/SLC_id/general_reference/Possible_species.txt', 'w') as f:
    for item in c:
        f.write("%s\n" % item)
        
subset=b[b['Species name'].isin(c)]
subset.to_csv('/home/shanedenecke/Documents/SLC_id/general_reference/subset_table.csv')

[x for x in mysp if x not in c]
[x for x in mysp if x not in rasp]
