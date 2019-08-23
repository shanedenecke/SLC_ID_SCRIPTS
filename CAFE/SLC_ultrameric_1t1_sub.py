#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 18 17:20:39 2019

@author: shanedenecke
"""

import os
import pandas as pd

os.chdir('/data2/shane/Documents/SLC_id/ultrametric_tree')

## import data
raw_orthodb=pd.read_csv('./OrthoDB10_clean_arthropod_suset.tsv',sep='\t',header=0)
taxid=pd.read_csv('/data2/shane/Documents/SLC_id/general_reference/non_model_proteomes/keys/Taxid_master_key_full.tsv',sep='\t',header=None)
taxid.columns=['taxid','abbreviation','Species_name']
raw_orthodb['taxid']=pd.to_numeric(raw_orthodb['taxid'])
all_taxids=list(set(raw_orthodb['taxid']))

    
## aggregate data frame based on how many counts each ortho group has
def unicount(series):
    return len(list(set(series)))
agg=raw_orthodb.groupby('OG')['taxid'].agg({'total_genes':'count','unique_taxids':unicount})


## filter for ortho groups with close to 1:1 orthology
final_og=[]
for i in range(agg.shape[0]):
    #if abs(agg.iloc[i]['total_genes']-agg.iloc[i]['unique_taxids'])>0 & abs(agg.iloc[i]['total_genes']-agg.iloc[i]['unique_taxids'])<3:
    if abs(agg.iloc[i]['total_genes']-agg.iloc[i]['unique_taxids'])<7:
   
        if agg.iloc[i]['unique_taxids']>160:
            final_og.append(agg.index[i])


## subset full arhtropod dicitonary with only those genes which have ortho groups with close to 1:1 orthology
subset=raw_orthodb[raw_orthodb.OG.isin(final_og)]
merged=pd.merge(subset,taxid)


unigene_codes=[str(merged['taxid'][x])+'_'+str(merged['gene'][x]) for x in range(0,merged.shape[0])]
final_names=[str(merged['Species_name'][x])+'_'+str(merged['abbreviation'][x])+'_'+str(merged['OG'][x])+'_'+str(merged['gene'][x]) for x in range(0,merged.shape[0])]


## write renaming dictionary
final_dict=pd.DataFrame(
        {'code':unigene_codes,
         'name':final_names
         })

final_dict.to_csv(r'./renaming_dictionary.csv',index=None)



## write test file with unicodes
with open('unicodes_fa_subset.txt', 'w') as f:
    for item in unigene_codes:
        f.write("%s\n" % item)



##write OGS
with open('oto_ish_arthropod_orthologues.txt', 'w') as f:
    for item in final_og:
        f.write("%s\n" % item)



#final_og = [line.strip() for line in open("1_to_1ish_arthropod_orthologues.txt", 'r')]