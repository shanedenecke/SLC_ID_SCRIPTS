#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 18 17:20:39 2019

@author: shanedenecke
"""

import os
import pandas as pd
import sys

os.path.dirname(os.path.abspath(__file__))
#os.chdir('/data2/shane/Documents/SLC_id/CAFE/Ultrametric_tree_CAFE')
sys.argv=['','./odb_files/odb10v0_OG2genes.33208.tab','Arthropod']
odb_input=sys.argv[1]
spec=[line.rstrip('\n') for line in open('./species_lists/'+sys.argv[2]+'_species_list.txt')]
#spec=[x+'_0' for x in spec]
## import data
raw_orthodb=pd.read_csv(odb_input,sep='\t',header=0)
raw_orthodb.columns=['og','junk']
raw_orthodb[['taxid','gene']]=raw_orthodb['junk'].str.split('_0:',expand=True)
taxid=pd.read_csv('/data2/shane/Documents/SLC_id/general_reference/non_model_proteomes/keys/Taxid_master_key_full.tsv',sep='\t',header=None)
taxid.columns=['taxid','abbreviation','Species_name']
raw_orthodb['taxid']=pd.to_numeric(raw_orthodb['taxid'])
all_taxids=list(set(raw_orthodb['taxid']))

#spec=[line.rstrip('\n') for line in open('/data2/shane/Documents/SLC_id/general_reference/CAFE/Lepidoptera_species_outgroups.txt')]
#species=[line.rstrip('\n') for line in open('/data2/shane/Documents/SLC_id/general_reference/CAFE/Hemiptera_species_outgroups.txt')]

spec_tax=list(taxid['taxid'][taxid.taxid.isin(spec)])

## aggregate data frame based on how many counts each ortho group has
def unicount(series):
    return len(list(set(series)))

spec_orthodb=raw_orthodb[raw_orthodb.taxid.isin(spec_tax)]

spec_agg=spec_orthodb.groupby('og')['taxid'].agg({'total_genes':'count','unique_taxids':unicount})
spec_agg['ortho_diff']=spec_agg['total_genes']-spec_agg['unique_taxids']
spec_agg=spec_agg.iloc[1:]
spec_final=spec_agg[(spec_agg['ortho_diff']==0) & (spec_agg['unique_taxids']==len(spec)-1)]

final_og=list(spec_final.index.values)[0:199]
final_og[]



######### TESTED UP TO HERE

## filter for ortho groups with close to 1:1 orthology
#final_og=[]
#for i in range(spec_final.shape[0]):
#    #if abs(spec_final.iloc[i]['total_genes']-spec_final.iloc[i]['unique_taxids'])>0 & abs(spec_final.iloc[i]['total_genes']-spec_final.iloc[i]['unique_taxids'])<3:
#    if abs(spec_final.iloc[i]['total_genes']-spec_final.iloc[i]['unique_taxids'])<2:
#   
#        if spec_final.iloc[i]['unique_taxids']<23:
#            final_og.append(spec_final.index[i])


## subset full arhtropod dicitonary with only those genes which have ortho groups with close to 1:1 orthology
subset=spec_orthodb[spec_orthodb.og.isin(final_og)]
merged=pd.merge(subset,taxid)


unigene_codes=[str(merged['taxid'][x])+'_'+str(merged['gene'][x]) for x in range(0,merged.shape[0])]
final_names=[str(merged['Species_name'][x])+'_'+str(merged['abbreviation'][x])+'_'+str(merged['og'][x])+'_'+str(merged['gene'][x]) for x in range(0,merged.shape[0])]


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
