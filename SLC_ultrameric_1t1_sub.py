#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 18 17:20:39 2019

@author: shanedenecke
"""

import os
import pandas as pd

os.chdir('/home/shanedenecke/Dropbox/wp7_prodrug/common_orthoDB')

pd.set_option('display.max_rows', 10)
pd.set_option('display.max_columns', 500)

#a=pd.read_csv('node_6656_taxid_85310_29058.og.eog',sep='\t',header=None)
a=pd.read_csv('/home/shanedenecke/Dropbox/wp7_prodrug/common_orthoDB/OrthoDB10_clean_arthropod_suset.tsv',sep='\t',header=0)
taxid=pd.read_csv('Taxid_OrthoDB_master_key.tsv',sep='\t',header=None)
taxid.columns=['taxid','abbreviation','Species_name']
#a=a.drop(a.columns[2:9],axis=1)
#a[['taxid','gene']]=a[1].str.split(':',expand=True)
#a.columns=['OG','junk','taxid','gene']

a['taxid']=pd.to_numeric(a['taxid'])


#b=list(set(a[0]))
all_taxids=list(set(a['taxid']))

#df = pd.DataFrame(columns=['OG','species_present','total'])
    
    
def unicount(series):
    return len(list(set(series)))

unicount(['A','A','B'])  


df=a.groupby('OG')['taxid'].agg({'total_genes':'count','unique_taxids':unicount})

d=[]
for i in range(df.shape[0]):
    #if abs(df.iloc[i]['total_genes']-df.iloc[i]['unique_taxids'])>0 & abs(df.iloc[i]['total_genes']-df.iloc[i]['unique_taxids'])<3:
    if abs(df.iloc[i]['total_genes']-df.iloc[i]['unique_taxids'])<7:
   
        if df.iloc[i]['unique_taxids']>160:
            d.append(df.index[i])
len(d)


final_og=d

with open('1_to_1ish_arthropod_orthologues', 'w') as f:
    for item in final_og:
        f.write("%s\n" % item)



final_og = [line.strip() for line in open("1_to_1ish_arthropod_orthologues.txt", 'r')]

subset=a[a.OG.isin(final_og)]

merged=pd.merge(subset,taxid)
unigene_codes=[str(merged['taxid'][x])+'_'+str(merged['gene'][x]) for x in range(0,merged.shape[0])]
final_names=[str(merged['Species_name'][x])+'_'+str(merged['abbreviation'][x])+'_'+str(merged['OG'][x])+'_'+str(merged['gene'][x]) for x in range(0,merged.shape[0])]

final_dict=pd.DataFrame(
        {'code':unigene_codes,
         'name':final_names
         })

final_dict.to_csv(r'./test.csv',index=None)




############## BIN
        
df=pd.DataFrame()
for i in range(len(b)):
    sub=a[a[0]==b[i]]
    sub_taxids=len(list(set(sub['taxid'])))
    dic={'OG':[b[1]], 'species_present':[sub_taxids], 'total':[sub.shape[0]]}
    df=df.append(pd.DataFrame.from_dict(data=dic))
 