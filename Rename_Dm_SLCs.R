#!/usr/bin/Rscript
library(tidyr)
library(data.table)
library(readr)
library(dplyr)


setwd('/home/shanedenecke/Documents/SLC_id')

key=fread('/home/shanedenecke/Documents/omics_data/Drosophila_melanogaster/keys/Dm_master_key_FB_fasta.csv',select=c(1,3,4))
colnames(key)=c('code','symbol','cg')
start.dict=fread('./Dm_Final_Database/SLC_dict.csv')


new.dict=merge(start.dict,key,by='code') %>% unique() %>% separate(col='name',into=c('slc','fam','junk')) %>% select(-junk) %>% unite(col='name',slc,fam,cg,symbol,sep="_")

cat(format_csv(new.dict))
