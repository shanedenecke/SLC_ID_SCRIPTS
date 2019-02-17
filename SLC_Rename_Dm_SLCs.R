#!/usr/bin/Rscript
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))


setwd('/home/shanedenecke/Documents/SLC_id')

key=fread('/home/shanedenecke/Documents/omics_data/Drosophila_melanogaster/keys/Dm_master_key_FB_fasta.csv',select=c(1,3,4))
colnames(key)=c('code','symbol','cg')
start.dict=fread('./Dm_Final_Database/SLC_source_dict.csv')


new.dict=merge(start.dict,key,by='code') %>% unique() %>% separate(col='name',into=c('slc','fam','junk')) %>% select(-junk) %>% unite(col='name',slc,fam,cg,symbol,sep="_")
new.dict=new.dict[!duplicated(new.dict$code),]
cat(format_csv(new.dict))
