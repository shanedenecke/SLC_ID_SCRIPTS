#!/usr/bin/Rscript
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))


setwd('~/Documents/SLC_id')

dm.key=fread('./general_reference/keys/Dm_master_key_FB_fasta.csv',select=c(1,3,4))
colnames(dm.key)=c('code','symbol','cg')
dm.start.dict=fread('./DroMel_Database/SLC_source_dict.csv')
dm.new.dict=merge(dm.start.dict,dm.key,by='code') %>% unique() %>% separate(col='name',into=c('slc','fam','junk')) %>% select(-junk) %>% unite(col='name',slc,fam,cg,symbol,sep="_")
dm.new.dict=dm.new.dict[!duplicated(dm.new.dict$code),]

cat(format_csv(dm.new.dict))



