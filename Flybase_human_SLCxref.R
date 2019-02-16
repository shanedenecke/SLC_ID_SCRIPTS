#!/usr/bin/Rscript

library(data.table)
library(dplyr)
library(tidyr)
library(readr)

key=fread('~/Documents/SLC_id/general_reference/SLC_info/Dm_master_key_FB_fasta.csv')
flybase.slc=fread('~/Documents/SLC_id/general_reference/SLC_info/DroMel_SLC_table_flybase.csv',select=c('FBgn','Family','ANNOTATION_SYMBOL')) %>% unite(col='name',Family,ANNOTATION_SYMBOL,sep='_')
colnames(flybase.slc)=c('code','name')
hs.slcs=fread('~/Documents/SLC_id/Drosophila_Database/Hs_to_Dm_Search/SLC_dict.csv')

fly.new=flybase.slc[!(code %in% hs.slcs$code)]

final_table=rbind(hs.slcs,fly.new)

cat(format_csv(final_table))
#hs.slcs[!(code %in% flybase.slc$code)]
