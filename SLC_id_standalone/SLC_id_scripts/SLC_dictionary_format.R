#!/usr/bin/Rscript

## SLC length analysis
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))
shhh(library(utils))

##setwd('/data2/shane/Documents/SLC_id/Dm_Database_Generate/Hs_to_DroMel_Search/')
setwd('./final_output')
raw.table=fread('total_slc_table.csv')
slc.fams=raw.table$name %>% unique()

l=list()
for(i in slc.fams){
  sub=subset(raw.table,name==i)
  for(j in 1:dim(sub)[1]){
    add=as.character(j)
    sub$name[j]=paste(sub$name[j],add,sep='_')
  }
  l[[i]]=sub
}
final=rbindlist(l)

cat(format_csv(final))















