#!/usr/bin/Rscript

## SLC length analysis
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))


##import data
#setwd('/home/shanedenecke/Documents/SLC_id/iterative_search/iterative_search_BomMor')
raw=fread('./final_output/total_slc_table.csv',select=c(1,2),sep=',',fill=T)


##get count series for each gene
count=c()
for(i in table(raw$family)){
  count=c(count,1:i)
}

raw$temp=count
final=unite(raw,col='name',family,temp)
colnames(final)=c('code','name')
final=final[!duplicated(final$code),]


cat(format_csv(final))
