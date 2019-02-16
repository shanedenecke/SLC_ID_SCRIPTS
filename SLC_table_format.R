#!/usr/bin/Rscript

## SLC length analysis
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))


##import data

raw=fread('./final_output/total_slc_table.csv',select=c(1,2),sep=',',fill=T)


## get count series for each gene
count=c()
for(i in table(raw$family)){
  count=c(count,1:i)
}

raw$temp=count
final=unite(raw,col='name',family,temp)
colnames(final)=c('code','name')


cat(format_csv(final))
