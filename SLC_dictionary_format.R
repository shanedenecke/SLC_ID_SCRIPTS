#!/usr/bin/Rscript

## SLC length analysis
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))
shhh(library(utils))

args = commandArgs(trailingOnly=TRUE)



##args[1]='/home/shanedenecke/Documents/SLC_id/Dm_Final_Database/SLC_source_dict.csv'
##setwd('/home/shanedenecke/Documents/SLC_id/Drosophila_Database/Hs_to_Dm_Search')
##args

a=readLines('dir.txt')
startdir= a %>% strsplit(.,split='_') %>% unlist()
start=startdir[grepl('[A-Z][a-z]{2}[A-Z][a-z]',startdir)]


if(length(start)==0){
  start='goodbye'
}

writeLines(start,'test.txt')

enddir=args[1] %>% strsplit(.,split=c('_|/')) %>% unlist() 
end=enddir[grepl('[A-Z][a-z]{2}[A-Z][a-z]',enddir)]

if(length(end)==0){
  end='hello'
}

##import data
raw=fread('./final_output/total_slc_table.csv',select=c(1,2),sep=',',fill=T)
source=fread(as.character(args[1]))

## get count series for each gene
count=c()
for(i in table(raw$family)){
  count=c(count,1:i)
}

raw$temp=count
final=unite(raw,col='name',family,temp)
colnames(final)=c('code','name')

if(start==end){
  missing=source[code %in% source$code[!(source$code %in% final$code)]] ## get values in source but not in final
  final=rbind(final,missing) %>% unique.data.frame()
}

cat(format_csv(final))
     