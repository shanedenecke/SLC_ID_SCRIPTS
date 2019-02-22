#!/usr/bin/Rscript

## SLC length analysis
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))
shhh(library(utils))

args = commandArgs(trailingOnly=TRUE)



##args[1]='/home/shanedenecke/Documents/SLC_id/iterative_database/iterative_database_BemTab/SLC_source_dict.csv'
##setwd('/home/shanedenecke/Documents/SLC_id/iterative_search/iterative_search_BemTab')
##args

a=readLines('dir.txt')
#a=getwd()
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


if(start==end){
  missing=source[code %in% source$code[!(source$code %in% raw$code)]] ## get values in source but not in final
  full=rbind(raw,missing) %>% unique.data.frame() %>% arrange(name)
}else{
  full=raw
}

colnames(full)=c('code','name')

## get count series for each gene
count=c()
for(i in table(full$name)){
  count=c(count,1:i)
}

full$temp=count
final=unite(full,col='name',name,temp)


## get count series for each gene



cat(format_csv(final))
     