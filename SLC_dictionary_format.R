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



















a=readLines('dir.txt')
#a=getwd()
startdir= a %>% strsplit(.,split='_') %>% unlist()
start=startdir[grepl('[A-Z][a-z]{2}[A-Z][a-z]',startdir)]
#

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
     