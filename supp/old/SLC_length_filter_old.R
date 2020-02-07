#!/usr/bin/Rscript

#rm(list=ls())
## SLC length analysis
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))

## already got lengths from expression analysis

#setwd('/home/shanedenecke/Documents/SLC_id/Drosophila_search/DROSOPHILA_DanPle')
setwd('./length_analysis')
## import all data calucated in bash SLC_HMM_Search script
gene_lengths=fread('./gene_lengths.txt',colClasses = 'character',col.names = c('gene','len','family'))
gene_lengths$family=gsub('(SLC_.+)_','\\1',gene_lengths$family,perl=T)

rows=dim(gene_lengths)[1]
slc.len=length(gene_lengths$family[grepl('SLC_',gene_lengths$family)])

if(slc.len<rows){
  gene_lengths[rows,]$family=gene_lengths[rows-1,]$family
}



## import keys and dictionaries for human SLC
hs.key=fread('./general_reference/SLC_info/Hs_master_key.csv')
hs.dict=fread('./general_reference/SLC_info/HomSap_SLC_dict_new.csv',col.names = c('sd_name','HUGO_name')) %>% 
  separate(sd_name,into=c('1','2','C'),sep="_") %>% unite("slc_family",'1','2') %>% select(-C) 
hs.dict$HUGO_name=gsub(" PE","",hs.dict$HUGO_name)

## get table of all SL lentths in humans
hs.slcs=merge(hs.dict,hs.key,by='HUGO_name') %>% group_by(slc_family) %>% summarise(avg.length=mean(len),maximum=max(len),minimum=min(len)) %>% mutate(sp='hs')
hs.slcs=merge(hs.dict,hs.key,by='HUGO_name') %>% group_by(slc_family) %>% 
  summarise(avg.length=mean(len),maximum=max(len),minimum=min(len),stdev=sd(len)) %>% mutate(sp='hs') 



l=list()

## for loop will iterate over potential SLCs and compare them to human SLC lengths
for(i in 1:nrow(gene_lengths)){
  
  ## get info for test gene
  test.len=gene_lengths[i,]$len %>% as.numeric()
  test.name=gene_lengths[i,]$gene %>% as.character()
  test.fam=gene_lengths[i,]$family %>% as.character()
  
  if(test.fam=='SLC_Unsorted' & test.len >100 & test.len<1500){
    l[[i]]=data.table(cbind(gene_lengths[i,],evaluation='UNK',hs.slcs[which(hs.slcs$slc_family=="SLC_X"),]))
    next
  }
  
  ## get spread for human SLC family. maximum of either mean+stdev or 1.5*maximum
  hs.stdev=hs.slcs[which(hs.slcs$slc_family==test.fam),'stdev'] %>% as.numeric()
  hs.max1=hs.slcs[which(hs.slcs$slc_family==test.fam),'maximum'] * 2
  hs.min1=hs.slcs[which(hs.slcs$slc_family==test.fam),'minimum'] * .5
 
  if(is.na(hs.stdev)){
    hs.max2=hs.slcs[which(hs.slcs$slc_family==test.fam),'maximum'] * 2
    hs.min2=hs.slcs[which(hs.slcs$slc_family==test.fam),'maximum'] * .5
  }else{
    hs.max2=hs.max1 + hs.stdev
    hs.min2=hs.min1 - hs.stdev
  }
  
  ############ EXCEPTION
  #if(dim(hs.max1)[1]==0 | dim(hs.min1)[1]==0){hs.max1=10000; hs.min1=0}
  #if(dim(hs.max2)[1]==0 | dim(hs.min2)[1]==0){hs.max2=10000; hs.min2=0}
  
  if(test.len>hs.max1 & test.len>hs.max2){
    l[[i]]=data.table(cbind(gene_lengths[i,],evaluation='LONG',hs.slcs[which(hs.slcs$slc_family==test.fam),]))
  }else if(test.len<hs.min1 & test.len<hs.min2){
    l[[i]]=data.table(cbind(gene_lengths[i,],evaluation='SHORT',hs.slcs[which(hs.slcs$slc_family==test.fam),]))
  }else{
    l[[i]]=data.table(cbind(gene_lengths[i,],evaluation='GOOD',hs.slcs[which(hs.slcs$slc_family==test.fam),]))
  }
}

total=rbindlist(l)
write.csv(total,'./Total_length_analysis.csv')
final=total %>% filter(evaluation=='GOOD' | evaluation=='UNK' ) %>% select(gene,family) %>% data.table()
colnames(final)=c('code','name')
cat(format_csv(final))

  