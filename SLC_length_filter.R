#!/usr/bin/Rscript


## SLC length analysis
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))

## already got lengths from expression analysis

#setwd('/home/shanedenecke/Documents/SLC_id/Drosophila_Database/Hs_to_Dm_Search/length_analysis')
setwd('./length_analysis')
## import all data calucated in bash SLC_HMM_Search script
gene_lengths=fread('./gene_lengths.txt',colClasses = 'character',col.names = c('gene','len','family'))
gene_lengths$family=gsub('(SLC_.+)_','\\1',gene_lengths$family,perl=T)

## import keys and dictionaries for human SLC
hs.key=fread('~/Documents/SLC_id/general_reference/SLC_info/Hs_master_key.csv')
hs.dict=fread('~/Documents/SLC_id/general_reference/SLC_info/HomSap_SLC_dict_new.csv',col.names = c('sd_name','HUGO_name')) %>% 
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
  hs.max1=hs.slcs[which(hs.slcs$slc_family==test.fam),'maximum'] * 1.5
  hs.min1=hs.slcs[which(hs.slcs$slc_family==test.fam),'minimum'] * .6667
  
  if(is.na(hs.stdev)){
    hs.max2=hs.slcs[which(hs.slcs$slc_family==test.fam),'maximum'] * 1.5
    hs.min2=hs.slcs[which(hs.slcs$slc_family==test.fam),'maximum'] * .6667
  }else{
    hs.max2=hs.max1 + hs.stdev
    hs.min2=hs.min1 - hs.stdev
  }
  
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
cat(format_csv(final))

  