#!/usr/bin/Rscript

#rm(list=ls())
## SLC length analysis
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))
shhh(library(stringr))

##already got lengths from expression analysis

args = commandArgs(trailingOnly=TRUE)
H=as.character(args[1])

#setwd('/mnt/disk/shane/Transporter_ID/SLC_id_pipeline/Dm_Database_Generate/Hs_to_DroMel_Search')
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
hs.key=fread(paste0(H,'/SLC_id_reference/HOMSAP_GENE_KEY.csv'))
hs.dict=fread(paste0(H,'/SLC_id_reference/HOMSAP_SLC_DICT.csv')) %>% rename(HUGO_name=code) %>% data.table()
hs.dict$family=gsub('(^SLC_[0-9|A-Z]+).+$','\\1',hs.dict$name)
hs.dict$HUGO_name=gsub(" PE","",hs.dict$HUGO_name)
dm.key=fread(paste0(H,'/SLC_id_reference/DROMEL_GENE_KEY.csv'))
dm.dict=fread(paste0(H,'/SLC_id_reference/DROMEL_SLC_LENGTH_KEY.csv'))
#dm.dict$family=dm.dict$Family


## get table of all SL lentths in humans
#hs.slcs=merge(hs.dict,hs.key,by='HUGO_name') %>% group_by(family) %>% summarise(avg.length=mean(len),maximum=max(len),minimum=min(len)) %>% mutate(sp='hs')
hs.slcs=merge(hs.dict,hs.key,by='HUGO_name') %>% group_by(family) %>% 
  summarise(avg.length=mean(len),maximum=max(len),minimum=min(len),stdev=sd(len)) %>% mutate(sp='hs') %>% data.table()
dm.slcs=merge(dm.dict,dm.key,by='Dm_FBgn') %>% group_by(family) %>% 
  summarise(avg.length=mean(len),maximum=max(len),minimum=min(len),stdev=sd(len)) %>% mutate(sp='hs') %>% data.table()

l2=list()
for(i in hs.slcs$family){
  dm=dm.slcs[family==i]
  hs=hs.slcs[family==i]
  min=min(dm$minimum,hs$minimum)
  max=max(dm$maximum,hs$maximum)
  stdev=max(dm$stdev,hs$stdev)
  l2[[i]]=data.table(family=i,minimum=min,maximum=max,stdev=stdev)
}
l2[['Unsorted']]=data.table(family='SLC_Unsorted',minimum=0,maximum=10000,stdev=100)
com.slcs=rbindlist(l2)

l=list()

## for loop will iterate over potential SLCs and compare them to human SLC lengths
for(i in 1:nrow(gene_lengths)){
  
  ## get info for test gene
  test.len=gene_lengths[i,]$len %>% as.numeric()
  test.name=gene_lengths[i,]$gene %>% as.character()
  test.fam=gene_lengths[i,]$family %>% as.character()
  
  if(test.fam=='SLC_Unsorted' & test.len >100 & test.len<1500){
    l[[i]]=data.table(cbind(gene_lengths[i,c('gene','len')],evaluation='UNK',com.slcs[family==test.fam]))
    next
  }
  
  ## get spread for human SLC family. maximum of either mean+stdev or 1.5*maximum
  hs.stdev=com.slcs[which(com.slcs$family==test.fam),'stdev'] %>% as.numeric()
  hs.max1=com.slcs[which(com.slcs$family==test.fam),'maximum'] * 1.5
  hs.min1=com.slcs[which(com.slcs$family==test.fam),'minimum'] * .667
 
  if(is.na(hs.stdev)){
    hs.max2=com.slcs[which(com.slcs$family==test.fam),'maximum'] * 1.5
    hs.min2=com.slcs[which(com.slcs$family==test.fam),'maximum'] * .667
  }else{
    hs.max2=hs.max1 + hs.stdev
    hs.min2=hs.min1 - hs.stdev
  }
  
  ############ EXCEPTION
  #if(dim(hs.max1)[1]==0 | dim(hs.min1)[1]==0){hs.max1=10000; hs.min1=0}
  #if(dim(hs.max2)[1]==0 | dim(hs.min2)[1]==0){hs.max2=10000; hs.min2=0}
  
  if(test.len>hs.max1 & test.len>hs.max2){
    l[[i]]=data.table(cbind(gene_lengths[i,c('gene','len')],evaluation='LONG',com.slcs[which(com.slcs$family==test.fam),]))
  }else if(test.len<hs.min1 & test.len<hs.min2){
    l[[i]]=data.table(cbind(gene_lengths[i,c('gene','len')],evaluation='SHORT',com.slcs[which(com.slcs$family==test.fam),]))
  }else{
    l[[i]]=data.table(cbind(gene_lengths[i,c('gene','len')],evaluation='GOOD',com.slcs[which(com.slcs$family==test.fam),]))
  }
}


total=rbindlist(l)

write.csv(total,'./Total_length_analysis.csv')
final=total %>% select(gene,family) %>% data.table()
colnames(final)=c('code','name')
cat(format_csv(final))

  
