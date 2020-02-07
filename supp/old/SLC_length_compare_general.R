## SLC length analysis
library(xlsx)
library(dplyr)
library(data.table)
library(tidyr)
library(readr)

## already got lengths from expression analysis

#setwd('/home/shanedenecke/Documents/SLC_id/Human_HMM_SLC/length_analysis_FIND')

gene_lengths=fread('./gene_lengths.txt',colClasses = 'character',col.names = c('gene','len','family'))
gene_lengths$family=gsub('(SLC_[0-9]+|X)_','\\1',gene_lengths$family,perl=T)
hs.key=fread('~/Documents/SLC_id/general_reference/SLC_info/Hs_master_key.csv')
hs.dict=fread('~/Documents/SLC_id/general_reference/SLC_info/HomSap_SLC_dict.csv',col.names = c('sd_name','HUGO_name')) %>% 
  separate(sd_name,into=c('1','2','C'),sep="_") %>% unite("slc_family",'1','2') %>% select(-C) 



hs.dict$HUGO_name=gsub(" PE","",hs.dict$HUGO_name)
hs.slcs=merge(hs.dict,hs.key,by='HUGO_name') %>% group_by(slc_family) %>% summarise(avg.length=mean(len),maximum=max(len),minimum=min(len)) %>% mutate(sp='hs')
hs.slcs=merge(hs.dict,hs.key,by='HUGO_name') %>% group_by(slc_family) %>% 
  summarise(avg.length=mean(len),maximum=max(len),minimum=min(len),stdev=sd(len)) %>% mutate(sp='hs')

#ha.slcs=ha.slcs.raw %>% select(ha_geneid,slc_family,len)
#ha.slcs$slc_family=
#colnames(ha.slcs)[1]='geneid'

#dm.slcs=merge(dm.dict,dm.key,by='Dm_FBgn') 
#colnames(dm.slcs)[1]='geneid'



l=list()
for(i in 1:nrow(gene_lengths)){
  test.len=gene_lengths[i,]$len %>% as.numeric()
  test.name=gene_lengths[i,]$gene %>% as.character()
  test.fam=gene_lengths[i,]$family %>% as.character()
  
  hs.stdev=hs.slcs[which(hs.slcs$slc_family==test.fam),'stdev']
  
  hs.max1=hs.slcs[which(hs.slcs$slc_family==test.fam),'maximum'] * 1.5
  hs.min1=hs.slcs[which(hs.slcs$slc_family==test.fam),'minimum'] * .6667
  
  if(is.na(hs.stdev)){
    hs.max2=hs.slcs[which(hs.slcs$slc_family==test.fam),'maximum'] * 1.5
    hs.min2=hs.slcs[which(hs.slcs$slc_family==test.fam),'maximum'] * .6667
  }else{
    hs.max2=hs.slcs[which(hs.slcs$slc_family==test.fam),'maximum'] + hs.stdev
    hs.min2=hs.slcs[which(hs.slcs$slc_family==test.fam),'maximum'] - hs.stdev
  }
  
  if(test.len>hs.max1 & test.len>hs.max2){
    l[[i]]=data.table(cbind(gene_lengths[i,],evaluation='LONG',hs.slcs[which(hs.slcs$slc_family==test.fam),]))
  }else if(test.len<hs.min1 & test.len<hs.min2){
    l[[i]]=data.table(cbind(gene_lengths[i,],evaluation='SHORT',hs.slcs[which(hs.slcs$slc_family==test.fam),]))
  }else{
    l[[i]]=data.table(cbind(gene_lengths[i,],evaluation='GOOD',hs.slcs[which(hs.slcs$slc_family==test.fam),]))
  }
}

final=rbindlist(l) %>% filter(evaluation=='GOOD') %>% select(gene,family) %>% data.table()
cat(format_csv(final))

  