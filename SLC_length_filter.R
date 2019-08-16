#!/usr/bin/Rscript

#rm(list=ls())
## SLC length analysis
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))
shhh(library(stringr))

## already got lengths from expression analysis

#setwd('/data2/shane/Documents/SLC_id/Dm_Database_Generate/Hs_to_DroMel_Search/')
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
hs.key=fread('/data2/shane/Documents/SLC_id/general_reference/keys/Hs_master_key.csv')
#hs.dict=fread('./general_reference/SLC_info/HomSap_SLC_dict_new.csv',col.names = c('sd_name','HUGO_name')) %>% 
  #separate(sd_name,into=c('1','2','C'),sep="_") %>% unite("slc_family",'1','2') %>% select(-C) 
hs.dict=fread('/data2/shane/Documents/SLC_id/general_reference/SLC_info/HomSap_SLC_dict.csv')
colnames(hs.dict)[2]='HUGO_name'
hs.dict$family=gsub('(SLC_.+)_.+$','\\1',hs.dict$name)
hs.dict$HUGO_name=gsub(" PE","",hs.dict$HUGO_name)
dm.key=fread('/data2/shane/Documents/SLC_id/general_reference/keys/Dm_master_key_by_gene.csv')
dm.dict=fread('/data2/shane/Documents/SLC_id/general_reference/SLC_info/Dm_SLC_length.csv')




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


comp.score=nrow(total[evaluation=='SHORT'])/nrow(total)
a=str_split(getwd(),c('_','/')) %>% unlist() 
score.name=a[length(a)-1]
fwrite(data.table(comp.score,score.name),'/data2/shane/Documents/SLC_id/genome_score/comp_score.txt',append=T)



write.csv(total,'./Total_length_analysis.csv')
final=total %>% filter(evaluation=='GOOD' | evaluation=='UNK' ) %>% select(gene,family) %>% data.table()
colnames(final)=c('code','name')
cat(format_csv(final))

  