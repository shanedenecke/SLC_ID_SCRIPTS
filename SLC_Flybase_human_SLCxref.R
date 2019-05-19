#!/usr/bin/Rscript

shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))
shhh(library(seqinr))

key=fread('/data2/shane/Documents/SLC_id/general_reference/keys/Dm_master_key_by_gene.csv') %>% select(Dm_FBgn)
  #mutate(full=paste(Dm_FBgn,name,CG,Dm_FBpp,sep='_')) %>% select(Dm_FBgn,full) %>% data.table() 

dm_unigene=read.fasta('/data2/shane/Documents/SLC_id/general_reference/model_proteomes/DroMel_unigene.faa',
                      set.attributes = F,as.string = T,forceDNAtolower = F) 
dm_unigene2=data.table(full=names(dm_unigene)) %>% mutate(FBgn=full) %>% separate(FBgn,into=c('FBgn','junk'),sep='_') %>% 
  select(-junk) %>% data.table()



fb.slc.gene=fread('/data2/shane/Documents/SLC_id/general_reference/SLC_info/DroMel_SLC_table_flybase.csv') %>% 
  select(FBgn,Family) %>%
  merge(dm_unigene2,by='FBgn') %>% select(-FBgn)
colnames(fb.slc.gene)=c('Family','code')
  
  


hs.slcs=fread('/data2/shane/Documents/SLC_id/Dm_Database_Generate/DroMel_iterative_search/final_output/SLC_final_output.csv')# %>%
  #mutate(FBgn=code) %>% separate(FBgn,into=c('FBgn','junk'),sep='_') %>% select(-junk)




a=merge(hs.slcs,fb.slc.gene,by='code',all=T) %>% data.table() %>% arrange(name) %>% data.table()
for(i in 1:nrow(a)){
  if(is.na(a[i]$name)){
    a[i,'name']=a[i]$Family
  }
}


num_strip=function(x){
  y=unlist(strsplit(x,split='_'))[1:2] %>% paste(collapse = '_') %>% as.character()
  names(y)=NULL
  return(y)
}

a$name=sapply(a$name,num_strip)
a=select(a,-Family)

l=list()
for(i in unique(a$name)){
  sub=a[name==i]
  series=1:dim(sub)[1]
  for(j in series){
    sub[j,'name']=paste(sub[j]$name,j,sep='_')
  }
  l[[i]]=sub
}

final=rbindlist(l)
cat(format_csv(final))











#hs.slcs[!(code %in% flybase.slc$code)]


#fly.new=flybase.slc[!(code %in% hs.slcs$code)]

#final_table=rbind(hs.slcs,fly.new)
#col.names = c('CG','refseq','name','family','Dm_FBgn','FBpp')) %>% select(CG,name,family,Dm_FBgn) %>%
#  merge(key,by='Dm_FBgn')

#colnames(flybase.slc)=c('code','name')




#                  ,select=c('FBgn','Family','SYMBOL','FBpp','ANNOTATION_SYMBOL')) #%>% 
#  unite(col='code',FBgn,SYMBOL,ANNOTATION_SYMBOL,FBpp,sep='_')
#colnames(flybase.slc)=c('Dm_FBgn','name')
