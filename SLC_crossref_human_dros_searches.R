#!/usr/bin/Rscript

shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))

## collate all SLC transporter lists from human and drosophila searches 
setwd('/data2/shane/Documents/SLC_id')

sp.list=c()
for(file in list.files('./proteomes')){
  if(grepl(".fa",file)){
    sp.list=c(sp.list,unlist(strsplit(file,'_'))[1])
  }
}
#sp.list=c('Mp','Tc','Ld','Lm','Hh','Bm','Ha','Bg','Ac','Am','Ag')

dirs=c(paste0('./Human_search/',list.files('./Human_search')),paste0('./Drosophila_search/',list.files('./Drosophila_search')))

l=list()
for(i in sp.list){
  par.match=dirs[sapply(i,grepl,dirs)]
  l2=list()
  for(j in par.match){
    
    ## catch any cases where file is missing
    slcs <- try(fread(paste0(j,'/final_output/SLC_final_output.csv')))
    if("try-error" %in% class(slcs)){
      next
    }
    if(dim(slcs)[1]==0){
      next
    }
    
    ## read in files which exist
    #slcs=fread(paste0(j,'/final_output/SLC_final_output.csv'))
    slcs$name=gsub("(SLC_.+_)[0-9]+$",'\\1',slcs$name)
    l2[[j]]=slcs
    
  }
  unique_hits=rbindlist(l2) %>% unique()
  #unique_hits$name=paste(unique_hits$name,seq(1:length(unique_hits$name)),sep="_")
  l3=list()
  for(k in unique(unique_hits$name)){
    sub=subset(unique_hits,name==k)
    mem=seq(1,nrow(sub),by=1)
    sub$name=paste(sub$nam,mem,sep='')
    l3[[k]]=sub
  }
  named.hits=rbindlist(l3)
  l[[i]]=named.hits
}
#

#### filter out duplicates. If there is a duplicated where 1 is unsorted then keep the named one. Otherwise just take the first
rem.dup=function(x){
  final=list()
  for (i in x$code){
    sub=subset(x,code==i)
    uns=table(sub$name %in% 'Unsorted')
    if(dim(sub)[1]>1 & length(uns[uns==T])==1){
      final[[i]]=sub %>% filter(grepl('Unsorted',name))
    }else if(dim(sub)[1]>1){
      final[[i]]=sub[1,]
    }else{
      final[[i]]=sub
    }
  }
  return(rbindlist(final))
}
l=lapply(l,rem.dup)   


l=lapply(l,function(x) x[!duplicated(x$code),])


dir.create('final_SLC_dicts')
setwd('final_SLC_dicts')

for (i in 1:length(l)){ 
  sub=l[[i]]
  n=names(l)[i]
  write.csv(sub,file=paste0(n,"Final_SLC_table.csv"),row.names = F)
} 


#l[['Dm']]=fread('/home/shanedenecke/Documents/SLC_id/Dm_Final_Database/SLC_dict.csv')
#l[['Hs']]=fread('/home/shanedenecke/Documents/SLC_id/Human_HMM_SLC/SLC_dict.csv')
   

