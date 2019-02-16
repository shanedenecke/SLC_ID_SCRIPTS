#
library(data.table)
library(dplyr)
library(tidyr)
library(readr)


## collate all SLC transporter lists from human and drosophila searches 
setwd('/home/shanedenecke/Documents/SLC_id')

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
    slcs=fread(paste0(j,'/SLC_dict.csv'))
    slcs$name=gsub("(SLC_[0-9]+)_[0-9]+$",'\\1',slcs$name)
    slcs$name=gsub("(SLC_X)_[0-9]+$",'\\1',slcs$name)
    l2[[j]]=slcs
  }
  unique_hits=rbindlist(l2) %>% unique()
  #unique_hits$name=paste(unique_hits$name,seq(1:length(unique_hits$name)),sep="_")
  l3=list()
  for(k in unique(unique_hits$name)){
    sub=subset(unique_hits,name==k)
    mem=seq(1,nrow(sub),by=1)
    sub$name=paste(sub$nam,mem,sep='_')
    l3[[k]]=sub
  }
  named.hits=rbindlist(l3)
  l[[i]]=named.hits
}





dir.create('Human_Drosophila_crossref')
setwd('Human_Drosophila_crossref')

for (i in 1:length(l)){ 
  sub=l[[i]]
  n=names(l)[i]
  write.csv(sub,file=paste0(n,"_SLC_dict.csv"),row.names = F)
} 


#l[['Dm']]=fread('/home/shanedenecke/Documents/SLC_id/Dm_Final_Database/SLC_dict.csv')
#l[['Hs']]=fread('/home/shanedenecke/Documents/SLC_id/Human_HMM_SLC/SLC_dict.csv')
   

