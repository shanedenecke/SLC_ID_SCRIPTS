#!/usr/bin/Rscript

shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))
shhh(library(seqinr))

### set working directory to source path
args <- commandArgs(trailingOnly = F) 
#args[1]='/mnt/disk/shane/Transporter_ID/SLC_id_pipeline' 
outdir='/mnt/disk/shane/Transporter_ID/SLC_id_pipeline'

scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
#scriptPath='/mnt/disk/shane/Transporter_ID/SLC_id_pipeline/SLC_id_standalone/SLC_id_scripts'
setwd(scriptPath)
setwd('..')


key=fread('./SLC_id_reference/Dm_master_key_by_gene.csv') %>% select(Dm_FBgn)

dm_unigene=read.fasta('./SLC_id_reference/DroMel_unigene.faa',
                      set.attributes = F,as.string = T,forceDNAtolower = F) 
dm_unigene2=data.table(full=names(dm_unigene)) %>% mutate(FBgn=full) %>% separate(FBgn,into=c('FBgn','junk'),sep='_') %>% 
  select(-junk) %>% data.table()



fb.slc.gene=fread('./SLC_id_reference/DroMel_SLC_table_flybase.csv') %>% 
  select(FBgn,Family) %>%
  merge(dm_unigene2,by='FBgn') %>% select(-FBgn)
colnames(fb.slc.gene)=c('Family','code')
  

hs.slcs=fread(paste0(outdir,'/Hs_to_DroMel_Search/final_output/SLC_final_output.csv'))


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


