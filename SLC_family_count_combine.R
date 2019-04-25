#!/usr/bin/Rscript

shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))

setwd('/data2/shane/Documents/SLC_id')

#### read all fastas into list
l=list()
for (i in list.files('./SLC_family_counts/',full.names = T)){
  a=fread(i)
  b=unlist(strsplit(i,split='/',fixed=T))
  c=b[length(b)]  
  d=unlist(strsplit(c,'_'))[1]
  e=gsub("Final","",d)
  colnames(a)[2]=e
  l[[i]]=a
}



g=Reduce(function(x, y) merge(x, y, ,by='family',all=TRUE), l)
k=fread('/data2/shane/Documents/SLC_id/general_reference/orthoDB_process/reference/Taxid_OrthoDB_master_key.tsv',
        col.names=c('taxid','short','long_name'))
for(i in 2:length(colnames(g))){
  s.name=colnames(g)[i]
  colnames(g)[i]=k[which(k$short==s.name)]$long_name
}
fwrite(g,'./SLC_family_counts/TOTAL_FAMILY_COUNTS.csv',row.names = F)



#h=dcast(melt(g, id.vars = "family"), variable ~ family) 
#colnames(h)[1]='SLC_family'