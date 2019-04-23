#!/usr/bin/Rscript

shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))

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
h=dcast(melt(g, id.vars = "family"), variable ~ family)
fwrite(h,'./SLC_family_counts/TOTAL_FAMILY_COUNTS.csv')