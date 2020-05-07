shhh <- suppressPackageStartupMessages
shhh(library(dplyr))
shhh(library(data.table))
shhh(library(stringr))
shhh(library(ape))
shhh(library(ggtree))
shhh(library(tidyr))

setwd('./phylogeny/phylip')

for(x in list.files()[grepl('SLC',list.files())]){
  a=fread(x)
  all=a$V1
  uni=unique(all)
  
  if(!is.null(uni)){
    new=c()
    for(i in uni){
      sub=all[all==i]
      if(length(sub)>1){
        for(j in 1:length(sub)){sub[j]=paste0(sub[j],'_',j)}
      }
      new=c(new,sub)
    }
  }
  b=cbind(new,a$V2)
  fwrite(b,x,sep=' ',col.names=F)
}