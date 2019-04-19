library(data.table)
library(tidyr)
library(dplyr)

a=fread('/data2/shane/Documents/SLC_id/general_reference/orthoDB_process/reference/Taxid_key.tsv')
b=a %>% separate(V2,into=c('Species','Genus'),sep='_')

firstup <- function(y) {
  substr(y, 1, 1) <- toupper(substr(y, 1, 1))
  y
}

fir3=function(x){
  vec=c()
  for(i in 1:length(x)){
    y=x[i]
    c=substr(y,1,3)
    d=firstup(c)
    vec[i]=d
  }
  return(vec)
}

b$Species=fir3(b$Species)
b$Genus=fir3(b$Genus)

f=unite(b,col=V2,Species,Genus,sep='')


fwrite(f,sep='\t',row.names = F,col.names = F,'/data2/shane/Documents/SLC_id/general_reference/orthoDB_process/reference/Taxid_key2.tsv')

cbind(b$V1,paste(fir3(b$Species),fir3(b$Genus),sep='_'))


setwd('/data2/shane/Documents/SLC_id/proteomes')
for(i in list.files()){
  a=unlist(strsplit(i,split='_'))
  b=sapply(a[1:2],fir3) %>% paste0(collapse = '') %>% paste(a[3],sep = '_')
  file.rename(i,b)
}