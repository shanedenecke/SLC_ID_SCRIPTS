library(data.table)
library(dplyr)
library(tidyr)
setwd('/data2/shane/Documents/OrthoDB_Mapping')

a=fread('/data2/shane/Documents/OrthoDB_Mapping/OrthoDB_10_family_counts/OrthoDB_tax_key_3col.tsv')
b=a$V1
c=paste(b,collapse='_|')

writeLines(c,'/data2/shane/Documents/OrthoDB_Mapping/OrthoDB_10_family_counts/Taxid_code_search.txt')

system('
cd /data2/shane/Documents/OrthoDB_Mapping/OrthoDB_10_family_counts/
grep -E (contents of Taxid_code_search.txt) odb10v0_OG2genes.tab > arth_og_reduced.tsv
  ')

d=fread('/data2/shane/Documents/OrthoDB_Mapping/OrthoDB_10_family_counts/arth_og_reduced2.tsv',header=F) %>% 
  separate(V2,into=c('taxid','geneid'),sep = '_')

e=d[taxid %in% a$V1]
e=e[grepl('at6656',V1)] ### extract only arthropoda node
#f1=e[taxid %in% a$V1] 
#test=e[taxid %in% c('6239','132113')]

#f= e %>% group_by(V1,taxid) %>% summarize(count=length(taxid))
f=e[ , .(count = .N), by = .(V1,taxid)]

g=f %>% data.table()

#h=g %>% spread(taxid,count) %>% data.table()
h=dcast(g,V1~taxid)


#h[is.na(h)]=0
k=h
taxcodes=as.character(b)
for(i in 2:length(colnames(k))){
  num=colnames(k)[i]
  if(num %in% taxcodes){
    colnames(k)[i]=a$V2[which(a$V1==num)]
  }
}
k[is.na(k)]=0
m=k[complete.cases(k)] 

fwrite(m,'./OrthoDB_10_family_counts/CAFE_orthoDB_RAW.tsv')

rs=m %>% select(-V1) %>% rowSums()
ra=c()
for(i in 1:nrow(m)){
  nums=as.numeric(m[i,2:(ncol(m)-2)])
  ra=c(ra,max(nums)-min(nums))
}

ra=m %>% select(-V1) %>% rowwise() %>% mutate(range=max(AcrEch:MyzCer)-min(AcrEch:MyzCer)) %>% data.table() %>% .[['range']]
rowm=m %>%  select(-V1) %>% rowMeans() %>% ceiling()
m$total=rs
m$range=ra
m$CaeEle=rowm


o=m[range<5 & range >3] %>% head(300)

colnames(o)[1]='Family ID'
o$Desc='(null)'
p=select(o,`Family ID`,Desc,everything()) %>% select(-total,-range)

fwrite(p,'./OrthoDB_10_family_counts/CAFE_orthoDB_supp.tsv',sep='\t')
p=fread('/data2/shane/Documents/OrthoDB_Mapping/OrthoDB_10_family_counts/CAFE_orthoDB_supp.tsv')
m[total<500] %>% dim()


count.summary=tmhmm.filtered.full %>% group_by(family,abbreviation) %>% summarize(count=length(family)) %>% spread(key=family,value=count) %>% data.table()


f=e[,.(count=length(taxid)),.(V1,taxid)];


e[,list(count=length(taxid)),by=taxid]


aql <- melt(airquality, id.vars = c("Month", "Day")) %>% data.table()
aqw <- dcast(aql, Month + Day ~ variable)
