#!/usr/bin/Rscript

shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))

setwd('/data2/shane/Documents/SLC_id')
#setwd('/home/shanedenecke/Dropbox/wp7_prodrug/SLC_id')
#co.variables=fread('./general_reference/Co_variables/co_variables.csv') 
co.variables=fread('./family_size_variation/Olympia_table_august_2019.csv',header=T) %>%
  select(Species_name,abbreviation,Order,Phagy,Vory,Host_use_category)
#### read all fastas into list

family.size=fread('/home/shanedenecke/Dropbox/wp7_prodrug/SLC_id/family_size_variation/TOTAL_FAMILY_COUNTS.csv')

shane.transpose=function(dt,newcol){
  new.rows=colnames(dt)[2:length(colnames(dt))]
  a=transpose(dt)
  colnames(a)=as.character(a[1,])
  a=a[-1,]
  b=mutate(a,newcol=new.rows) %>% select(newcol,everything())
  return(b)
}

family.size.new=shane.transpose(family.size) %>% data.table()
dtnew <- family.size.new[, lapply(.SD, as.numeric), by=newcol]
colnames(dtnew)[1]='Family'


comps=c('Order','Phagy','Vory')

for(i in 1:length(dtnew$Family)){
  for(j in comps){
    slc=dtnew$Family[i]
    sub=dtnew[Family==slc]
    subcomp=select(co.variables,abbreviation,j) %>% split(subcomp$Order)
    fil.vec=lapply(subcomp,function(x) dim(x)[1])>8
    subcomp[fil.vec]
    
l=list()
for(i in slcs){
  s=annot[[i]] %>% as.character() %>% as.numeric()
  reg=lm(formula=s~annot$Host_use_category)
  coef=summary(reg)$coefficients 
  comparison=rownames(coef)
  
  l[[i]]= coef %>% data.table(comp=comparison,family=i) 
}




library(dplyr)
PATH <- "https://raw.githubusercontent.com/guru99-edu/R-Programming/master/poisons.csv"
df <- read.csv(PATH) %>%
  select(-X) %>% 
  mutate(poison = factor(poison, ordered = TRUE))
glimpse(df)

lm(formula=time~poison,data=df) %>% summary()
anova_one_way <- aov(time~poison, data = df)
TukeyHSD(anova_one_way)






apply(select(family.size.new,-newcol),1,as.numeric)





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
k=fread('/data2/shane/Documents/SLC_id/general_reference/non_model_proteomes/keys/Taxid_master_key_full.tsv',
        col.names=c('taxid','short','long_name'))
for(i in 2:length(colnames(g))){
  s.name=colnames(g)[i]
  colnames(g)[i]=k[which(k$short==s.name)]$long_name
}



shane.transpose=function(dt,newcol){
  new.rows=colnames(dt)[2:length(colnames(dt))]
  a=transpose(dt)
  colnames(a)=as.character(a[1,])
  a=a[-1,]
  b=mutate(a,newcol=new.rows) %>% select(newcol,everything())
  return(b)
}

trans=shane.transpose(g,'SLC_family') %>% data.table()
colnames(trans)[1]='Species_name'

annot=merge(trans,co.variables,by='Species_name')

cols=colnames(annot)
expl=cols[68:73]
slcs=cols[2:66]

for(i in expl){annot[[i]]=as.factor(annot[[i]])}
annot=annot[!is.na(Host_use_category)]
l=list()
for(i in slcs){
  s=annot[[i]] %>% as.character() %>% as.numeric()
  reg=lm(formula=s~annot$Host_use_category)
  coef=summary(reg)$coefficients 
  comparison=rownames(coef)
  
  l[[i]]= coef %>% data.table(comp=comparison,family=i) 
}
rbindlist(l) %>% filter(!grepl("Intercept",comp)) %>% filter(`Pr(>|t|)`<.05)%>% arrange(`Pr(>|t|)`) %>% data.table()




for(i in expl){
  for(j in slcs){
    e=annot[[i]]
    s=annot[[j]] %>% as.numeric()
    
    cor(e,s)
  }
}












fwrite(g,'./SLC_family_counts/TOTAL_FAMILY_COUNTS.csv',row.names = F)



#h=dcast(melt(g, id.vars = "family"), variable ~ family) 
#colnames(h)[1]='SLC_family'