#!/usr/bin/Rscript

shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))
dir.create('SLC_family_counts')
setwd('/data2/shane/Documents/SLC_id')
co.variables=fread('./general_reference/Co_variables/Olympia_table_august_2019_order_condense.csv',header=T)

file.copy('./DroMel_Database/SLC_source_dict.csv','/data2/shane/Documents/SLC_id/final_SLC_dicts/DroMelFinal_SLC_table.csv')
file.copy('./HomSap_Database/SLC_source_dict.csv','/data2/shane/Documents/SLC_id/final_SLC_dicts/HomSapFinal_SLC_table.csv')


######################## functions
dash.remove=function(x){
  split=unlist(strsplit(x,'_'))
  toge=paste(split[1],split[2],sep='_')
  final=toge
  return(final)
}

var.me=function(x){
  return(var(x)/mean(x))
}
shane.transpose=function(dt,newcol){
  new.rows=colnames(dt)[2:length(colnames(dt))]
  a=transpose(dt)
  colnames(a)=as.character(a[1,])
  a=a[-1,]
  b=mutate(a,newcol=new.rows) %>% select(newcol,everything())
  return(b)
}
########################################3

slc_fams=readLines('/data2/shane/Documents/SLC_id/general_reference/SLC_info/SLC_families.txt')
slc_fams=sapply(slc_fams,dash.remove)
names(slc_fams)=NULL

l=list()
##read dataframes into python list
for (i in list.files('./final_SLC_dicts/')){
  slc_table=fread(paste0('./final_SLC_dicts/',i))
  abbreviation=gsub('Final','',unlist(strsplit(i,'_'))[1])
  slc_table$family=sapply(slc_table$name,dash.remove)
  with.count=slc_table %>% group_by(family) %>% summarise(count=length(family))
  colnames(with.count)[2]=abbreviation
  l[[i]]=with.count
}

count.summary=Reduce(function(x, y) merge(x, y, ,by='family',all=TRUE), l) %>% data.table()
count.summary[is.na(count.summary)]=0
count.summary=shane.transpose(count.summary) 
nams=count.summary$newcol
count.summary$newcol=NULL
totals.matrix=apply(count.summary,2,as.numeric)

varrow=apply(totals.matrix,2,var.me)
var.table=t(data.frame(varrow)) %>% data.frame()
var.table$abbrevaition='blank'
var.table$slc_total='blank'
#var.table=select(var.table,abbreviation,slc_total,everything())

totals=rowSums(totals.matrix)

count.summary=data.table(count.summary)
count.summary$abbreviation=nams
count.summary$slc_total=totals
count.summary=rbind(count.summary,var.table,use.names=F)
 


count.summary=select(count.summary,abbreviation,slc_total,everything())
#rbindlist(count.summary,varrow2,use.names = T,fill=T)
fwrite(count.summary,'./SLC_family_counts/count_summary.csv')


hist(count.summary$slc_total)

## merge to olympia's data
full=merge(count.summary,co.variables,by="abbreviation")

comps=c('Order','Phagy','Vory','Pre-adult_feeding_category','Adult_feeding_category')
l=list()
for(i in slc_fams){
  for(j in comps){
    
    ##Filter for infrequent categories
    good=names(which(table(full[[j]])>8))
    sub=full[which(full[[j]] %in% good)]
    counts=sub[[i]] %>% as.character() %>% as.numeric()
    
    model=lm(formula=counts~sub[[j]])
    coef=summary(model)$coefficients 
    comparison=rownames(coef)
  l[[paste(i,j,sep='_')]]= coef %>% data.table(comp=comparison,family=i,co_variable=j) 
}}

final=rbindlist(l) %>% filter(!grepl("Intercept",comp)) %>% data.table()
colnames(final)=gsub('Pr(>|t|)','p.value',colnames(final),fixed=T)
final$comp=gsub('sub[[j]]','',final$comp,fixed=T)
final$fdr=p.adjust(final$p.value,method='fdr')
final$bonf=p.adjust(final$p.value,method='bonferroni')
final=arrange(final,bonf) %>% filter(bonf<.05) %>% data.table()
fwrite(final,'./SLC_family_counts/correlation_anlaysis.csv')












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