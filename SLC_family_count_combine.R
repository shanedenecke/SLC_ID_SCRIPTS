#!/usr/bin/Rscript

shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))
#setwd('/home/shanedenecke/Dropbox/wp7_prodrug/SLC_id') #local
setwd('/data2/shane/Documents/SLC_id')
dir.create('SLC_family_counts')
co.variables=fread('./general_reference/Co_variables/Olympia_table_august_2019_order_condense.csv',header=T)
#co.variables=fread('/home/shanedenecke/Dropbox/wp7_prodrug/SLC_id/family_size_variation/Olympia_table_august_2019.csv',header=T)



#file.copy('./DroMel_Database/SLC_source_dict.csv','/data2/shane/Documents/SLC_id/final_SLC_dicts/DroMelFinal_SLC_table.csv')
#file.copy('./HomSap_Database/SLC_source_dict.csv','/data2/shane/Documents/SLC_id/final_SLC_dicts/HomSapFinal_SLC_table.csv')


######################## functions
dash.remove=function(x){
  split=unlist(strsplit(x,'_'))
  toge=paste(split[1],split[2],sep='_')
  final=toge
  return(final)
}

var.me=function(x){
  range(x)[2]-range(x)[1]
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
#slc_fams=readLines('/home/shanedenecke/Dropbox/wp7_prodrug/SLC_id/old/general_reference/SLC_info/SLC_families.txt')
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
#good.slcs=names(which(colSums(totals.matrix)>(3*140)))

count.summary=data.table(count.summary)
count.summary$abbreviation=nams
count.summary$slc_total=totals
draft.sum=rbind(count.summary,var.table,use.names=F)


############### SLC53 listed twice
for(i in slc_fams){
  sub=as.numeric(draft.sum[[i]])
  
  #if(max(sub)<10){# & sub[173]<.5){
  if(sub[173]<10){
    print(i)
    draft.sum[[i]]=NULL
  }
}



draft.sum=select(draft.sum,abbreviation,slc_total,everything())
#rbindlist(draft.sum,varrow2,use.names = T,fill=T)
fwrite(draft.sum,'./SLC_family_counts/count_summary.csv')

draft.sum=draft.sum[1:172]
hist(as.numeric(draft.sum$slc_total))

test=count.summary %>%
  mutate_at(vars(matches("SLC")), funs(as.numeric(as.character(.)))) %>%
  data.table()
apply(select(test,matches('SLC')),2,mean)

mns=apply(select(test,matches('SLC')),2,mean)
size.adjust=data.table(fam=names(mns),avg=mns)


## merge to olympia's data
full=merge(draft.sum,co.variables,by="abbreviation")

comps=c('overall_phylo','Diet_category')
l=list()
for(i in slc_fams[slc_fams %in% colnames(full)]){
  for(j in comps){
    
    ##Filter for infrequent categories
    good=names(which(table(full[[j]])>8))
    sub=full[which(full[[j]] %in% good)]
    counts=sub[[i]] %>% as.character() %>% as.numeric()
    #model=aov(formula=counts~sub[[j]])
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

for(i in 1:nrow(final)){final$adjusted_Estimate[i]=final[i]$Estimate/(size.adjust[fam==final[i]$family])$avg}



final=arrange(final,bonf) %>% filter(bonf<1e-10)  %>% filter(abs(adjusted_Estimate)>1) %>% data.table()
fwrite(final,'./SLC_family_counts/correlation_anlaysis.csv')


dir.create('comparison_plots')

for(i in 1:nrow(final)){
  row=final[i]
  co.var=row$co_variable
  fam=row$family
  plot=select(full,fam,co.var)
  plot[[fam]]=as.numeric(plot[[fam]])
  
  plot=plot[plot[[co.var]] %in% names(which(table(plot[[co.var]])>5))]
  
  plot[duplicated(plot)]
  
  pdf(paste0('./comparison_plots/',co.var,'_',fam,'.pdf'))
  gp=ggplot(plot,aes_string(x=co.var,y=fam,fill=co.var))
  gp=gp+geom_boxplot()
  gp=gp+labs(x=paste0('\n',co.var),y=paste0(fam,'\n'))
  gp=gp+theme_bw()
  gp=gp+theme(text=element_text(face="bold",family="serif"),panel.grid=element_blank(),
              axis.ticks.x=element_line(),panel.border=element_rect(colour="black",fill=NA),
              strip.text=element_text(size=20),strip.background=element_rect("white"),
              axis.title=element_text(size=17),axis.text.x=element_text(angle=45,hjust=1))
  
  print(gp)
  dev.off()
}
  


############################################## JUNK







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