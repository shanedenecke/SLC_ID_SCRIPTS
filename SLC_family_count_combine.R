#!/usr/bin/Rscript

shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))
shhh(library(ggplot2))
shhh(library(ggsci))
shhh(library(stringi))


#setwd('/home/shanedenecke/Dropbox/wp7_prodrug/SLC_id') #local
setwd('/data2/shane/Documents/SLC_id')
dir.create('SLC_family_counts')
co.variables=fread('./general_reference/Co_variables/Olympia_table_august_2019_shaemod.csv',header=T)
slc.function=fread('/data2/shane/Documents/SLC_id/general_reference/SLC_info/SLC_function_groups.csv')
#co.variables=fread('/home/shanedenecke/Dropbox/wp7_prodrug/SLC_id/family_size_variation/Olympia_table_august_2019.csv',header=T) ##local


### copy DroMel and HomSap databases 
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
#slc_fams=readLines('/home/shanedenecke/Dropbox/wp7_prodrug/SLC_id/old/general_reference/SLC_info/SLC_families.txt') ### local
slc_fams=sapply(slc_fams,dash.remove)
names(slc_fams)=NULL

l=list()

##Get all counts data into master table. Rows SLC families. Columns species
for (i in list.files('./final_SLC_dicts/')){
  slc_table=fread(paste0('./final_SLC_dicts/',i))
  #slc_table=fread('/home/shanedenecke/Dropbox/wp7_prodrug/SLC_id/HomSap_SLC_dict.csv')
  abbreviation=gsub('Final','',unlist(strsplit(i,'_'))[1])
  slc_table$family=sapply(slc_table$name,dash.remove)
  with.count=slc_table %>% group_by(family) %>% summarise(count=length(family))
  colnames(with.count)[2]=abbreviation
  l[[i]]=with.count
}
## combine list into count.summary variable and transpose to revers rows/cols
count.summary=Reduce(function(x, y) merge(x, y, ,by='family',all=TRUE), l)  %>% 
  shane.transpose() %>%
  mutate_at(vars(matches("SLC")), funs(as.numeric(as.character(.)))) %>%
  data.table()
count.summary[is.na(count.summary)]=0

## calculate totals
count.summary$slc_total=rowSums(count.summary[,2:68])
colnames(count.summary)[1]='abbreviation'
count.summary=select(count.summary,abbreviation,everything())

##Produce histogram of total SLC famiy size
hist(as.numeric(count.summary$slc_total))


##write totals to file
fwrite(count.summary,'./SLC_family_counts/count_summary.csv')


################################# COUNTS ANALYSIS

aa=slc.function$AA
sugar=stri_remove_empty(slc.function$Sugar)
drug=stri_remove_empty(slc.function$Drug)
ion=stri_remove_empty(slc.function$Ion)

anal.counts=count.summary %>% select(-slc_total,-SLC_X,-SLC_Unsorted)

anal.counts$SLC_AA=anal.counts %>% select(aa) %>% rowSums()
anal.counts$SLC_sugar=anal.counts %>% select(sugar) %>% rowSums()
anal.counts$SLC_drug=anal.counts %>% select(drug) %>% rowSums()
anal.counts$SLC_ion=anal.counts %>% select(ion) %>% rowSums()




############### filter out SLC families which have <100 total members in dataset 
for(i in colnames(anal.counts)[2:70]){
  sub=as.numeric(anal.counts[[i]])
  #if(max(sub)<10){# & sub[173]<.5){
  if(sum(sub)<100 & max(sub)<10){
    print(i)
    anal.counts[[i]]=NULL
  }
}

mean.fam.size=apply(select(anal.counts,matches('SLC')),2,mean)
size.adjust=data.table(fam=names(mean.fam.size),avg=mean.fam.size)


## merge to olympia's data
full.counts=merge(anal.counts,co.variables,by="abbreviation")


## set varaibles for loop
comps=c('Taxanomic_Classification','Diet_category')
l=list()

### Perform loop which goes through each SL family and comparison and builds linear model for each
for(i in colnames(full.counts)[grep('SLC_',colnames(full.counts))]){
  for(j in comps){
    
    ##Filter for infrequent categories
    good=names(which(table(full.counts[[j]])>5)) ## filter for occurances which don't occur at least 5 times
    sub=full.counts[which(full.counts[[j]] %in% good)] ## subset for comparisons which have at least 5 occurances
    #model=aov(formula=counts~sub[[j]])
    model=lm(formula=sub[[i]]~sub[[j]])
    coef=summary(model)$coefficients 
    comparison=rownames(coef)
  l[[paste(i,j,sep='_')]]= coef %>% data.table(comp=comparison,family=i,co_variable=j) 
}}

## merge all lists into large data table with all associations
raw.model=rbindlist(l) %>% select(Estimate,`Pr(>|t|)`,comp,family,co_variable) %>% filter(!grepl("Intercept",comp)) %>% data.table() 
colnames(raw.model)=gsub('Pr(>|t|)','p.value',colnames(raw.model),fixed=T) 
raw.model$comp=gsub('sub[[j]]','',raw.model$comp,fixed=T)

## Add adjusted estimate and p value corrections
for(i in 1:nrow(raw.model)){raw.model$adjusted_Estimate[i]=raw.model[i]$Estimate/(size.adjust[fam==raw.model[i]$family])$avg}
raw.model$fdr=p.adjust(raw.model$p.value,method='fdr') 
raw.model$bonf=p.adjust(raw.model$p.value,method='bonferroni') 


final=arrange(raw.model,bonf) %>% filter(bonf<1e-10)  %>% filter(abs(adjusted_Estimate)>1) %>% filter(Estimate>5) %>% data.table()
fwrite(raw.model,'./SLC_family_counts/full_correlation_anlaysis.csv')
fwrite(final,'./SLC_family_counts/significant_correlation_anlaysis.csv')


dir.create('comparison_plots')

for(i in 1:nrow(final)){
  row=final[i]
  co.var=row$co_variable
  fam=row$family
  red=select(full.counts,fam,co.var)
  red[[fam]]=as.numeric(red[[fam]])
  
  plot=red[red[[co.var]] %in% names(which(table(red[[co.var]])>5))] ### remove elements that aren't present at least 5 times
  
  ylab=paste0('Number of ',gsub('_',' ',fam),' Members','\n')
  tit=paste0(gsub('_',' ',fam),' Family Size vrs. ',gsub('_',' ',co.var))
  xlab=paste0('\n',gsub('_',' ',co.var))
  
  
  pdf(paste0('./comparison_plots/',co.var,'_',fam,'.pdf'))
  
  gp=ggplot(plot,aes_string(co.var,y=fam,fill=co.var))
  gp=gp+geom_boxplot()
  gp=gp+labs(x=xlab,y=ylab)
  gp=gp+scale_fill_rickandmorty()
  gp=gp+ggtitle(tit)
  gp=gp+theme_bw()
  gp=gp+theme(text=element_text(face="bold",family="serif"),panel.grid=element_blank(),
              axis.ticks.x=element_line(),panel.border=element_rect(colour="black",fill=NA),
              strip.text=element_text(size=20),strip.background=element_rect("white"),
              axis.title=element_text(size=17),axis.text.x=element_text(angle=45,hjust=1),
              legend.position = 'none',plot.title = element_text(hjust = 0.5))
  
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









fwrite(g,'./SLC_family_counts/TOTAL_FAMILY_COUNTS.csv',row.names = F)



#h=dcast(melt(g, id.vars = "family"), variable ~ family) 
#colnames(h)[1]='SLC_family'