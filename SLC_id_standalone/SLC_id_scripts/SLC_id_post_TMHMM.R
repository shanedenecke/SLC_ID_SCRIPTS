#!/usr/bin/Rscript

### Import packages
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))
shhh(library(ggplot2))
shhh(library(ggsci))
shhh(library(stringi))
shhh(library(seqinr))



### functions
shane.transpose=function(dt,newcol){
  new.rows=colnames(dt)[2:length(colnames(dt))]
  a=transpose(dt)
  colnames(a)=as.character(a[1,])
  a=a[-1,]
  b=mutate(a,newcol=new.rows) %>% select(newcol,everything())
  return(b)
}


### functions
dash.remove=function(x){
  split=unlist(strsplit(x,'_'))
  toge=paste(split[1],split[2],sep='_')
  final=toge
  return(final)
}


### Read in HMM files and family names
human.hmm=fread(paste0(sourcePath,'/SLC_id_reference/Human_SLC_HMM.txt'))
colnames(human.hmm)=c('code','tm_domains')

dros.hmm=fread(paste0(sourcePath,'/SLC_id_reference/Drosophila_Flybase_SLC_TMHMM.csv'))
colnames(dros.hmm)=c('code','tm_domains','family')

slc_fams=readLines(paste0(sourcePath,'/SLC_id_reference/SLC_Families.txt'))
slc_fams=sapply(slc_fams,dash.remove)
names(slc_fams)=NULL




###Generate HMM filter tables
human.hmm$family=sapply(human.hmm$code,dash.remove)
model.hmm=rbindlist(list(human.hmm,dros.hmm))
model.hmm.key=human.hmm %>% group_by(family) %>% summarize(minimum=max(0,min((min(tm_domains)-2),(min(tm_domains)/2)))) %>% data.table()
model.hmm.key=rbindlist(list(model.hmm.key,data.table(family='SLC_Unsorted',minimum=0)),use.names = T)



### Import argument for metadata file and BUSCO files
args <- commandArgs(trailingOnly = F) 
#args[1]=paste0(sourcePath,'/SLC_id_reference/Arthropod_species_metadata.tsv') 
meta.data=fread(args[1],header=T,sep='\t')
busco.data=fread('./BUSCO/BUSCO_final_summary.tsv') %>% rename(abbreviation=Species)
meta.full=merge(meta.data,busco.data,by='abbreviation')
meta.full=rbind(meta.full,meta.data[abbreviation=='DroMel'],fill=T) %>% rbind(meta.data[abbreviation=='HomSap'],fill=T)

### Set working directories and source directory
setwd('/mnt/disk/shane/Transporter_ID/SLC_id_pipeline')
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
scriptPath='/mnt/disk/shane/Transporter_ID/SLC_id_pipeline/SLC_id_standalone/SLC_id_scripts'
sourcePath=dirname(scriptPath)

### Read in ulfiltered dictionary from SLC_id_Pre_THMHM.R
unfiltered_dict=fread('./TMHMM_filter/Renamed_unfiltered_SLC_dict.csv')
arth.hmm=fread('./TMHMM_filter/SLC_TMHMM_scores.txt',header=F,sep='\t') 
colnames(arth.hmm)=c('name','tm_count')
arth.hmm=separate(data=arth.hmm,col=name,sep='___',into=c('Species_name','abbreviation','code','family'))
arth.hmm.annot=merge(arth.hmm,meta.full,by=c('Species_name','abbreviation'))
arth.hmm.annot$family=sapply(arth.hmm.annot$family,dash.remove)

#### Filter for SLCs meeting the TM minimum threshold
g.l=list()
b.l=list()
for(i in 1:nrow(model.hmm.key)){
  row=model.hmm.key[i]
  fam=row$family
  mini=row$minimum
  
  sub=arth.hmm.annot[family==fam]
  
  good=sub[tm_count>=mini]
  #good=sub[tm_count>=mini | tm_count>4]
  g.l[[i]]=good
  
  bad=sub[tm_count<mini & tm_count<=4]
  b.l[[i]]=bad
}
tmhmm.filtered.full=rbindlist(g.l)
tmhmm.filtered.full=tmhmm.filtered.full %>% mutate(name=paste0(Species_name,'___',abbreviation,'___',code,'___',family,'_')) %>% data.table()
tmhmm.removed=rbindlist(b.l)

### Do some additional post processing
#dros.addback=tmhmm.removed[abbreviation=='DroMel'] ## add back the one thing that was filtered out in DroMel
#tmhmm.filtered.full=rbindlist(list(tmhmm.filtered.full,dros.addback))
count.summary=tmhmm.filtered.full %>% group_by(family,abbreviation) %>% summarize(count=length(family)) %>% spread(key=family,value=count) %>% data.table()
count.summary[is.na(count.summary)]=0
#select(count.summary,matches('SLC')) %>% as.matrix()  
count.summary$SLC_total=rowSums(select(count.summary,matches('SLC')))


### Write out a fasta file with all of the SLCs found in all of your proteomes
writeLines(tmhmm.filtered.full$name,'./TMHMM_filter/Filtered_SLC_codes.txt')
system(paste0(scriptPath,'/unigene_fa_sub.sh ./TMHMM_filter/Renamed_unfiltered_SLC.faa ./TMHMM_filter/Filtered_SLC_codes.txt > ./Final_raw_outputs/Total_Fasta_sequences.faa'))

##write totals to files
fwrite(count.summary,'./Final_raw_outputs/Total_count_summary.csv')
fwrite(tmhmm.removed,'./Final_raw_outputs/TMM_filtered_out.csv')
fwrite(tmhmm.filtered.full,'./Final_raw_outputs/Total_dictionary.csv')
#save(fa.final,file='./Final_raw_outputs//All_SLCs_fasta.Robj')
  

### output final SLC dicts for each species
out.base=tmhmm.filtered.full %>% select(abbreviation,code,family) %>% mutate(name=paste0(abbreviation,'__',code,'__',family,'_')) %>% 
  select(-family) %>% unique.data.frame() %>% data.table()  

dir.create('real_final_SLC_dicts')
lapply(split(out.base,out.base$abbreviation),function(x) fwrite(select(x,-abbreviation),file=paste0('./real_final_SLC_dicts/',x$abbreviation[1],'_final_SLC_table.csv')))


