#!/usr/bin/Rscript


#################### SETUP

### Import packages
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(seqinr))
shhh(library(argparser))


### Set Source directory
getScriptPath <- function(){
    cmd.args <- commandArgs()
    m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
    script.dir <- dirname(regmatches(cmd.args, m))
    if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
    if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
    return(script.dir)
}
scriptPath=getScriptPath()
sourcePath=dirname(scriptPath)

### Get arguments. Only argument is the metadata file with columns "Species_name" and  "abbreviation"
p=arg_parser('SLC_pre_TMHMM')
p <- add_argument(p, "--meta", help="path to metadat")
argv=parse_args(p)


#### Create new directories
dir.create('Final_outputs')

### functions
shane.transpose=function(dt,newcol){
  new.rows=colnames(dt)[2:length(colnames(dt))]
  a=transpose(dt)
  colnames(a)=as.character(a[1,])
  a=a[-1,]
  b=mutate(a,newcol=new.rows) %>% select(newcol,everything())
  return(b)
}

dash.remove=function(x){
  split=unlist(strsplit(x,'_'))
  toge=paste(split[1],split[2],sep='_')
  final=toge
  return(final)
}

last.remove=function(x){
  sp=unlist(strsplit(x,'__'))
  sp2=sp[1:length(sp)-1]
  f=paste(sp2,collapse='___')
  return(f)
}


### Debugging
#setwd('/mnt/disk/shane/Transporter_ID/SLC_id_pipeline')
#scriptPath='/mnt/disk/shane/Transporter_ID/SLC_id_pipeline/SLC_ID_SCRIPTS/SLC_id_standalone/SLC_id_scripts'
#sourcePath=dirname(scriptPath)
#argv$meta='/mnt/disk/shane/Transporter_ID/SLC_id_pipeline/GENERAL_REFERENCE/keys/Arthropod_species_metadata.tsv'


################### IMPORT DATA AND DO SOME BASIC PROCESSING

### Read in HMM files and family names
human.hmm=fread(paste0(sourcePath,'/SLC_id_reference/HOMSAP_SLC_TMHMM.txt'))
colnames(human.hmm)=c('code','tm_domains')

dros.hmm=fread(paste0(sourcePath,'/SLC_id_reference/DROMEL_SLC_TMHMM.csv'))
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
meta.data=fread(argv$meta,header=T,sep='\t')
busco.data=fread('./BUSCO/BUSCO_final_summary.tsv') %>% rename(abbreviation=Species)
meta.full=merge(meta.data,busco.data,by='abbreviation')
meta.full=rbind(meta.full,meta.data[abbreviation=='DroMel'],fill=T) %>% rbind(meta.data[abbreviation=='HomSap'],fill=T)



### Read in ulfiltered dictionary from SLC_id_Pre_THMHM.R
unfiltered_dict=fread('./TMHMM_filter/Renamed_unfiltered_SLC_dict.csv')
arth.hmm=fread('./TMHMM_filter/SLC_TMHMM_scores.txt',header=F,sep='\t',fill=T) 
colnames(arth.hmm)=c('name','tm_count')
arth.hmm=arth.hmm[grepl('__',name)]
arth.hmm=separate(data=arth.hmm,col=name,sep='___',into=c('Species_name','abbreviation','code','family'))
arth.hmm.annot=merge(arth.hmm,meta.full,by=c('Species_name','abbreviation'))
arth.hmm.annot$family=sapply(arth.hmm.annot$family,dash.remove)



#### PERFORM TM FILTERING ON SLCS
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


#### FORMAT AND OUTPUT DICTIONARY
tmhmm.filtered.full=rbindlist(g.l)
tmhmm.filtered.full=tmhmm.filtered.full %>% mutate(name=paste0(Species_name,'___',abbreviation,'___',code,'___',family,'_')) %>% data.table()
tmhmm.removed=rbindlist(b.l)
fwrite(tmhmm.filtered.full,'./Final_outputs/Total_dictionary.csv')
fwrite(tmhmm.removed,'./Final_outputs/TMM_filtered_out.csv')


### OUTPUT COUNTS SUMMARY
count.summary=tmhmm.filtered.full %>% group_by(family,abbreviation) %>% summarize(count=length(family)) %>% spread(key=family,value=count) %>% data.table()
count.summary[is.na(count.summary)]=0
count.summary$SLC_total=rowSums(select(count.summary,matches('SLC')))
transposed.counts=shane.transpose(count.summary,abbreviation)
full.meta.counts=merge(meta.full,count.summary,by='abbreviation')


fwrite(count.summary,'./Final_outputs/Total_count_summary.csv')
fwrite(transposed.counts,'./Final_outputs/Total_count_summary_transpose.csv')
fwrite(full.meta.counts,'./Final_outputs/Full_Metadata_summary.csv')



### WRITE FASTA FILES
writeLines(tmhmm.filtered.full$name,'./TMHMM_filter/Filtered_SLC_codes.txt')
raw.fasta=read.fasta('./TMHMM_filter/Renamed_unfiltered_SLC.faa',as.string = T,forceDNAtolower = F,set.attributes = F)
table.names=as.character(sapply(tmhmm.filtered.full$name,last.remove)) ### remove end bit of names from dicitonary to make match
tmhmm.filtered.full$fasta_name=table.names
filtered.fasta=raw.fasta[as.character(sapply(names(raw.fasta),last.remove)) %in% table.names]
write.fasta(filtered.fasta,names=names(filtered.fasta),nbchar=10000,file.out='./Final_outputs/Total_raw_fasta.faa')

### OUTPUT SPECIES SPECIFIC DICTIONARIES
#out.base=tmhmm.filtered.full %>% select(abbreviation,code,family) %>% mutate(name=paste0(abbreviation,'__',code,'__',family,'_')) %>% 
#  select(-family) %>% unique.data.frame() %>% data.table()  
out.base=tmhmm.filtered.full %>% select(abbreviation,code,family)

dir.create('./Final_outputs/Final_SLC_dicts')
lapply(split(out.base,out.base$abbreviation),function(x) fwrite(select(x,-abbreviation),file=paste0('./Final_outputs/Final_SLC_dicts/',x$abbreviation[1],'_final_SLC_table.csv')))

### OUTPUT SPECIES SPECIFIC FASTA
dir.create('./Final_outputs/Final_SLC_fasta')
for(i in unique(tmhmm.filtered.full$abbreviation)){
  system2(command="grep",args=c(i,'./Final_outputs/Total_raw_fasta.faa'),stdout=paste0('./Final_outputs/Final_SLC_fasta/',i,'SLCs.faa'))
}