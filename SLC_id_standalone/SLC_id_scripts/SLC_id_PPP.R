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
dash.remove=function(x){
  split=unlist(strsplit(x,'_'))
  toge=paste(split[1],split[2],sep='_')
  final=toge
  return(final)
}

shane.transpose=function(dt,newcol){
  new.rows=colnames(dt)[2:length(colnames(dt))]
  a=transpose(dt)
  colnames(a)=as.character(a[1,])
  a=a[-1,]
  b=mutate(a,newcol=new.rows) %>% select(newcol,everything())
  return(b)
}

### Set working directories and source directory
setwd('/mnt/disk/shane/Transporter_ID/SLC_id_pipeline')
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
scriptPath='/mnt/disk/shane/Transporter_ID/SLC_id_pipeline/SLC_id_standalone/SLC_id_scripts'
sourcePath=dirname(scriptPath)

### Import argument for metadata file
args <- commandArgs(trailingOnly = F) 
args[1]=paste0(sourcePath,'/SLC_id_reference/Arthropod_species_metadata.tsv') 

#### Create new directories for outputs
dir.create('TMHMM_filter')
dir.create('Final_raw_outputs')
dir.create('Figures')

## read in metadata
meta.data=fread(args[1],header=T,sep='\t')
busco.data=fread('./BUSCO/BUSCO_final_summary.tsv') %>% rename(abbreviation=Species)
meta.full=merge(meta.data,busco.data,by='abbreviation') %>% filter(Completeness>75) %>% data.table()
meta.full=rbind(meta.full,meta.data[abbreviation=='DroMel'],fill=T) %>% rbind(meta.data[abbreviation=='HomSap'],fill=T)


### Read in HMM files and family names
human.hmm=fread(paste0(sourcePath,'/SLC_id_reference/Human_SLC_HMM.txt'))
colnames(human.hmm)=c('code','tm_domains')

dros.hmm=fread(paste0(sourcePath,'/SLC_id_reference/Drosophila_Flybase_SLC_TMHMM.csv'))
colnames(dros.hmm)=c('code','tm_domains','family')

slc_fams=readLines(paste0(sourcePath,'/SLC_id_reference/SLC_Families.txt'))
slc_fams=sapply(slc_fams,dash.remove)
names(slc_fams)=NULL


#### Copy Drosphila and Human dictionaries and proteomes into directory with the other targets 
file.remove('./preliminary_SLC_dicts/DroMelPreliminary_SLC_table.csv')
file.remove('./preliminary_SLC_dicts/HomSapPreliminary_SLC_table.csv')
file.copy('./DroMel_Database/SLC_source_dict.csv','./preliminary_SLC_dicts/DroMelPreliminary_SLC_table.csv')
file.copy('./HomSap_Database/SLC_source_dict.csv','./preliminary_SLC_dicts/HomSapPreliminary_SLC_table.csv')


file.remove('./filtered_proteomes/DroMel_unigene.faa')
file.remove('./filtered_proteomes/HomSap_unigene.faa')
file.copy(paste0(sourcePath,'/SLC_id_reference/DroMel_unigene.faa'),'./filtered_proteomes/')
file.copy(paste0(sourcePath,'/SLC_id_reference/HomSap_unigene.faa'),'./filtered_proteomes/')



######################################### START THE ANALYSIS PROPER


###Generate HMM filter tables
human.hmm$family=sapply(human.hmm$code,dash.remove)
model.hmm=rbindlist(list(human.hmm,dros.hmm))
model.hmm.key=human.hmm %>% group_by(family) %>% summarize(minimum=max(0,min((min(tm_domains)-2),(min(tm_domains)/2)))) %>% data.table()
model.hmm.key=rbindlist(list(model.hmm.key,data.table(family='SLC_Unsorted',minimum=0)),use.names = T)


########## GET PRELIMINARY LIST AND SUBSET FASTA
l=list()
for(i in list.files('./preliminary_SLC_dicts',full.names = T)){
  dict=fread(i)
  abbrev=gsub('./preliminary_SLC_dicts/','',i,fixed=T)
  abbrev=gsub('Preliminary_SLC_table.csv','',abbrev,fixed=T)
  fam=sapply(dict$name,dash.remove)
  dict$abbreviation=abbrev
  dict$family=fam
  dict2=merge(dict,select(meta.data,Species_name,abbreviation),by='abbreviation',use.names=T) %>% mutate(name)
  dict2$name=paste(dict2$Species_name,dict2$abbreviation,dict2$code,dict2$name,sep='___')
  dict2$name=gsub(' ','_',dict2$name)
  dict2$code=as.character(dict2$code)
  l[[i]]=dict2
  ind.rename.dict= dict2 %>% select(name,code)
  writeLines(ind.rename.dict$code,'./TMHMM_filter/slc_unfiltered_codes.txt')
  fwrite(ind.rename.dict,'./TMHMM_filter/temp_rename_dict.csv')
  file.copy(paste0('./filtered_proteomes/',abbrev,'_unigene.faa'),'./TMHMM_filter/temp_proteome.faa',overwrite = T)
  
  system(paste0(scriptPath,'/unigene_fa_sub.sh ./TMHMM_filter/temp_proteome.faa ./TMHMM_filter/slc_unfiltered_codes.txt > ./TMHMM_filter/SLC_unfiltered_all_raw.faa'))
  system(paste0(scriptPath,'/fasta_rename.py ./TMHMM_filter/SLC_unfiltered_all_raw.faa ./TMHMM_filter/temp_rename_dict.csv >> ./TMHMM_filter/Renamed_unfiltered_SLC.faa'))
}

### remove any fasta files that didn't get renamed properly
unnamed.fasta=read.fasta('./TMHMM_filter/Renamed_unfiltered_SLC.faa',seqtype='AA',as.string=T,set.attributes = F,strip.desc=T)
unnamed.fasta2=unnamed.fasta[grepl('__',names(unnamed.fasta))]
write.fasta(unnamed.fasta2,names(unnamed.fasta2),file.out='./TMHMM_filter/Renamed_unfiltered_SLC.faa')

system("tmhmm  ./TMHMM_filter/Renamed_unfiltered_SLC.faa > ./TMHMM_filter/SLC_TMHMM_full.txt")
system(" cat ./TMHMM_filter/SLC_TMHMM_full.txt | grep 'Number of predicted' | perl -pe 's/..(.+) Number of predicted TMHs:\\s+(\\S+)/$1\t$2/g' > ./TMHMM_filter/SLC_TMHMM_scores.txt")


all=rbindlist(l,use.names = T)

arth.hmm=fread('./TMHMM_filter/SLC_TMHMM_scores.txt',header=F,sep='\t') 
colnames(arth.hmm)=c('name','tm_count')
#arth.hmm=arth.hmm[grepl('___',name)]

m=merge(arth.hmm,all,by='name') %>% merge(meta.full,by=c('Species_name','abbreviation'))
rbind(m,)

#m=merge(arth.hmm,all,by='name') 


g.l=list()
b.l=list()
for(i in 1:nrow(model.hmm.key)){
  row=model.hmm.key[i]
  fam=row$family
  mini=row$minimum
  
  sub=m[family==fam]
  
  good=sub[tm_count>=mini]
  #good=sub[tm_count>=mini | tm_count>4]
  g.l[[i]]=good
  
  bad=sub[tm_count<mini & tm_count<=4]
  b.l[[i]]=bad
}
tmhmm.filtered.full=rbindlist(g.l)
tmhmm.removed=rbindlist(b.l)
dros.addback=tmhmm.removed[abbreviation=='DroMel'] ## add back the one thing that was filtered out in DroMel
tmhmm.filtered.full=rbindlist(list(tmhmm.filtered.full,dros.addback))
count.summary=tmhmm.filtered.full %>% group_by(family,abbreviation) %>% summarize(count=length(family)) %>% spread(key=family,value=count) %>% data.table()
count.summary[is.na(count.summary)]=0
#select(count.summary,matches('SLC')) %>% as.matrix()  

count.summary$SLC_total=rowSums(select(count.summary,matches('SLC')))

writeLines(tmhmm.filtered.full$name,'./TMHMM_filter/Filtered_SLC_codes.txt')

system('~/Applications/Custom_Applications/unigene_fa_sub.sh ./TMHMM_filter/Renamed_unfiltered_SLC.faa ./TMHMM_filter/Filtered_SLC_codes.txt > ./Final_raw_outputs/FileS1_All_SLC_final.faa')

##write totals to file
fwrite(count.summary,'./Final_raw_outputs/TableS7_count_summary.csv')
fwrite(tmhmm.removed,'./Final_raw_outputs/TMM_filtered_out.csv')
fwrite(tmhmm.filtered.full,'./Final_raw_outputs/TableS6_Full_dict_table.csv')
#save(fa.final,file='./Final_raw_outputs//All_SLCs_fasta.Robj')
  

a=fread('./Final_raw_outputs/TableS6_Full_dict_table.csv') %>% select(abbreviation,code,family) %>% mutate(name=paste0(abbreviation,'__',code,'__',family,'_')) %>% 
  select(-family) %>% unique.data.frame() %>% data.table()  ##### don't know why removing the family is here
  #unique.data.frame() %>% data.table()

dir.create('real_final_SLC_dicts')
lapply(split(a,a$abbreviation),function(x) fwrite(select(x,-abbreviation),file=paste0('./real_final_SLC_dicts/',x$abbreviation[1],'_final_SLC_table.csv')))


