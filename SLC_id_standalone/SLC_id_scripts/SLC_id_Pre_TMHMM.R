#!/usr/bin/Rscript

### Import packages
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))

### functions
dash.remove=function(x){
  split=unlist(strsplit(x,'_'))
  toge=paste(split[1],split[2],sep='_')
  final=toge
  return(final)
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
unfilterd_dict=rbindlist(l,use.names = T)
fwrite(unfilterd_dict,'./TMHMM_filter/Renamed_unfiltered_SLC_dict.csv')

### remove any fasta files that didn't get renamed properly
unnamed.fasta=read.fasta('./TMHMM_filter/Renamed_unfiltered_SLC.faa',seqtype='AA',as.string=T,set.attributes = F,strip.desc=T)
unnamed.fasta2=unnamed.fasta[grepl('__',names(unnamed.fasta))]
write.fasta(unnamed.fasta2,names(unnamed.fasta2),file.out='./TMHMM_filter/Renamed_unfiltered_SLC.faa')

