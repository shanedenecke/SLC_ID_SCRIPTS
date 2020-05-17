#!/usr/bin/Rscript


######################### SETUP

### Import packages
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(seqinr))
shhh(library(argparser))

### functions
dash.remove=function(x){
  split=unlist(strsplit(x,'_'))
  toge=paste(split[1],split[2],sep='_')
  final=toge
  return(final)
}


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




### Debugging
#setwd('/mnt/disk/shane/Transporter_ID/SLC_id_pipeline')
#scriptPath='/mnt/disk/shane/Transporter_ID/SLC_id_pipeline/SLC_ID_SCRIPTS/SLC_id_standalone/SLC_id_scripts'
#sourcePath=dirname(scriptPath)
#argv$meta='./GENERAL_REFERENCE/keys/Arthropod_species_metadata.tsv'

print(argv$meta)
print(sourcePath)
print(scriptPath)

#### Create new directories for outputs
dir.create('TMHMM_filter')

## read in metadata
meta.data=fread(argv$meta,header=T,sep='\t')


#### Copy Drosphila and Human dictionaries and proteomes into directory with the other targets 
file.remove('./preliminary_SLC_dicts/DroMelPreliminary_SLC_table.csv')
file.remove('./preliminary_SLC_dicts/HomSapPreliminary_SLC_table.csv')
file.copy('./DroMel_Database/SLC_source_dict.csv','./preliminary_SLC_dicts/DroMelPreliminary_SLC_table.csv')
file.copy('./HomSap_Database/SLC_source_dict.csv','./preliminary_SLC_dicts/HomSapPreliminary_SLC_table.csv')


file.remove('./filtered_proteomes/DroMel_unigene.faa')
file.remove('./filtered_proteomes/HomSap_unigene.faa')
file.copy(paste0(sourcePath,'/SLC_id_reference/DroMel_unigene.faa'),'./filtered_proteomes/')
file.copy(paste0(sourcePath,'/SLC_id_reference/HomSap_unigene.faa'),'./filtered_proteomes/')


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
  system(paste0(scriptPath,'/fasta_rename.py ./TMHMM_filter/SLC_unfiltered_all_raw.faa ./TMHMM_filter/temp_rename_dict.csv >> ./TMHMM_filter/Renamed_unfiltered_SLC.faa 2> TM_errors.txt'))
}
unfilterd_dict=rbindlist(l,use.names = T)
fwrite(unfilterd_dict,'./TMHMM_filter/Renamed_unfiltered_SLC_dict.csv')


