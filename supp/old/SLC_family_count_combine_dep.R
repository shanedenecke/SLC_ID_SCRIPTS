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
file.remove('/data2/shane/Documents/SLC_id/final_SLC_dicts/DroMelFinal_SLC_table.csv')
file.remove('/data2/shane/Documents/SLC_id/final_SLC_dicts/HomSapFinal_SLC_table.csv')
file.copy('./Dm_Database_Generate/SLC_source_dict_flybaseXref.csv','/data2/shane/Documents/SLC_id/final_SLC_dicts/DroMelFinal_SLC_table.csv')
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


#### NEED to rewrite based on TMHMM values
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
count.summary$slc_total=rowSums(count.summary[,2:67])
colnames(count.summary)[1]='abbreviation'
count.summary=select(count.summary,abbreviation,everything())

##write totals to file
fwrite(count.summary,'./SLC_family_counts/count_summary.csv')

