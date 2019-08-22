library(data.table)
library(dplyr)
library(ggplot2)
library(seqinr)


## Set Working directory
setwd('/data2/shane/Documents/SLC_id')


#### Create dictionary 

## table with some info on each species
co.var=fread('./general_reference/Co_variables/Olympia_table_august_2019_shaemod.csv',header=T) %>%
  select(Species_name,abbreviation,Taxanomic_Classification,Phagy,Vory,Diet_category)

######################## functions
dash.remove=function(x){
  split=unlist(strsplit(x,'_'))
  toge=paste(split[1],split[2],sep='_')
  final=toge
  return(final)
}
#########################


l=list()
for(i in list.files('./final_SLC_dicts',full.names = T)){
  dict=fread(i)
  
  abbrev=gsub('./final_SLC_dicts/','',i,fixed=T)
  abbrev=gsub('Final_SLC_table.csv','',abbrev,fixed=T)
  
  fam=sapply(dict$name,dash.remove)
  
  dict$abbreviation=abbrev
  dict$family=fam
  
  l[[i]]=dict
}

all=rbindlist(l) %>% merge(co.var,by='abbreviation')
fasta.name=paste(all$Species_name,all$abbreviation,all$code,all$name,sep='_')
rename.dict=data.table(name=fasta.name,code=all$code)
#rename.dict$name=paste(rename.dict$name,rena)
all$Species_name=gsub(' ','_',all$Species_name)

fwrite(all,'./shiny_prep/Reference_csv_dictionary.csv')
fwrite(rename.dict,'./shiny_prep/Rename_SLC_dict.csv')
writeLines(rename.dict$code,'shiny_prep/slc_codes.txt')


##################### FASTA

system('
cd /data2/shane/Documents/SLC_id
cat ./proteomes/* ./general_reference/model_proteomes/*.faa > ./shiny_prep/all_proteomes.faa
/data2/shane/Applications/custom/unigene_fa_sub.sh ./shiny_prep/all_proteomes.faa ./shiny_prep/slc_codes.txt > ./shiny_prep/SLC_all_raw.faa
/data2/shane/Applications/custom/fasta_rename.py ./shiny_prep/SLC_all_raw.faa ./shiny_prep/Rename_SLC_dict.csv > ./shiny_prep/Renamed_SLC.faa
')

raw.fa=read.fasta('./shiny_prep/Renamed_SLC.faa',seqtype='AA',as.string=T,set.attributes = F,strip.desc=T)



