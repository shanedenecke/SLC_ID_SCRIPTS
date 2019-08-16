library(data.table)
library(dplyr)
library(tidyr)

dros.slc.dict=fread('~/Documents/SLC_id/DroMel_Database/SLC_source_dict.csv',col.names = c('name','code'))
dros.key=fread('~/Documents/SLC_id/general_reference/keys/Dm_master_key_FB_fasta.csv') %>% select(Dm_FBgn,name) %>% unique.data.frame()


colnames(dros.slc.dict)=c('Dm_FBgn','SLC')
dros.slc.dict=dros.slc.dict %>% separate(SLC,into=c('slc','fam','no'),sep='_') %>% unite(SLC,slc,fam) %>% select(-no)


new.dict=merge(dros.slc.dict,dros.key,by='Dm_FBgn') %>% unite(full_name,SLC,name)
colnames(new.dict)=c('code','name')
fwrite(new.dict,'~/Documents/SLC_id/SLC_align/updated_Dros_names.csv',row.names = F)



