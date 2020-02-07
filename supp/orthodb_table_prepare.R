library(dplyr)
library(data.table)

a=list.files('/data2/shane/Documents/SLC_id/general_reference/non_model_proteomes/raw_fasta/')
orthodb170=gsub('_0.fs','',a)

b=fread('/data2/shane/Documents/SLC_id/general_reference/non_model_proteomes/keys/Taxid_master_key_full.tsv')
mycodes=b$V1 %>% as.character()



d=fread('/data2/shane/Documents/SLC_id/general_reference/non_model_proteomes/176_insect_species_genes.csv')
e=readLines('/data2/shane/Documents/SLC_id/general_reference/non_model_proteomes/keys/OrthoDB_tax_key_taxid_only.tsv')

f=d[V2 %in% e]

fwrite(f,'/data2/shane/Documents/SLC_id/general_reference/non_model_proteomes/168_insect_orthoDB.tsv',sep='\t',col.names = F)