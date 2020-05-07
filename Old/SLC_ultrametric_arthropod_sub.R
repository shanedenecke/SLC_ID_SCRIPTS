library(data.table)
library(dplyr)

#setwd('/data2/shane/Documents/SLC_id/ultrametric_tree')

args = commandArgs(trailingOnly=TRUE)
H=as.character(args[1])


odb10.met=fread(paste0(H,'/GENERAL_REFERENCE/CAFE/odb10v0_OG2genes.33208.tab')) %>%
colnames()

tax=fread('./odb10_taxid.txt')
genes=fread('./odb10_genes.txt',header=F)
groups=fread('./odb10_groups.txt',header=F)
oly=fread('/data2/shane/Documents/SLC_id/general_reference/Co_variables/Olympia_table_august_2019_shaemod.csv',header=T)

tax.key=fread('/data2/shane/Documents/SLC_id/general_reference/non_model_proteomes/keys/Taxid_master_key_full.tsv')

ids=tax.key[V2 %in% oly$abbreviation]$V1

index=tax$V1 %in% ids

ids.sub=tax[index]
gene.sub=genes[index]
group.sub=groups[index]

colnames(ids.sub)='taxid'
colnames(gene.sub)='gene'
colnames(group.sub)='OG'

final=data.table(ids.sub,gene.sub,group.sub)

fwrite(final,'./OrthoDB10_clean_arthropod_suset.tsv',sep='\t')
