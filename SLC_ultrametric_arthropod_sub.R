library(data.table)
library(dplyr)

tax=fread('/home/shanedenecke/Dropbox/wp7_prodrug/common_orthoDB/odb10_taxid.txt')
genes=fread('/home/shanedenecke/Dropbox/wp7_prodrug/common_orthoDB/odb10_genes.txt',header=F)
groups=fread('/home/shanedenecke/Dropbox/wp7_prodrug/common_orthoDB/odb10_groups.txt',header=F)
oly=fread('/home/shanedenecke/Dropbox/wp7_prodrug/SLC_id/family_size_variation/Olympia_table_august_2019.csv',header=T)

tax.key=fread('/home/shanedenecke/Dropbox/wp7_prodrug/common_orthoDB/Taxid_OrthoDB_master_key.tsv')

ids=tax.key[V2 %in% oly$abbreviation]$V1

index=tax$V1 %in% ids

ids.sub=tax[index]
gene.sub=genes[index]
group.sub=groups[index]

colnames(ids.sub)='taxid'
colnames(gene.sub)='gene'
colnames(group.sub)='OG'

final=data.table(ids.sub,gene.sub,group.sub)

fwrite(final,'/home/shanedenecke/Dropbox/wp7_prodrug/common_orthoDB/OrthoDB10_clean_arthropod_suset.tsv',sep='\t')
