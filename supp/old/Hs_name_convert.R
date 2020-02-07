library(dplyr)
library(data.table)
library(biomaRt)
library(tidyr)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
hs_genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), mart = ensembl)


slc.names=readLines('~/Documents/SLC_id/general_reference/HGNC_names.txt')

convert=hs_genes %>% filter(hgnc_symbol %in% slc.names) %>% data.table()


ortho.table=fread('~/Dropbox/wp4_drug_target_search/target_identify/helicoverpa/ref_files/Ha_orthology.txt',
                  sep='\t',header=F,stringsAsFactors = F, colClasses = "character",col.names = c('group','id')) %>% 
                  separate(id,into=c('taxid','geneid'),sep="_")


slc.groups=ortho.table[geneid %in% convert$ensembl_gene_id]$group %>% unique()
sp=c('9606','29058','7227')

b=ortho.table %>% filter(group %in% slc.groups) %>% filter(taxid %in% sp) %>% data.table()


for(i in unique(b$group)){
  sub=subset(b,group==i)


View(b)
