library(dplyr)
library(data.table)
library(stringr)
library(ape)

setwd('/data2/shane/Documents/SLC_id')
counts=fread('./Final_raw_outputs/TableS4_count_summary.csv')
oly=fread('./general_reference/Co_variables/Arthropod_species_metadata.csv',header=T) %>% 
  select(Species_name,abbreviation)
taxid.codes=fread('./general_reference/non_model_proteomes/keys/OrthoDB_tax_key_3col.tsv',
                  col.names = c('taxid_code','abbreviation','species_name'))
 


abbrev=function(x){
  a=unlist(strsplit(x,'_'))
  both=substr(a,1,3)
  first=both[1]
  sec=toupper(substr(both[2],1,1))
  thir=substr(both[2],2,3)
  final=paste0(first,sec,thir)
  return(final)
}


##### junk
v=c()

for(i in colnames(counts)[2:length(counts)]){
  sub=as.numeric(counts[[i]])
  #if(max(sub)<10){# & sub[173]<.5){
  if(sum(sub)<100 | max(sub)<10){
    v=c(v,i)
    counts[[i]]=NULL
  }
}
counts[['SLC_total']]=NULL
### 

#hem=readLines('./general_reference/CAFE/Hemiptera_species_outgroups.txt')
#lep=readLines('./general_reference/CAFE/Lepidoptera_species_outgroups.txt')
#lep=lep[-grep("Heliconiuserato_lativitta",lep)]
#rep=readLines('./general_reference/CAFE/Arthropod_species.txt')

hem=sapply(readLines('./general_reference/CAFE/Hemiptera_species_outgroups.txt'),abbrev)
lep=sapply(readLines('./general_reference/CAFE/Lepidoptera_species_outgroups.txt'),abbrev)
lep=lep[-grep("HelLat",lep)]
rep=sapply(readLines('./general_reference/CAFE/Arthropod_species.txt'),abbrev)



shane.transpose=function(dt,newcol){
  new.rows=colnames(dt)[2:length(colnames(dt))]
  a=transpose(dt)
  colnames(a)=as.character(a[1,])
  a=a[-1,]
  b=mutate(a,newcol=new.rows) %>% select(newcol,everything())
  return(b)
}

trans.counts=shane.transpose(counts) %>% data.table()

#for(i in 2:length(colnames(trans.counts))){
#  c=colnames(trans.counts)[i]
#  colnames(trans.counts)[i]=oly[abbreviation==c]$Species_name
#}

colnames(trans.counts)[1]='Family ID'
trans.counts$`Family ID`=gsub('_','',trans.counts$`Family ID`)
trans.counts$Desc='(null)'

full.counts=select(trans.counts,Desc,everything())
#colnames(full.counts)[3:length(colnames(full.counts))]=sapply(colnames(full.counts)[3:length(colnames(full.counts))],abbrev)
lepi=full.counts %>% select(c('Desc','Family ID',lep)) #%>% filter(`Family ID`=='SLC_22')
hemi=full.counts %>% select(c('Desc','Family ID',hem)) #%>% filter(`Family ID`=='SLC_33')
repi=full.counts %>% select(c('Desc','Family ID',rep)) #%>% filter(`Family ID`=='SLC_22')

fwrite(full.counts,'./CAFE/CAFE_FULL.tsv',sep='\t')
    
fwrite(lepi,'./CAFE/Lepidotperan_SLC_CAFE_table.tsv',sep='\t')
fwrite(hemi,'./CAFE/Hemipteran_SLC_CAFE_table.tsv',sep='\t')
fwrite(repi,'./CAFE/Arthropod_SLC_CAFE_table.tsv',sep='\t')
  


######################### CLEAN UP TRE FILES


l.tree=readLines('./general_reference/CAFE/Lepi_RAxML_bipartitions.concat_genes.fs.aln.trimmed.phy.tre')
h.tree=readLines('./general_reference/CAFE/Hemi_RAxML_bipartitions.concat_genes.fs.aln.trimmed.phy.tre')
a.tree=readLines('./general_reference/CAFE/Arth_RAxML_bipartitions.concat_genes.fs.aln.trimmed.phy.tre')
full_list=c(l.tree,h.tree,a.tree)
names(full_list)=c('Lepi','Hemi','Arth')


taxid.codes$cafe_code=paste(taxid.codes$taxid_code,'0',sep='_')
tax_key_sub=taxid.codes[abbreviation %in% unique(c(lep,hem,rep))]

l=list()
for(i in 1:length(full_list)){
  tr=full_list[i]
  #sec=tr
  for(j in taxid.codes$cafe_code){
    if(grepl(j,tr)){
      tr=gsub(j,taxid.codes$abbreviation[which(taxid.codes$cafe_code==j)],tr)
    }
  }
  tr=gsub('7227_0','DroMel',tr)
  l[[i]]=tr
  writeLines(tr,paste0('/data2/shane/Documents/SLC_id/CAFE/raxml_tree',names(full_list)[i],'.tre'))
} 

l.tree=read.tree('./CAFE/raxml_treeLepi.tre')
l.tree.ch=chronopl(l.tree, lambda=0.1)
write.tree(l.tree.ch, file="./CAFE/Lepi_tree_ultrametric.tre")

l.tree=read.tree('./CAFE/raxml_treeHemi.tre')
l.tree.ch=chronopl(l.tree, lambda=0.1)
write.tree(l.tree.ch, file="./CAFE/Hemi_tree_ultrametric.tre")


l.tree=read.tree('./CAFE/raxml_treeArth.tre')
l.tree.ch=chronopl(l.tree, lambda=0.1)
write.tree(l.tree.ch, file="./CAFE/Arth_tree_ultrametric.tre")

