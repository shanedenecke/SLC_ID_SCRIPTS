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
ce.counts=fread('/data2/shane/Documents/SLC_id/general_reference/CAFE/Ce_counts.csv')
sup.counts=fread('/data2/shane/Documents/OrthoDB_Mapping/OrthoDB_10_family_counts/CAFE_orthoDB_supp.tsv') %>% 
  select(-MenMol,-LimLun,-ChiSup,-HelLat,-PacVen,-EurMay,-CerMar)


abbrev=function(x){
  a=unlist(strsplit(x,'_'))
  both=substr(a,1,3)
  first=both[1]
  sec=toupper(substr(both[2],1,1))
  thir=substr(both[2],2,3)
  final=paste0(first,sec,thir)
  return(final)
}


shane.transpose=function(dt,newcol){
  new.rows=colnames(dt)[2:length(colnames(dt))]
  a=transpose(dt)
  colnames(a)=as.character(a[1,])
  a=a[-1,]
  b=mutate(a,newcol=new.rows) %>% select(newcol,everything())
  return(b)
}

lambda.convert=function(x){
  a=gsub("[A-z]+:[0-9]+.[0-9]+",1,x)
  b=gsub("[0-9]{2,}:0.[0-9]+",1,a)
  c=gsub(":[0-9]+.[0-9]+",1,b)
  d=gsub('[A-z]+:[0-9]+',1,c)
  return(d)
}



#### clean up coutns
######## CLEAN UP CAFE TABLES
counts2=rbindlist(list(counts,ce.counts),use.names = T,fill=T)


##################################### ONLY RUN FOR FILTER
counts2[is.na(counts2)]=0
v=c()

for(i in colnames(counts2)[2:length(counts2)]){
  sub=as.numeric(counts2[[i]])
  if(max(sub)-min(sub)>300 | sum(sub)<60){
    #if(sum(sub)<100 | max(sub)<10){
    v=c(v,i)
    counts2[[i]]=NULL
  }
}
counts2[['SLC_total']]=NULL

trans.counts=shane.transpose(counts2) %>% data.table()

#for(i in 2:length(colnames(trans.counts))){
#  c=colnames(trans.counts)[i]
#  colnames(trans.counts)[i]=oly[abbreviation==c]$Species_name
#}

colnames(trans.counts)[1]='Family ID'
trans.counts$`Family ID`=gsub('_','',trans.counts$`Family ID`)
trans.counts$Desc='(null)'
#trans.counts=rbindlist(list(trans.counts,sup.counts),use.names = T,fill=T)
full.counts=select(trans.counts,Desc,everything())

fwrite(full.counts,'./CAFE/CAFE_tables/CAFE_FULL.tsv',sep='\t')

#### run loop

######################### CLEAN UP TRE FILES
iter=list.files('/data2/shane/Documents/SLC_id/CAFE/trees')[grepl('raxml_tree',list.files('/data2/shane/Documents/SLC_id/CAFE/trees'))] %>% 
  {.[!str_detect(., "named_")]} %>%
  str_remove('raxml_tree_') %>% str_remove('.tre')
taxid.codes$cafe_code=paste(taxid.codes$taxid_code,'0',sep='_')
taxid.codes=rbind(taxid.codes,data.table(taxid_code=6239,abbreviation='CaeEle',species_name='Caenorhabditis elegans',cafe_code="6239_0"))
taxid.codes=rbind(taxid.codes,data.table(taxid_code=29058,abbreviation='HelArm',species_name='Helicoverpa armigera',cafe_code="29058_0"))
taxid.codes=rbind(taxid.codes,data.table(taxid_code=7165,abbreviation='AnoGam',species_name='Anopheles gambiae',cafe_code="7165_0"))



for (i in iter){
  
  ### rename tree
  rax.tree=readLines(paste0('/data2/shane/Documents/SLC_id/CAFE/trees/raxml_tree_',i,'.tre'))
  for(j in taxid.codes$cafe_code){
    if(grepl(j,rax.tree)){
      rax.tree=gsub(j,taxid.codes$abbreviation[which(taxid.codes$cafe_code==j)],rax.tree)
    }
  }
  rax.tree.name=gsub('7227_0','DroMel',rax.tree)
  writeLines(rax.tree.name,paste0('/data2/shane/Documents/SLC_id/CAFE/trees/raxml_tree_named_',i,'.tre'))
  
 
  ## create ultrametric tree
  l.tree.ch=chronopl(read.tree(paste0('/data2/shane/Documents/SLC_id/CAFE/trees/raxml_tree_named_',i,'.tre')), lambda=0.1)
  l.tree.ch$edge.length=l.tree.ch$edge.length*1000
  l.tree.ch$node.label=NULL
  write.tree(l.tree.ch, file=paste0("./CAFE/trees/",i,'_tree_ultrametric.tre'))
  
  ## create lambda file
  writeLines(lambda.convert(readLines(paste0("./CAFE/trees/",i,'_tree_ultrametric.tre'))),paste0("./CAFE/trees/",i,'_tree_lambda.txt')) 
  
  
  sp=str_extract_all(readLines(paste0("./CAFE/trees/",i,'_tree_ultrametric.tre')),pattern = "[A-z]+",simplify = T)
  #sp_clean=sp[sp='_']
  
  sp.counts=full.counts %>% select(c('Desc','Family ID',sp))
  fwrite(sp.counts,paste0('./CAFE/CAFE_tables/',i,'_SLC_CAFE_table.tsv'),sep='\t')
}
  

