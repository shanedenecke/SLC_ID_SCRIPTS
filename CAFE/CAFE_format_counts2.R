shhh <- suppressPackageStartupMessages
shhh(library(dplyr))
shhh(library(data.table))
shhh(library(stringr))
shhh(library(ape))
shhh(library(ggtree))
shhh(library(tidyr))

setwd('/data2/shane/Documents/SLC_id')
counts=fread('./Final_raw_outputs/TableS7_count_summary.csv')
oly=fread('./general_reference/Co_variables/Arthropod_species_metadata.csv',header=T) %>% 
  select(Species_name,abbreviation)
taxid.codes=fread('./general_reference/non_model_proteomes/keys/OrthoDB_tax_key_3col.tsv',
                  col.names = c('taxid_code','abbreviation','species_name'))
taxid.codes$cafe_code=paste(taxid.codes$taxid_code,'0',sep='_')
taxid.codes=rbind(taxid.codes,data.table(taxid_code=6239,abbreviation='CaeEle',species_name='Caenorhabditis elegans',cafe_code="6239_0"))
taxid.codes=rbind(taxid.codes,data.table(taxid_code=29058,abbreviation='HelArm',species_name='Helicoverpa armigera',cafe_code="29058_0"))
taxid.codes=rbind(taxid.codes,data.table(taxid_code=7165,abbreviation='AnoGam',species_name='Anopheles gambiae',cafe_code="7165_0"))
taxid.codes=rbind(taxid.codes,data.table(taxid_code=7227,abbreviation='DroMel',species_name='Drosophila_melanogaster',cafe_code="7227_0"))


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

#x=readLines(paste0("./CAFE/trees/",i,'_tree_ultrametric.tre'))
lambda.convert=function(x){
  a=gsub("[A-z]+:[0-9]+.[0-9]+",1,x)
  b=gsub("[0-9]{2,}:[0-9]+.[0-9]+",1,a)
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
dir.create('CAFE')
dir.create('./CAFE/CAFE_tables')
fwrite(full.counts,'./CAFE/CAFE_tables/CAFE_FULL.tsv',sep='\t')



##### OrthoDB CAFE tables 

#arth.orthodb=fread('/data2/shane/Documents/SLC_id/CAFE/Ultrametric_tree_CAFE/odb_files/odb10v0_OG2genes.33208.tab') %>% 
#  separate(V2,into=c('taxid','geneid'),sep = ':')
arth.orthodb=fread('/data2/shane/Documents/SLC_id/CAFE/Ultrametric_tree_CAFE/odb_files/odb10v0_OG2genes.33208.tab') %>% 
  separate(V2,into=c('taxid','geneid'),sep = ':')

arth.orthodb.sel=arth.orthodb[taxid %in% taxid.codes$cafe_code]
arth.orthodb.sel.counts=arth.orthodb.sel[ , .(count = .N), by = .(V1,taxid)] %>% dcast(V1~taxid)

taxcodes=as.character(taxid.codes$cafe_code)
for(i in 2:length(colnames(arth.orthodb.sel.counts))){
  num=colnames(arth.orthodb.sel.counts)[i]
  if(num %in% taxcodes){
    colnames(arth.orthodb.sel.counts)[i]=taxid.codes$abbreviation[which(taxid.codes$cafe_code==num)]
  }
}
na.count=apply(arth.orthodb.sel.counts, 1, function(x) sum(is.na(x)))
arth.orthodb.sel.counts[is.na(arth.orthodb.sel.counts)]=0
ra=arth.orthodb.sel.counts %>% select(-V1) %>% rowwise() %>% mutate(range=max(AcrEch:MyzCer)-min(AcrEch:MyzCer)) %>% data.table() %>% .[['range']]
fwrite(arth.orthodb.sel.counts,'/data2/shane/Documents/SLC_id/CAFE/CAFE_tables/CAFE_orthoDB_RAW.tsv')


#rs=arth.orthodb.sel.counts %>% select(-V1) %>% rowSums()
arth.orthodb.sel.counts$missing=na.count
arth.orthodb.sel.counts$range=ra
orthodb.final=arth.orthodb.sel.counts[missing<10]
orthodb.final$missing=NULL
orthodb.final$range=NULL
colnames(orthodb.final)[1]='Family ID'
orthodb.final$Desc='(null)'
orthodb.final$`Family ID`=as.character(orthodb.final$`Family ID`)
orthodb.final=select(orthodb.final,Desc,everything())

#ra=orthodb.final %>% select(-`Family ID`,-Desc) %>% rowwise() %>% mutate(range=max(AcrEch:MyzCer)-min(AcrEch:MyzCer)) %>% data.table() %>% .[['range']]



######################### CLEAN UP TRE FILES
iter=list.files('/data2/shane/Documents/SLC_id/CAFE/trees')[grepl('raxml_tree',list.files('/data2/shane/Documents/SLC_id/CAFE/trees'))] %>% 
  {.[!str_detect(., "named_")]} %>%
  str_remove('raxml_tree_') %>% str_remove('.tre')


###### RUN LOOP
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
  
  ## read in named raxmL tree
  tr=read.tree(paste0('/data2/shane/Documents/SLC_id/CAFE/trees/raxml_tree_named_',i,'.tre'))
  
  
  nodes <- c(); maxes=c()
  maxes=c()
  mins=c()
  if(("CaeEle" %in% tr$tip.label) & ("DroMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("CaeEle","DroMel")));maxes=c(maxes,1000);mins=c(800)}
  if(("AcyPis" %in% tr$tip.label) & ("DroMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("AcyPis","DroMel")));maxes=c(maxes,401);mins=c(mins,345)}
  if(("ApiMel" %in% tr$tip.label) & ("DroMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("ApiMel","DroMel")));maxes=c(maxes,372);mins=c(mins,317)} ## has fossil
  if(("AedAeg" %in% tr$tip.label) & ("DroMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("AedAeg","DroMel")));maxes=c(maxes,206);mins=c(mins,107)} ## has fossil
  if(("NilLug" %in% tr$tip.label) & ("AcyPis" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("AcyPis","NilLug")));maxes=c(maxes,346);mins=c(mins,232)}
  if(("TetUrt" %in% tr$tip.label) & ("DroMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("TetUrt","DroMel")));maxes=c(maxes,579);mins=c(mins,539)}
  if(("PluXyl" %in% tr$tip.label) & ("BomMor" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("PluXyl","BomMor")));maxes=c(maxes,178);mins=c(mins,116)} ## has fossil

 
  ## create ultrametric tree
  ### Credit to Alex SL for format https://phylobotanist.blogspot.com/2019/
  mycalibration <- makeChronosCalib(tr, node=c(nodes), age.min=mins,age.max=maxes)
  mytimetree <- chronos(tr, lambda = 1, model = "correlated", calibration = mycalibration)
  num=mytimetree$node.label %>% as.numeric()
  mytimetree$node.label=NULL
  #mytimetree$edge.length=round(mytimetree$edge.length)
  write.tree(mytimetree, file=paste0("./CAFE/trees/",i,'_tree_ultrametric.tre'))

  
  #### Plot tree
  plot.tree=read.tree(paste0("./CAFE/trees/",i,'_tree_ultrametric.tre'))
  cols=c()
  for(j in num){
    if(is.na(j)){cols=c(cols,'black')
    }else if(j>90){cols=c(cols,'green4')
    }else if(j>70){cols=c(cols,'blue4')
    }else{cols=c(cols,'red')}
  }
  
  ma=max(mytimetree$edge.length)
  xma=ma+100
  ma.r=seq(0,round(ma,-2),by=100)
  
  diff=ma-round(ma,-2)
  
  pdf(paste0("./CAFE/trees/",i,'_tree_ultrametric.pdf'),width=14,height=7)
    gp=ggtree(plot.tree)#, mrsd = "2010-01-01")
    gp=gp+geom_tiplab(size=6)
    gp=gp+geom_nodepoint(size=5,col=cols)
    gp=gp+xlab('Millions of years ago (Mya)')
    gp=gp+theme_tree2()
    gp=gp+theme(axis.text.x=element_text(size=20,face='bold',color = 'black'))
    gp=gp+scale_x_continuous(breaks=diff+ma.r,labels=as.character(rev(ma.r)),limits=c(0,xma))
    print(gp)
  dev.off()
  
  #l.tree.ch=chronopl(read.tree(paste0('/data2/shane/Documents/SLC_id/CAFE/trees/raxml_tree_named_',i,'.tre')), lambda=0.1)
  #l.tree.ch$edge.length=l.tree.ch$edge.length*1000
  #l.tree.ch$node.label=NULL
  #write.tree(l.tree.ch, file=paste0("./CAFE/trees/",i,'_tree_ultrametric.tre'))
  
  ## create lambda file
  writeLines(lambda.convert(readLines(paste0("./CAFE/trees/",i,'_tree_ultrametric.tre'))),paste0("./CAFE/trees/",i,'_tree_lambda.txt')) 
  
  
  sp=str_extract_all(readLines(paste0("./CAFE/trees/",i,'_tree_ultrametric.tre')),pattern = "[A-z]+",simplify = T)
  #sp_clean=sp[sp='_']
  
  sp.counts=full.counts %>% select(c('Desc','Family ID',sp))
  fwrite(sp.counts,paste0('./CAFE/CAFE_tables/',i,'_SLC_CAFE_table.tsv'),sep='\t')
  
  orthodb.counts=orthodb.final %>% select(c('Desc','Family ID',sp))
  fwrite(orthodb.counts,paste0('./CAFE/CAFE_tables/',i,'_OrthoDB_CAFE_table.tsv'),sep='\t')
}
  

