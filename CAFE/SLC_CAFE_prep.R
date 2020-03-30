shhh <- suppressPackageStartupMessages
shhh(library(dplyr))
shhh(library(data.table))
shhh(library(stringr))
shhh(library(ape))
shhh(library(ggtree))
shhh(library(tidyr))
shhh(library(ggplot2))


#args = commandArgs(trailingOnly=TRUE)
#H=as.character(args[1])
#setwd('/data2/shane/Transporter_ID/SLC_id')
used.species=readLines('./Final_raw_outputs/Good_quality_species.txt')
dir.create('./CAFE',showWarnings = F)
dir.create('./CAFE/CAFE_tables',showWarnings = F)


### Import data
slc.counts=fread('./Final_raw_outputs/TableS7_count_summary.csv')
metadata=fread('./GENERAL_REFERENCE/Co_variables/Arthropod_species_metadata.tsv',header=T) %>% 
  select(Species_name,abbreviation,taxid_code) %>% filter(taxid_code %in% used.species) %>% data.table()
#taxid_key=fread('./GENERAL_REFERENCE/species_metadata/taxid_key.tsv',col.names = c('taxid_code','abbreviation','species_name'))
#system('cp ./GENERAL_REFERENCE/ultrametric_tree_backup/* ./CAFE/clean_raxml_trees')

### Import functions
lambda.convert=function(x){
  a=gsub("[A-z]+:[0-9]+.[0-9]+",1,x)
  b=gsub("[0-9]{2,}:[0-9]+.[0-9]+",1,a)
  c=gsub(":[0-9]+.[0-9]+",1,b)
  d=gsub('[A-z]+:[0-9]+',1,c)
  return(d)
}


shane.transpose=function(dt,newcol){
  new.rows=colnames(dt)[2:length(colnames(dt))]
  a=transpose(dt)
  colnames(a)=as.character(a[1,])
  a=a[-1,]
  b=mutate(a,newcol=new.rows) %>% select(newcol,everything())
  return(b)
}


### format table for CAFE

slc.counts=slc.counts %>% select(-SLC_total,-SLC_Unsorted) %>% shane.transpose()
colnames(slc.counts)[1]='Family ID'
#slc.counts$`Family ID`=gsub('_','',slc.counts$`Family ID`)
slc.counts$Desc='(null)'
#trans.counts=rbindlist(list(trans.counts,sup.counts),use.names = T,fill=T)
slc.counts=select(slc.counts,Desc,'Family ID',metadata$abbreviation)

fwrite(slc.counts,'./CAFE/SLC_COUNTS_CAFE_FULL.tsv',sep='\t')


### create list of tree base names
iter=list.files('./CAFE/clean_raxml_trees')[grepl('RAxML_bipartitions.',list.files('./CAFE/clean_raxml_trees'))] %>%
  str_remove('RAxML_bipartitions.') %>% str_remove('.tre')

###### RUN LOOP
for (i in iter){
  
  ### rename tree
  rax.tree=readLines(paste0('./CAFE/clean_raxml_trees/RAxML_bipartitions.',i,'.tre'))
  for(j in metadata$taxid_code){
    if(grepl(j,rax.tree)){
      rax.tree=gsub(j,metadata[taxid_code==j]$abbreviation,rax.tree)
    }
  }
  rax.tree.name=gsub('7227_0','DroMel',rax.tree)
  writeLines(rax.tree.name,paste0('./CAFE/clean_raxml_trees/raxml_tree_named_',i,'.tre'))
  
  ## read in named raxmL tree
  tr=read.tree(paste0('./CAFE/clean_raxml_trees/raxml_tree_named_',i,'.tre'))
  
  
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
  write.tree(mytimetree, file=paste0("./CAFE/clean_raxml_trees/",i,'_tree_ultrametric.tre'))
  
  
  #### Plot tree
  plot.tree=read.tree(paste0("./CAFE/clean_raxml_trees/",i,'_tree_ultrametric.tre'))
  cols=c()
  for(j in num){
    if(is.na(j)){cols=c(cols,'black')
    }else if(j>80){cols=c(cols,'green4')
    }else if(j>50){cols=c(cols,'blue4')
    }else{cols=c(cols,'red')}
  }
  
  ma=max(mytimetree$edge.length)
  xma=ma+100
  ma.r=seq(0,round(ma,-2),by=100)
  
  diff=ma-round(ma,-2)
  
  #pdf(paste0("./CAFE/clean_raxml_trees/",i,'_tree_ultrametric.pdf'),width=14,height=7)
  gp=ggtree(plot.tree)#, mrsd = "2010-01-01")
  gp=gp+geom_tiplab(size=6)
  gp=gp+geom_nodepoint(size=5,col=cols)
  gp=gp+xlab('Millions of years ago (Mya)')
  gp=gp+theme_tree2()
  gp=gp+theme(axis.text.x=element_text(size=20,face='bold',color = 'black'))
  gp=gp+scale_x_continuous(breaks=diff+ma.r,labels=as.character(rev(ma.r)),limits=c(0,xma))
  #print(gp)
  #dev.off()
  
  ggsave(paste0("./CAFE/clean_raxml_trees/",i,'_tree_ultrametric.pdf'),plot=gp,width=14,height=10)
  
  #l.tree.ch=chronopl(read.tree(paste0(paste0(H,'CAFE/trees/raxml_tree_named_',i,'.tre')), lambda=0.1)
  #l.tree.ch$edge.length=l.tree.ch$edge.length*1000
  #l.tree.ch$node.label=NULL
  #write.tree(l.tree.ch, file=paste0("./CAFE/clean_raxml_trees/",i,'_tree_ultrametric.tre'))
  
  ## create lambda file
  writeLines(lambda.convert(readLines(paste0("./CAFE/clean_raxml_trees/",i,'_tree_ultrametric.tre'))),paste0("./CAFE/clean_raxml_trees/",i,'_tree_lambda.txt')) 
  
  
  sp=str_extract_all(readLines(paste0("./CAFE/clean_raxml_trees/",i,'_tree_ultrametric.tre')),pattern = "[A-z]+",simplify = T)
  sp=as.character(sp)
  
  sp.counts=slc.counts %>% select(c('Desc','Family ID',all_of(sp)))
  fwrite(sp.counts,paste0('./CAFE/CAFE_tables/',i,'_SLC_CAFE_table.tsv'),sep='\t')
  #orthodb.counts=orthodb.final %>% select(c('Desc','Family ID',sp))
  #fwrite(orthodb.counts,paste0('./CAFE/CAFE_tables/',i,'_OrthoDB_CAFE_table.tsv'),sep='\t')
}





