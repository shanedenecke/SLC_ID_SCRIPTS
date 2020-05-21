#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(dplyr))
shhh(library(data.table))
shhh(library(stringr))
shhh(library(ape))
shhh(library(ggtree))
shhh(library(tidyr))
shhh(library(ggplot2))
shhh(library(treeio))
setwd('/mnt/disk/shane/Transporter_ID/SLC_id_pipeline')
dir.create('./CAFE/CAFE_figures')

iter=list.files('./CAFE/CAFE_tables/') %>% gsub('_SLC_CAFE_table.tsv','',.)
full.metadata=fread('./Final_outputs/Full_Metadata_summary.csv')
#group='Hemiptera_species'
#family='SLC_33'
#node.annot=''
#label.annot=''


xma.calc=function(tree){
  table=as_tibble(tree) %>% data.table()
  sum1=table[node %in% tidytree::ancestor(tree,1)]$branch.length %>% sum(na.rm=T)
  sum2=table[node==1]$branch.length
  return(sum1+sum2)
}

### tree.fig function                 
tree.fig=function(group,family,node.annot='',label.annot=''){
  
  ##import ultrametric tree
  ultra.tree=read.tree(paste0('./CAFE/clean_raxml_trees/',group,'_tree_ultrametric.nwk'))
  tbl=as_tibble(ultra.tree) %>% data.table()
  
  #### set colors
  node.tree=read.tree(paste0('./CAFE/clean_raxml_trees/RAxML_bipartitions.',group,'_subset.nwk'))
  node.scores=as.numeric(node.tree$node.label)
  cols=c()
  for(j in node.scores){
    if(is.na(j)){cols=c(cols,'grey50')
    }else if(j>80){cols=c(cols,'green4')
    }else if(j>50){cols=c(cols,'gold4')
    }else{cols=c(cols,'red')}
  }
  
  
  ### Import node labels
  lab.text=readLines(paste0('./CAFE/outputs/',group,'/','Base_asr.tre'))[3] %>% 
    gsub("^.+ = ",'',.) %>% gsub('>_[0-9|\\:|\\.]+','>',.) %>%
    gsub(';','',.)
  #lab.text=gsub('# The labeled CAFE tree:\t','',readLines(paste0('./CAFE/outputs/',group,'_SLC_summary.txt_fams.txt'))[1])
  lab.tree=read.tree(text=paste0('(',lab.text,')',';'))
  
  ## import counts
  #count.table=fread(paste0('./CAFE/outputs/',group,'_SLC_summary.txt_anc.txt'))[`Family ID`==family]
  count.table=fread(paste0('./CAFE/outputs/',group,'/','Base_count.tab'))[`FamilyID`==family]
  if(nrow(count.table)==0){
    break
  }
  count.reduce=count.table %>% select(-matches('[A-z]'))
  count.term=count.table %>% select(matches('[A-z]')) %>% select(-`FamilyID`)
  colnames(count.term)=gsub('<[0-9]+>','',colnames(count.term))
  
  ## add node labels to ultrametric
  sorted=c()
  for(i in lab.tree$node.label){
    temp=count.reduce[[i]]
    sorted=c(sorted,temp)
  }
  final=c('',as.character(sorted))
  final2=final[-1]
  ultra.tree$node.label=final2
  
  ## add tip labels to ultrametric with numbers
  for(i in ultra.tree$tip.label){
    temp=count.term[[i]]
    ultra.tree$tip.label[which(ultra.tree$tip.label==i)]=paste0(i,' (',temp,')')
  }
  
  ##Set scaling factors
  #ma=max(sapply(1:ultra.tree$Nnode,function(x) tbl[node %in% ancestor(ultra.tree,x)]$branch.length %>% sum(na.rm = T)))
  #xma=ma*1.4
  #ma.r=seq(0,round(ma,-2),by=50)
  #diff=ma-round(ma,-2)
  
  #ma=max(ultra.tree$edge.length)
  #if(grepl('ArachIn',group)){ma=(254.90342+254.90367+45.39890)}
  ma=xma.calc(ultra.tree)
  xma=ma*1.3
  ma.r=seq(0,round(ma,-2),by=100)
  diff=ma-round(ma,-2)
  
  
  
  #### Make plot
  
  gp=ggtree(ultra.tree,size=2)
  gp=gp+geom_tiplab(size=10,fontface='bold')#,aes(label=paste0('bold(', label, ')')), parse=TRUE)
  gp=gp+geom_nodepoint(size=21,color=cols)
  gp=gp+geom_nodelab(hjust=.6,size=12,fontface='bold',color='white')
  
  if(is.list(node.annot)){
    names(node.annot)=label.annot
    plot.annot=vector("list",length=length(label.annot))
    names(plot.annot)=label.annot
    for(i in names(plot.annot)){
      plot.annot[[i]][1]=tbl[grepl(node.annot[[i]][1],label)]$node %>% as.numeric()
      plot.annot[[i]][2]=tbl[grepl(node.annot[[i]][2],label)]$node %>% as.numeric()
      gp=gp+geom_strip(taxa1=plot.annot[[i]][1],taxa2=plot.annot[[i]][2],offset.text=4,fontsize=10,
                       barsize=5,color='black',label=i,offset=ma/6.5)
    }
  }
  #gp=gp+geom_hilight(node=list(node1,55),fill='darkgreen',alpha=.3)
  gp=gp+theme_tree2()
  gp=gp+theme(axis.text.x=element_text(size=20,face='bold',color = 'black'),axis.line.x=element_line(size=3),
              axis.title.x=element_text(size=20))
  gp=gp+scale_x_continuous(breaks=diff+ma.r,labels=as.character(rev(ma.r)),limits=c(0,xma))
  #print(gp)
}

a=tree.fig(group='Hemiptera_species',family='SLC_33')
ggsave(a,filename='test.pdf',device='pdf',width=15,height=10)
#fams=counts[['Family ID']][!grepl('Unsorted',counts[['Family ID']])]
fams=colnames(full.metadata)[grepl('SLC_',colnames(full.metadata))]
fams=fams[fams!='SLC_14' & fams!='SLC_Unsorted' & fams!='SLC_total']


iter=iter[c(1,3,5)]
iter=iter[2]
for (i in iter){ 
  counts=fread(paste0('./CAFE/CAFE_tables/',i,'_SLC_CAFE_table.tsv'))
  for(j in fams){ 
    res=try(tree.fig(group = i,family = j),silent=T)
    if(class(res)!='try-error'){
      temp=tree.fig(group = i,family = j)
      ggsave(plot=temp,filename=paste0('./CAFE/CAFE_figures/',i,'_',j,'.pdf'),device='pdf',width=20,height=15)
    }
  }
}
