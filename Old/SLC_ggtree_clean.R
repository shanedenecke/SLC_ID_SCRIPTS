shhh <- suppressPackageStartupMessages
shhh(library(dplyr))
shhh(library(data.table))
shhh(library(ape))
shhh(library(ggtree))
shhh(library(treeio))
shhh(library(ggplot2))


#setwd('/data2/shane/Transporter_ID/SLC_id')
dir.create('./SLC_phylogeny/clean_ggtrees')


############### FUNCTION
di2multi4node <- function (phy, tol = 0.5) {
  # Adapted di2multi function from the ape package to plot polytomies
  # based on numeric node support values
  # (di2multi does this based on edge lengths)
  # Needs adjustment for unrooted trees as currently skips the first edge
  if (is.null(phy$edge.length)) 
    stop("the tree has no branch length")
  if (is.na(as.numeric(phy$node.label[2])))
    stop("node labels can't be converted to numeric values")
  if (is.null(phy$node.label))
    stop("the tree has no node labels")
  ind <- which(phy$edge[, 2] > length(phy$tip.label))[as.numeric(phy$node.label[2:length(phy$node.label)]) < tol]
  n <- length(ind)
  if (!n) 
    return(phy)
  foo <- function(ancestor, des2del) {
    wh <- which(phy$edge[, 1] == des2del)
    for (k in wh) {
      if (phy$edge[k, 2] %in% node2del) 
        foo(ancestor, phy$edge[k, 2])
      else phy$edge[k, 1] <<- ancestor
    }
  }
  node2del <- phy$edge[ind, 2]
  anc <- phy$edge[ind, 1]
  for (i in 1:n) {
    if (anc[i] %in% node2del) 
      next
    foo(anc[i], node2del[i])
  }
  phy$edge <- phy$edge[-ind, ]
  phy$edge.length <- phy$edge.length[-ind]
  phy$Nnode <- phy$Nnode - n
  sel <- phy$edge > min(node2del)
  for (i in which(sel)) phy$edge[i] <- phy$edge[i] - sum(node2del < 
                                                           phy$edge[i])
  if (!is.null(phy$node.label)) 
    phy$node.label <- phy$node.label[-(node2del - length(phy$tip.label))]
  phy
}

max.branch=function(tree){
  tb=as_tibble(tree)
  unitip=unique(tb$label)[grepl('[A-z]',unique(tb$label))]
  tot=c()
  for(i in unitip){
    anc=sum(ancestor(tb,i)$branch.length,na.rm = T)
    last=tbl %>% filter(label==i) %>% .[['branch.length']]
    tot=c(tot,anc+last) #### calculate all total branch lengths 
  }
  final=max(tot)
  return(final)
}

dromel.short=function(x){
  short=gsub('CG[0-9]+_','',x) %>% gsub('FBpp[0-9]+_','',.)
}


unispec=c('TriCas','AcyPis','ApiMel','BomMor','DroMel','HomSap','SpoFru','HelArm','NezVir')
names(unispec)=c('blue4','gold4','red4','purple','black','grey50','blue4','gold4','red4')

dims=list()
for(i in list.files('./SLC_phylogeny/raxml_trees/',full.names = T,include.dirs = F)){
  #import tree
  bname=gsub('RAxML_bipartitions.(.+).tre$','\\1',basename(i))
  base.tree=read.tree(i)
  base.tree$tip.label=dromel.short(base.tree$tip.label)
  
  ### collapse nodes with poor bootstrap values
  base.tree$node.label=as.numeric(base.tree$node.label)
  base.tree$node.label[is.na(base.tree$node.label)]=99.99
  tree.collapse=di2multi4node(base.tree,30)
  tree.collapse$node.label=as.character(tree.collapse$node.label)
  tree.collapse$node.label[tree.collapse$node.label=='99.99']=''

  ##set colors
  cols=c()
  for(j in tree.collapse$tip.label){cols=c(cols,names(unispec)[sapply(unispec,grepl,j)])}  
  
  ### root tree
  tbl=as_tibble(tree.collapse)
  out=tbl$label[1]
  tree.root=root(tree.collapse,outgroup=out,edgelabel=T)
  
  ### set maximum value
  xma=max.branch(tree.root)+4
  yma=length(tree.root$tip.label)*.3
  
  gp=ggtree(tree.root,size=2)
  gp=gp+geom_tiplab(color=cols,size=4,fontface='bold',align=T)
  gp=gp+geom_nodepoint(size=4,col='black')
  gp=gp+geom_nodelab(hjust=1,vjust=.3,size=1.5,fontface='bold',col='white')
  gp=gp+lims(x=c(0,xma))
  gp=gp+theme(title = element_text(size=12))
  #print(gp)
  ggsave(plot=gp,filename = paste0('./SLC_phylogeny/clean_ggtrees/',bname,'_preliminary_phylogeny.pdf'),device='pdf',height=yma,width=xma*1.5,limitsize = F)
  dims[[i]]=data.table(xma,yma,bname)
}
#rbindlist(dims)