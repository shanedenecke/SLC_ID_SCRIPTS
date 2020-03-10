shhh <- suppressPackageStartupMessages
shhh(library(dplyr))
shhh(library(data.table))
shhh(library(ape))
shhh(library(ggtree))
shhh(library(phytools))


#setwd('/data2/shane/Transporter_ID/SLC_id')
dir.create('./SLC_phylogeny/clean_ggtrees')
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

unispec=c('TriCas','AcyPis','ApiMel','BomMor','DroMel','HomSap','SpoFru','HelArm','NezVir')
names(unispec)=c('blue4','gold4','red4','purple','black','grey50','blue4','gold4','red4')

for(i in list.files('./SLC_phylogeny/raxml_trees/',full.names = T,include.dirs = F)){
  bname=gsub('RAxML_bipartitions.(.+).tre$','\\1',basename(i))
  base.tree=read.tree(i)
  
  ### root tree
  base.tree=root(base.tree,outgroup='HomSap_SLC9A4',edgelabel=T)
  
  
  
  ##set colors
  cols=c()
  for(j in base.tree$tip.label){cols=c(cols,names(unispec)[sapply(unispec,grepl,j)])}  
  
  
  ### collapse nodes with poor bootstrap values
  base.tree$node.label=as.numeric(base.tree$node.label)
  base.tree$node.label[is.na(base.tree$node.label)]=99.99
  tree2=di2multi4node(base.tree,30)
  tree2$node.label=as.character(tree2$node.label)
  tree2$node.label[tree2$node.label=='99.99']=''
  
  gp=ggtree(base.tree,size=2)
  gp=gp+geom_tiplab(color=cols,size=4,fontface='bold')
  gp=gp+geom_nodepoint(size=4,col='black')
  gp=gp+geom_nodelab(hjust=1,vjust=.3,size=1.5,fontface='bold',col='white')
  gp=gp+lims(x=c(0,10),y=c(0,length(base.tree$tip.label)))
  gp=gp+theme(title = element_text(size=12))
  print(gp)
  ggsave(plot=gp,filename = paste0('./SLC_phylogeny/clean_ggtrees',bname,'.pdf'),device='pdf')
  
}