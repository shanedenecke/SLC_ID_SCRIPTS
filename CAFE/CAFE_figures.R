shhh <- suppressPackageStartupMessages
shhh(library(dplyr))
shhh(library(data.table))
shhh(library(stringr))
shhh(library(ape))
shhh(library(ggtree))
shhh(library(tidyr))

setwd('/data2/shane/Documents/SLC_id/Figures')
dir.create('CAFE_figures')

arth.ids=gsub("# Output format for: ' Average Expansion', 'Expansions', 'No Change', 'Contractions', and 'Branch-specific P-values' = (node ID, node ID): ",
              "",readLines('/data2/shane/Documents/SLC_id/CAFE/outputs/Arthropod_SLC_cafe_output.cafe')[5],fixed=T)
arth.ids2=unlist(strsplit(arth.ids,split=' '))
sp=gsub('# IDs of nodes:','',readLines('/data2/shane/Documents/SLC_id/CAFE/outputs/Arthropod_SLC_cafe_output.cafe')[4])
sp2=unlist(str_extract_all(sp,"[A-z]{6,}<[0-9]+>")) %>% str_remove_all('<|>')
arth.cafe=fread('/data2/shane/Documents/SLC_id/CAFE/outputs/Arthropod_SLC_cafe_output.cafe',skip=11,sep='\t') %>% select(V1,V3,V4) %>% 
  filter(V3<.05) %>% separate(col=V4,into=arth.ids2,sep='\\),\\(') %>%
    data.table()
#colnames(arth.cafe)[3]=arth.ids2
arth.cafe %>% fwrite('./CAFE_figures/Arth_CAFE.csv')


### arachnid
arac.ids=gsub("# Output format for: ' Average Expansion', 'Expansions', 'No Change', 'Contractions', and 'Branch-specific P-values' = (node ID, node ID): ",
              "",readLines('/data2/shane/Documents/SLC_id/CAFE/outputs/Arachnid_SLC_cafe_output.cafe')[5],fixed=T)
arac.ids2=unlist(strsplit(arac.ids,split=' '))
sp=gsub('# IDs of nodes:','',readLines('/data2/shane/Documents/SLC_id/CAFE/outputs/Arachnid_SLC_cafe_output.cafe')[4])
sp2=unlist(str_extract_all(sp,"[A-z]{6,}<[0-9]+>")) %>% str_remove_all('<|>')
arac.cafe=fread('/data2/shane/Documents/SLC_id/CAFE/outputs/Arachnid_SLC_cafe_output.cafe',skip=11,sep='\t') %>% select(V1,V3,V4) %>% 
  filter(V3<.05) %>% separate(col=V4,into=arac.ids2,sep='\\),\\(') %>%
  data.table()
#colnames(arac.cafe)[3]=arac.ids2
arac.cafe %>% fwrite('./CAFE_figures/Arac_CAFE.csv')



### hemi
hemi.ids=gsub("# Output format for: ' Average Expansion', 'Expansions', 'No Change', 'Contractions', and 'Branch-specific P-values' = (node ID, node ID): ",
              "",readLines('/data2/shane/Documents/SLC_id/CAFE/outputs/Hemipteran_SLC_cafe_output.cafe')[5],fixed=T)
hemi.ids2=unlist(strsplit(hemi.ids,split=' '))
sp=gsub('# IDs of nodes:','',readLines('/data2/shane/Documents/SLC_id/CAFE/outputs/Hemipteran_SLC_cafe_output.cafe')[4])
sp2=unlist(str_extract_all(sp,"[A-z]{6,}<[0-9]+>")) %>% str_remove_all('<|>')
hemi.cafe=fread('/data2/shane/Documents/SLC_id/CAFE/outputs/Hemipteran_SLC_cafe_output.cafe',skip=11,sep='\t') %>% select(V1,V3,V4) %>% 
  filter(V3<.05) %>% separate(col=V4,into=hemi.ids2,sep='\\),\\(') %>%
  data.table()
#colnames(hemi.cafe)[3]=hemi.ids2
hemi.cafe %>% fwrite('./CAFE_figures/hemi_CAFE.csv')



#lepi.cafe=fread('/data2/shane/Documents/SLC_id/CAFE/outputs/Lepidopteran_SLC_cafe_output.cafe',skip=11,sep='\t') %>% select(V1,V3,V4) %>% 
#  filter(V3<.011) %>% data.table()

group='Arthropod'
family='SLC36'

#ree.fig=function(raxtree,ultratree,fams.summary,counts,family){
tree.fig=function(group,family){
  
  ##import ultrametric tree
  base.tree=read.tree(paste0('/data2/shane/Documents/SLC_id/CAFE/trees/',group,'_tree_ultrametric.tre'))
  
  ##import node labels
  lab.text=gsub('# The labeled CAFE tree:\t','',readLines(paste0('/data2/shane/Documents/SLC_id/CAFE/outputs/',group,'_SLC_summary.txt_fams.txt'))[1])
  lab.tree=read.tree(text=paste0('(',lab.text,')',';'))
  
  ## import counts
  count.table=fread(paste0('/data2/shane/Documents/SLC_id/CAFE/outputs/',group,'_SLC_summary.txt_anc.txt'))[`Family ID`==family]
  count.reduce=count.table %>% select(-matches('[A-z]'))
  count.term=count.table %>% select(matches('[A-z]')) %>% select(-`Family ID`)
  colnames(count.term)=gsub('<[0-9]+>','',colnames(count.term))
  
  ## add node labels to ultrametric
  sorted=c()
  for(i in lab.tree$node.label){
    temp=count.reduce[[i]]
    sorted=c(sorted,temp)
  }
  final=c('',as.character(sorted))
  base.tree$node.label=final
  
  ## add tip labels to ultrametric with numbers
  for(i in base.tree$tip.label){
    temp=count.term[[i]]
    base.tree$tip.label[which(base.tree$tip.label==i)]=paste0(i,' (',temp,')')
  }
  
  ## create color scheme from ultrametric
  node.scores=read.tree(paste0('/data2/shane/Documents/SLC_id/CAFE/trees/raxml_tree_named_',group,'.tre'))$node.label %>% as.numeric()
  cols=c()
  for(j in node.scores){
    if(is.na(j)){cols=c(cols,'grey50')
    }else if(j>90){cols=c(cols,'green3')
    }else if(j>70){cols=c(cols,'#e4d948ff')
    }else{cols=c(cols,'red')}
  }
  
  ##Set scaling factors
  ma=max(base.tree$edge.length)
  xma=ma+100
  ma.r=seq(0,round(ma,-2),by=100)
  diff=ma-round(ma,-2)
  #gsub("(^[0-9]+)","\\(\\1\\)",base.tree$node.label)
  
  gp=ggtree(base.tree,size=2)
  gp=gp+geom_tiplab(size=8,fontface='bold')#,aes(label=paste0('bold(', label, ')')), parse=TRUE)
  gp=gp+geom_nodepoint(size=16,col=cols)
  gp=gp+geom_nodelab(hjust=.75,size=8,fontface='bold')
  #gp=gp+geom_nodelab(hjust=1.9,vjust=-.4,size=8,fontface='bold')
  gp=gp+theme_tree2()
  gp=gp+theme(axis.text.x=element_text(size=20,face='bold',color = 'black'),axis.line.x=element_line(size=3),
              axis.title.x=element_text(size=20))
  gp=gp+scale_x_continuous(breaks=diff+ma.r,labels=as.character(rev(ma.r)),limits=c(0,xma))
  print(gp)
}

#dir.create('CAFE_figures')
arth.21=tree.fig(group='Arthropod','SLC21')
ggsave(filename='./CAFE_figures/Arth_21.pdf',plot=arth.21,device='pdf',width=10,height=8)

hemi.33=tree.fig(group='Hemipteran','SLC33')
ggsave(filename='./CAFE_figures/Hemi_33.pdf',plot=hemi.33,device='pdf',width=15,height=10)

arac.2=tree.fig(group='Arachnid','SLC2')
ggsave(filename='./CAFE_figures/Arac_2.pdf',plot=arac.2,device='pdf',width=15,height=10)
#lepi.22=tree.fig(group='Lepidopteran','SLC22')
#ggsave(filename='./CAFE_figures/Lepi.22.pdf',plot=lepi.22,device='pdf',width=15,height=10)

hemi.36=tree.fig(group='Arthropod','SLC36')

arac.60=tree.fig(group='Arachnid','SLC60')
ggsave(filename='./CAFE_figures/Arac.60.pdf',plot=arac.60,device='pdf',width=15,height=10)

arac.35=tree.fig(group='Arachnid','SLC35')
ggsave(filename='./CAFE_figures/Arac.35.pdf',plot=arac.35,device='pdf',width=15,height=10)

arth.2=tree.fig(group='Arthropod','SLC2')
ggsave(filename='./CAFE_figures/Arth.2.pdf',plot=arth.2,device='pdf',width=15,height=10)