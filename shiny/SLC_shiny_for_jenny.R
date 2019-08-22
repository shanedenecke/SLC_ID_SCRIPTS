library(data.table)
library(dplyr)
library(ggplot2)
library(seqinr)
library(ggsci)


## Set Working directory
setwd('/data2/shane/Documents/SLC_id/shiny_prep/')





#### Import Raw base data
raw.fa=read.fasta('./Renamed_SLC.faa',seqtype='AA',as.string=T,set.attributes = F,strip.desc=T)
full.table=fread('./Reference_csv_dictionary.csv')
species.meta=fread('./Olympia_table_august_2019_shaemod.csv',header=T) 

##### SPECIES ####################

## CHOOSE INPUT
sp.input='Tribolium castaneum'
sp.input.mod=gsub(' ','_',sp.input)

### SUBSET SPECIES
species.sub=full.table[Species_name==sp.input.mod]

## CREATE FULL DICTIONARY/TABLE FOR SHINY
sp.dict=select(species.sub,Species_name,code,name,family)
head(sp.dict)

## COUNT DATA
sp.count.table=sp.dict %>% group_by(family) %>% summarize(count=length(family)) %>% data.table()
sp.count.table$species=sp.input.mod
head(sp.count.table)


## GENERATE HISTOGRAM
total.counts=full.table %>% group_by(Species_name) %>% summarize(fam_size=length(Species_name)) %>% data.table()
sub.size=total.counts[Species_name==sp.input.mod]$fam_size

gp=ggplot(total.counts,aes(x=fam_size))
gp=gp+geom_histogram(colour="black", fill="white",binwidth=30)
gp=gp+geom_vline(xintercept=sub.size,color="red", linetype="solid", size=1)
gp=gp+geom_density(alpha=.2, fill="#FF6666")
gp=gp+labs(x='Family Size',y='Number of Occurances')
gp=gp+theme_bw()
gp=gp+theme(text=element_text(face="bold",family="serif"),panel.grid=element_blank(),
            axis.ticks.x=element_line(),panel.border=element_rect(colour="black",fill=NA),
            strip.text=element_text(size=20),strip.background=element_rect("white"),
            axis.title=element_text(size=25),axis.text.x=element_text(size=20),
            legend.position = 'none',plot.title = element_text(hjust = 0.5))
print(gp)


## GENERATE FASTA
ind=which(grepl(sp.input.mod,names(raw.fa)))
sp.fasta.sub=raw.fa[ind]

a=c()
for(i in names(sp.fasta.sub)){
  n=paste0('>',i)
  a=c(a,n)
  a=c(a,sp.fasta.sub[[i]])
}
writeLines(a,'test.fa')

head(sp.fasta.sub)


write.fasta(sp.fasta.sub,names=names(sp.fasta.sub),file.out=paste0(sp.input.mod,'_Fasta_subset.faa'),nbchar=10000)



##### BASED ON SLC FAMILY


## GIVE INPUT
fam.input='SLC_9'


## SUBSET FAMILY
fam.sub=full.table[family==fam.input]


##CREATE BASE DICTIONARY
fam.dict=select(fam.sub,Species_name,code,name,family)
head(fam.dict)

## DIVIDE BASED ON COUNTS
fam.count.table=fam.dict %>% group_by(Species_name) %>% summarize(count=length(Species_name)) %>% data.table()
fam.count.table$family=fam.input
head(fam.count.table)


### SUBSET FASTA
fam.input.mod=paste(fam.input,'_',sep='')
index=which(grepl(fam.input.mod,names(raw.fa)))
fam.fasta.sub=raw.fa[index]
head(fam.fasta.sub)
write.fasta(fam.fasta.sub,names=names(fam.fasta.sub),file.out=paste0(sp.input.mod,'_Fasta_subset.faa'),nbchar=10000)


## Tree

### Pre processing  
load('tree_list')
tree.groups=gsub("RAxML_bipartitions.","",names(l))


## make searchable form of input variable
tree.search=paste0(fam,'.tre')

if(tree.search %in% tree.groups){ ## check to see if there are any arhtropod SLCs in fam variable
  
  tree=l[grepl(tree.search,names(l))][[1]]  ## subset tree
  
  ## scale size to number of observations
  no.obs=length(tree$tip.label)
  text.size=min(250/no.obs,3)
  
  
  ### plot graph
  gp=ggtree(tree,branch.length=.1)
  gp=gp+geom_tiplab(size=text.size)
  gp=gp+xlim(0, 8)
  gp=gp+theme_tree()
  print(gp)
  
  ggsave(paste0('./test_trees/',fam,'.pdf'),gp,device='pdf',width=20,height=10)  ## allow users to download pdf or png 
  write.tree(tree) ## write raw tree output 
}else{ ## if no SLCs detected print error message
  print(paste0('Error: No Tree can be generated for the ',fam,' family. Probably no arthropod SLCs were detected'))
}


#### Boxplots
## can be either "Phagy", "Taxanomic_Classification", "Vory", "Diet_category"
co.var='Taxanomic_Classification'


merged=merge(fam.count.table,species.meta,by='Species_name')


red=select(merged,count,co.var)
red[[fam.input]]=as.numeric(red[['count']])

plot=red[red[[co.var]] %in% names(which(table(red[[co.var]])>5))] ### remove elements that aren't present at least 5 times

ylab=paste0('Number of ',gsub('_',' ',fam.input),' Members','\n')
tit=paste0(gsub('_',' ',fam.input),' Family Size vrs. ',gsub('_',' ',co.var))
xlab=paste0('\n',gsub('_',' ',co.var))

gp=ggplot(plot,aes_string(co.var,y=fam.input,fill=co.var))
gp=gp+geom_boxplot()
gp=gp+labs(x=xlab,y=ylab)
gp=gp+scale_fill_rickandmorty()
gp=gp+ggtitle(tit)
gp=gp+theme_bw()
gp=gp+theme(text=element_text(face="bold",family="serif"),panel.grid=element_blank(),
            axis.ticks.x=element_line(),panel.border=element_rect(colour="black",fill=NA),
            strip.text=element_text(size=20),strip.background=element_rect("white"),
            axis.title=element_text(size=17),axis.text.x=element_text(angle=45,hjust=1),
            legend.position = 'none',plot.title = element_text(hjust = 0.5))

print(gp)



