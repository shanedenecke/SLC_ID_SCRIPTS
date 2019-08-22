library(data.table)
library(dplyr)
library(ggplot2)
library(seqinr)


## Set Working directory
setwd('/data2/shane/Documents/SLC_id')


#### Create dictionary 

## table with some info on each species
co.var=fread('./general_reference/Co_variables/Olympia_table_august_2019_shaemod.csv',header=T) %>%
  select(Species_name,abbreviation,Taxanomic_Classification,Phagy,Vory,Diet_category)

######################## functions
dash.remove=function(x){
  split=unlist(strsplit(x,'_'))
  toge=paste(split[1],split[2],sep='_')
  final=toge
  return(final)
}
#########################


l=list()
for(i in list.files('./final_SLC_dicts',full.names = T)){
  dict=fread(i)
  
  abbrev=gsub('./final_SLC_dicts/','',i,fixed=T)
  abbrev=gsub('Final_SLC_table.csv','',abbrev,fixed=T)
  
  fam=sapply(dict$name,dash.remove)
  
  dict$abbreviation=abbrev
  dict$family=fam
  
  l[[i]]=dict
}

all=rbindlist(l) %>% merge(co.var,by='abbreviation')
fasta.name=paste(all$Species_name,all$abbreviation,all$code,all$name,sep='_')
rename.dict=data.table(name=fasta.name,code=all$code)
#rename.dict$name=paste(rename.dict$name,rena)
all$Species_name=gsub(' ','_',all$Species_name)

fwrite(all,'./shiny_prep/Reference_csv_dictionary.csv')
fwrite(rename.dict,'./shiny_prep/Rename_SLC_dict.csv')
writeLines(rename.dict$code,'shiny_prep/slc_codes.txt')
raw.fa=read.fasta('./shiny_prep/Renamed_SLC.faa',seqtype='AA',as.string=T,set.attributes = F,strip.desc=T)


##################### FASTA

system('
cd /data2/shane/Documents/SLC_id
cat ./proteomes/* ./general_reference/model_proteomes/*.faa > ./shiny_prep/all_proteomes.faa
/data2/shane/Applications/custom/unigene_fa_sub.sh ./shiny_prep/all_proteomes.faa ./shiny_prep/slc_codes.txt > ./shiny_prep/SLC_all_raw.faa
/data2/shane/Applications/custom/fasta_rename.py ./shiny_prep/SLC_all_raw.faa ./shiny_prep/Rename_SLC_dict.csv > ./shiny_prep/Renamed_SLC.faa
')






##### Everything based on species
sp='Tribolium castaneum'
input=gsub(' ','_',sp)

species.sub=all[Species_name==input]

## create SLC dictionary
dict=select(species.sub,Species_name,code,name,family)

## SLC counts data
count.table=dict %>% group_by(family) %>% summarize(count=length(family)) %>% data.table()
count.table$species=input

## Histogram

total.counts=all %>% group_by(Species_name) %>% summarize(fam_size=length(Species_name)) %>% data.table()
sub.size=total.counts[Species_name==input]$fam_size


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


## fasta
ind=which(grepl(input,names(raw.fa)))
#names(ind)=NULL
fasta.sub=raw.fa[ind]




##### BASED ON FAMILY

fam='SLC_9'


fam.sub=all[family==fam]


##dictionary
dict=select(fam.sub,Species_name,code,name,family)

## SLC counts data
count.table=dict %>% group_by(Species_name) %>% summarize(count=length(Species_name)) %>% data.table()
count.table$family=fam

fa.search=paste(fam,'_',sep='')
##fasta
ind=which(grepl(fa.search,names(raw.fa)))
#names(ind)=NULL
fasta.sub=raw.fa[ind]


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


