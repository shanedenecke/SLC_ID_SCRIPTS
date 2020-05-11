#!/usr/bin/Rscript

shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))
shhh(library(seqinr))
shhh(library(argparser))


### set working directory to source path
p=arg_parser('Flybase Xref')
p <- add_argument(p, "--out", help="path to outdir")
argv=parse_args(p)
outdir=argv$out

getScriptPath <- function(){
  cmd.args <- commandArgs()
  m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
  script.dir <- dirname(regmatches(cmd.args, m))
  if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
  if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
  return(script.dir)
}
scriptPath=getScriptPath()
sourcePath=dirname(scriptPath)

setwd(scriptPath)
setwd('..')


key=fread('./SLC_id_reference/DROMEL_GENE_KEY.csv') %>% select(Dm_FBgn)

dm_unigene=read.fasta('./SLC_id_reference/DroMel_unigene.faa',
                      set.attributes = F,as.string = T,forceDNAtolower = F) 
dm_unigene2=data.table(full=names(dm_unigene)) %>% mutate(FBgn=full) %>% separate(FBgn,into=c('FBgn','junk'),sep='_') %>% 
  select(-junk) %>% data.table()



fb.slc.gene=fread('./SLC_id_reference/DROMEL_SLC_TABLE_FLYBASE.csv') %>% 
  select(FBgn,Family) %>%
  merge(dm_unigene2,by='FBgn') %>% select(-FBgn)
colnames(fb.slc.gene)=c('Family','code')
  

hs.slcs=fread(paste0(outdir,'/Hs_to_DroMel_Search/final_output/SLC_final_output.csv'))


a=merge(hs.slcs,fb.slc.gene,by='code',all=T) %>% data.table() %>% arrange(name) %>% data.table()
for(i in 1:nrow(a)){
  if(is.na(a[i]$name)){
    a[i,'name']=a[i]$Family
  }
}


num_strip=function(x){
  y=unlist(strsplit(x,split='_'))[1:2] %>% paste(collapse = '_') %>% as.character()
  names(y)=NULL
  return(y)
}

a$name=sapply(a$name,num_strip)
a=select(a,-Family)

l=list()
for(i in unique(a$name)){
  sub=a[name==i]
  series=1:dim(sub)[1]
  for(j in series){
    sub[j,'name']=paste(sub[j]$name,j,sep='_')
  }
  l[[i]]=sub
}

final=rbindlist(l)
cat(format_csv(final))


