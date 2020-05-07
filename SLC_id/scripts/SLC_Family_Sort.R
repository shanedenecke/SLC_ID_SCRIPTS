#!/usr/bin/Rscript

## R script to sort SLC hits. The input is a formatted blast output table. 
### Each gene is sorted into family or is discarded based on criteria similar to Hoglund et. al 2011

## 1st argument is an SLC dictionary (usually form the source database) which will be used against the reciprocal blast results

#rm(list=ls())
#setwd('/data2/shane/Transporter_ID/SLC_id/Dm_Database_Generate/Hs_to_DroMel_Search')
##Import libraries
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))

############################ Set WD and arguments
#
args = commandArgs(trailingOnly=TRUE)

##debug
#setwd('/home/shanedenecke/Documents/SLC_id/Human_search/HUMAN_AttCep_unigene.faa')
#args[1]="/data2/shane/Transporter_ID/SLC_id/HomSap_Database/SLC_source_dict.csv"

setwd('./recip_blast')

## read in source table
source.table=fread(as.character(args[1]))


##summarize number of SLCs in each family from the soure table
source.count.summary=source.table %>% separate(col=name,sep='_',remove=T,into=c('slc','fam','subfam')) %>% 
  unite(col='slc_fam',slc,fam,sep="_") %>% group_by(slc_fam) %>% summarize(fam_count=length(slc_fam))
source.count.summary$slc_fam=sapply(source.count.summary$slc_fam, function(x) paste0(x,"_")) ## Add extra _ in each


## functions
dash.remove=function(x){
  split=unlist(strsplit(x,'_'))
  toge=paste(split[1],split[2],sep='_')
  final=paste0(toge,'_')
  return(final)
}

## initialize some empty lists
slc.total=list()
filter.list=list()
used.list=c()

all.blast=list()
for (i in list.files()){
  all.blast[[i]]=fread(i,colClasses = c('character','character',rep('numeric',3)),sep='\t',header=F)
}
all.blast=rbindlist(all.blast) %>% unique.data.frame()
colnames(all.blast)=c('query','subject','pident','evalue','querycov')
all.blast$subject=sapply(all.blast$subject,dash.remove)
all.blast$subject=gsub('_NA','',all.blast$subject)
unigenes=all.blast$query %>% unique()



for (x in 1:length(unigenes)){
  
  ####### Set variables
  
  i=unigenes[x]
  sub.blast=subset(all.blast,query==i)
  slc_fams=sub.blast$subject
  slc_tab=table(slc_fams) %>% sort(decreasing = T)
  #top_hit_fam=names(slc_tab)[1]
  top_hit_fam=slc_fams[1] 
  top_hit_gene=sub.blast$query[1] 
  pident=sub.blast$pident[2] 
  if(is.na(pident)){pident=100}
  evalues=sub.blast$evalue
  if(length(evalues)==1){evalues=c(evalues,rep(100,4))}
  fam_count=as.numeric(source.count.summary[grepl(top_hit_fam,source.count.summary$slc_fam),"fam_count"])
  
  ##############################Iterative
  if(100 %in% sub.blast$pident){
    if(pident<20){
      next
    }
    else if(grepl('SLC_Unsorted',sub.blast$subject[1])){ ### if the gene is previously unsorted
      if(slc_tab[1]>=length(slc_tab)-1){ ## all hits but one are in concordance
        slc.total[[i]]=data.table(geneid=top_hit_gene,family=top_hit_fam)
       next
      }else{
        slc.total[[i]]=data.table(geneid=top_hit_gene,family='SLC_Unsorted_')
      }
    }else if(grepl('SLC_',top_hit_fam)){ ## if the gene was previously an SLC (not unsorted)
      slc.total[[i]]=data.table(geneid=top_hit_gene,family=top_hit_fam)
      next
    }else{ ## if the gene was not previously called
      sub.blast=sub.blast[-1,]
      if(dim(sub.blast)[1]==0){next} ## if we have eliminated the rest of the data frame
      slc_fams=sub.blast$subject
      slc_tab=table(slc_fams) %>% sort(decreasing = T)
      top_hit_fam=names(slc_tab)[1]
      top_hit_gene=sub.blast$query[1]
      evalues=sub.blast$evalue
      if(length(evalues)==1){evalues=c(evalues,rep(100,4))}
      fam_count=as.numeric(source.count.summary[grepl(top_hit_fam,source.count.summary$slc_fam),"fam_count"])
    }
  }
  
  
  
  ###exceptions
  if(is.nan(evalues[2]/evalues[1])){evalues[2]=1e-100;evalues[1]=1e-100} ## when both evalues are 0
  if(is.na(fam_count)){fam_count=100}

  ####################### Non iterative
  
  if(!grepl('SLC_',sub.blast$subject[1])){ ## filter out any that have top blast hit not an SLC
    filter.list[[i]]=data.table(geneid=top_hit_gene,family=top_hit_fam)
  }else if(length(which(slc_fams %in% top_hit_fam))==fam_count){ ## where all members of a family are present in the list
    slc.total[[i]]=data.table(geneid=top_hit_gene,family=top_hit_fam)
  #}else if(slc_tab[1]>=length(slc_tab)-1){ ## include any where all are from from one SLC family except 1
  }else if(slc_tab[1]>=nrow(sub.blast)-1){
    slc.total[[i]]=data.table(geneid=top_hit_gene,family=top_hit_fam)
  }else if((evalues[2]/evalues[1]> 1e30) & grepl('SLC_',top_hit_fam)){ ## Keep where top hit is SLC and overwhelmingly significant
    slc.total[[i]]=data.table(geneid=top_hit_gene,family=top_hit_fam)
  }else if(length(grep('SLC_',slc_fams))>=length(slc_tab)-1){ ## keep where 4 out of 5 are SLCs but can be any family
    slc.total[[i]]=data.table(geneid=top_hit_gene,family='SLC_Unsorted_')
  }else{
    filter.list[[i]]=data.table(geneid=top_hit_gene,family=top_hit_fam)
  }
}

filter.table=rbindlist(filter.list) 
slc.table=rbindlist(slc.total) %>% arrange(family) %>% data.table()
if(nrow(filter.table)>0){
  colnames(filter.table)=c('code','name') 
  filter.output=filter.table %>% merge(source.table,by='code',all=T)
  write.csv(filter.output,'../prelim_summary/filtered_out.csv')
}

if(nrow(slc.table)>0){
  colnames(slc.table)=c('code','name')
  cat(format_csv(slc.table))
}



