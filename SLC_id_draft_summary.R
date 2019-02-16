#!/usr/bin/Rscript

library(data.table)
library(dplyr)
library(tidyr)
library(readr)

#setwd('/home/shanedenecke/Documents/SLC_id/Drosophila_search/DROSOPHILA_Ha')
setwd('./recip_blast_FIND')
args = commandArgs(trailingOnly=TRUE)
#cat(args[1])
#args[1]="/home/shanedenecke/Documents/SLC_id/Drosophila_search/DROSOPHILA_Ha/SLC_dict.csv"

dict=fread(as.character(args[1]))
dict.summary=dict %>% separate(col=name,sep='_',remove=T,into=c('slc','fam','subfam')) %>% unite(col='slc_fam',slc,fam,sep="_") %>% group_by(slc_fam) %>% summarize(fam_count=length(slc_fam))
dict.summary$slc_fam=sapply(dict.summary$slc_fam, function(x) paste0(x,"_"))
slc.total=list()
not.found=c()
used.list=c()
for (i in list.files()){
  blast=fread(i,colClasses = c('character','character',rep('numeric',3)))
  target.family=paste0(unlist(strsplit(i,split='.',fixed=T))[1],"_")
  target.fam.size=dict.summary[dict.summary$slc_fam==target.family,'fam_count']

  for (j in unique(blast$V1)){
    
    if(j %in% used.list){
      next
    }
    
    sub=subset(blast,V1==j)[1:5,] %>% na.omit
    slc_fams=sapply(sub$V2, function(x) strsplit(x,split='_') %>% unlist %>% .[c(1,2)] %>% paste0(collapse = "_") %>% paste0("_"))
    slc.tab=table(slc_fams)
    com.name=names(sort(table(slc_fams),decreasing=TRUE)[1])
    
    if(length(unique(slc_fams))==1 & com.name==target.family){ ## all hits the same
      slc.total[[j]]=data.table(geneid=sub[1,'V1'],family=unique(slc_fams))
      used.list=c(used.list,j)
    }
    
    else if(!(FALSE %in% (unique(slc_fams[1:3]) %in% target.family))){ ## only two families detected
         ## three top evalue hits are target family
          slc.total[[j]]=data.table(geneid=sub[1,'V1'],family=unique(slc_fams[1])) 
          used.list=c(used.list,j)
    }
    
    else if(TRUE %in% grepl(target.family,slc_fams[1])){ ## if the target family is the most significant hit
      if(slc_fams[grepl(target.family,slc_fams)]==target.fam.size){ ## all members of the family are detected
        slc.total[[j]]=data.table(geneid=sub[1,'V1'],family=target.family)
        used.list=c(used.list,j)
      }else if(sub$V4[1]<1e-60){ # & (sub$V4[2]/sub$V4[1]>1e20) It is overwhelmingly significant
        slc.total[[j]]=data.table(geneid=sub[1,'V1'],family=target.family)
        used.list=c(used.list,j)
      }
      else if(slc_fams[sapply(slc_fams,function(x) grepl('SLC',x))] %>% length()==length(slc_fams)){ ##if all 5 members are SLCs but none of the above conditions are met
        slc.total[[j]]=data.table(geneid=sub[1,'V1'],family='SLC_X')
        used.list=c(used.list,j)
      }
    }
  }
}

a=rbindlist(slc.total) %>% arrange(family)
cat(format_csv(a))
#write.csv(a,file='./HMMsearch_SLC_table.csv',row.names = F)
