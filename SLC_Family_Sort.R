#!/usr/bin/Rscript

## R script to sort SLC hits. The input is a formatted blast output table. 
### Each gene is sorted into family or is discarded based on criteria similar to Hoglund et. al 2011

## 1st argument is an SLC dictionary (usually form the source database) which will be used against the reciprocal blast results

#rm(list=ls())

##Import libraries
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))

############################ Set WD and arguments

args = commandArgs(trailingOnly=TRUE)
##debug
#setwd(/home/shanedenecke/Documents/SLC_id/iterative_search/iterative_search_AcyPis')
#args[1]="/home/shanedenecke/Documents/SLC_id/iterative_database/iterative_database_AcyPis/SLC_source_dict.csv"

setwd('./recip_blast')


divNA=function(x,y){
  if(is.na(x) | is.na(y)){
    return(1)
  }else if(is.numeric(x) & is.numeric(y)){
    return(x/y)
  }
}

##format some tables. dict will be your starting SLC dict. dict.summary a counts table
dict=fread(as.character(args[1]))
dict.summary=dict %>% separate(col=name,sep='_',remove=T,into=c('slc','fam','subfam')) %>% unite(col='slc_fam',slc,fam,sep="_") %>% group_by(slc_fam) %>% summarize(fam_count=length(slc_fam))
dict.summary$slc_fam=sapply(dict.summary$slc_fam, function(x) paste0(x,"_"))


## initialize some empty lists
slc.total=list()
filter.list=list()
used.list=c()

## Now start this moster for loop that will iterate through all blast files

for (i in list.files()){ ### iterate through each blast output file
  
  ## import blast database and set target family
  
  blast=fread(i,colClasses = c('character','character',rep('numeric',3)),sep='\t',col.names = c('query','subject','pident','evalue','querycov')) ## raw blast output table
  target.family=unlist(strsplit(i,split='.',fixed=T))[1] ### target family== what family is to be searched
  target.fam.size=dict.summary[dict.summary$slc_fam==target.family,'fam_count'] ### How many members in the target family are in the soruce genome

  for (j in unique(blast$query)){ ## iterate over genes in blast output

    if(j %in% used.list){ ## if the gene has already been asigned you can skip it
      next
    }
    
    ## calculate some useful variables for the coming for loop
    gene_subset=subset(blast,query==j)### subset blast results for gene
    slc_fams=sapply(gene_subset$subject, function(x) strsplit(x,split='_') %>% unlist %>% .[c(1,2)] %>% paste0(collapse = "_") %>% paste0("_")) ## extract families and format with extra "_"
    if(length(slc_fams)==1){
      slc_fams=c(slc_fams,rep(slc_fams,4))
    }
    slc.tab=table(slc_fams) ## returns table of which SLC families are there
    com.name=names(sort(slc.tab,decreasing=TRUE)[1]) ## gets the most common SLC family
    pident=gene_subset$pident
    if(length(pident)==1){
      pident=c(pident,rep(100,4))
    }
    
    ####################### ITERATIVE SEARCH
    
    if(pident[1]==100 & slc_fams[2]==target.family & pident[2]<20){ ## for recursive search if not 20% identity
      filter.list[[j]]=data.table(geneid=gene_subset[1,'query'],family=target.family)
      next
    }

    
    if(pident[1]==100 & grepl(target.family,slc_fams[1])){
      slc_fams=sapply(gene_subset$subject, function(x) strsplit(x,split='_') %>% unlist %>% .[c(1,2)] %>% paste0(collapse = "_") %>% paste0("_"))[1] ## extract families and format with extra "_"
      slc.total[[j]]=data.table(geneid=gene_subset[1,'query'],family=target.family)
      used.list=c(used.list,j)
      next
    }else if(gene_subset$pident[1]==100 & grepl('Unsorted',gene_subset$subject[1])){ ## if there is a perfect blast hit but all other blast hits map to other family
      gene_subset=gene_subset %>% filter(pident!=100)
      slc_fams=sapply(gene_subset$subject, function(x) strsplit(x,split='_') %>% unlist %>% .[c(1,2)] %>% paste0(collapse = "_") %>% paste0("_"))[1] ## extract families and format with extra "_"
      if(length(slc_fams)==1){ ### if all 4 remaining hits map to the same family
        slc.total[[j]]=data.table(geneid=gene_subset[1,'query'],family=target.family)
        used.list=c(used.list,j)
        next
      }
    }else if(gene_subset$pident[1]==100){
      gene_subset=gene_subset %>% filter(pident!=100) %>% head(5)
      slc_fams=sapply(gene_subset$subject, function(x) strsplit(x,split='_') %>% unlist %>% .[c(1,2)] %>% paste0(collapse = "_") %>% paste0("_"))[1] ## extract families and format with extra "_"
      slc_fams=slc_fams %>% na.omit()
    }
      
    ###################### NORMAL SEARCH
    if (!grepl('_X|Unsorted',i)){
      gene_subset=gene_subset %>% filter(!grepl('_X|Unsorted',subject)) ### if it is not the SLCX case then filter out all the SLCX hits
      if(dim(gene_subset)[1]==0){
        next
      }
    }else if(grepl('_X',i)){
      gene_subset=gene_subset
    }
    
    if(dim(gene_subset)[1]==0){ ## catch any instances where we might have deleted the entire dataframe
      next
    }
    
   
    ##filter candidates out immediately
    
    if(slc_fams[1]!=target.family){ ## if the most significant hit is not from the target family
      next
    }
  
    #### Now sort remaining genes
    
    if(length(unique(slc_fams))==1 & com.name==target.family){ ## all hits the same and this corresponds to the target family
      slc.total[[j]]=data.table(geneid=gene_subset[1,'query'],family=target.family)
      used.list=c(used.list,j)
    }else if(!(FALSE %in% (unique(slc_fams[1:3]) %in% target.family))){ ## Three top evalue hits are in the target family
          slc.total[[j]]=data.table(geneid=gene_subset[1,'query'],family=target.family) 
          used.list=c(used.list,j)
    }else if(slc_fams[grepl(target.family,slc_fams)]==target.fam.size){ ## all members of the family are detected in the top5
      slc.total[[j]]=data.table(geneid=gene_subset[1,'query'],family=target.family)
      used.list=c(used.list,j)
      
    }else if(gene_subset$evalue[1]<1e-60 | divNA(gene_subset$evalue[2],gene_subset$evalue[1])>1e10){ # It is overwhelmingly significant
      slc.total[[j]]=data.table(geneid=gene_subset[1,'query'],family=target.family)
      used.list=c(used.list,j)
    }else if(length(slc.tab[names(slc.tab)==target.family][slc.tab>3])>0){ ## 4 out of the top 5 blast hits are the target family
      slc.total[[j]]=data.table(geneid=gene_subset[1,'query'],family=target.family)
      used.list=c(used.list,j)
    }else if(slc_fams[sapply(slc_fams,function(x) grepl('SLC',x))] %>% length()==length(slc_fams)){ ##if all 5 members are SLCs but none of the above conditions are met
        slc.total[[j]]=data.table(geneid=gene_subset[1,'query'],family='SLC_Unsorted')
        used.list=c(used.list,j)
    }
  }
}

a=rbindlist(slc.total) %>% data.table() 
colnames(a)=c('geneid','family')

a$family=unlist(a$family)
final=a %>% arrange(family) %>% unique.data.frame()
b=rbindlist(filter.list) %>% arrange(family) %>% unique.data.frame()
write.csv(b,'../prelim_summary/filtered_out.csv')
## sometimes things are detected in initial search but not in iterative (e.g. SLC9 tr|H9JQM4|H9JQM4_BOMMO in bombyx mori just below threshold'
## add in any mismatches from dictionary if recursive


#filter(a,family=='SLC_2_')
cat(format_csv(final))
#write.csv(a,file='./HMMsearch_SLC_table.csv',row.names = F)

