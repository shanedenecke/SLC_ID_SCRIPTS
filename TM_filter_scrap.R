

##########################################################################



noname=arth.hmm[!grep('___',name)]
noname.index=which(!grepl('___',arth.hmm$name))






fasta.name=paste(all$Species_name,all$abbreviation,all$code,all$name,sep='___')
rename.dict=data.table(name=fasta.name,code=all$code)
all$Species_name=gsub(' ','_',all$Species_name)

fwrite(rename.dict,'./TMHMM_filter/Rename_SLC_dict.csv')
fwrite(all,'./TMHMM_filter/Full_unfiltered_SLC_table.csv')
writeLines(rename.dict$code,'./TMHMM_filter/slc_unfiltered_codes.txt')
##################### FASTA

system('
       cd /data2/shane/Documents/SLC_id
       cat ./proteomes/* ./general_reference/model_proteomes/*.faa > ./TMHMM_filter/all_proteomes.faa
       /data2/shane/Applications/custom/unigene_fa_sub.sh ./TMHMM_filter/all_proteomes.faa ./TMHMM_filter/slc_unfiltered_codes.txt > ./TMHMM_filter/SLC_unfiltered_all_raw.faa
       /data2/shane/Applications/custom/fasta_rename.py ./TMHMM_filter/SLC_unfiltered_all_raw.faa ./TMHMM_filter/Rename_SLC_dict.csv > ./TMHMM_filter/Renamed_unfiltered_SLC.faa
       
       /home/pioannidis/Programs/tmhmm-2.0c/bin/tmhmm  $1 | grep "Number of predicted" |\
       perl -pe "s/..(.+) Number of predicted TMHs:\s+(\S+)/$1\t$2/g"" | awk "{ if ($2 >="$2") { print } }"
       
       /data2/shane/Applications/custom/tmhmm_filter.sh ./TMHMM_filter/Renamed_unfiltered_SLC.faa 0 > ./TMHMM_filter/SLC_TMHMM_scores.txt
       ')

SLC_named_fasta=read.fasta('./TMHMM_filter/Renamed_unfiltered_SLC.faa',seqtype='AA',as.string=T,set.attributes = F,strip.desc=T)


####################FILTER BASED ON TMHMM VALUES
arth.hmm=fread('./TMHMM_filter/SLC_TMHMM_scores.txt') 
colnames(arth.hmm)=c('name','tm_count')

noname=arth.hmm[!grep('___',name)]
noname.index=which(!grepl('___',arth.hmm$name))

clean.arth.hmm=arth.hmm %>% separate(sep = '___',col=name,into=c('Species_name','abbreviation','code','slc_code')) %>%
  mutate(family=dash.remove(slc_code)) %>% select(code,tm_count) %>% data.table()

merge(all,clean.arth.hmm,by='code')


for(i in )
  
  
  
  
  sep_sub=function(x,se='_',collap='_',start=1,end=1){
    x1=unlist(strsplit(x,se,fixed=T))
    paste(x1[start:end],collapse = collap)
  }



remove_sep=function(x){
  az=unlist(strsplit(x,'|',fixed=T))
  bz=az[length(az)]
  return(bz)
}


full.l=list()
count.l=list()
non.match=list()
iterframe=distinct(all,family,abbreviation)

for(i in 1:nrow(iterframe)){
  f=paste0(iterframe[i]$family)
  a=iterframe[i]$abbreviation
  
  tm.min=model.hmm.key[family==f]$minimum
  
  sub.raw=all[family==f & abbreviation==a]
  
  search.code=sapply(sub.raw$code,remove_sep) %>% paste(collapse = '|')
  
  sub.tm=arth.hmm[grepl(search.code, arth.hmm$V1) & grepl(a,arth.hmm$V1)]
  
  
  for(j in nrow(sub.raw)_{
    cod=sub.raw[j]$code
    tms=sub.tm[grepl(code,V1)]
    
  }
  
  
  
  if(nrow(sub.raw!=nrow(sub.tm))){
    non.match[[i]]=sub.raw
  }else{
    
    fil=sub.raw[sub.tm$V2>tm.min]
    
    full.l[[i]]=fil
    cou=fil %>% group_by(family) %>% summarise(count=length(family)) 
    colnames(cou)[2]=a
    count.l[[paste0(a,'_',f,'_')]]=cou
  }
}
full.dict.table=rbindlist(full.l)

##### CALCULATE COUNTS TABLE
count.l2=list()
for(i in unique(iterframe$abbreviation)){
  spec=count.l[which(grepl(i,names(count.l)))]
  count.l2[[i]]=rbindlist(spec)
}

## combine list into count.summary variable and transpose to revers rows/cols
count.summary=Reduce(function(x, y) merge(x, y, ,by='family',all=TRUE), count.l2)  %>% 
  shane.transpose() %>%
  mutate_at(vars(matches("SLC")), funs(as.numeric(as.character(.)))) %>%
  data.table()
count.summary[is.na(count.summary)]=0

## calculate totals
count.summary$slc_total=rowSums(count.summary[,2:67])
colnames(count.summary)[1]='abbreviation'
count.summary=select(count.summary,abbreviation,everything())



##write totals to file
fwrite(count.summary,'./Final_raw_outputs/count_summary.csv')
fwrite(full.dict.table,'./Final_raw_outputs/Full_dict_table.csv')


##################### FASTA

system('
       cd /data2/shane/Documents/SLC_id
       /data2/shane/Applications/custom/unigene_fa_sub.sh ./TMHMM_filter/all_proteomes.faa ./TMHMM_filter/slc_unfiltered_codes.txt > ./TMHMM_filter/SLC_unfiltered_all_raw.faa
       /data2/shane/Applications/custom/fasta_rename.py ./TMHMM_filter/SLC_unfiltered_all_raw.faa ./TMHMM_filter/Rename_SLC_dict.csv > ./TMHMM_filter/Renamed_unfiltered_SLC.faa
       
       /data2/shane/Applications/custom/tmhmm_filter.sh ./TMHMM_filter/Renamed_unfiltered_SLC.faa 0 > ./TMHMM_filter/SLC_TMHMM_scores.txt
       ')



















fwrite(tmhmm.filtered.out,'./SLC_family_counts/TMM_filtered_out.csv')








count.l=list()
full.l=list()
remove.l=list()
### FILTER TM values 
##Get all counts data into master table. Rows SLC families. Columns species
for (i in list.files('./final_SLC_dicts/')){
  slc_table=fread(paste0('./final_SLC_dicts/',i))
  #slc_table=fread('/home/shanedenecke/Dropbox/wp7_prodrug/SLC_id/HomSap_SLC_dict.csv')
  abbreviation=gsub('Final','',unlist(strsplit(i,'_'))[1])
  slc_table$family=sapply(slc_table$name,dash.remove)
  
  tm.l=list()
  tm.filter.out=list()
  arth.index=arth.hmm[grepl(abbreviation,arth.hmm$V1)]
  #TMHMM filter
  for(j in 1:nrow(slc_table)){
    code=slc_table[j]$code
    tfam=slc_table[j]$family
    human.min=model.hmm.key[family==tfam]$minimum
    test.tm=arth.index[grepl(code,arth.index$V1,fixed=T)]$V2
    if(length(test.tm)==0){
      tm.l[[j]]=slc_table[j]
    }else if(test.tm>=human.min){
      tm.l[[j]]=slc_table[j]
    }else{
      tm.filter.out[[j]]=slc_table[j]
    }
  }
  tm.filter=rbindlist(tm.l)
  
  remove.l[[i]]=rbindlist(tm.filter.out)
  full.l[[i]]=tm.filter
  
  with.count=tm.filter %>% group_by(family) %>% summarise(count=length(family))
  colnames(with.count)[2]=abbreviation
  count.l[[i]]=with.count
}

tmhmm.filtered.out=rbindlist(remove.l)
full.dict.table=rbindlist(full.l)  

## combine list into count.summary variable and transpose to revers rows/cols
count.summary=Reduce(function(x, y) merge(x, y, ,by='family',all=TRUE), count.l)  %>% 
  shane.transpose() %>%
  mutate_at(vars(matches("SLC")), funs(as.numeric(as.character(.)))) %>%
  data.table()
count.summary[is.na(count.summary)]=0

## calculate totals
count.summary$slc_total=rowSums(count.summary[,2:67])
colnames(count.summary)[1]='abbreviation'
count.summary=select(count.summary,abbreviation,everything())

##write totals to file
fwrite(count.summary,'./SLC_family_counts/count_summary.csv')
fwrite(tmhmm.filtered.out,'./SLC_family_counts/TMM_filtered_out.csv')
fwrite(full.dict.table,'./SLC_family_counts/Full_dict_table.csv')







############### Drosophila verify

### Drosophila verify
#dros=fread('/data2/shane/Documents/SLC_id/dros_tmhmm.txt')
#dros$family=sapply(dros$V1,dash.remove)
#slc_table=dros
#tm.l=list()
#tm.filter.out=list()
#for(j in 1:nrow(slc_table)){
#  test.tm=slc_table[j]$V2
#  tfam=slc_table[j]$family
#  human.min=human.hmm.key[family==tfam]$minimum

#  if(length(test.tm)==0){
#    tm.l[[j]]=slc_table[j]
#  }else if(test.tm>=human.min){
#    tm.l[[j]]=slc_table[j]
#  }else{
#    tm.filter.out[[j]]=slc_table[j]
#  }
#}
#rbindlist(tm.filter.out)